from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic import HttpUrl
from enum import Enum, unique
from typing import NewType, NamedTuple, Any, Callable, TypeVar
from snakemake.io import expand, InputFiles  # type: ignore

X = TypeVar("X")
Y = TypeVar("Y")


def fmap_maybe(f: Callable[[X], Y], x: X | None) -> None | Y:
    pass


BuildKey = NewType("BuildKey", str)
RefKey = NewType("RefKey", str)


@unique
class ChrIndex(Enum):
    _ignore_ = "ChrIndex i"
    ChrIndex = vars()
    for i in range(23):
        ChrIndex[f"CHR{i}"] = i
    CHRX = 23
    CHRY = 24

    def __init__(self, i: int) -> None:
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def chr_name_full(self, prefix: str) -> str:
        return f"{prefix}{self.chr_name}"


class StratOutputs(NamedTuple):
    low_complexity: InputFiles
    xy_sex: InputFiles
    xy_auto: InputFiles
    map: InputFiles


class BaseModel(BaseModel_):
    class Config:
        frozen = True
        extra = "forbid"


class Paths(BaseModel):
    resources: Path
    results: Path


class Tools(BaseModel):
    repseq: HttpUrl
    gemlib: HttpUrl


class BedFile(BaseModel):
    url: HttpUrl
    chr_prefix: str = "chr"


class LowComplexity(BaseModel):
    rmsk: BedFile
    simreps: BedFile


class XY(BaseModel):
    x_features: BedFile
    y_features: BedFile
    x_par: BedFile


class SegDups(BaseModel):
    self_chain: BedFile
    self_chain_link: BedFile


class Include(BaseModel):
    low_complexity: bool
    xy: bool
    map: bool


class Build(BaseModel):
    chr_filter: set[ChrIndex]
    include: Include


class RefFile(BaseModel):
    url: HttpUrl
    chr_prefix: str


class Stratification(BaseModel):
    ref: RefFile
    gap: BedFile | None
    low_complexity: LowComplexity | None
    xy: XY | None
    segdups: SegDups | None
    builds: dict[BuildKey, Build]


class GiabStrats(BaseModel):
    paths: Paths
    tools: Tools
    stratifications: dict[str, Stratification]

    # hack to make rmd scripts work with this (note this will totally kill
    # the config as it passes into an rmd script)
    def items(self) -> Any:
        return {}.items()

    def refkey_to_strat(self, k: RefKey) -> Stratification:
        return self.stratifications[k]

    def refkey_to_ref_url(self, k: RefKey) -> str:
        return self.stratifications[k].ref.url

    def refkey_to_gap_url(self, k: RefKey) -> str | None:
        return fmap_maybe(lambda x: x.url, self.stratifications[k].gap)

    def refkey_to_x_features_url(self, k: RefKey) -> str | None:
        return fmap_maybe(lambda x: x.x_features.url, self.stratifications[k].xy)

    def refkey_to_y_features_url(self, k: RefKey) -> str | None:
        return fmap_maybe(lambda x: x.y_features.url, self.stratifications[k].xy)

    def refkey_to_x_par_url(self, k: RefKey) -> str | None:
        return fmap_maybe(lambda x: x.x_par.url, self.stratifications[k].xy)

    def refkey_to_simreps_url(self, k: RefKey) -> str | None:
        return fmap_maybe(
            lambda x: x.simreps.url, self.stratifications[k].low_complexity
        )

    def refkey_to_rmsk_url(self, k: RefKey) -> str | None:
        return fmap_maybe(lambda x: x.rmsk.url, self.stratifications[k].low_complexity)

    def refkey_to_self_chain_url(self, k: RefKey) -> str | None:
        return fmap_maybe(lambda x: x.self_chain.url, self.stratifications[k].segdups)

    def refkey_to_self_chain_link_url(self, k: RefKey) -> str | None:
        return fmap_maybe(
            lambda x: x.self_chain_link.url, self.stratifications[k].segdups
        )

    def buildkey_to_build(self, rk: RefKey, bk: BuildKey) -> Build:
        return self.stratifications[rk].builds[bk]

    def buildkey_to_chr_indices(self, rk: RefKey, bk: BuildKey) -> set[ChrIndex]:
        cs = self.stratifications[rk].builds[bk].chr_filter
        return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    def buildkey_to_chr_names(self, rk: RefKey, bk: BuildKey) -> list[str]:
        # TODO don't hardcode this in the future
        prefix = "chr"
        return [i.chr_name_full(prefix) for i in self.buildkey_to_chr_indices(rk, bk)]

    def buildkey_to_chr_pattern(self, rk: RefKey, bk: BuildKey) -> str:
        return "\\|".join(self.buildkey_to_chr_names(rk, bk))

    def strat_targets(
        self, rk: RefKey, bk: BuildKey, s: StratOutputs
    ) -> list[InputFiles]:
        inc = self.buildkey_to_build(rk, bk).include
        cis = self.buildkey_to_chr_indices(rk, bk)

        # XY is special since we should not include the XY strats if we don't
        # also have X and Y in the filter. Likewise we should not include an
        # autosomes if they are not in the filter
        want_xy = ChrIndex.CHRX in cis and ChrIndex.CHRY in cis and inc.xy
        want_autosomes = len(cis - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

        all_targets = [
            (s.low_complexity, inc.low_complexity),
            (s.xy_sex, want_xy),
            (s.xy_auto, want_autosomes),
            (s.map, inc.map),
        ]

        return [
            p
            for ps in [
                expand(target, allow_missing=True, ref_key=rk, build_key=bk)
                for target, wants in all_targets
                if wants
            ]
            for p in ps
        ]
