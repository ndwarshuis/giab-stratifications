from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic import HttpUrl, FilePath
from enum import Enum, unique
from typing import NewType, NamedTuple, Any, Callable, TypeVar
from snakemake.io import expand, InputFiles  # type: ignore
from more_itertools import flatten

X = TypeVar("X")
Y = TypeVar("Y")


BuildKey = NewType("BuildKey", str)
RefKey = NewType("RefKey", str)


def fmap_maybe(f: Callable[[X], Y], x: X | None) -> None | Y:
    return None if x is None else f(x)


def _flatten_targets(
    all_targets: list[tuple[InputFiles, bool]],
    rk: RefKey,
    bk: BuildKey,
) -> InputFiles:
    return [
        *flatten(
            expand(target, allow_missing=True, ref_key=rk, build_key=bk)
            for target, wants in all_targets
            if wants
        )
    ]


class ZipFmt(Enum):
    NOZIP = "nozip"
    GZIP = "gzip"
    BGZIP = "bgzip"


class XYFeature(Enum):
    XTR = "XTR"
    Ampliconic = "Ampliconic"


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


class SatelliteOutputs(NamedTuple):
    rmsk: InputFiles
    censat: InputFiles


class LowComplexityOutputs(NamedTuple):
    uniform_repeats: InputFiles
    all_repeats: InputFiles


class StratOutputs(NamedTuple):
    low_complexity: LowComplexityOutputs
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


# TODO non-negative ints which cannot equal each other
class BedColumns(BaseModel):
    chr: int = 1
    start: int = 2
    end: int = 3


class FileSrc(BaseModel):
    filepath: FilePath


class HttpSrc(BaseModel):
    url: HttpUrl
    zipfmt: ZipFmt = ZipFmt.BGZIP


class BedFile(BaseModel):
    src: FileSrc | HttpSrc
    chr_prefix: str = "chr"
    bed_cols: BedColumns = BedColumns()


class RMSKFile(BedFile):
    class_col: int


class LowComplexity(BaseModel):
    rmsk: RMSKFile
    simreps: BedFile
    satellites: BedFile | None


class XY(BaseModel):
    x_features: BedFile
    y_features: BedFile
    features: set[XYFeature]
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
    low_complexity: LowComplexity
    xy: XY
    segdups: SegDups
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

    # def refkey_to_gap_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(lambda x: x.url, self.stratifications[k].gap)

    # def refkey_to_x_features_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(lambda x: x.x_features.url, self.stratifications[k].xy)

    # def refkey_to_y_features_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(lambda x: x.y_features.url, self.stratifications[k].xy)

    # def refkey_to_x_par_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(lambda x: x.x_par.url, self.stratifications[k].xy)

    # def refkey_to_simreps_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(
    #         lambda x: x.simreps.url, self.stratifications[k].low_complexity
    #     )

    # def refkey_to_rmsk_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(lambda x: x.rmsk.url, self.stratifications[k].low_complexity)

    # def refkey_to_self_chain_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(lambda x: x.self_chain.url, self.stratifications[k].segdups)

    # def refkey_to_self_chain_link_url(self, k: RefKey) -> str | None:
    #     return fmap_maybe(
    #         lambda x: x.self_chain_link.url, self.stratifications[k].segdups
    #     )

    def refkey_to_final_chr_prefix(self, k: RefKey) -> str:
        return self.stratifications[k].ref.chr_prefix

    def refkey_to_input_chr_prefix(
        self,
        f: Callable[[Stratification], BedFile],
        k: RefKey,
    ) -> str:
        return f(self.stratifications[k]).chr_prefix

    def buildkey_to_init_chr_mapping(
        self,
        f: Callable[[Stratification], BedFile],
        rk: RefKey,
        bk: BuildKey,
    ) -> dict[str, int]:
        p = self.refkey_to_input_chr_prefix(f, rk)
        cs = self.buildkey_to_chr_indices(rk, bk)
        return {c.chr_name_full(p): c.value for c in cs}

    def buildkey_to_final_chr_mapping(self, rk: RefKey, bk: BuildKey) -> dict[int, str]:
        p = self.refkey_to_final_chr_prefix(rk)
        cs = self.buildkey_to_chr_indices(rk, bk)
        return {c.value: c.chr_name_full(p) for c in cs}

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

    def satellite_targets(self, rk: RefKey, s: SatelliteOutputs) -> InputFiles:
        return (
            s.censat
            if self.stratifications[rk].low_complexity.satellites is not None
            else s.rmsk
        )

    def low_complexity_targets(
        self, rk: RefKey, bk: BuildKey, s: LowComplexityOutputs
    ) -> InputFiles:
        if self.buildkey_to_build(rk, bk).include.low_complexity:
            r = self.stratifications[rk].low_complexity
            all_tgts = [
                (s.uniform_repeats, True),
                (s.all_repeats, r.rmsk is not None and r.simreps is not None),
            ]
            return _flatten_targets(all_tgts, rk, bk)
        else:
            return []

    def strat_targets(self, rk: RefKey, bk: BuildKey, s: StratOutputs) -> InputFiles:
        inc = self.buildkey_to_build(rk, bk).include
        cis = self.buildkey_to_chr_indices(rk, bk)

        # XY is special since we should not include the XY strats if we don't
        # also have X and Y in the filter. Likewise we should not include an
        # autosomes if they are not in the filter
        want_xy = ChrIndex.CHRX in cis and ChrIndex.CHRY in cis and inc.xy
        want_autosomes = len(cis - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

        all_targets = [
            (s.xy_sex, want_xy),
            (s.xy_auto, want_autosomes),
            (s.map, inc.map),
        ]

        return _flatten_targets(all_targets, rk, bk) + self.low_complexity_targets(
            rk, bk, s.low_complexity
        )
