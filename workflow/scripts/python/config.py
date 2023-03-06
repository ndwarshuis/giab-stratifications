from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic import HttpUrl
from enum import Enum, auto
from typing import NewType, NamedTuple, Any
from snakemake.io import expand, InputFiles  # type: ignore

BuildKey = NewType("BuildKey", str)
RefKey = NewType("RefKey", str)


class ChrIndex(Enum):
    CHR1 = auto()
    CHR2 = auto()
    CHR3 = auto()
    CHR4 = auto()
    CHR5 = auto()
    CHR6 = auto()
    CHR7 = auto()
    CHR8 = auto()
    CHR9 = auto()
    CHR10 = auto()
    CHR11 = auto()
    CHR12 = auto()
    CHR13 = auto()
    CHR14 = auto()
    CHR15 = auto()
    CHR16 = auto()
    CHR17 = auto()
    CHR18 = auto()
    CHR19 = auto()
    CHR20 = auto()
    CHR21 = auto()
    CHR22 = auto()
    CHRX = auto()
    CHRY = auto()

    @property
    def chr_name(self) -> str:
        if self.value == self.CHRY:
            return "Y"
        elif self.value == self.CHRX:
            return "X"
        return str(self.value)

    def chr_name_full(self, prefix: str) -> str:
        return f"{prefix}{self.chr_name}"


class StratOutputs(NamedTuple):
    low_complexity: InputFiles
    xy_sex: InputFiles
    xy_auto: InputFiles
    map: InputFiles


class BaseModel(BaseModel_):
    pass


class Paths(BaseModel):
    resources: Path
    results: Path


class Tools(BaseModel):
    repseq: HttpUrl
    gemlib: HttpUrl


class LowComplexity(BaseModel):
    rmsk_url: HttpUrl
    simreps_url: HttpUrl


class XY(BaseModel):
    x_features: HttpUrl
    y_features: HttpUrl
    x_par: HttpUrl


class SegDups(BaseModel):
    self_chain: HttpUrl
    self_chain_link: HttpUrl


class Include(BaseModel):
    low_complexity: bool
    xy: bool
    map: bool


class Build(BaseModel):
    chr_filter: set[ChrIndex]
    include: Include


class Stratification(BaseModel):
    ref: HttpUrl
    gap: HttpUrl | None
    chr_prefix: str
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
        return self.stratifications[k].ref

    def refkey_to_gap_url(self, k: RefKey) -> str | None:
        return self.stratifications[k].gap

    def refkey_to_x_features_url(self, k: RefKey) -> str:
        return self.stratifications[k].xy.x_features

    def refkey_to_y_features_url(self, k: RefKey) -> str:
        return self.stratifications[k].xy.y_features

    def refkey_to_x_par_url(self, k: RefKey) -> str:
        return self.stratifications[k].xy.x_par

    def refkey_to_simreps_url(self, k: RefKey) -> str:
        return self.stratifications[k].low_complexity.simreps_url

    def refkey_to_rmsk_url(self, k: RefKey) -> str:
        return self.stratifications[k].low_complexity.rmsk_url

    def refkey_to_self_chain_url(self, k: RefKey) -> str:
        return self.stratifications[k].segdups.self_chain

    def refkey_to_self_chain_link_url(self, k: RefKey) -> str:
        return self.stratifications[k].segdups.self_chain_link

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
