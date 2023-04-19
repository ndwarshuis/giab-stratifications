from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic import validator, HttpUrl, FilePath, NonNegativeInt
from enum import Enum, unique
from typing import NewType, Any, Callable, TypeVar, Type, NamedTuple
from typing_extensions import Self
from more_itertools import unzip, unique_everseen
from common.io import is_gzip, is_bgzip

X = TypeVar("X")
Y = TypeVar("Y")


BuildKey = NewType("BuildKey", str)
RefKey = NewType("RefKey", str)


def fmap_maybe(f: Callable[[X], Y], x: X | None) -> None | Y:
    return None if x is None else f(x)


@unique
class ChrIndex(Enum):
    """Represents a valid chromosome index.

    Chromosomes are numbered by integers 1-24 (23 and 24 being X and Y
    respectively). These integers reflect the sort order in output bed files.
    """

    _ignore_ = "ChrIndex i"
    ChrIndex = vars()
    for i in range(1, 23):
        ChrIndex[f"CHR{i}"] = i
    CHRX = 23
    CHRY = 24

    @classmethod
    def from_name(cls, n: str) -> Self:
        try:
            return next(i for i in cls if i.chr_name == n)
        except StopIteration:
            raise ValueError(f"could make chr index from name '{n}'")

    def __init__(self, i: int) -> None:
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def chr_name_full(self, prefix: str) -> str:
        return f"{prefix}{self.chr_name}"

    def matches(self, prefix: str, query: str) -> bool:
        return query.startswith(self.chr_name_full(prefix))


def chr_is_unknown(prefix: str, query: str) -> bool:
    return query.startswith(f"{prefix}Un")


class ChrConversion(NamedTuple):
    """Data to filter, sort, and standardize chromosome names.

    Members:
    fromPrefix - chr prefix of input file (usually a bed file)
    toPrefix - chr prefix of output stratification file (should match ref)
    indices - the desired chromosome indices to keep
    """

    fromPrefix: str
    toPrefix: str
    indices: set[ChrIndex]


# For instances where we simply need to sort all the chromosomes and they are
# already named appropriately
def fullset_conv(prefix: str) -> ChrConversion:
    return ChrConversion(prefix, prefix, set([i for i in ChrIndex]))


def conversion_to_init_mapper(c: ChrConversion) -> dict[str, int]:
    return {i.chr_name_full(c.fromPrefix): i.value for i in c.indices}


def conversion_to_final_mapper(c: ChrConversion) -> dict[int, str]:
    return {i.value: i.chr_name_full(c.toPrefix) for i in c.indices}


class BaseModel(BaseModel_):
    class Config:
        frozen = True
        extra = "forbid"


class Paths(BaseModel):
    """Local build paths for snakemake."""

    resources: Path
    results: Path


class Tools(BaseModel):
    """Urls for tools to download/build/use in the pipeline."""

    repseq: HttpUrl
    gemlib: HttpUrl


# TODO non-negative ints which cannot equal each other
class BedColumns(BaseModel):
    """Denotes coordinate columns in a bed file (0-indexed)."""

    chr: NonNegativeInt = 0
    start: NonNegativeInt = 1
    end: NonNegativeInt = 2

    @validator("start")
    def start_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, int],
    ) -> NonNegativeInt:
        try:
            assert values["chr"] != v, "Bed columns must be different"
        except KeyError:
            pass
        return v

    @validator("end")
    def end_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, int],
    ) -> NonNegativeInt:
        try:
            assert (
                values["chr"] != v and values["start"] != v
            ), "Bed columns must be different"
        except KeyError:
            pass
        return v

    def assert_different(self, x: int) -> None:
        assert (
            self.chr != x and self.start != x and self.end != x
        ), "Column must be different index"

    @property
    def columns(self) -> dict[int, Type[int | str]]:
        return {self.chr: str, self.start: int, self.end: int}


class HashedSrc_(BaseModel):
    md5: str | None = None


class FileSrc_(HashedSrc_):
    """Base class for local src files."""

    filepath: FilePath


class BedFileSrc(FileSrc_):
    """Filepath for bedfile."""

    @validator("filepath")
    def is_gzip(cls, v: FilePath) -> FilePath:
        assert is_gzip(v), "not in gzip format"
        return v


class RefFileSrc(FileSrc_):
    """Filepath for reference."""

    @validator("filepath")
    def path_is_bgzip(cls, v: FilePath) -> FilePath:
        assert is_bgzip(v), "not in bgzip format"
        return v


class HttpSrc_(HashedSrc_):
    """Base class for downloaded src files."""

    url: HttpUrl


class BedHttpSrc(HttpSrc_):
    """Url for bed file"""

    pass


class RefHttpSrc(HttpSrc_):
    """Url for reference"""

    pass


RefSrc = RefFileSrc | RefHttpSrc

# TODO this is for more than just "bed files" (right now it basically means "a
# file that is either not zipped or gzipped but not bgzipped")
BedSrc = BedFileSrc | BedHttpSrc


class BedFileParams(BaseModel):
    """Parameters decribing how to parse a bed-like file.

    Members:
    chr_prefix - the prefix on the chromosomes
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    """

    chr_prefix: str = "chr"
    bed_cols: BedColumns = BedColumns()
    skip_lines: NonNegativeInt = 0
    sep: str = "\t"


class BedFile(BaseModel):
    """Inport specs for a bed-like file.

    Members:
    src - how to get the file
    chr_prefix - the prefix on the chromosomes
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    """

    src: BedSrc
    params: BedFileParams = BedFileParams()


class RMSKFile(BedFile):
    """Input file for repeat masker stratification."""

    class_col: NonNegativeInt

    @validator("class_col")
    def end_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, BedColumns],
    ) -> NonNegativeInt:
        try:
            values["bed_cols"].assert_different(v)
        except KeyError:
            pass
        return v


class LowComplexity(BaseModel):
    """Configuration for low complexity stratification."""

    rmsk: RMSKFile
    simreps: BedFile
    satellites: BedFile | None


class XYFile(BedFile):
    """Bed file input for XY features."""

    level_col: NonNegativeInt

    @validator("level_col")
    def level_different(
        cls,
        v: NonNegativeInt,
        values: dict[str, BedColumns],
    ) -> NonNegativeInt:
        try:
            values["bed_cols"].assert_different(v)
        except KeyError:
            pass
        return v


class XYFeatures(BaseModel):
    """Configuration for XY features stratifications."""

    x_bed: XYFile
    y_bed: XYFile
    ampliconic: bool
    xtr: bool


class XYPar(BaseModel):
    """Regions for the PARs on the X/Y chromosomes."""

    start: tuple[NonNegativeInt, NonNegativeInt]
    end: tuple[NonNegativeInt, NonNegativeInt]

    @validator("start", "end")
    def positive_region(cls, v: tuple[int, int]) -> tuple[int, int]:
        assert v[1] > v[0], "End must be greater than start"
        return v

    def fmt(self, i: ChrIndex, prefix: str) -> str:
        # TODO this smells like something I'll be doing alot
        c = i.chr_name_full(prefix)
        return "\n".join(
            [
                f"{c}\t{self.start[0]}\t{self.start[1]}",
                f"{c}\t{self.end[0]}\t{self.end[1]}",
            ]
        )


class XY(BaseModel):
    """Configuration for the XY stratification."""

    features: XYFeatures
    x_par: XYPar
    y_par: XYPar

    def fmt_x_par(self, prefix: str) -> str:
        return self.x_par.fmt(ChrIndex.CHRX, prefix)

    def fmt_y_par(self, prefix: str) -> str:
        return self.y_par.fmt(ChrIndex.CHRY, prefix)


class SegDups(BaseModel):
    """Configuration for Segdup stratifications."""

    superdups: BedFile | None


class LowMapParams(BaseModel):
    """Parameters for a single mappability bed file."""

    length: NonNegativeInt
    mismatches: NonNegativeInt
    indels: NonNegativeInt


class Include(BaseModel):
    """Flags to control which stratification levels are included."""

    low_complexity: bool
    xy: bool
    gc: bool
    functional: bool
    segdups: bool
    union: bool
    mappability: set[LowMapParams]


class Build(BaseModel):
    """Spec for a stratification build."""

    chr_filter: set[ChrIndex]
    include: Include


class RefFile(BaseModel):
    """Specification for reference file."""

    src: RefSrc
    chr_prefix: str


class Functional(BaseModel):
    """Configuration for Functional stratifications."""

    ftbl_src: BedSrc
    gff_src: BedSrc


class Stratification(BaseModel):
    """Configuration for stratifications for a given reference."""

    ref: RefFile
    gap: BedFile | None
    low_complexity: LowComplexity
    xy: XY
    segdups: SegDups
    functional: Functional
    builds: dict[BuildKey, Build]


class GiabStrats(BaseModel):
    """Top level stratification object."""

    paths: Paths
    tools: Tools
    stratifications: dict[RefKey, Stratification]

    # hack to make rmd scripts work with this (note this will totally kill
    # the config as it passes into an rmd script)
    def items(self) -> Any:
        return {}.items()

    # file paths

    @property
    def resources_dir(self) -> Path:
        return self.paths.resources

    @property
    def tools_src_dir(self) -> Path:
        return self.resources_dir / "tools"

    @property
    def ref_src_dir(self) -> Path:
        return self.paths.resources / "{ref_key}"

    @property
    def results_dir(self) -> Path:
        return self.paths.results

    @property
    def tools_make_dir(self) -> Path:
        return self.results_dir / "tools" / "make"

    @property
    def tools_bin_dir(self) -> Path:
        return self.results_dir / "tools" / "bin"

    @property
    def final_root_dir(self) -> Path:
        return self.results_dir / "final"

    @property
    def final_build_dir(self) -> Path:
        return self.results_dir / "final" / "{ref_key}@{build_key}"

    @property
    def intermediate_root_dir(self) -> Path:
        return self.results_dir / "intermediates"

    @property
    def intermediate_build_dir(self) -> Path:
        return self.intermediate_root_dir / "{ref_key}@{build_key}"

    @property
    def intermediate_ref_dir(self) -> Path:
        return self.intermediate_root_dir / "ref"

    @property
    def log_root_dir(self) -> Path:
        return self.results_dir / "log"

    @property
    def bench_root_dir(self) -> Path:
        return self.results_dir / "bench"

    @property
    def log_src_dir(self) -> Path:
        return self.log_root_dir / "{ref_key}"

    @property
    def log_build_dir(self) -> Path:
        return self.log_root_dir / "{ref_key}@{build_key}"

    @property
    def bench_build_dir(self) -> Path:
        return self.bench_root_dir / "{ref_key}@{build_key}"

    def build_strat_path(self, level: str, name: str) -> Path:
        return self.final_build_dir / level / f"{{ref_key}}_{name}.bed.gz"

    # because smk doesn't check these for existence yet:
    # https://github.com/snakemake/snakemake/issues/1657
    def _workflow_path(self, components: list[str]) -> Path:
        p = Path("workflow", *components)
        # except that it doesn't work too well in subworkflows...
        # assert p.exists(), f"{p} does not exist"
        return p

    def env_path(self, envname: str) -> Path:
        return self._workflow_path(["envs", f"{envname}.yml"])

    def _scripts_dir(self, rest: list[str]) -> Path:
        return self._workflow_path(["scripts", *rest])

    def python_script(self, basename: str) -> Path:
        return self._scripts_dir(["python", basename])

    def rmd_script(self, basename: str) -> Path:
        return self._scripts_dir(["rmarkdown", basename])

    # general accessors

    def refkey_to_strat(self, k: RefKey) -> Stratification:
        return self.stratifications[k]

    def buildkey_to_build(self, rk: RefKey, bk: BuildKey) -> Build:
        return self.stratifications[rk].builds[bk]

    def buildkey_to_include(self, rk: RefKey, bk: BuildKey) -> Include:
        return self.buildkey_to_build(rk, bk).include

    def buildkey_to_chr_indices(self, rk: RefKey, bk: BuildKey) -> set[ChrIndex]:
        cs = self.stratifications[rk].builds[bk].chr_filter
        return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    # src getters (for use in downloading inputs)

    def refkey_to_ref_src(self, k: RefKey) -> RefSrc:
        return self.refkey_to_strat(k).ref.src

    def refkey_to_gap_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(lambda x: x.src, self.stratifications[k].gap)

    def refkey_to_x_features_src(self, k: RefKey) -> BedSrc:
        return self.stratifications[k].xy.features.x_bed.src

    def refkey_to_y_features_src(self, k: RefKey) -> BedSrc:
        return self.stratifications[k].xy.features.y_bed.src

    def refkey_to_simreps_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(
            lambda x: x.simreps.src, self.stratifications[k].low_complexity
        )

    def refkey_to_rmsk_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(lambda x: x.rmsk.src, self.stratifications[k].low_complexity)

    def refkey_to_satellite_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(
            lambda x: fmap_maybe(lambda x: x.src, x.satellites),
            self.stratifications[k].low_complexity,
        )

    def refkey_to_superdups_src(self, k: RefKey) -> BedSrc | None:
        return fmap_maybe(
            lambda x: fmap_maybe(lambda x: x.src, x.superdups),
            self.stratifications[k].segdups,
        )

    def refkey_to_functional_ftbl_src(self, k: RefKey) -> BedSrc:
        return self.stratifications[k].functional.ftbl_src

    def refkey_to_functional_gff_src(self, k: RefKey) -> BedSrc:
        return self.stratifications[k].functional.gff_src

    # checksum getters (for use in ensuring inputs are correct)

    def refkey_to_ref_checksum(self, k: RefKey) -> str | None:
        return self.refkey_to_ref_src(k).md5

    # chromosome standardization

    def refkey_to_final_chr_prefix(self, k: RefKey) -> str:
        return self.stratifications[k].ref.chr_prefix

    def buildkey_to_final_chr_mapping(self, rk: RefKey, bk: BuildKey) -> dict[int, str]:
        p = self.refkey_to_final_chr_prefix(rk)
        cs = self.buildkey_to_chr_indices(rk, bk)
        return {c.value: c.chr_name_full(p) for c in cs}

    def buildkey_to_chr_conversion(
        self,
        rk: RefKey,
        bk: BuildKey,
        fromChr: str,
    ) -> ChrConversion:
        toChr = self.refkey_to_final_chr_prefix(rk)
        cis = self.buildkey_to_chr_indices(rk, bk)
        return ChrConversion(fromChr, toChr, cis)

    def buildkey_to_mappability(
        self,
        rk: RefKey,
        bk: BuildKey,
    ) -> tuple[list[int], list[int], list[int]]:
        ms = self.buildkey_to_include(rk, bk).mappability
        l, m, e = unzip((m.length, m.mismatches, m.indels) for m in ms)
        return ([*l], [*m], [*e])

    # include switches (for controlling which snakemake rules to activate)

    def want_low_complexity_censat(self, rk: RefKey) -> bool:
        return self.stratifications[rk].low_complexity.satellites is not None

    def _want_chr_index(self, rk: RefKey, bk: BuildKey, i: ChrIndex) -> bool:
        cis = self.buildkey_to_chr_indices(rk, bk)
        return i in cis

    def want_xy_x(self, rk: RefKey, bk: BuildKey) -> bool:
        return (
            self._want_chr_index(rk, bk, ChrIndex.CHRX)
            and self.buildkey_to_include(rk, bk).xy
        )

    def want_xy_y(self, rk: RefKey, bk: BuildKey) -> bool:
        return (
            self._want_chr_index(rk, bk, ChrIndex.CHRY)
            and self.buildkey_to_include(rk, bk).xy
        )

    def wanted_xy_chr_names(self, rk: RefKey, bk: BuildKey) -> list[str]:
        return [
            i.chr_name
            for i in [ChrIndex.CHRX, ChrIndex.CHRY]
            if self._want_chr_index(rk, bk, i)
        ]

    def want_xy_sex(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.want_xy_x(rk, bk) and self.want_xy_y(rk, bk)

    def want_xy_auto(self, rk: RefKey, bk: BuildKey) -> bool:
        cis = self.buildkey_to_chr_indices(rk, bk)
        return len(cis - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

    def want_xy_XTR(self, rk: RefKey) -> bool:
        return self.refkey_to_strat(rk).xy.features.xtr

    def want_xy_ampliconic(self, rk: RefKey) -> bool:
        return self.refkey_to_strat(rk).xy.features.ampliconic

    def want_low_complexity(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_include(rk, bk).low_complexity

    def want_gc(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_include(rk, bk).gc

    def want_functional(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_include(rk, bk).functional

    def want_segdups(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_include(rk, bk).segdups

    def _want_union(self, rk: RefKey, bk: BuildKey) -> bool:
        return self.buildkey_to_include(rk, bk).union

    def want_mappability(self, rk: RefKey, bk: BuildKey) -> bool:
        return len(self.buildkey_to_include(rk, bk).mappability) > 0

    def want_segdup_and_map(self, rk: RefKey, bk: BuildKey) -> bool:
        return (
            self.buildkey_to_include(rk, bk).union
            and self.want_segdups(rk, bk)
            and self.want_mappability(rk, bk)
        )

    def want_alldifficult(self, rk: RefKey, bk: BuildKey) -> bool:
        return (
            self.want_segdup_and_map(rk, bk)
            and self.want_low_complexity(rk, bk)
            and self.want_gc(rk, bk)
        )

    # key lists for downloading resources

    @property
    def _all_builds(self) -> list[tuple[RefKey, BuildKey]]:
        return [(rk, bk) for rk, s in self.stratifications.items() for bk in s.builds]

    @property
    def all_refkeys(self) -> list[RefKey]:
        return [*self.stratifications]

    def _all_refkey_from_want(
        self,
        f: Callable[[RefKey, BuildKey], bool],
    ) -> list[RefKey]:
        return [*unique_everseen(x[0] for x in self._all_builds if f(*x))]

    @property
    def all_refkey_gap(self) -> list[RefKey]:
        return [k for k in self.all_refkeys if self.refkey_to_gap_src(k) is not None]

    @property
    def all_refkey_rmsk_trf(self) -> list[RefKey]:
        return self._all_refkey_from_want(self.want_low_complexity)

    @property
    def all_refkey_censat(self) -> list[RefKey]:
        return self._all_refkey_from_want(
            lambda r, b: self.want_low_complexity(r, b)
            and self.want_low_complexity_censat(r)
        )

    @property
    def all_refkey_functional(self) -> list[RefKey]:
        return self._all_refkey_from_want(self.want_functional)

    @property
    def all_refkey_segdups(self) -> list[RefKey]:
        return self._all_refkey_from_want(self.want_segdups)
