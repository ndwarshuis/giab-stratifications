"""
Conventions:
* Functions ending in "_unsafe" should never throw errors; if they do then the
  code is incorrect. This is in contrast with other errors which may happen due
  to network issues, invalid inputs, etc.
"""
import sys
import json
import pandas as pd
import re
from pathlib import Path
from pydantic import BaseModel as BaseModel_
from pydantic.generics import GenericModel as GenericModel_
from pydantic.generics import GenericModelT
from pydantic import validator, HttpUrl, FilePath, NonNegativeInt, Field
from dataclasses import dataclass
from enum import Enum, unique
from typing import (
    Union,
    NewType,
    Any,
    Callable,
    TypeVar,
    Type,
    NamedTuple,
    cast,
    Annotated,
    Generic,
    TypeGuard,
    Protocol,
)
from typing_extensions import Self, assert_never
from common.functional import (
    fmap_maybe,
    fmap_maybe_def,
    both,
    with_first,
    DesignError,
    match1_unsafe,
    match2_unsafe,
    not_none_unsafe,
    none_unsafe,
)
from common.io import is_gzip, is_bgzip
import common.bed as bed
from more_itertools import unzip


W = TypeVar("W")
X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


Percent = Annotated[int, Field(ge=0, le=100)]


class RefKey_(str):
    pass


class HapRefKey_(RefKey_):
    pass


class Dip1RefKey_(RefKey_):
    pass


class Dip2RefKey_(RefKey_):
    pass


RefKeyT = TypeVar("RefKeyT", HapRefKey_, Dip1RefKey_, Dip2RefKey_)


HapRefKey = HapRefKey_
Dip1RefKey = Dip1RefKey_
Dip2RefKey = Dip2RefKey_
AnyRefKey = HapRefKey | Dip1RefKey | Dip2RefKey


@dataclass(frozen=True)
class RefKeyFull(Generic[RefKeyT]):
    key: RefKeyT
    hap: "Haplotype | None"

    @property
    def strip(self) -> RefKey_:
        return RefKey_(self.key)

    @property
    def as_tuple(self) -> "tuple[str, Haplotype | None]":
        return (str(self.key), self.hap)

    @property
    def name(self) -> str:
        k, h = self.as_tuple
        return f"{k}.{h.name}" if h is not None else k


class RefFinalKey(RefKeyFull[RefKeyT]):
    pass


class RefSrcKey(RefKeyFull[RefKeyT]):
    pass


AnyRefFinalKey = (
    RefFinalKey[HapRefKey] | RefFinalKey[Dip1RefKey] | RefFinalKey[Dip2RefKey]
)


# class wrapper so I can pattern match on them (which newtype won't allow)
class HapBuildKey(str):
    pass


class Dip1BuildKey(str):
    pass


class Dip2BuildKey(str):
    pass


BuildKey = HapBuildKey | Dip1BuildKey | Dip2BuildKey
RefKey = HapRefKey | Dip1RefKey | Dip2RefKey


BuildKeyT = TypeVar("BuildKeyT", HapBuildKey, Dip1BuildKey, Dip2BuildKey)

CompareKey = NewType("CompareKey", str)
OtherLevelKey = NewType("OtherLevelKey", str)
OtherStratKey = NewType("OtherStratKey", str)
HaplotypeName = NewType("HaplotypeName", str)


CHR_INDEX_PLACEHOLDER = "%i"
CHR_HAP_PLACEHOLDER = "%h"


class BaseModel(BaseModel_):
    class Config:
        frozen = True
        extra = "forbid"


class GenericModel(GenericModel_):
    class Config:
        frozen = True
        extra = "forbid"

    # dirty hack to get pickling to work for generic model types; see
    # https://github.com/pydantic/pydantic/issues/1667
    #
    # it seems this was fixed, but there might be some issue with getting
    # snakemake to recognize to get its paths correct
    def __class_getitem__(
        cls: Type[GenericModelT], params: Union[Type[Any], tuple[Type[Any], ...]]
    ) -> Type[Any]:
        created_class = super().__class_getitem__(params)
        setattr(
            sys.modules[created_class.__module__], created_class.__name__, created_class
        )
        return created_class


class Haplotype(Enum):
    HAP1: int = 0
    HAP2: int = 1

    @classmethod
    def from_name(cls, n: str) -> Self:
        try:
            return next(i for i in cls if i.name == n)
        except StopIteration:
            raise ValueError(f"could make haplotype from name '{n}'")

    @classmethod
    def with_haps(cls, f: "Callable[[Haplotype], X]") -> tuple[X, X]:
        x0 = f(cls.HAP1)
        x1 = f(cls.HAP2)
        return (x0, x1)

    @property
    def name(self) -> HaplotypeName:
        return HaplotypeName(f"hap{self.value + 1}")

    def from_either(self, left: X, right: X) -> X:
        if self is Haplotype.HAP1:
            return left
        elif self is Haplotype.HAP2:
            return right
        else:
            assert_never(self)


class _Src:
    def to_refkeys(self, rk: RefKeyT) -> list[RefKeyFull[RefKeyT]]:
        return NotImplemented

    def to_str_refkeys(self, rk: RefKeyT) -> list[str]:
        return [k.name for k in self.to_refkeys(rk)]


# TODO this is silly, I don't want to be required to put "hap: blablab"
# whenever I just have "one thing"
class Haploid_(GenericModel, Generic[X], _Src):
    hap: X

    def to_refkeys(self, rk: RefKeyT) -> list[RefKeyFull[RefKeyT]]:
        return [self.key(rk)]

    def key(self, rk: RefKeyT) -> RefKeyFull[RefKeyT]:
        return RefKeyFull(rk, None)


class Diploid_(GenericModel, Generic[X], _Src):
    hap1: X
    hap2: X

    def from_either(self, hap: Haplotype) -> X:
        return hap.from_either(self.hap1, self.hap2)

    def key1(self, rk: RefKeyT) -> RefKeyFull[RefKeyT]:
        return RefKeyFull(rk, Haplotype.HAP1)

    def key2(self, rk: RefKeyT) -> RefKeyFull[RefKeyT]:
        return RefKeyFull(rk, Haplotype.HAP2)

    def keys(self, rk: RefKeyT) -> tuple[RefKeyFull[RefKeyT], RefKeyFull[RefKeyT]]:
        return (self.key1(rk), self.key2(rk))

    def to_refkeys(self, rk: RefKeyT) -> list[RefKeyFull[RefKeyT]]:
        return list(self.keys(rk))


def parse_final_refkey(s: str) -> tuple[str, Haplotype | None]:
    m = re.match("(.+)\\.(hap[12])", s)
    # ASSUME this will never fail due to the hap1/2 permitted match pattern
    return (s, None) if m is None else (m[1], Haplotype.from_name(m[2]))


# dummy Identity type to make higher-order types more consistent
# HapOnly = Union[Haploid_[X]]
HapOrDip = Haploid_[X] | Diploid_[X]


def choose_xy_unsafe(c: "ChrIndex", x_res: X, y_res: X) -> X:
    if c is ChrIndex.CHRX:
        return x_res
    elif c is ChrIndex.CHRY:
        return y_res
    else:
        raise DesignError(f"I am not an X or Y, I am a {c}")


@unique
class ChrIndex(Enum):
    """Represents a valid chromosome index.

    Chromosomes are numbered by integers 1-24 (23 and 24 being X and Y
    respectively). These integers reflect the sort order in output bed files.
    """

    # NOTE: these start at 1 not 0 to conincide with the names of (most)
    # the chromosomes
    CHR1: int = 1
    CHR2: int = 2
    CHR3: int = 3
    CHR4: int = 4
    CHR5: int = 5
    CHR6: int = 6
    CHR7: int = 7
    CHR8: int = 8
    CHR9: int = 9
    CHR10: int = 10
    CHR11: int = 11
    CHR12: int = 12
    CHR13: int = 13
    CHR14: int = 14
    CHR15: int = 15
    CHR16: int = 16
    CHR17: int = 17
    CHR18: int = 18
    CHR19: int = 19
    CHR20: int = 20
    CHR21: int = 21
    CHR22: int = 22
    CHRX: int = 23
    CHRY: int = 24

    @classmethod
    def from_name(cls, n: str) -> Self:
        try:
            return next(i for i in cls if i.chr_name == n)
        except StopIteration:
            raise ValueError(f"could make chr index from name '{n}'")

    @classmethod
    def from_name_unsafe(cls, n: str) -> Self:
        try:
            return cls.from_name(n)
        except ValueError as e:
            raise DesignError(e)

    def __init__(self, i: int) -> None:
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def chr_name_full(self, p: "HapChrPattern") -> str | None:
        return p.to_chr_name(self)

    def chr_name_full_dip(self, p: "DipChrPattern", hap: Haplotype) -> str | None:
        return p.to_chr_name(self, hap)

    def to_internal_index(self, hap: Haplotype) -> bed.InternalChrIndex:
        return bed.InternalChrIndex(hap.value * 24 + self.value - 1)

    @property
    def xy_to_hap_unsafe(self) -> Haplotype:
        return choose_xy_unsafe(self, Haplotype.HAP2, Haplotype.HAP1)


class ChrPattern:
    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        return NotImplemented


class HapChrPattern(BaseModel, ChrPattern):
    template: str = "chr%i"
    special: dict[ChrIndex, str] = {}
    exclusions: list[ChrIndex] = []

    @validator("template")
    def is_valid_template(cls, v: str) -> str:
        assert v.count(CHR_INDEX_PLACEHOLDER) == 1, "chr template must have '%i' in it"
        return v

    def to_chr_name(self, i: ChrIndex) -> str | None:
        if i in self.exclusions:
            return None
        elif i in self.special:
            return self.special[i]
        else:
            return self.template.replace(CHR_INDEX_PLACEHOLDER, i.chr_name)

    def to_pairs(
        self,
        cs: set[ChrIndex],
        h: Haplotype,
    ) -> list[tuple[bed.InternalChrIndex, str]]:
        return [
            (c.to_internal_index(h), n)
            for c in cs
            if (n := self.to_chr_name(c)) is not None
        ]

    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        # NOTE: the haplotype argument is doing nothing since it is only
        # used to make the index which I remove before returning here
        return [x[1] for x in self.to_pairs(cs, Haplotype.HAP1)]

    def init_mapper(self, cs: set[ChrIndex], hap: Haplotype) -> bed.InitMapper:
        return {n: i for i, n in self.to_pairs(cs, hap)}

    def final_mapper(self, cs: set[ChrIndex], hap: Haplotype) -> bed.FinalMapper:
        return {i: n for i, n in self.to_pairs(cs, hap)}


class DipChrPattern(BaseModel, ChrPattern):
    template: str = "chr%i_%h"
    special: dict[ChrIndex, str] = {}
    hapnames: Diploid_[HaplotypeName] = Diploid_(
        hap1=HaplotypeName("PATERNAL"),
        hap2=HaplotypeName("MATERNAL"),
    )
    exclusions: Diploid_[list[ChrIndex]] = Diploid_(
        hap1=[ChrIndex.CHRX],
        hap2=[ChrIndex.CHRY],
    )

    @validator("template")
    def is_valid_template(cls, v: str) -> str:
        assert (
            v.count(CHR_INDEX_PLACEHOLDER) == 1 and v.count(CHR_HAP_PLACEHOLDER) == 1
        ), "chr template must have '%i' and '%h' in it"
        return v

    def to_chr_name(self, i: ChrIndex, h: Haplotype) -> str | None:
        exc = self.exclusions.from_either(h)
        name = self.hapnames.from_either(h)
        if i in exc:
            return None
        elif i in self.special:
            return self.special[i]
        else:
            return self.template.replace(CHR_INDEX_PLACEHOLDER, i.chr_name).replace(
                CHR_HAP_PLACEHOLDER, name
            )

    def to_pairs(self, cs: set[ChrIndex]) -> list[tuple[bed.InternalChrIndex, str]]:
        return [
            (c.to_internal_index(h), n)
            for c in cs
            for h in Haplotype
            if (n := self.to_chr_name(c, h)) is not None
        ]

    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        return [x[1] for x in self.to_pairs(cs)]

    def init_mapper(self, cs: set[ChrIndex]) -> bed.InitMapper:
        return {n: i for i, n in self.to_pairs(cs)}

    def final_mapper(self, cs: set[ChrIndex]) -> bed.FinalMapper:
        return {i: n for i, n in self.to_pairs(cs)}

    def to_hap_pattern(self, hap: Haplotype) -> HapChrPattern:
        hs = self.hapnames.from_either(hap)
        xs = self.exclusions.from_either(hap)
        return HapChrPattern(
            template=self.template.replace(CHR_HAP_PLACEHOLDER, hs),
            special=self.special,
            exclusions=xs,
        )


AnyChrPatternT = TypeVar("AnyChrPatternT", HapChrPattern, DipChrPattern)


class HapChrSource(GenericModel, Generic[X]):
    """Specification for a haploid source file."""

    src: Haploid_[X]
    chr_pattern: HapChrPattern = HapChrPattern()


class DipChrSource1(GenericModel, Generic[X]):
    """Specification for a combined diploid source file.

    The 'src' is assumed to have all chromosomes for both haplotypes in one
    file, which implies they are labeled so as to distinguish the haps. The
    pattern will match both the chromosome number and the haplotype within the
    chromosome name.
    """

    src: Haploid_[X]
    chr_pattern: DipChrPattern = DipChrPattern()


class DipChrSource2(GenericModel, Generic[X]):
    """Specification for split diploid source files.

    Each source may or may not have each haplotype labeled; the identity of each
    haplotype in either source file is determined based on the configuration key
    under which it appears (hap1 or hap2) and the chromosome names for each are
    matched according to its corresponding entry in `chr_pattern`.
    """

    src: Diploid_[X]
    chr_pattern: Diploid_[HapChrPattern] = Diploid_(
        hap1=HapChrPattern(),
        hap2=HapChrPattern(),
    )


# DipChrSource: TypeAlias = DipChrSource1[X] | DipChrSource2[X]
# AnyChrSource = Union[HapChrSource[X], DipChrSource[X]]


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

AnySrcT = TypeVar("AnySrcT", Haploid_[BedSrc], Diploid_[BedSrc])


# TODO clean this up with real polymorphism when mypy catches up with Haskell
# 98, see https://github.com/python/typing/issues/548
AnyBedT = TypeVar(
    "AnyBedT",
    HapChrSource[BedSrc],
    DipChrSource1[BedSrc] | DipChrSource2[BedSrc],
    # DipChrSource[BedSrc],
)
AnyBedT_ = TypeVar(
    "AnyBedT_",
    HapChrSource[BedSrc],
    DipChrSource1[BedSrc],
    DipChrSource2[BedSrc],
)


@unique
class CoreLevel(Enum):
    """A stratification level (eg "GCcontent" or "mappability")

    These are the only "built-in" levels contained within the pipeline.
    Users may add other levels if they wish to include other source
    files, but these must be specified manually (see below).
    """

    FUNCTIONAL = "Functional"
    LOWCOMPLEXITY = "LowComplexity"
    GC = "GCcontent"
    MAPPABILITY = "Mappability"
    SEGDUPS = "SegmentalDuplications"
    UNION = "Union"
    TELOMERES = "Telomere"
    XY = "XY"
    # overlaps with "other" strat categories, needed because this is where
    # the gaps strat will go
    OtherDifficult = "OtherDifficult"


@dataclass
class HapToHapChrConversion:
    fromPattern: HapChrPattern
    toPattern: HapChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices, Haplotype.HAP1)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices, Haplotype.HAP1)


@dataclass
class DipToDipChrConversion:
    fromPattern: DipChrPattern
    toPattern: DipChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass
class HapToDipChrConversion:
    fromPattern: Diploid_[HapChrPattern]
    toPattern: DipChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.InitMapper]:
        p = self.fromPattern
        i = self.indices
        return (
            p.hap1.init_mapper(i, Haplotype.HAP1),
            p.hap2.init_mapper(i, Haplotype.HAP2),
        )

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass
class DipToHapChrConversion:
    fromPattern: DipChrPattern
    toPattern: Diploid_[HapChrPattern]
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.SplitMapper]:
        im = self.fromPattern.init_mapper(self.indices)
        fm0 = self.toPattern.hap1.final_mapper(self.indices, Haplotype.HAP1)
        return (im, bed.make_split_mapper(im, fm0))

    @property
    def final_mapper(self) -> tuple[bed.FinalMapper, bed.FinalMapper]:
        p = self.toPattern
        i = self.indices
        return (
            p.hap1.final_mapper(i, Haplotype.HAP1),
            p.hap2.final_mapper(i, Haplotype.HAP2),
        )


class Paths(BaseModel):
    """Local build paths for snakemake."""

    resources: Path = Path("resources")
    results: Path = Path("results")


class Tools(BaseModel):
    """Urls for tools to download/build/use in the pipeline."""

    repseq: HttpUrl = "https://github.com/ndwarshuis/repseq/archive/refs/tags/v1.1.0.tar.gz"  # type: ignore
    gemlib: HttpUrl = "https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download"  # type: ignore


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


class BedFileParams(BaseModel):
    """Parameters decribing how to parse a bed-like file.

    Members:
    chr_pattern - the pattern on the chromosomes; must include the special
      directive '%i' which will denote a standardized name (eg 1, 2 ...X, Y);
      the pattern is assumed to match the whole chromosome name.
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    """

    # chr_pattern: HapChrPattern = HapChrPattern()
    bed_cols: BedColumns = BedColumns()
    skip_lines: NonNegativeInt = 0
    sep: str = "\t"


class BedFile(GenericModel, Generic[X]):
    """Inport specs for a bed-like file."""

    data: X
    params: BedFileParams = BedFileParams()

    # @property
    # def src_list(self) -> list[BedSrc]:
    #     return diploid_to_list(self.data)

    def read(self, path: Path) -> pd.DataFrame:
        return self._read(path, [])

    def _read(self, path: Path, more: list[int] = []) -> pd.DataFrame:
        p = self.params
        return bed.read_bed(path, p.bed_cols.columns, p.skip_lines, p.sep, more)

    # def read_filter_sort_bed(
    #     self,
    #     ipath: Path,
    #     opath: Path,
    #     conv:
    #     more: list[int] = [],
    # ) -> None:
    #     """Read a haploid bed file, sort it, and write it in bgzip format."""
    #     # conv = sconf.haploid_stratifications.to_build_data_unsafe(
    #     #     rk, bk
    #     # ).chr_conversion(pat)
    #     im = None
    #     fm = None
    #     df = self._read(ipath, more)
    #     df_ = bed.filter_sort_bed(im, fm, df)
    #     bed.write_bed(opath, df_)


HapBedSrc = HapChrSource[BedSrc]
DipBedSrc = DipChrSource1[BedSrc] | DipChrSource2[BedSrc]
Dip1BedSrc = DipChrSource1[BedSrc]
Dip2BedSrc = DipChrSource2[BedSrc]

HapBedFile = BedFile[HapBedSrc]
DipBedFile = BedFile[DipBedSrc]
Dip1BedFile = BedFile[Dip1BedSrc]
Dip2BedFile = BedFile[Dip2BedSrc]
AnyBedFileT = BedFile[AnyBedT]


# TODO mypy for some reason doesn't understand how to narrow a
# Something[Union[X, Y]] to a Something[X] using 'isinstance'
def is_dip1_bed(
    x: BedFile[DipChrSource1[X] | DipChrSource2[X]],
) -> TypeGuard[BedFile[DipChrSource1[X]]]:
    return isinstance(x.data, DipChrSource1)


def is_dip2_bed(
    x: BedFile[DipChrSource1[X] | DipChrSource2[X]],
) -> TypeGuard[BedFile[DipChrSource2[X]]]:
    return isinstance(x.data, DipChrSource2)


def with_dip_bedfile(
    bf: BedFile[DipChrSource1[X] | DipChrSource2[X]],
    dip1: Callable[[BedFile[DipChrSource1[X]]], Y],
    dip2: Callable[[BedFile[DipChrSource2[X]]], Y],
) -> Y:
    if is_dip1_bed(bf):
        return dip1(bf)
    elif is_dip2_bed(bf):
        return dip2(bf)
    else:
        assert False, "TODO this is a mypy FP"
        # assert_never(bf)


def with_hap_or_dip(
    x: Haploid_[X] | Diploid_[X],
    dip1: Callable[[Haploid_[X]], Y],
    dip2: Callable[[Diploid_[X]], Y],
) -> Y:
    if isinstance(x, Diploid_):
        return dip2(x)
    if isinstance(x, Haploid_):
        return dip1(x)
    else:
        assert_never(x)


def from_hap_or_dip(x: Haploid_[X] | Diploid_[X], hap: Haplotype | None) -> X:
    return with_hap_or_dip(
        x,
        lambda x: none_unsafe(hap, x.hap),
        lambda x: not_none_unsafe(hap, lambda h: x.from_either(h)),
    )


class VCFFile(BedFile[X], Generic[X]):
    """Inport specs for a vcf file."""

    data: X  # type narrowing won't work without this redfinition


class RMSKFile(BedFile[X], Generic[X]):
    """Input file for repeat masker stratification."""

    data: X  # type narrowing won't work without this redfinition
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

    def read(self, path: Path) -> pd.DataFrame:
        return super()._read(path, [self.class_col])


class SatFile(BedFile[X], Generic[X]):
    """Configuration for a satellites file."""

    sat_col: NonNegativeInt

    def read(self, path: Path) -> pd.DataFrame:
        return super()._read(path, [self.sat_col])


class LowComplexity(GenericModel, Generic[X]):
    """Configuration for low complexity stratification."""

    rmsk: RMSKFile[X] | None
    simreps: BedFile[X] | None
    satellites: SatFile[X] | None


class XYFile(HapBedFile):
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

    def read(self, path: Path) -> pd.DataFrame:
        return super()._read(path, [self.level_col])


# TODO what if the reference is XX?
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

    def fmt(self, i: ChrIndex, pattern: HapChrPattern) -> str:
        # TODO this smells like something I'll be doing alot
        c = i.chr_name_full(pattern)
        return "\n".join(
            [
                f"{c}\t{self.start[0]}\t{self.start[1]}",
                f"{c}\t{self.end[0]}\t{self.end[1]}",
            ]
        )


class XY(BaseModel):
    """Configuration for the XY stratification."""

    features: XYFeatures | None
    x_par: XYPar | None
    y_par: XYPar | None

    def fmt_x_par(self, pattern: HapChrPattern) -> str | None:
        return fmap_maybe(lambda x: x.fmt(ChrIndex.CHRX, pattern), self.x_par)

    def fmt_x_par_unsafe(self, pattern: HapChrPattern) -> str:
        s = self.fmt_x_par(pattern)
        if s is None:
            raise DesignError("X PAR does not exist")
        return s

    def fmt_y_par(self, pattern: HapChrPattern) -> str | None:
        return fmap_maybe(lambda x: x.fmt(ChrIndex.CHRY, pattern), self.y_par)

    def fmt_y_par_unsafe(self, pattern: HapChrPattern) -> str:
        s = self.fmt_y_par(pattern)
        if s is None:
            raise DesignError("Y PAR does not exist")
        return s


class Mappability(BaseModel):
    """Configuration for Mappability stratification.

    members:
    - unplaced_chr_patterns: a list of regexps that will be used to identify
      non-primary chromosomes in the reference to be included in mappability
      evaluation.
    """

    unplaced_chr_patterns: list[str]


class SegDups(GenericModel, Generic[X]):
    """Configuration for Segdup stratifications."""

    superdups: BedFile[X] | None


class LowMapParams(BaseModel):
    """Parameters for a single mappability bed file."""

    length: NonNegativeInt
    mismatches: NonNegativeInt
    indels: NonNegativeInt


GCBound = tuple[Percent, bool]


class GCParams(BaseModel):
    low: list[GCBound] = [
        (15, False),
        (20, False),
        (25, True),
        (30, True),
    ]
    high: list[GCBound] = [
        (55, True),
        (60, False),
        (65, True),
        (70, False),
        (75, False),
        (80, False),
        (85, False),
    ]

    @validator("high")
    def has_balanced_ranges(
        cls,
        high: list[GCBound],
        values: dict[str, Any],
    ) -> list[tuple[Percent, bool]]:
        try:
            low = cast(list[GCBound], values["low"])
            assert len([x for x in low if x[1]]) == len(
                [x for x in high if x[1]]
            ), "GC low/high must have same number of range boundaries"
        except KeyError:
            pass
        return high

    @property
    def low_sorted(self) -> list[GCBound]:
        return sorted(self.low, key=lambda x: x[0])

    @property
    def high_sorted(self) -> list[GCBound]:
        return sorted(self.high, key=lambda x: x[0])

    @property
    def low_fractions(self) -> list[int]:
        return [x[0] for x in self.low_sorted]

    @property
    def high_fractions(self) -> list[int]:
        return [x[0] for x in self.high_sorted]

    @property
    def low_bounds(self) -> tuple[int, list[int]]:
        bounds = self.low_fractions
        return (bounds[0], bounds[1:])

    @property
    def high_bounds(self) -> tuple[int, list[int]]:
        bounds = self.high_fractions
        return (bounds[-1], bounds[:-1])


class Include(BaseModel):
    """Flags to control which stratification levels are included."""

    low_complexity: bool = True
    xy: bool = True
    functional: bool = True
    segdups: bool = True
    union: bool = True
    telomeres: bool = True
    # TODO also add KIR and MHC since these should be derivable from refseq
    vdj: bool = True
    mappability: set[LowMapParams] = {
        LowMapParams(length=250, mismatches=0, indels=0),
        LowMapParams(length=100, mismatches=2, indels=1),
    }
    gc: GCParams | None = GCParams()
    # NOTE: This is crude but it should a) work, b) provide a decent user xp
    # and c) typecheck nicely without requiring me to use a zillionth typevar
    # in all my signatures.
    #
    # This parameter is only used for diploid configurations. While it will be
    # present for all configurations, it will only be read when needed.
    # Defaulting to True means that the user never needs to specify it if they
    # want this for diploid (which they probably do) and they don't need to care
    # in haploid cases. The only issue would be if the user specified this in
    # the haploid case; it technically should be a validation error since it
    # makes no sense in the case of haploid, but here it is setup to not hurt
    # anything.
    hets: bool = True


class OtherBedFile(BedFile[AnyBedT], Generic[AnyBedT]):
    remove_gaps: bool = False


class Bench(GenericModel, Generic[AnyBedT, AnyBedT_]):
    """Configuration for benchmark to use when validating stratifications."""

    bench_vcf: VCFFile[AnyBedT_]
    query_vcf: VCFFile[AnyBedT_]
    bench_bed: BedFile[AnyBedT]


class BuildCompare(BaseModel):
    """Configuration for comparing generated strats to previous versions."""

    other: CompareKey
    path_mapper: dict[Path, Path] = {}
    replacements: list[tuple[str, str]] = []
    ignore_other: list[str] = []
    ignore_generated: list[str] = []


class Build_(GenericModel, Generic[AnyBedT, AnyBedT_]):
    chr_filter: set[ChrIndex]
    comparison: BuildCompare | None = None
    bench: Bench[AnyBedT, AnyBedT_] | None = None
    other_strats: dict[OtherLevelKey, dict[OtherStratKey, OtherBedFile[AnyBedT]]] = {}
    # TODO if I really want I could validate this such that the user would be
    # politely alerted in case they specify any diploid params for a haploid
    # config.
    include: Include = Include()

    @property
    def compare_key(self) -> CompareKey | None:
        return fmap_maybe(lambda x: x.other, self.comparison)


# # need class overrides here to get default values for include, which are
# # different b/t hap and dip cases
# class HaploidBuild(Build_[AnyBedT, AnyBedT_]):
#     """Spec for a stratification build."""

#     include: HapInclude = IncludeHaploid()


# # # AnyDipBedSrcT = TypeVar("AnyDipBedSrcT", DipChrSource1[BedSrc], DipChrSource2[BedSrc])


# class DiploidBuild(Build_[AnyBedT, AnyBedT_]):
#     """Spec for a stratification build."""

#     include: DipInclude = DipInclude()


# BenchT = TypeVar(
#     "BenchT",
#     Bench[HapChrSource[BedSrc]],
#     Bench[DipChrSource1[BedSrc]],
#     Bench[DipChrSource2[BedSrc]],
# )

# RefFile = AnyChrSource[RefSrc]


# class RefFile(BaseModel, Generic[X, Y]):
#     """Specification for a reference file."""

#     src: X
#     chr_pattern: Y


# class RefFileHaploid(BaseModel):
#     """Specification for a reference file."""

#     src: RefSrc
#     chr_pattern: HapChrPattern = HapChrPattern()


# RefFileDiploid = Diploid_[RefFileHaploid]


class Functional(GenericModel, Generic[X]):
    """Configuration for Functional stratifications."""

    ftbl_src: X
    gff_src: X


class StratInputs_(GenericModel, Generic[AnyBedT, AnySrcT]):
    gap: BedFile[AnyBedT] | None
    low_complexity: LowComplexity[AnyBedT]
    xy: XY
    mappability: Mappability | None
    segdups: SegDups[AnyBedT]
    functional: Functional[AnySrcT] | None

    @property
    def xy_features_unsafe(self) -> XYFeatures:
        f = self.xy.features
        if f is None:
            raise DesignError("XY features does not exist")
        return f

    def xy_feature_bed_unsafe(self, i: ChrIndex) -> XYFile:
        f = self.xy_features_unsafe
        return choose_xy_unsafe(i, f.x_bed, f.y_bed)


class StratInputToBed(Protocol):
    A = TypeVar(
        "A", HapChrSource[BedSrc], DipChrSource1[BedSrc] | DipChrSource2[BedSrc]
    )

    def __call__(self, __x: StratInputs_[A, AnySrcT]) -> BedFile[A] | None:
        pass


class StratInputToSrc(Protocol):
    A = TypeVar("A", Haploid_[BedSrc], Diploid_[BedSrc])

    def __call__(self, __x: StratInputs_[AnyBedT, A]) -> A | None:
        pass


HaploidStratInputs = StratInputs_[HapChrSource[BedSrc], Haploid_[BedSrc]]
DiploidStratInputs = StratInputs_[
    DipChrSource1[BedSrc] | DipChrSource2[BedSrc], Diploid_[BedSrc]
]
AnyStratInputs = HaploidStratInputs | DiploidStratInputs

StratInputT = TypeVar("StratInputT", HaploidStratInputs, DiploidStratInputs)

RefSourceT = TypeVar(
    "RefSourceT",
    HapChrSource[RefSrc],
    DipChrSource1[RefSrc],
    DipChrSource2[RefSrc],
)


@dataclass(frozen=True)
class RefData_(Generic[RefKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT]):
    refkey: RefKeyT
    ref: RefSourceT
    strat_inputs: StratInputs_[AnyBedT, AnySrcT]
    builds: dict[BuildKeyT, Build_[AnyBedT, AnyBedT_]]

    @property
    def ref_refkeys(self) -> list[RefKeyFull[RefKeyT]]:
        return self.ref.src.to_refkeys(self.refkey)

    @property
    def ref_str_refkeys(self) -> list[str]:
        return self.ref.src.to_str_refkeys(self.refkey)

    @property
    def mappability_patterns(self) -> list[str]:
        return fmap_maybe_def(
            [],
            lambda m: m.unplaced_chr_patterns,
            self.strat_inputs.mappability,
        )

    def to_bed_refsrckeys(
        self,
        f: Callable[[StratInputs_[AnyBedT, AnySrcT]], BedFile[AnyBedT] | None],
    ) -> list[str]:
        return fmap_maybe_def(
            [],
            lambda b: b.data.src.to_str_refkeys(self.refkey),
            f(self.strat_inputs),
        )

    def to_build_data_unsafe(
        self,
        bk: BuildKeyT,
    ) -> "BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT]":
        bd = self.to_build_data(bk)
        if bd is None:
            raise DesignError(f"Could not create build data from key '{bk}'")
        return bd

    def to_build_data(
        self,
        bk: BuildKeyT,
    ) -> (
        "BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT] | None"
    ):
        try:
            return BuildData_(self, bk, self.builds[bk])
        except KeyError:
            return None

    def get_refkeys_unsafe_(self, f: "RefDataToSrc") -> list[str]:
        return not_none_unsafe(
            f(self),
            lambda s: s.to_str_refkeys(self.refkey),
        )

    def get_refkeys_unsafe(self, f: StratInputToSrc) -> list[str]:
        return not_none_unsafe(
            f(self.strat_inputs),
            lambda s: s.to_str_refkeys(self.refkey),
        )

    @property
    def has_low_complexity_rmsk(self) -> bool:
        return self.strat_inputs.low_complexity.rmsk is not None

    @property
    def has_low_complexity_simreps(self) -> bool:
        return self.strat_inputs.low_complexity.simreps is not None

    @property
    def has_low_complexity_censat(self) -> bool:
        return self.strat_inputs.low_complexity.satellites is not None


class RefDataToBed(Protocol):
    A = TypeVar(
        "A", HapChrSource[BedSrc], DipChrSource1[BedSrc] | DipChrSource2[BedSrc]
    )

    def __call__(
        self,
        __x: RefData_[RefKeyT, RefSourceT, A, AnyBedT_, AnySrcT, BuildKeyT],
    ) -> BedFile[A] | None:
        pass


class HapRefData(
    RefData_[
        HapRefKey_,
        HapChrSource[RefSrc],
        HapBedSrc,
        HapBedSrc,
        Haploid_[BedSrc],
        HapBuildKey,
    ]
):
    pass


class Dip1RefData(
    RefData_[
        Dip1RefKey_,
        DipChrSource1[RefSrc],
        DipBedSrc,
        Dip1BedSrc,
        Diploid_[BedSrc],
        Dip1BuildKey,
    ]
):
    pass


class Dip2RefData(
    RefData_[
        Dip2RefKey_,
        DipChrSource2[RefSrc],
        DipBedSrc,
        Dip2BedSrc,
        Diploid_[BedSrc],
        Dip2BuildKey,
    ]
):
    pass


AnyRefData = HapRefData | Dip1RefData | Dip2RefData


@dataclass
class BuildData_(Generic[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT]):
    refdata: RefData_[RefKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT]
    buildkey: BuildKeyT
    build: Build_[AnyBedT, AnyBedT_]

    @property
    def chr_indices(self) -> set[ChrIndex]:
        cs = self.build.chr_filter
        return set([x for x in ChrIndex]) if len(cs) == 0 else cs

    @property
    def want_xy_x(self) -> bool:
        return ChrIndex.CHRX in self.chr_indices and self.build.include.xy

    @property
    def want_xy_y(self) -> bool:
        return ChrIndex.CHRY in self.chr_indices and self.build.include.xy

    @property
    def wanted_xy_chr_names(self) -> list[str]:
        return [
            i.chr_name for i in [ChrIndex.CHRX, ChrIndex.CHRY] if i in self.chr_indices
        ]

    # def want_xy_sex(self, rk: RefKey, bk: BuildKey) -> bool:
    #     return self.want_xy_x(rk, bk) and self.want_xy_y(rk, bk)

    @property
    def want_x_PAR(self) -> bool:
        return self.want_xy_x and self.refdata.strat_inputs.xy.x_par is not None

    @property
    def want_y_PAR(self) -> bool:
        return self.want_xy_y and self.refdata.strat_inputs.xy.y_par is not None

    @property
    def want_xy_auto(self) -> bool:
        return len(self.chr_indices - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

    @property
    def want_xy_XTR(self) -> bool:
        return fmap_maybe_def(
            False, lambda x: x.xtr, self.refdata.strat_inputs.xy.features
        )

    @property
    def want_xy_ampliconic(self) -> bool:
        return fmap_maybe_def(
            False, lambda x: x.ampliconic, self.refdata.strat_inputs.xy.features
        )

    @property
    def want_low_complexity(self) -> bool:
        return self.build.include.low_complexity

    @property
    def want_gc(self) -> bool:
        return self.build.include.gc is not None

    @property
    def want_functional(self) -> bool:
        return (
            self.build.include.functional
            and self.refdata.strat_inputs.functional is not None
        )

    @property
    def want_telomeres(self) -> bool:
        return self.build.include.telomeres

    @property
    def want_segdups(self) -> bool:
        return (
            self.refdata.strat_inputs.segdups.superdups is not None
            and self.build.include.segdups
        )

    @property
    def _want_union(self) -> bool:
        return self.build.include.union

    @property
    def want_mappability(self) -> bool:
        return (
            self.refdata.strat_inputs.mappability is not None
            and len(self.build.include.mappability) > 0
        )

    @property
    def want_segdup_and_map(self) -> bool:
        return self.build.include.union and self.want_segdups and self.want_mappability

    @property
    def want_alldifficult(self) -> bool:
        return self.want_segdup_and_map and self.want_low_complexity and self.want_gc

    @property
    def want_benchmark(self) -> bool:
        return self.build.bench is not None

    @property
    def want_gaps(self) -> bool:
        return self.refdata.strat_inputs.gap is not None

    @property
    def want_vdj(self) -> bool:
        vdj_chrs = {ChrIndex(i) for i in [2, 7, 14, 22]}
        return (
            self.build.include.vdj
            and self.refdata.strat_inputs.functional is not None
            and len(self.chr_indices & vdj_chrs) > 0
        )

    @property
    def mappability_params(
        self,
    ) -> tuple[list[int], list[int], list[int]]:
        ms = self.build.include.mappability
        xs = [(m.length, m.mismatches, m.indels) for m in ms]
        # TODO mypy doesn't like unzip here for some reason :(
        return ([x[0] for x in xs], [x[1] for x in xs], [x[2] for x in xs])


def with_ref_data(
    rd: AnyRefData,
    hap_f: Callable[[HapRefData], X],
    dip1_f: Callable[[Dip1RefData], X],
    dip2_f: Callable[[Dip2RefData], X],
) -> X:
    if isinstance(rd.ref, HapChrSource):
        return hap_f(rd)
    elif isinstance(rd.ref, DipChrSource1):
        return dip1_f(rd)
    elif isinstance(rd.ref, DipChrSource2):
        return dip2_f(rd)
    else:
        assert_never(rd)


class RefDataToSrc(Protocol):
    A = TypeVar("A", Haploid_[BedSrc], Diploid_[BedSrc])

    def __call__(
        self,
        __x: RefData_[RefKeyT, RefSourceT, AnyBedT, AnyBedT_, A, BuildKeyT],
    ) -> A | None:
        pass


def ref_data_to_src_(rd: AnyRefData, hap: Haplotype | None, f: RefDataToSrc) -> BedSrc:
    src = with_ref_data(rd, lambda rd: f(rd), lambda rd: f(rd), lambda rd: f(rd))
    # TODO mypy doens't like me using my 'maybe' functional functions
    if src is None:
        raise DesignError()
    return from_hap_or_dip(src, hap)


def ref_data_to_src(
    rd: AnyRefData, hap: Haplotype | None, f: StratInputToSrc
) -> BedSrc:
    return ref_data_to_src_(rd, hap, lambda rd: f(rd.strat_inputs))


class BuildDataToBed(Protocol):
    A = TypeVar("A", HapBedSrc, DipBedSrc)

    def __call__(
        self,
        __x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, A, AnyBedT_, AnySrcT],
    ) -> BedFile[A] | None:
        pass


class BuildDataToVCF(Protocol):
    A = TypeVar("A", HapBedSrc, Dip1BedSrc, Dip2BedSrc)

    def __call__(
        self,
        __x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, A, AnySrcT],
    ) -> BedFile[A] | None:
        pass


class BuildDataToSrc(Protocol):
    A = TypeVar("A", Haploid_[BedSrc], Diploid_[BedSrc])

    def __call__(
        self,
        __x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, A],
    ) -> A | None:
        pass


Xcov = TypeVar("Xcov", covariant=True)


class OutputPattern(Protocol, Generic[Xcov]):
    A = TypeVar("A", HapRefKey_, Dip1RefKey_, Dip2RefKey_)

    def __call__(self, __x: RefKeyFull[A]) -> Xcov:
        pass


class Stratification(
    GenericModel, Generic[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT]
):
    """Configuration for stratifications for a given reference."""

    ref: RefSourceT
    strat_inputs: StratInputs_[AnyBedT, AnySrcT]
    builds: dict[BuildKeyT, Build_[AnyBedT, AnyBedT_]]


@dataclass
class HapBuildData(
    BuildData_[
        HapRefKey_,
        HapBuildKey,
        HapChrSource[RefSrc],
        HapBedSrc,
        HapBedSrc,
        Haploid_[BedSrc],
    ]
):
    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.refdata.ref.chr_pattern.final_mapper(
            self.chr_indices, Haplotype.HAP1
        )

    @property
    def ref_chr_conversion(self) -> HapToHapChrConversion:
        p = self.refdata.ref.chr_pattern
        return HapToHapChrConversion(p, p, self.chr_indices)

    def chr_conversion(self, fromChr: HapChrPattern) -> HapToHapChrConversion:
        return HapToHapChrConversion(
            fromChr,
            self.refdata.ref.chr_pattern,
            self.chr_indices,
        )

    def read_filter_sort_hap_bed(self, bf: HapBedFile, ipath: Path) -> pd.DataFrame:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        conv = self.chr_conversion(bf.data.chr_pattern)
        df = bf.read(ipath)
        return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)

    def read_write_filter_sort_hap_bed(
        self,
        bf: HapBedFile,
        ipath: Path,
        opath: Path,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df = self.read_filter_sort_hap_bed(bf, ipath)
        bed.write_bed(opath, g(df))


@dataclass
class Dip1BuildData(
    BuildData_[
        Dip1RefKey_,
        Dip1BuildKey,
        DipChrSource1[RefSrc],
        DipBedSrc,
        Dip1BedSrc,
        Diploid_[BedSrc],
    ]
):
    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.refdata.ref.chr_pattern.final_mapper(self.chr_indices)

    @property
    def ref_chr_conversion(self) -> DipToDipChrConversion:
        p = self.refdata.ref.chr_pattern
        return DipToDipChrConversion(p, p, self.chr_indices)

    def hap_chr_conversion(
        self,
        fromChr: Diploid_[HapChrPattern],
    ) -> HapToDipChrConversion:
        return HapToDipChrConversion(
            fromChr,
            self.refdata.ref.chr_pattern,
            self.chr_indices,
        )

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
    ) -> DipToDipChrConversion:
        return DipToDipChrConversion(
            fromChr,
            self.refdata.ref.chr_pattern,
            self.chr_indices,
        )

    def read_filter_sort_dip1_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
    ) -> pd.DataFrame:
        """Read a diploid bed file, sort it, and write it in bgzip format."""
        conv = self.dip_chr_conversion(bf.data.chr_pattern)
        df = bf.read(ipath)
        return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)

    def read_write_filter_sort_dip1_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
        opath: Path,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df = self.read_filter_sort_dip1_bed(bf, ipath)
        bed.write_bed(opath, g(df))

    def read_filter_sort_dip2_bed(
        self,
        bf: Dip2BedFile,
        ipath: tuple[Path, Path],
    ) -> pd.DataFrame:
        """Read two haploid bed files, combine and sort them as diploid, and write
        it in bgzip format.
        """

        def go(b: Dip2BedFile, i: Path, imap: bed.InitMapper) -> pd.DataFrame:
            df = b.read(i)
            return bed.filter_sort_bed(imap, fmap, df)

        conv = self.hap_chr_conversion(bf.data.chr_pattern)
        imap1, imap2 = conv.init_mapper
        fmap = conv.final_mapper

        return pd.concat(
            [
                go(bf, *x)
                for x in [
                    (ipath[0], imap1),
                    (ipath[1], imap2),
                ]
            ]
        )

    def read_write_filter_sort_dip2_bed(
        self,
        bf: Dip2BedFile,
        ipath: tuple[Path, Path],
        opath: Path,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df = self.read_filter_sort_dip2_bed(bf, ipath)
        bed.write_bed(opath, g(df))


@dataclass
class Dip2BuildData(
    BuildData_[
        Dip2RefKey_,
        Dip2BuildKey,
        DipChrSource2[RefSrc],
        DipBedSrc,
        Dip2BedSrc,
        Diploid_[BedSrc],
    ]
):
    @property
    def final_mapper(self) -> tuple[bed.FinalMapper, bed.FinalMapper]:
        p = self.refdata.ref.chr_pattern
        i = self.chr_indices
        return (
            p.hap1.final_mapper(i, Haplotype.HAP1),
            p.hap2.final_mapper(i, Haplotype.HAP2),
        )

    @property
    def ref_chr_conversion(self) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
        def go(h: HapChrPattern) -> HapToHapChrConversion:
            return HapToHapChrConversion(h, h, self.chr_indices)

        p = self.refdata.ref.chr_pattern
        return (go(p.hap1), go(p.hap2))

    def hap_chr_conversion(
        self,
        fromChr: Diploid_[HapChrPattern],
    ) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
        toChr = self.refdata.ref.chr_pattern
        cis = self.chr_indices
        return (
            HapToHapChrConversion(fromChr.hap1, toChr.hap1, cis),
            HapToHapChrConversion(fromChr.hap2, toChr.hap2, cis),
        )

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
    ) -> DipToHapChrConversion:
        return DipToHapChrConversion(
            fromChr,
            self.refdata.ref.chr_pattern,
            self.chr_indices,
        )

    def read_filter_sort_dip1_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        conv = self.dip_chr_conversion(bf.data.chr_pattern)
        imap, splitter = conv.init_mapper
        fmap0, fmap1 = conv.final_mapper

        def go(df: pd.DataFrame, fmap: bed.FinalMapper) -> pd.DataFrame:
            return bed.filter_sort_bed(imap, fmap, df)

        df = bf.read(ipath)
        df0, df1 = bed.split_bed(splitter, df)
        return (go(df0, fmap0), go(df1, fmap1))

    def read_write_filter_sort_dip1_bed(
        self,
        bf: Dip1BedFile,
        ipath: Path,
        opath: tuple[Path, Path],
        g0: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
        g1: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        """Read a haploid bed file, sort it, and write it in bgzip format."""
        df0, df1 = self.read_filter_sort_dip1_bed(bf, ipath)
        bed.write_bed(opath[0], g0(df0))
        bed.write_bed(opath[1], g1(df1))

    def read_filter_sort_dip2_bed(
        self,
        bf: Dip2BedFile,
        ipath: Path,
        hap: Haplotype,
    ) -> pd.DataFrame:
        conv = self.hap_chr_conversion(bf.data.chr_pattern)
        df = bf.read(ipath)
        conv_ = hap.from_either(conv[0], conv[1])
        return bed.filter_sort_bed(conv_.init_mapper, conv_.final_mapper, df)

    def read_write_filter_sort_dip2_bed(
        self,
        bf: Dip2BedFile,
        ipath: Path,
        opath: Path,
        hap: Haplotype,
        g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    ) -> None:
        df = self.read_filter_sort_dip2_bed(bf, ipath, hap)
        bed.write_bed(opath, g(df))


AnyBuildData = HapBuildData | Dip1BuildData | Dip2BuildData


def with_build_data(
    bd: AnyBuildData,
    hap_f: Callable[[HapBuildData], X],
    dip1_f: Callable[[Dip1BuildData], X],
    dip2_f: Callable[[Dip2BuildData], X],
) -> X:
    if isinstance(bd, HapBuildData):
        return hap_f(bd)
    elif isinstance(bd, Dip1BuildData):
        return dip1_f(bd)
    elif isinstance(bd, Dip2BuildData):
        return dip2_f(bd)
    else:
        assert_never(bd)


def build_data_to_src(
    bd: AnyBuildData,
    hap: Haplotype | None,
    f: BuildDataToSrc,
) -> BedSrc:
    src = with_build_data(bd, lambda bd: f(bd), lambda bd: f(bd), lambda bd: f(bd))
    # TODO mypy doens't like me using my 'maybe' functional functions
    if src is None:
        raise DesignError()
    return from_hap_or_dip(src, hap)


def to_ref_data_unsafe(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
    rk: RefKeyT,
) -> RefData_[RefKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT]:
    rd = to_ref_data(xs, rk)
    if rd is None:
        raise DesignError(f"Could not get ref data for key '{rk}'")
    return rd


def to_ref_data(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
    rk: RefKeyT,
) -> RefData_[RefKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT] | None:
    try:
        s = xs[rk]
        return RefData_(rk, s.ref, s.strat_inputs, s.builds)
    except KeyError:
        return None


def all_ref_data(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
) -> list[RefData_[RefKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT]]:
    return [to_ref_data_unsafe(xs, rk) for rk in xs]


def all_ref_keys(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
) -> list[RefKeyT]:
    return [r.refkey for r in all_ref_data(xs)]


def all_ref_refsrckeys(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
) -> list[str]:
    return [s for k, v in xs.items() for s in v.ref.src.to_str_refkeys(k)]


def all_build_data(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
) -> list[BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT]]:
    return [r.to_build_data_unsafe(b) for r in all_ref_data(xs) for b in r.builds]


def all_bed_build_and_refsrckeys(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
    f: BuildDataToSrc,
) -> list[tuple[str, str]]:
    return [
        (rk, str(b.buildkey))
        for b in all_build_data(xs)
        if (src := f(b)) is not None
        for rk in src.to_str_refkeys(b.refdata.refkey)
    ]


def all_bed_refsrckeys(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
    f: BuildDataToSrc,
) -> list[str]:
    return [rk for rk, _ in all_bed_build_and_refsrckeys(xs, f)]


def all_build_keys(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
) -> list[tuple[RefKeyT, BuildKeyT]]:
    return [(r.refdata.refkey, r.buildkey) for r in all_build_data(xs)]


def all_ref_build_keys(
    xs: dict[
        RefKeyT,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT, BuildKeyT],
    ],
) -> list[tuple[str, str]]:
    return [
        (rk, str(r.buildkey))
        for r in all_build_data(xs)
        for rk in r.refdata.ref.src.to_str_refkeys(r.refdata.refkey)
    ]


HapStrat = Stratification[
    HapChrSource[RefSrc],
    HapBedSrc,
    HapBedSrc,
    Haploid_[BedSrc],
    HapBuildKey,
]
Dip1Strat = Stratification[
    DipChrSource1[RefSrc],
    DipBedSrc,
    Dip1BedSrc,
    Diploid_[BedSrc],
    Dip1BuildKey,
]
Dip2Strat = Stratification[
    DipChrSource2[RefSrc],
    DipBedSrc,
    Dip2BedSrc,
    Diploid_[BedSrc],
    Dip2BuildKey,
]
AnyStratification = HapStrat | Dip1Strat | Dip2Strat

# HaploidStratDict = StratDict_[
#     HapRefKey_,
#     HapChrSource[RefSrc],
#     HapBedSrc,
#     HapBedSrc,
#     Haploid_[BedSrc],
#     HapBuildKey,
#     HapInclude,
# ]

# Diploid1StratDict = StratDict_[
#     Dip1RefKey_,
#     DipChrSource1[RefSrc],
#     DipBedSrc,
#     Dip1BedSrc,
#     Diploid_[BedSrc],
#     Dip1BuildKey,
#     DipInclude,
# ]

# Diploid2StratDict = StratDict_[
#     Dip2RefKey_,
#     DipChrSource2[RefSrc],
#     DipBedSrc,
#     Dip2BedSrc,
#     Diploid_[BedSrc],
#     Dip2BuildKey,
#     DipInclude,
# ]


def sub_output_path(pat: str, rk: RefKeyFull[RefKeyT]) -> Path:
    if "{" in pat or "}" in pat:
        raise DesignError(f"not all wildcards replaced in pattern {pat}")
    return Path(pat.replace("%s", rk.name))


def prepare_output_path(path: Path) -> Path:
    return Path(str(path).replace("{ref_key}", "%s"))


class DataLogDirs(NamedTuple):
    data: Path
    log: Path


class DataLogBenchDirs(NamedTuple):
    data: Path
    log: Path
    bench: Path


class FilterSortDirs(NamedTuple):
    data: Path
    bench: Path
    log: Path
    subbed: Path


class BedInterDirs(NamedTuple):
    filtersort: FilterSortDirs
    postsort: DataLogBenchDirs


class BedDirs(NamedTuple):
    src: DataLogDirs
    inter: BedInterDirs
    final: Callable[[str], Path]


class RefInterDirs(NamedTuple):
    prebuild: DataLogBenchDirs
    filtersort: FilterSortDirs
    build: DataLogBenchDirs


class RefSrcDirs(NamedTuple):
    benchmark: DataLogDirs
    reference: DataLogDirs


class RefDirs(NamedTuple):
    src: RefSrcDirs
    inter: RefInterDirs


# TODO not sure if we need to break the functions apart this finely
def bd_to_si(
    f: StratInputToBed,
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return f(x.refdata.strat_inputs)


# def bed_to_src(
#     x: BedFile[AnyBedT] | None,
# ) -> Haploid_[BedSrc] | Diploid_[BedSrc] | None:
#     return fmap_maybe(lambda x: x.data.src, x)


def si_to_trf(x: StratInputs_[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.low_complexity.simreps


def bd_to_trf(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_trf(x.refdata.strat_inputs) if x.want_low_complexity else None


def si_to_rmsk(x: StratInputs_[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.low_complexity.rmsk


def bd_to_rmsk(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_rmsk(x.refdata.strat_inputs) if x.want_low_complexity else None


def si_to_satellites(x: StratInputs_[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.low_complexity.satellites


def bd_to_satellites(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_satellites(x.refdata.strat_inputs) if x.want_low_complexity else None


def si_to_superdups(x: StratInputs_[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.segdups.superdups


def bd_to_superdups(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_superdups(x.refdata.strat_inputs) if x.want_segdups else None


def si_to_gaps(x: StratInputs_[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.gap


def bd_to_gaps(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_gaps(x.refdata.strat_inputs) if x.want_gaps else None


def bd_to_other(
    lk: OtherLevelKey,
    sk: OtherStratKey,
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> OtherBedFile[AnyBedT] | None:
    return x.build.other_strats[lk][sk]


def bd_to_bench_bed(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return fmap_maybe(lambda y: y.bench_bed, x.build.bench)


def bd_to_bench_vcf(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> VCFFile[AnyBedT_] | None:
    return fmap_maybe(lambda y: y.bench_vcf, x.build.bench)


def bd_to_query_vcf(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> VCFFile[AnyBedT_] | None:
    return fmap_maybe(lambda y: y.query_vcf, x.build.bench)


def si_to_ftbl(x: StratInputs_[AnyBedT, AnySrcT]) -> AnySrcT | None:
    return fmap_maybe(lambda y: y.ftbl_src, x.functional)


def si_to_gff(x: StratInputs_[AnyBedT, AnySrcT]) -> AnySrcT | None:
    return fmap_maybe(lambda y: y.gff_src, x.functional)


def bd_to_ftbl(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> AnySrcT | None:
    return si_to_ftbl(x.refdata.strat_inputs) if x.want_functional else None


def bd_to_gff(
    x: BuildData_[RefKeyT, BuildKeyT, RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> AnySrcT | None:
    return si_to_gff(x.refdata.strat_inputs) if x.want_functional else None


# TODO add validator to ensure none of the keys in the strat/build dicts overlap
class GiabStrats(BaseModel):
    """Top level stratification object."""

    other_levels: list[OtherLevelKey] = [
        OtherLevelKey("Ancestry"),
        OtherLevelKey("FunctionalTechnicallyDifficult"),
        OtherLevelKey("GenomeSpecific"),
        OtherLevelKey("OtherDifficult"),
    ]
    paths: Paths = Paths()
    tools: Tools = Tools()
    # haploid_stratifications: HaploidStratDict = HaploidStratDict()
    # diploid1_stratifications: Diploid1StratDict = Diploid1StratDict()
    # diploid2_stratifications: Diploid2StratDict = Diploid2StratDict()
    haploid_stratifications: dict[HapRefKey_, HapStrat] = {}
    diploid1_stratifications: dict[Dip1RefKey_, Dip1Strat] = {}
    diploid2_stratifications: dict[Dip2RefKey_, Dip2Strat] = {}
    comparison_strats: dict[CompareKey, HttpUrl] = {}
    benchmark_subsets: list[str] = [
        "AllAutosomes",
        "AllTandemRepeats",
        "AllHomopolymers_ge7bp_imperfectge11bp_slop5",
        "SimpleRepeat_diTR_10to49_slop5",
        "SimpleRepeat_homopolymer_7to11_slop5",
        "SimpleRepeat_homopolymer_ge21_slop5",
        "SimpleRepeat_imperfecthomopolge11_slop5",
        "SimpleRepeat_imperfecthomopolge21_slop5",
        "SimpleRepeat_homopolymer_7to11_AT_slop5",
        "SimpleRepeat_homopolymer_ge21_AT_slop5",
        "SimpleRepeat_imperfecthomopolge11_AT_slop5",
        "SimpleRepeat_imperfecthomopolge21_AT_slop5",
        "SimpleRepeat_homopolymer_7to11_GC_slop5",
        "SimpleRepeat_homopolymer_ge21_GC_slop5",
        "SimpleRepeat_imperfecthomopolge11_GC_slop5",
        "SimpleRepeat_imperfecthomopolge21_GC_slop5",
        "alldifficultregions",
        "alllowmapandsegdupregions",
        "chrX_PAR",
        "chrX_XTR",
        "chrY_XTR",
        "notinalldifficultregions",
        "notinAllHomopolymers_ge7bp_imperfectge11bp_slop5",
        "notinAllTandemRepeatsandHomopolymers_slop5",
        "segdups",
    ]

    # @validator("stratifications", each_item=True)
    # def builds_have_valid_existing(
    #     cls,
    #     v: HapStrat,
    #     values: dict[str, Any],
    # ) -> HapStrat:
    #     try:
    #         levels = cast(list[OtherLevelKey], values["other_levels"])
    #         bad = [
    #             f"level='{lk}'; build='{bk}'"
    #             for bk, b in v.builds.items()
    #             for lk in b.other_strats
    #             if lk not in levels
    #         ]
    #         if len(bad) > 0:
    #             assert (
    #                 False
    #             ), f"builds referencing invalid strat categories: {', '.join(bad)}"
    #     except KeyError:
    #         pass
    #     return v

    # TODO extend this to the other types
    @validator("haploid_stratifications", each_item=True)
    def builds_have_valid_old_version(
        cls,
        v: HapStrat,
        values: dict[str, Any],
    ) -> HapStrat:
        try:
            prev = cast(dict[CompareKey, HttpUrl], values["previous"])
            bad = [
                f"version='{pk}'; build='{bk}'"
                for bk, b in v.builds.items()
                if b.comparison is not None
                for pk in b.comparison.other
                if pk not in prev
            ]
            if len(bad) > 0:
                assert (
                    False
                ), f"builds referencing invalid previous version keys: {', '.join(bad)}"
        except KeyError:
            pass
        return v

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
        return self.paths.resources / "{ref_src_key}"

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
        return self.final_root_dir / "{ref_final_key}@{build_key}"

    @property
    def intermediate_root_dir(self) -> Path:
        return self.results_dir / "intermediates"

    @property
    def intermediate_build_hapless_dir(self) -> Path:
        return self.intermediate_root_dir / "{ref_key}@{build_key}"

    @property
    def intermediate_build_dir(self) -> Path:
        return self.intermediate_root_dir / "{ref_final_key}@{build_key}"

    @property
    def log_root_dir(self) -> Path:
        return self.results_dir / "log"

    @property
    def bench_root_dir(self) -> Path:
        return self.results_dir / "bench"

    @property
    def log_src_dir(self) -> Path:
        return self.log_root_dir / "resources" / "{ref_src_key}"

    @property
    def log_results_dir(self) -> Path:
        return self.log_root_dir / "results"

    @property
    def log_build_dir(self) -> Path:
        return self.log_results_dir / "builds" / "{ref_final_key}@{build_key}"

    @property
    def log_build_hapless_dir(self) -> Path:
        return self.log_results_dir / "builds" / "{ref_key}@{build_key}"

    @property
    def bench_build_dir(self) -> Path:
        return self.bench_root_dir / "{ref_final_key}@{build_key}"

    @property
    def bench_build_hapless_dir(self) -> Path:
        return self.bench_root_dir / "{ref_key}@{build_key}"

    def build_final_strat_path(self, level: str, name: str) -> Path:
        return self.final_build_dir / level / f"{{ref_final_key}}_{name}.bed.gz"

    def build_strat_path(self, level: CoreLevel, name: str) -> Path:
        return self.build_final_strat_path(level.value, name)

    @property
    def ref_dirs(self) -> RefDirs:
        return RefDirs(
            src=RefSrcDirs(
                reference=DataLogDirs(
                    data=self.ref_src_dir / "reference",
                    log=self.log_src_dir / "reference",
                ),
                benchmark=DataLogDirs(
                    data=self.ref_src_dir / "benchmark" / "{build_key}",
                    log=self.log_src_dir / "benchmark" / "{build_key}",
                ),
            ),
            inter=RefInterDirs(
                prebuild=DataLogBenchDirs(
                    data=self.intermediate_root_dir / "{ref_final_key}",
                    log=self.log_results_dir / "{ref_final_key}",
                    bench=self.bench_root_dir / "{ref_final_key}",
                ),
                filtersort=FilterSortDirs(
                    data=self.intermediate_build_hapless_dir / "ref",
                    subbed=prepare_output_path(
                        self.intermediate_build_hapless_dir / "ref"
                    ),
                    log=self.log_build_hapless_dir / "ref",
                    bench=self.bench_build_hapless_dir / "ref",
                ),
                build=DataLogBenchDirs(
                    data=self.intermediate_root_dir / "{ref_final_key}@{build_key}",
                    log=self.log_results_dir / "{ref_final_key}@{build_key}",
                    bench=self.bench_root_dir / "{ref_final_key}@{build_key}",
                ),
            ),
        )

    def to_bed_dirs(self, level: CoreLevel) -> BedDirs:
        v = level.value
        return BedDirs(
            src=DataLogDirs(
                self.ref_src_dir / v,
                self.log_src_dir / v,
            ),
            inter=BedInterDirs(
                filtersort=FilterSortDirs(
                    data=self.intermediate_build_hapless_dir / v,
                    log=self.log_build_hapless_dir / v,
                    bench=self.bench_build_hapless_dir / v,
                    subbed=prepare_output_path(self.intermediate_build_hapless_dir / v),
                ),
                postsort=DataLogBenchDirs(
                    data=self.intermediate_build_dir / v,
                    log=self.log_build_dir / v,
                    bench=self.bench_build_dir / v,
                ),
            ),
            final=lambda name: self.build_strat_path(level, name),
        )

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

    def refsrckey_to_ref_src(self, rsk: str) -> RefSrc | None:
        # TODO misleading name
        rk, hap = parse_final_refkey(rsk)
        src = self.to_ref_data(rk).ref.src
        return from_hap_or_dip(src, hap)

    def refkey_to_bed_refsrckeys(self, f: StratInputToBed, rk: str) -> list[str]:
        # TODO this seems like a useful glue function (the labmda that is)
        return self.to_ref_data(rk).get_refkeys_unsafe(
            lambda si: fmap_maybe(lambda x: x.data.src, f(si))
        )

    def refsrckey_to_bed_src(self, f: StratInputToBed, rsk: str) -> BedSrc:
        rk, hap = parse_final_refkey(rsk)
        rd = self.to_ref_data(rk)
        return ref_data_to_src(
            rd, hap, lambda rd: fmap_maybe(lambda x: x.data.src, f(rd))
        )

    def refsrckey_to_x_features_src(self, rsk: str) -> BedSrc:
        # TODO this is confusing
        rk, _ = parse_final_refkey(rsk)
        return not_none_unsafe(
            self.to_ref_data(rk).strat_inputs.xy.features,
            lambda x: x.x_bed.data.src.hap,
        )

    def refsrckey_to_y_features_src(self, rsk: str) -> BedSrc:
        rk, _ = parse_final_refkey(rsk)
        return not_none_unsafe(
            self.to_ref_data(rk).strat_inputs.xy.features,
            lambda x: x.y_bed.data.src.hap,
        )

    def buildkey_to_bed_refsrckeys(
        self, f: BuildDataToBed, rk: str, bk: str
    ) -> list[str]:
        # TODO this "update" function is not DRY
        # return self.refkey_to_bed_refsrckeys(lambda rd: f(rd.to_build_data(bk)), rk)
        return self.to_ref_data(rk).get_refkeys_unsafe_(
            lambda rd: fmap_maybe(lambda x: x.data.src, f(rd.to_build_data(bk)))
        )

    def buildkey_to_bed_src(self, f: BuildDataToBed, rsk: str, bk: str) -> BedSrc:
        # return self.refsrckey_to_bed_src(lambda rd: f(rd.to_build_data(bk)), rsk)
        rk, hap = parse_final_refkey(rsk)
        rd = self.to_ref_data(rk)
        return ref_data_to_src_(
            rd,
            hap,
            lambda rd: fmap_maybe(lambda x: x.data.src, f(rd.to_build_data(bk))),
        )

    def buildkey_to_vcf_src(self, f: BuildDataToVCF, rsk: str, bk: str) -> BedSrc:
        rk, hap = parse_final_refkey(rsk)
        bd = self.to_build_data(rk, bk)
        return build_data_to_src(
            bd,
            hap,
            lambda bd: fmap_maybe(lambda x: x.data.src, f(bd)),
        )

    def refkey_to_functional_refsrckeys(self, f: StratInputToSrc, rk: str) -> list[str]:
        return self.to_ref_data(rk).get_refkeys_unsafe(f)

    def refsrckey_to_functional_src(self, f: StratInputToSrc, rsk: str) -> BedSrc:
        rk, hap = parse_final_refkey(rsk)
        rd = self.to_ref_data(rk)
        return ref_data_to_src(rd, hap, f)

    # # include switches (for controlling which snakemake rules to activate)

    # def has_low_complexity_rmsk(self, rk: RefKey) -> bool:
    #     return self.refkey_to_strat(rk).strat_inputs.low_complexity.rmsk is not None

    # def has_low_complexity_simreps(self, rk: RefKey) -> bool:
    #     return self.refkey_to_strat(rk).strat_inputs.low_complexity.simreps is not None

    # def has_low_complexity_censat(self, rk: RefKey) -> bool:
    #     return (
    #         self.refkey_to_strat(rk).strat_inputs.low_complexity.satellites is not None
    #     )

    # def _want_chr_index(self, p: BuildPair, i: ChrIndex) -> bool:
    #     cis = self.buildkey_to_chr_indices(p)
    #     return i in cis

    # def want_xy_x(self, p: BuildPair) -> bool:
    #     return self._want_chr_index(p, ChrIndex.CHRX) and self.buildkey_to_include(p).xy

    # def want_xy_y(self, p: BuildPair) -> bool:
    #     return self._want_chr_index(p, ChrIndex.CHRY) and self.buildkey_to_include(p).xy

    # def wanted_xy_chr_names(self, p: BuildPair) -> list[str]:
    #     return [
    #         i.chr_name
    #         for i in [ChrIndex.CHRX, ChrIndex.CHRY]
    #         if self._want_chr_index(p, i)
    #     ]

    # # def want_xy_sex(self, rk: RefKey, bk: BuildKey) -> bool:
    # #     return self.want_xy_x(rk, bk) and self.want_xy_y(rk, bk)

    # def want_x_PAR(self, p: BuildPair) -> bool:
    #     return self.want_xy_x(p) and self.refkey_to_x_PAR(p.ref) is not None

    # def want_y_PAR(self, p: BuildPair) -> bool:
    #     return self.want_xy_y(p) and self.refkey_to_y_PAR(p.ref) is not None

    # def want_xy_auto(self, p: BuildPair) -> bool:
    #     cis = self.buildkey_to_chr_indices(p)
    #     return len(cis - set([ChrIndex.CHRX, ChrIndex.CHRY])) > 0

    def want_xy_XTR(self, rfk: str) -> bool:
        rk, _ = parse_final_refkey(rfk)
        f = self.to_ref_data(rk).strat_inputs.xy.features
        return f is not None and f.xtr

    def want_xy_ampliconic(self, rfk: str) -> bool:
        rk, _ = parse_final_refkey(rfk)
        f = self.to_ref_data(rk).strat_inputs.xy.features
        return f is not None and f.ampliconic

    def want_low_complexity(self, rfk: str, bk: str) -> bool:
        rk, _ = parse_final_refkey(rfk)
        return self.to_build_data(rk, bk).build.include.low_complexity

    # def want_gc(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).gc is not None

    # def want_functional(self, p: BuildPair) -> bool:
    #     return (
    #         self.buildkey_to_include(p).functional
    #         and self.refkey_to_functional_ftbl_src(p.ref) is not None
    #         and self.refkey_to_functional_gff_src(p.ref) is not None
    #     )

    # def want_telomeres(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).telomeres

    # def want_segdups(self, p: BuildPair) -> bool:
    #     return (
    #         self.refkey_to_superdups_src(p.ref) is not None
    #         and self.buildkey_to_include(p).segdups
    #     )

    # def _want_union(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_include(p).union

    # def want_mappability(self, p: BuildPair) -> bool:
    #     return (
    #         self.refkey_to_strat(p.ref).strat_inputs.mappability is not None
    #         and len(self.buildkey_to_include(p).mappability) > 0
    #     )

    # def want_segdup_and_map(self, p: BuildPair) -> bool:
    #     return (
    #         self.buildkey_to_include(p).union
    #         and self.want_segdups(p)
    #         and self.want_mappability(p)
    #     )

    # def want_alldifficult(self, p: BuildPair) -> bool:
    #     return (
    #         self.want_segdup_and_map(p)
    #         and self.want_low_complexity(p)
    #         and self.want_gc(p)
    #     )

    # def want_benchmark(self, p: BuildPair) -> bool:
    #     return self.buildkey_to_build(p).bench is not None

    # def want_gaps(self, rk: RefKey) -> bool:
    #     return self.refkey_to_strat(rk).strat_inputs.gap is not None

    # def want_vdj(self, p: BuildPair) -> bool:
    #     cis = self.buildkey_to_chr_indices(p)
    #     vdj_chrs = {ChrIndex(i) for i in [2, 7, 14, 22]}
    #     return (
    #         self.buildkey_to_include(p).vdj
    #         and self.refkey_to_functional_ftbl_src(p.ref) is not None
    #         and self.refkey_to_functional_gff_src(p.ref) is not None
    #         and len(cis & vdj_chrs) > 0
    #     )

    # key lists for downloading resources

    # @property
    # def _all_haploid_refdata(self) -> list[HapRefData]:
    #     return self.haploid_stratifications.all_ref_data

    @property
    def _all_haploid_builds(self) -> list[tuple[HapRefKey, HapBuildKey]]:
        return all_build_keys(self.haploid_stratifications)

    @property
    def _all_diploid1_builds(self) -> list[tuple[Dip1RefKey, Dip1BuildKey]]:
        return all_build_keys(self.diploid1_stratifications)

    @property
    def _all_diploid2_builds(self) -> list[tuple[Dip2RefKey, Dip2BuildKey]]:
        return all_build_keys(self.diploid2_stratifications)

    @property
    def all_build_keys(self) -> tuple[list[str], list[str]]:
        rs, bs = unzip(
            [(str(x), str(y)) for x, y in self._all_haploid_builds]
            + [(str(x), str(y)) for x, y in self._all_diploid1_builds]
            + [(str(x), str(y)) for x, y in self._all_diploid2_builds]
        )
        return (list(rs), list(bs))

    @property
    def all_full_build_keys(self) -> tuple[list[str], list[str]]:
        rs, bs = unzip(
            all_ref_build_keys(self.haploid_stratifications)
            + all_ref_build_keys(self.diploid1_stratifications)
            + all_ref_build_keys(self.diploid2_stratifications)
        )
        return (list(rs), list(bs))

    # @property
    # def all_refdata(self) -> list[AnyRefData]:
    #     # TODO lame....
    #     return (
    #         [HapRefData(**asdict(x)) for x in self.haploid_stratifications.all_ref_data]
    #         + [
    #             Dip1RefData(**asdict(x))
    #             for x in self.diploid1_stratifications.all_ref_data
    #         ]
    #         + [
    #             Dip2RefData(**asdict(x))
    #             for x in self.diploid2_stratifications.all_ref_data
    #         ]
    #     )

    @property
    def all_ref_refsrckeys(self) -> list[str]:
        return (
            all_ref_refsrckeys(self.haploid_stratifications)
            + all_ref_refsrckeys(self.diploid1_stratifications)
            + all_ref_refsrckeys(self.diploid2_stratifications)
        )

    # TODO this seems too messy
    def parse_refkey(self, rk: str) -> AnyRefKey:
        if HapRefKey_(rk) in self.haploid_stratifications:
            return HapRefKey_(rk)
        elif Dip1RefKey_(rk) in self.diploid1_stratifications:
            return Dip1RefKey_(rk)
        elif Dip2RefKey_(rk) in self.diploid2_stratifications:
            return Dip2RefKey_(rk)
        else:
            raise DesignError(f"invalid ref key: '{rk}'")

    def to_ref_data(self, rk: str) -> AnyRefData:
        k = self.parse_refkey(rk)
        if isinstance(k, HapRefKey_):
            return HapRefData(
                (rd := to_ref_data_unsafe(self.haploid_stratifications, k)).refkey,
                rd.ref,
                rd.strat_inputs,
                rd.builds,
            )
        elif isinstance(k, Dip1RefKey_):
            return Dip1RefData(
                (rd0 := to_ref_data_unsafe(self.diploid1_stratifications, k)).refkey,
                rd0.ref,
                rd0.strat_inputs,
                rd0.builds,
            )
        elif isinstance(k, Dip2RefKey_):
            return Dip2RefData(
                (rd1 := to_ref_data_unsafe(self.diploid2_stratifications, k)).refkey,
                rd1.ref,
                rd1.strat_inputs,
                rd1.builds,
            )
        else:
            assert_never(k)

    def with_ref_data(
        self,
        rk: str,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Dip1RefData], X],
        dip2_f: Callable[[Dip2RefData], X],
    ) -> X:
        return with_ref_data(self.to_ref_data(rk), hap_f, dip1_f, dip2_f)

    def to_build_data(self, rk: str, bk: str) -> AnyBuildData:
        # TODO gross...
        def hap(rd: HapRefData) -> AnyBuildData:
            return HapBuildData(
                (bd := rd.to_build_data_unsafe(HapBuildKey(bk))).refdata,
                bd.buildkey,
                bd.build,
            )

        def dip1(rd: Dip1RefData) -> AnyBuildData:
            return Dip1BuildData(
                (bd := rd.to_build_data_unsafe(Dip1BuildKey(bk))).refdata,
                bd.buildkey,
                bd.build,
            )

        def dip2(rd: Dip2RefData) -> AnyBuildData:
            return Dip2BuildData(
                (bd := rd.to_build_data_unsafe(Dip2BuildKey(bk))).refdata,
                bd.buildkey,
                bd.build,
            )

        return with_ref_data(self.to_ref_data(rk), hap, dip1, dip2)

    def with_build_data(
        self,
        rk: str,
        bk: str,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Dip2BuildData], X],
    ) -> X:
        return with_build_data(self.to_build_data(rk, bk), hap_f, dip1_f, dip2_f)

    def with_build_data_final(
        self,
        rfk: str,
        bk: str,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        rk_, hap = parse_final_refkey(rfk)
        return self.with_build_data(
            rk_,
            bk,
            lambda bd: none_unsafe(hap, hap_f(bd)),
            lambda bd: none_unsafe(hap, dip1_f(bd)),
            lambda bd: not_none_unsafe(hap, lambda hap: dip2_f(hap, bd)),
        )

    def with_ref_data_and_bed(
        self,
        rk: str,
        get_bed_f: RefDataToBed,
        hap_f: Callable[[HapRefData, HapBedFile], X],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFile], X],
        dip_1to2_f: Callable[[Dip2RefData, Dip1BedFile], X],
        dip_2to1_f: Callable[[Dip1RefData, Dip2BedFile], X],
        dip_2to2_f: Callable[[Dip2RefData, Dip2BedFile], X],
    ) -> X:
        return self.with_ref_data(
            rk,
            lambda rd: not_none_unsafe(get_bed_f(rd), lambda bd: hap_f(rd, bd)),
            lambda rd: with_dip_bedfile(
                not_none_unsafe(get_bed_f(rd), lambda bd: bd),
                lambda bf: dip_1to1_f(rd, bf),
                lambda bf: dip_2to1_f(rd, bf),
            ),
            lambda rd: with_dip_bedfile(
                not_none_unsafe(get_bed_f(rd), lambda bd: bd),
                lambda bf: dip_1to2_f(rd, bf),
                lambda bf: dip_2to2_f(rd, bf),
            ),
        )

    def with_ref_data_and_bed_hap(
        self,
        rfk: str,
        get_bed_f: RefDataToBed,
        hap_f: Callable[[HapRefData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2RefData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1RefData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2RefData, Dip2BedFile], Z],
    ) -> Z:
        rk, hap = parse_final_refkey(rfk)
        return self.with_ref_data_and_bed(
            rk,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_1to2_f(h, rd, bf)),
            dip_2to1_f,
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_2to2_f(h, rd, bf)),
        )

    def with_ref_data_bed_unsafe(
        self,
        rsk: str,
        get_bed_f: StratInputToBed,
        hap_f: Callable[[HapRefData, HapBedFile], X],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFile], X],
        dip_1to2_f: Callable[[Dip2RefData, Dip1BedFile], X],
        dip_2to1_f: Callable[[Haplotype, Dip1RefData, Dip2BedFile], X],
        dip_2to2_f: Callable[[Haplotype, Dip2RefData, Dip2BedFile], X],
    ) -> X:
        rk, hap = parse_final_refkey(rsk)
        return self.with_ref_data_and_bed(
            rsk,
            lambda rd: get_bed_f(rd.strat_inputs),
            lambda rd, bd: none_unsafe(hap, hap_f(rd, bd)),
            lambda rd, bd: none_unsafe(hap, dip_1to1_f(rd, bd)),
            lambda rd, bd: none_unsafe(hap, dip_1to2_f(rd, bd)),
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_2to1_f(h, rd, bf)),
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_2to2_f(h, rd, bf)),
        )

    def with_build_data_and_bed(
        self,
        rk: str,
        bk: str,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        # TODO gross...
        return self.with_ref_data_and_bed(
            rk,
            lambda rd: get_bed_f(rd.to_build_data_unsafe(bk)),
            lambda rd, bf: hap_f(
                HapBuildData(
                    (bd := rd.to_build_data_unsafe(HapBuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda rd, bf: dip_1to1_f(
                Dip1BuildData(
                    (bd := rd.to_build_data_unsafe(Dip1BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda rd, bf: dip_1to2_f(
                Dip2BuildData(
                    (bd := rd.to_build_data_unsafe(Dip2BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda rd, bf: dip_2to1_f(
                Dip1BuildData(
                    (bd := rd.to_build_data_unsafe(Dip1BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda rd, bf: dip_2to2_f(
                Dip2BuildData(
                    (bd := rd.to_build_data_unsafe(Dip2BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
        )

    def with_build_data_and_bed_hap(
        self,
        rfk: str,
        bk: str,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        return self.with_ref_data_and_bed_hap(
            rfk,
            lambda rd: get_bed_f(rd.to_build_data_unsafe(bk)),
            lambda rd, bf: hap_f(
                HapBuildData(
                    (bd := rd.to_build_data_unsafe(HapBuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda rd, bf: dip_1to1_f(
                Dip1BuildData(
                    (bd := rd.to_build_data_unsafe(Dip1BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda hap, rd, bf: dip_1to2_f(
                hap,
                Dip2BuildData(
                    (bd := rd.to_build_data_unsafe(Dip2BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda rd, bf: dip_2to1_f(
                Dip1BuildData(
                    (bd := rd.to_build_data_unsafe(Dip1BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
            lambda hap, rd, bf: dip_2to2_f(
                hap,
                Dip2BuildData(
                    (bd := rd.to_build_data_unsafe(Dip2BuildKey(bk))).refdata,
                    bd.buildkey,
                    bd.build,
                ),
                bf,
            ),
        )

    def with_build_data_and_bed_i(
        self,
        rk: str,
        bk: str,
        inputs: list[X],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[tuple[X, X], Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[tuple[X, X], Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        return self.with_build_data_and_bed(
            rk,
            bk,
            get_bed_f,
            lambda bd, bf: match1_unsafe(inputs, lambda i: hap_f(i, bd, bf)),
            lambda bd, bf: match1_unsafe(inputs, lambda i: dip_1to1_f(i, bd, bf)),
            lambda bd, bf: match1_unsafe(inputs, lambda i: dip_1to2_f(i, bd, bf)),
            lambda bd, bf: match2_unsafe(
                inputs, lambda i0, i1: dip_2to1_f((i0, i1), bd, bf)
            ),
            lambda bd, bf: match2_unsafe(
                inputs, lambda i0, i1: dip_2to2_f((i0, i1), bd, bf)
            ),
        )

    def with_build_data_and_bed_io(
        self,
        rk: str,
        bk: str,
        inputs: list[X],
        output_f: OutputPattern[Y],
        write_outputs: Callable[[list[Y]], None],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Y, HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Y, Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, tuple[Y, Y], Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[tuple[X, X], Y, Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[tuple[X, X], tuple[Y, Y], Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        def out1(src: Haploid_[RefSrc], rk: RefKeyT) -> Y:
            return with_first(output_f(src.key(rk)), lambda o: write_outputs([o]))

        def out2(src: Diploid_[RefSrc], rk: RefKeyT) -> tuple[Y, Y]:
            return with_first(
                both(output_f, src.keys(rk)), lambda o: write_outputs([*o])
            )

        return self.with_build_data_and_bed_i(
            rk,
            bk,
            inputs,
            get_bed_f,
            lambda i, bd, bf: hap_f(
                i, out1(bd.refdata.ref.src, bd.refdata.refkey), bd, bf
            ),
            lambda i, bd, bf: dip_1to1_f(
                i, out1(bd.refdata.ref.src, bd.refdata.refkey), bd, bf
            ),
            lambda i, bd, bf: dip_1to2_f(
                i, out2(bd.refdata.ref.src, bd.refdata.refkey), bd, bf
            ),
            lambda i, bd, bf: dip_2to1_f(
                i, out1(bd.refdata.ref.src, bd.refdata.refkey), bd, bf
            ),
            lambda i, bd, bf: dip_2to2_f(
                i, out2(bd.refdata.ref.src, bd.refdata.refkey), bd, bf
            ),
        )

    # TODO refactor this so it takes only a builddata and not these weird
    # strings, which will make this function independent of the class
    def with_build_data_and_bed_io_(
        self,
        rk: str,
        bk: str,
        inputs: list[X],
        output: Path,
        output_pattern: str,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Path, HapBuildData, HapBedFile], None],
        dip_1to1_f: Callable[[X, Path, Dip1BuildData, Dip1BedFile], None],
        dip_1to2_f: Callable[[X, tuple[Path, Path], Dip2BuildData, Dip1BedFile], None],
        dip_2to1_f: Callable[[tuple[X, X], Path, Dip1BuildData, Dip2BedFile], None],
        dip_2to2_f: Callable[[X, Path, Haplotype, Dip2BuildData, Dip2BedFile], None],
    ) -> None:
        def _dip_2to2_f(
            i: tuple[X, X],
            o: tuple[Path, Path],
            bd: Dip2BuildData,
            bf: Dip2BedFile,
        ) -> None:
            dip_2to2_f(i[0], o[0], Haplotype.HAP1, bd, bf)
            dip_2to2_f(i[1], o[1], Haplotype.HAP2, bd, bf)

        def write_output(ps: list[Path]) -> None:
            with open(output, "w") as f:
                json.dump([str(p) for p in ps], f)

        self.with_build_data_and_bed_io(
            rk,
            bk,
            inputs,
            lambda rk: sub_output_path(output_pattern, rk),
            write_output,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            dip_1to2_f,
            dip_2to1_f,
            _dip_2to2_f,
        )

    def _all_bed_build_and_refsrckeys(self, f: BuildDataToSrc) -> list[tuple[str, str]]:
        return (
            all_bed_build_and_refsrckeys(self.haploid_stratifications, f)
            + all_bed_build_and_refsrckeys(self.diploid1_stratifications, f)
            + all_bed_build_and_refsrckeys(self.diploid2_stratifications, f)
        )

    def _all_bed_refsrckeys(self, f: BuildDataToSrc) -> list[str]:
        return [rk for rk, _ in self._all_bed_build_and_refsrckeys(f)]

    @property
    def all_refkey_gap(self) -> list[str]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_gaps(bd))
        )

    @property
    def all_refkey_rmsk(self) -> list[str]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_rmsk(bd))
        )

    @property
    def all_refkey_trf(self) -> list[str]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_trf(bd))
        )

    @property
    def all_refkey_censat(self) -> list[str]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_satellites(bd))
        )

    @property
    def all_refkey_segdups(self) -> list[str]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_superdups(bd))
        )

    @property
    def all_buildkey_bench_bed(self) -> list[tuple[str, str]]:
        return self._all_bed_build_and_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_bench_bed(bd))
        )

    @property
    def all_buildkey_bench_vcf(self) -> list[tuple[str, str]]:
        return self._all_bed_build_and_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_bench_vcf(bd))
        )

    @property
    def all_buildkey_query_vcf(self) -> list[tuple[str, str]]:
        return self._all_bed_build_and_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_query_vcf(bd))
        )

    @property
    def all_refkey_ftbl(self) -> list[str]:
        return self._all_bed_refsrckeys(bd_to_ftbl)

    @property
    def all_refkey_gff(self) -> list[str]:
        return self._all_bed_refsrckeys(bd_to_gff)


def filter_sort_bed_main(
    f: BuildDataToBed,
    smk: Any,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a bed and filter/sort it appropriately.

    This is meant to be called in snakemake scripts, as this operations is
    very common. 'smk' is the snakemake object and 'f' is a function to
    retrieve the bed configuration from the config instance (which will be
    obtained from the snakemake object).
    """
    sconf: GiabStrats = smk.config
    ws: dict[str, str] = smk.wildcards

    if not isinstance((ins := smk.input), list) and not all(
        isinstance(x, str) for x in ins
    ):
        raise DesignError(f"Inputs must be a list of strings, got {ins}")

    if not isinstance(output_pattern := smk.params["output_pattern"], str):
        raise DesignError(f"Output pattern must be a string, got {output_pattern}")

    sconf.with_build_data_and_bed_io_(
        ws["ref_key"],
        ws["build_key"],
        [Path(i) for i in ins],
        smk.output[0],
        output_pattern,
        f,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, g),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip1_bed(b, i, o, g),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip1_bed(b, i, o, g, g),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip2_bed(b, i, o, g),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_dip2_bed(b, i, o, hap, g),
    )
