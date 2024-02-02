"""Master configuration definition for the entire pipeline :)

Overview:

The entire pipeline configuration is defined here as a pydantic model, which
validates and types all data prior to downstream processing in snakemake rules
or in scripts called by snakemake. Since pydantic interfaces with mypy, and
since mypy can be used to lint python code for correctness, the logic of this
pipeline is heavily biased toward python code rather than snakemake.

In general, snakemake is used to handle build dependencies, and almost all other
logic is handled via python scripts which are type-checked "at the boundaries.

Data hierarchy and flow:

Each target reference is index by a key called a ref key ('ref_key'). Each
reference can have one of more builds indexed by a build key ('build_key') which
describes what should be included (chromosome numbers and various stratification
levels) within a given stratification 'build' for the target reference under
which it is located. The number of builds is found by the Cartesian product of
all ref keys and build keys.

Additionally, ref/build keys are grouped under one of three categories
corresponding to the way its haplotypes are configured: haploid, diploid1, and
diploid2 (see next section). In the case of diploid2, each haplotype is spread
across two files. For this reason (and others), we further distinguish the ref
key into a "full ref key" which may have a haplotype associated with it. A full
ref key associated with the reference is also called a "final ref key" and one
associated with an input is called a "source ref key."

Any given reference (indexed by the ref key) can have input files associated
with it; these will be used to make some of the stratifications. In general
these are bed files or files that can be coerced into being bed files by
selecting certain columns. Note that in the case of diploid assemblies, these
inputs may have to be merged or split depending on if the target reference and
the input file have each haplotype in one file or two files (see next section).

Haploid vs Diploid:

This pipeline is meant to generate stratifications for human genomes which may
either be haploid or diploid. In the case of diploid references, there are two
cases to consider: dip1 or dip2 (see below for terminology). In the case of dip1
we only output one final directory since the two haplotypes are in one file (and
presumably will be consumed as such). In the case of dip2, we make two output
directories (one per haplotype) since the haplotypes are spread across two
reference files. Thus from the reference perspective (ie the fasta for which the
stratifications must be built) we have 3 cases to consider: dip1, dip2, or hap.

To further complicate matters, the inputs for the dip1 and dip2 cases may
themselves be some combination of dip1 or dip2 (ie a reference might be dip1 but
a bed file for the reference might be dip2). Thus for input files, we have 5
cases to consider:
* hap -> hap: input and reference are both haploid
* dip1 -> dip1: input and reference are both dip1
* dip1 -> dip2: input is dip1, reference is dip2
* dip2 -> dip2: input is dip2, reference is dip1
* dip2 -> dip2: input and reference are both dip2

In the case of dip1 -> dip2, the input must be split. In the case of
dip2 -> dip1, the inputs (two) must be combined. For the other three, each input
can be used more or less as-is since the cardinality of the inputs and reference
matches. To make the dip2 -> dip1 case even more complex, the chromosome names
might be identical across the two input files, but will need to be distinguished
when combined by adding a suffix to them.

To make this process as simple as possible, input bed files (and related) will
be downloaded and processed immediately into bed output that correspond directly
to their target reference. This post-download processing step is called
"normalization" and generally consists of the following:
* filters desired chromosomes
* sorts all chromosomes by numeric index (with X and Y being 23 and 24, and the
  paternal haplotype sorted prior to the maternal haplotype when applicable)
* splits or combines the input as described above when applicable
* renames all chromosomes to match the target reference
* rearranges all columns into proper bed format

Each normalization is also generally a checkpoint, since it is not known until
runtime how many files it needs to consume or produce. (Note that this is a
design choice to reduce repeated code. It is technically feasible to produce
snakemake rules and scripts that individual handle each of the 5 cases in a way
that doesn't require checkpoints, but these 5 rules would need to be repeated
for each input file which would be more error prone than the solution chosen
here)

Terminology:
* haploid (or "hap"): half a diploid dataset, ie one haplotype
* diploid1 (or "dip1"): a diploid dataset that is all in one file
* diploid2 (or "dip2"): a diploid dataset that is in two files

Conventions:
* 'DesignError' exceptions are those that should not happen; if they do the code
  is incorrect
* bed files are gzip'ed, fasta files are bgzip'ed, and all files are compressed
  with one of these after downloading
* paternal sorts before maternal
* chromosomes are numbered and sorted 1-24 where X and Y are 23/24 respectively
"""

from __future__ import annotations
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
from more_itertools import duplicates_everseen
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
    unzip2,
    unzip3,
    noop,
)
from common.io import is_gzip, is_bgzip
import common.bed as bed

################################################################################
# Constants

CHR_INDEX_PLACEHOLDER = "%i"
CHR_HAP_PLACEHOLDER = "%h"


################################################################################
# Type aliases

Percent = Annotated[int, Field(ge=0, le=100)]

RefKey = NewType("RefKey", str)
# full refkey represented as a string (see below for class)
RefKeyFullS = NewType("RefKeyFullS", str)
BuildKey = NewType("BuildKey", str)
CompareKey = NewType("CompareKey", str)
OtherLevelKey = NewType("OtherLevelKey", str)
OtherStratKey = NewType("OtherStratKey", str)
HaplotypeName = NewType("HaplotypeName", str)

# helper type representing snakemake wildcards
SmkWildcards = dict[str, Any]

GCBound = tuple[Percent, bool]

################################################################################
# Type variables

W = TypeVar("W")
X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


################################################################################
# Helper functions


def flip_hap(h: Haplotype) -> Haplotype:
    return h.from_either(Haplotype.HAP2, Haplotype.HAP1)


def parse_full_refkey_class(s: RefKeyFullS) -> RefKeyFull:
    m = re.match("(.+)\\.(hap[12])", s)
    # ASSUME this will never fail due to the hap1/2 permitted match pattern
    rk, hap = (s, None) if m is None else (m[1], Haplotype.from_name(m[2]))
    return RefKeyFull(RefKey(rk), hap)


def parse_full_refkey(s: RefKeyFullS) -> tuple[RefKey, Haplotype | None]:
    return parse_full_refkey_class(s).as_tuple


def strip_full_refkey(s: RefKeyFullS) -> RefKey:
    return parse_full_refkey(s)[0]


def flip_full_refkey_class(r: RefKeyFull) -> RefKeyFull:
    return RefKeyFull(r.key, fmap_maybe(flip_hap, r.hap))


def flip_full_refkey(s: RefKeyFullS) -> RefKeyFullS:
    return flip_full_refkey_class(parse_full_refkey_class(s)).name


def choose_xy_unsafe(c: "ChrIndex", x_res: X, y_res: X) -> X:
    if c is ChrIndex.CHRX:
        return x_res
    elif c is ChrIndex.CHRY:
        return y_res
    else:
        raise DesignError(f"I am not an X or Y, I am a {c}")


# type helpers


# TODO mypy for some reason doesn't understand how to narrow a
# Something[Union[X, Y]] to a Something[X] using 'isinstance'
def is_dip1_bed(
    x: BedFile[Dip1ChrSrc[X] | Dip2ChrSrc[X]],
) -> TypeGuard[BedFile[Dip1ChrSrc[X]]]:
    return isinstance(x.data, Dip1ChrSrc)


def is_dip2_bed(
    x: BedFile[Dip1ChrSrc[X] | Dip2ChrSrc[X]],
) -> TypeGuard[BedFile[Dip2ChrSrc[X]]]:
    return isinstance(x.data, Dip2ChrSrc)


# union detanglers


def with_dip_bedfile(
    bf: BedFile[Dip1ChrSrc[X] | Dip2ChrSrc[X]],
    dip1: Callable[[BedFile[Dip1ChrSrc[X]]], Y],
    dip2: Callable[[BedFile[Dip2ChrSrc[X]]], Y],
) -> Y:
    if is_dip1_bed(bf):
        return dip1(bf)
    elif is_dip2_bed(bf):
        return dip2(bf)
    else:
        # TODO this is a mypy bug, I should be able to use assert_never here
        raise DesignError("not a dip1 or dip2")
        # assert_never(bf)


def with_hap_or_dip(
    x: Haploid[X] | Diploid[X],
    dip1: Callable[[Haploid[X]], Y],
    dip2: Callable[[Diploid[X]], Y],
) -> Y:
    if isinstance(x, Diploid):
        return dip2(x)
    if isinstance(x, Haploid):
        return dip1(x)
    else:
        assert_never(x)


def from_hap_or_dip(x: Haploid[X] | Diploid[X], hap: Haplotype | None) -> X:
    return with_hap_or_dip(
        x,
        lambda x: none_unsafe(hap, x.hap),
        lambda x: not_none_unsafe(hap, lambda h: x.from_either(h)),
    )


def with_ref_data(
    rd: AnyRefData,
    hap_f: Callable[[HapRefData], X],
    dip1_f: Callable[[Dip1RefData], X],
    dip2_f: Callable[[Dip2RefData], X],
) -> X:
    if isinstance(rd.ref, HapChrSrc):
        return hap_f(rd)
    elif isinstance(rd.ref, Dip1ChrSrc):
        return dip1_f(rd)
    elif isinstance(rd.ref, Dip2ChrSrc):
        return dip2_f(rd)
    else:
        assert_never(rd)


def with_build_data(
    bd: AnyBuildData,
    hap_f: Callable[[HapBuildData], X],
    dip1_f: Callable[[Dip1BuildData], X],
    dip2_f: Callable[[Dip2BuildData], X],
) -> X:
    if isinstance(bd.refdata.ref, HapChrSrc):
        return hap_f(bd)
    elif isinstance(bd.refdata.ref, Dip1ChrSrc):
        return dip1_f(bd)
    elif isinstance(bd.refdata.ref, Dip2ChrSrc):
        return dip2_f(bd)
    else:
        assert_never(bd)


# functions for dealing with 'dict[RefKey, X' type things


def to_ref_data_unsafe(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
    rk: RefKey,
) -> RefData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]:
    try:
        s = xs[rk]
        return RefData_(rk, s.ref, s.strat_inputs, s.builds)
    except KeyError:
        raise DesignError(f"Could not get ref data for key '{rk}'")


def all_ref_data(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
) -> list[RefData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]]:
    return [to_ref_data_unsafe(xs, rk) for rk in xs]


def all_ref_refsrckeys(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
) -> list[RefKeyFullS]:
    return [s for k, v in xs.items() for s in v.ref.src.to_str_refkeys(k)]


def all_build_data(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
) -> list[BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]]:
    return [r.to_build_data_unsafe(b) for r in all_ref_data(xs) for b in r.builds]


def all_bed_build_and_refsrckeys(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
    f: BuildDataToSrc,
) -> list[tuple[RefKeyFullS, BuildKey]]:
    return [
        (rk, b.buildkey)
        for b in all_build_data(xs)
        if (src := f(b)) is not None
        for rk in src.to_str_refkeys(b.refdata.refkey)
    ]


def all_bed_refsrckeys(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
    f: BuildDataToSrc,
) -> list[RefKeyFullS]:
    return [rk for rk, _ in all_bed_build_and_refsrckeys(xs, f)]


def all_build_keys(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
) -> list[tuple[RefKey, BuildKey]]:
    return [(r.refdata.refkey, r.buildkey) for r in all_build_data(xs)]


def all_ref_build_keys(
    xs: dict[
        RefKey,
        Stratification[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
    ],
) -> list[tuple[RefKeyFullS, BuildKey]]:
    return [
        (rk, r.buildkey)
        for r in all_build_data(xs)
        for rk in r.refdata.ref.src.to_str_refkeys(r.refdata.refkey)
    ]


# path formatters


def sub_output_path(pat: str, rk: RefKeyFull) -> Path:
    if "{" in pat or "}" in pat:
        raise DesignError(f"not all wildcards replaced in pattern {pat}")
    return Path(pat.replace("%s", rk.name))


def prepare_output_path(path: Path) -> Path:
    return Path(str(path).replace("{ref_key}", "%s"))


# itty bitty accessor functions


def bd_to_si(
    f: StratInputToBed,
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return f(x.refdata.strat_inputs)


def si_to_simreps(x: StratInputs[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.low_complexity.simreps


def bd_to_simreps(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_simreps(x.refdata.strat_inputs) if x.want_low_complexity else None


def si_to_rmsk(x: StratInputs[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.low_complexity.rmsk


def bd_to_rmsk(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_rmsk(x.refdata.strat_inputs) if x.want_low_complexity else None


def si_to_satellites(x: StratInputs[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.low_complexity.satellites


def bd_to_satellites(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_satellites(x.refdata.strat_inputs) if x.want_low_complexity else None


def si_to_superdups(x: StratInputs[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.segdups.superdups


def bd_to_superdups(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_superdups(x.refdata.strat_inputs) if x.want_segdups else None


def si_to_gaps(x: StratInputs[AnyBedT, AnySrcT]) -> BedFile[AnyBedT] | None:
    return x.gap


def bd_to_gaps(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return si_to_gaps(x.refdata.strat_inputs) if x.want_gaps else None


def bd_to_other(
    lk: OtherLevelKey,
    sk: OtherStratKey,
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> OtherBedFile[AnyBedT] | None:
    return x.build.other_strats[lk][sk]


def bd_to_bench_bed(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> BedFile[AnyBedT] | None:
    return fmap_maybe(lambda y: y.bench_bed, x.build.bench)


def bd_to_bench_vcf(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> VCFFile[AnyBedT_] | None:
    return fmap_maybe(lambda y: y.bench_vcf, x.build.bench)


def bd_to_query_vcf(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> VCFFile[AnyBedT_] | None:
    return fmap_maybe(lambda y: y.query_vcf, x.build.bench)


def si_to_ftbl(x: StratInputs[AnyBedT, AnySrcT]) -> AnySrcT | None:
    return fmap_maybe(lambda y: y.ftbl_src, x.functional)


def si_to_gff(x: StratInputs[AnyBedT, AnySrcT]) -> AnySrcT | None:
    return fmap_maybe(lambda y: y.gff_src, x.functional)


def bd_to_ftbl(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> AnySrcT | None:
    return si_to_ftbl(x.refdata.strat_inputs) if x.want_functional else None


def bd_to_gff(
    x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT],
) -> AnySrcT | None:
    return si_to_gff(x.refdata.strat_inputs) if x.want_functional else None


# snakemake wildcard helpers


def wc_lookup(ws: SmkWildcards, k: str) -> Any:
    try:
        return ws[k]
    except KeyError:
        raise DesignError(f"Could not find {k} in wildcards")


def wc_to_refkey(ws: SmkWildcards) -> RefKey:
    return RefKey(wc_lookup(ws, "ref_key"))


def wc_to_buildkey(ws: SmkWildcards) -> BuildKey:
    return BuildKey(wc_lookup(ws, "build_key"))


def wc_to_reffinalkey(ws: SmkWildcards) -> RefKeyFullS:
    return RefKeyFullS(wc_lookup(ws, "ref_final_key"))


# IO functions for processing bed files of various flavors


def read_filter_sort_hap_bed(
    bd: HapBuildData, bf: HapBedFile, ipath: Path
) -> pd.DataFrame:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    conv = bd.refdata.ref.chr_conversion(bf.data.chr_pattern, bd.chr_indices)
    df = bf.read(ipath)
    return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)


def read_write_filter_sort_hap_bed(
    ipath: Path,
    opath: Path,
    bd: HapBuildData,
    bf: HapBedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df = read_filter_sort_hap_bed(bd, bf, ipath)
    bed.write_bed(opath, g(df))


def read_filter_sort_dip1to1_bed(
    bd: Dip1BuildData,
    bf: Dip1BedFile,
    ipath: Path,
) -> pd.DataFrame:
    """Read a diploid bed file, sort it, and write it in bgzip format."""
    conv = bd.refdata.ref.dip_chr_conversion(bf.data.chr_pattern, bd.chr_indices)
    df = bf.read(ipath)
    return bed.filter_sort_bed(conv.init_mapper, conv.final_mapper, df)


def read_write_filter_sort_dip1to1_bed(
    ipath: Path,
    opath: Path,
    bd: Dip1BuildData,
    bf: Dip1BedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df = read_filter_sort_dip1to1_bed(bd, bf, ipath)
    bed.write_bed(opath, g(df))


def read_filter_sort_dip2to1_bed(
    bd: Dip1BuildData,
    bf: Dip2BedFile,
    ipath: tuple[Path, Path],
) -> pd.DataFrame:
    """Read two haploid bed files, combine and sort them as diploid, and write
    it in bgzip format.
    """

    def go(b: Dip2BedFile, i: Path, imap: bed.InitMapper) -> pd.DataFrame:
        df = b.read(i)
        return bed.filter_sort_bed(imap, fmap, df)

    conv = bd.refdata.ref.hap_chr_conversion(bf.data.chr_pattern, bd.chr_indices)
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


def read_write_filter_sort_dip2to1_bed(
    ipath: tuple[Path, Path],
    opath: Path,
    bd: Dip1BuildData,
    bf: Dip2BedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df = read_filter_sort_dip2to1_bed(bd, bf, ipath)
    bed.write_bed(opath, g(df))


def read_filter_sort_dip1to2_bed(
    bd: Dip2BuildData,
    bf: Dip1BedFile,
    ipath: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    conv = bd.refdata.ref.dip_chr_conversion(bf.data.chr_pattern, bd.chr_indices)
    imap, splitter = conv.init_mapper
    fmap0, fmap1 = conv.final_mapper

    def go(df: pd.DataFrame, fmap: bed.FinalMapper) -> pd.DataFrame:
        return bed.filter_sort_bed(imap, fmap, df)

    df = bf.read(ipath)
    df0, df1 = bed.split_bed(splitter, df)
    return (go(df0, fmap0), go(df1, fmap1))


def read_write_filter_sort_dip1to2_bed(
    ipath: Path,
    opath: tuple[Path, Path],
    bd: Dip2BuildData,
    bf: Dip1BedFile,
    g0: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
    g1: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    """Read a haploid bed file, sort it, and write it in bgzip format."""
    df0, df1 = read_filter_sort_dip1to2_bed(bd, bf, ipath)
    bed.write_bed(opath[0], g0(df0))
    bed.write_bed(opath[1], g1(df1))


def read_filter_sort_dip2to2_bed(
    bd: Dip2BuildData,
    bf: Dip2BedFile,
    ipath: Path,
    hap: Haplotype,
) -> pd.DataFrame:
    conv = bd.refdata.ref.hap_chr_conversion(bf.data.chr_pattern, bd.chr_indices)
    df = bf.read(ipath)
    conv_ = hap.from_either(conv[0], conv[1])
    return bed.filter_sort_bed(conv_.init_mapper, conv_.final_mapper, df)


def read_write_filter_sort_dip2to2_bed(
    ipath: Path,
    opath: Path,
    hap: Haplotype,
    bd: Dip2BuildData,
    bf: Dip2BedFile,
    g: Callable[[pd.DataFrame], pd.DataFrame] = lambda x: x,
) -> None:
    df = read_filter_sort_dip2to2_bed(bd, bf, ipath, hap)
    bed.write_bed(opath, g(df))


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
    ws: SmkWildcards = smk.wildcards

    if not isinstance((ins := smk.input), list) and not all(
        isinstance(x, str) for x in ins
    ):
        raise DesignError(f"Inputs must be a list of strings, got {ins}")

    if not isinstance(output_pattern := smk.params["output_pattern"], str):
        raise DesignError(f"Output pattern must be a string, got {output_pattern}")

    sconf.with_build_data_and_bed_io(
        wc_to_refkey(ws),
        wc_to_buildkey(ws),
        [Path(i) for i in ins],
        smk.output[0],
        output_pattern,
        f,
        read_write_filter_sort_hap_bed,
        read_write_filter_sort_dip1to1_bed,
        read_write_filter_sort_dip1to2_bed,
        read_write_filter_sort_dip2to1_bed,
        read_write_filter_sort_dip2to2_bed,
    )


################################################################################
# Helper classes


@dataclass(frozen=True)
class RefKeyFull:
    """Ref key which may or may not have a haplotype appended to it."""

    key: RefKey
    hap: Haplotype | None

    @property
    def strip(self) -> RefKey:
        return RefKey(self.key)

    @property
    def has_hap(self) -> bool:
        return self.hap is not None

    @property
    def as_tuple(self) -> tuple[RefKey, Haplotype | None]:
        return (self.key, self.hap)

    @property
    def name(self) -> RefKeyFullS:
        k, h = self.as_tuple
        return RefKeyFullS(f"{k}.{h.name}" if h is not None else k)


class Haplotype(Enum):
    "One of the human diploid haplotypes. 0 = Paternal, 1 = Maternal"
    HAP1: int = 0
    HAP2: int = 1

    @classmethod
    def from_name(cls, n: str) -> Self:
        "Build haplotype from a string. Must be exactly 'hap1' or 'hap2'."
        try:
            return next(i for i in cls if i.name == n)
        except StopIteration:
            raise ValueError(f"could make haplotype from name '{n}'")

    @property
    def name(self) -> HaplotypeName:
        return HaplotypeName(f"hap{self.value + 1}")

    def from_either(self, left: X, right: X) -> X:
        "Do either left (pat) or right (mat) depending on the haplotype."
        if self is Haplotype.HAP1:
            return left
        elif self is Haplotype.HAP2:
            return right
        else:
            assert_never(self)


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
        "Build chr index from a string. Must be a valid digit or 'X' or 'Y'"
        try:
            return next(i for i in cls if i.chr_name == n)
        except StopIteration:
            raise ValueError(f"could make chr index from name '{n}'")

    @classmethod
    def from_name_unsafe(cls, n: str) -> Self:
        "Like 'from_name' but raises DesignError"
        try:
            return cls.from_name(n)
        except ValueError as e:
            raise DesignError(e)

    def __init__(self, i: int) -> None:
        "Build chr index from an integer (which must be in [1,24])"
        self.chr_name: str = "X" if i == 23 else ("Y" if i == 24 else str(i))

    def to_internal_index(self, hap: Haplotype) -> bed.InternalChrIndex:
        "Convert this index into an integer corresponding to sort order"
        return bed.InternalChrIndex(hap.value * 24 + self.value - 1)

    # TODO this obviously only makes sense for males
    @property
    def xy_to_hap_unsafe(self) -> Haplotype:
        """
        Convert this index to a haplotype given it is either X or Y.

        Throw DesignError if not X or Y.
        """
        return choose_xy_unsafe(self, Haplotype.HAP2, Haplotype.HAP1)


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
    OTHER_DIFFICULT = "OtherDifficult"
    DIPLOID = "Diploid"


# chromosome name conversions


class _NonDivergentConversion:
    """A chromosome name conversion that doesn't involve a split or merge"""

    @property
    def init_mapper(self) -> bed.InitMapper:
        return NotImplemented

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return NotImplemented


@dataclass(frozen=True)
class HapToHapChrConversion(_NonDivergentConversion):
    fromPattern: HapChrPattern
    toPattern: HapChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices, Haplotype.HAP1)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices, Haplotype.HAP1)


@dataclass(frozen=True)
class DipToDipChrConversion(_NonDivergentConversion):
    fromPattern: DipChrPattern
    toPattern: DipChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> bed.InitMapper:
        return self.fromPattern.init_mapper(self.indices)

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass(frozen=True)
class HapToDipChrConversion:
    fromPattern: Diploid[HapChrPattern]
    toPattern: DipChrPattern
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.InitMapper]:
        return self.fromPattern.both(lambda p, h: p.init_mapper(self.indices, h))

    @property
    def final_mapper(self) -> bed.FinalMapper:
        return self.toPattern.final_mapper(self.indices)


@dataclass(frozen=True)
class DipToHapChrConversion:
    fromPattern: DipChrPattern
    toPattern: Diploid[HapChrPattern]
    indices: set[ChrIndex]

    @property
    def init_mapper(self) -> tuple[bed.InitMapper, bed.SplitMapper]:
        im = self.fromPattern.init_mapper(self.indices)
        fm0 = self.toPattern.hap1.final_mapper(self.indices, Haplotype.HAP1)
        return (im, bed.make_split_mapper(im, fm0))

    @property
    def final_mapper(self) -> tuple[bed.FinalMapper, bed.FinalMapper]:
        return self.toPattern.both(lambda p, h: p.final_mapper(self.indices, h))


# tuples representing file paths for the pipeline


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


################################################################################
# Snakemake configuration model


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


class _Src:
    """Helper class providing means to convert a refkey to a full refkey
    depending one subclass-specific implementation.

    """

    def to_refkeys(self, rk: RefKey) -> list[RefKeyFull]:
        return NotImplemented

    def to_str_refkeys(self, rk: RefKey) -> list[RefKeyFullS]:
        return [k.name for k in self.to_refkeys(rk)]


# TODO this is silly, I don't want to be required to put "hap: blablab"
# whenever I just have "one thing"
class Haploid(GenericModel, Generic[X], _Src):
    """A haploid thing"""

    hap: X

    def to_refkeys(self, rk: RefKey) -> list[RefKeyFull]:
        return [self.key(rk)]

    def key(self, rk: RefKey) -> RefKeyFull:
        return RefKeyFull(rk, None)


class Diploid(GenericModel, Generic[X], _Src):
    """A diploid thing"""

    hap1: X
    hap2: X

    def from_either(self, hap: Haplotype) -> X:
        return hap.from_either(self.hap1, self.hap2)

    def both(self, f: Callable[[X, Haplotype], Y]) -> tuple[Y, Y]:
        return (f(self.hap1, Haplotype.HAP1), f(self.hap2, Haplotype.HAP2))

    def key1(self, rk: RefKey) -> RefKeyFull:
        return RefKeyFull(rk, Haplotype.HAP1)

    def key2(self, rk: RefKey) -> RefKeyFull:
        return RefKeyFull(rk, Haplotype.HAP2)

    def keys(self, rk: RefKey) -> tuple[RefKeyFull, RefKeyFull]:
        return (self.key1(rk), self.key2(rk))

    def to_refkeys(self, rk: RefKey) -> list[RefKeyFull]:
        return list(self.keys(rk))


class ChrPattern:
    """A general chromosome pattern providing interface to convert indices to
    names."""

    def to_names(self, cs: set[ChrIndex]) -> list[str]:
        return NotImplemented


class HapChrPattern(BaseModel, ChrPattern):
    """Chromosome pattern for a haploid file.

    'template' contains a placeholder for the chromosome index.
    """

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
    """Chromosome pattern for a haploid file.

    'template' contains placeholders for both the bare chromosome name/index
    and the haplotype, which maps to a specific haplotype index via 'hapnames'.
    """

    template: str = "chr%i_%h"
    special: dict[ChrIndex, str] = {}
    hapnames: Diploid[HaplotypeName] = Diploid(
        hap1=HaplotypeName("PATERNAL"),
        hap2=HaplotypeName("MATERNAL"),
    )
    # By default, paternal doesn't have X and maternal doesn't have Y
    exclusions: Diploid[list[ChrIndex]] = Diploid(
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


class HapChrSrc(GenericModel, Generic[X]):
    """Specification for a haploid source file."""

    src: Haploid[X]
    chr_pattern: HapChrPattern = HapChrPattern()

    def chr_conversion(
        self, fromChr: HapChrPattern, cis: set[ChrIndex]
    ) -> HapToHapChrConversion:
        """Create a chromosome names conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return HapToHapChrConversion(fromChr, self.chr_pattern, cis)

    def noop_conversion(self, cis: set[ChrIndex]) -> HapToHapChrConversion:
        """Create a chromosome conversion for this source itself.

        Useful for anything generated via the reference itself, since those
        will have chromosome names corresponding directly to the reference.
        """
        return self.chr_conversion(self.chr_pattern, cis)


class Dip1ChrSrc(GenericModel, Generic[X]):
    """Specification for a combined diploid source file.

    The 'src' is assumed to have all chromosomes for both haplotypes in one
    file, which implies they are labeled so as to distinguish the haps. The
    pattern will match both the chromosome number and the haplotype within the
    chromosome name.
    """

    src: Haploid[X]
    chr_pattern: DipChrPattern = DipChrPattern()

    def hap_chr_conversion(
        self,
        fromChr: Diploid[HapChrPattern],
        cis: set[ChrIndex],
    ) -> HapToDipChrConversion:
        """Create a dip2->dip1 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return HapToDipChrConversion(fromChr, self.chr_pattern, cis)

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
        cis: set[ChrIndex],
    ) -> DipToDipChrConversion:
        """Create a dip1->dip1 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return DipToDipChrConversion(fromChr, self.chr_pattern, cis)

    def noop_conversion(self, cis: set[ChrIndex]) -> DipToDipChrConversion:
        """Create a chromosome conversion for this source itself.

        Useful for anything generated via the reference itself, since those
        will have chromosome names corresponding directly to the reference.
        """
        return self.dip_chr_conversion(self.chr_pattern, cis)


class Dip2ChrSrc(GenericModel, Generic[X]):
    """Specification for split diploid source files.

    Each source may or may not have each haplotype labeled; the identity of each
    haplotype in either source file is determined based on the configuration key
    under which it appears (hap1 or hap2) and the chromosome names for each are
    matched according to its corresponding entry in `chr_pattern`.
    """

    src: Diploid[X]
    chr_pattern: Diploid[HapChrPattern] = Diploid(
        hap1=HapChrPattern(),
        hap2=HapChrPattern(),
    )

    def hap_chr_conversion(
        self,
        fromChr: Diploid[HapChrPattern],
        cis: set[ChrIndex],
    ) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
        """Create a dip2->dip2 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        toChr = self.chr_pattern
        return (
            HapToHapChrConversion(fromChr.hap1, toChr.hap1, cis),
            HapToHapChrConversion(fromChr.hap2, toChr.hap2, cis),
        )

    def dip_chr_conversion(
        self,
        fromChr: DipChrPattern,
        cis: set[ChrIndex],
    ) -> DipToHapChrConversion:
        """Create a dip1->dip2 conversion corresponding to 'fromChr'.

        'cis' is the list of chromosome indices that the conversion will
        consider not excluded.
        """
        return DipToHapChrConversion(fromChr, self.chr_pattern, cis)

    def noop_conversion(
        self, cis: set[ChrIndex]
    ) -> tuple[HapToHapChrConversion, HapToHapChrConversion]:
        """Create a chromosome conversion for this source itself.

        Useful for anything generated via the reference itself, since those
        will have chromosome names corresponding directly to the reference.
        """
        return self.hap_chr_conversion(self.chr_pattern, cis)


class HashedSrc_(BaseModel):
    """A source that may be hashed to verify its integrity"""

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


class Paths(BaseModel):
    """Local build paths for snakemake."""

    resources: Path = Path("resources")
    results: Path = Path("results")


class Tools(BaseModel):
    """Urls for tools to download/build/use in the pipeline."""

    repseq: HttpUrl = "https://github.com/ndwarshuis/repseq/archive/refs/tags/v1.1.0.tar.gz"  # type: ignore
    paftools: HttpUrl = "https://raw.githubusercontent.com/lh3/minimap2/e28a55be86b298708a2a67c924d665a00b8d829c/misc/paftools.js"  # type: ignore
    dipcall_aux: HttpUrl = "https://raw.githubusercontent.com/lh3/dipcall/6bd5d7724699491f215aeb5fb628490ebf2cc3ae/dipcall-aux.js"  # type: ignore
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
    bed_cols - the columns for the bed coordinates
    skip_lines - how many input lines to skip
    sep - column separator regexp (for "beds" with spaces instead of tabs)
    """

    bed_cols: BedColumns = BedColumns()
    skip_lines: NonNegativeInt = 0
    sep: str = "\t"


class BedFile(GenericModel, Generic[X]):
    """Inport specs for a bed-like file."""

    data: X
    params: BedFileParams = BedFileParams()

    def read(self, path: Path) -> pd.DataFrame:
        """Read bed file with params from Path."""
        return self._read(path, [])

    def _read(self, path: Path, more: list[int] = []) -> pd.DataFrame:
        "Read bed file with params from Path, optionally with 'more' columns."
        p = self.params
        return bed.read_bed(path, p.bed_cols.columns, p.skip_lines, p.sep, more)


RefSrc = RefFileSrc | RefHttpSrc
HapRefSrc = Haploid[RefSrc]
DipRefSrc = Diploid[RefSrc]

BedSrc = BedFileSrc | BedHttpSrc

HapBedSrc = HapChrSrc[BedSrc]
DipBedSrc = Dip1ChrSrc[BedSrc] | Dip2ChrSrc[BedSrc]
Dip1BedSrc = Dip1ChrSrc[BedSrc]
Dip2BedSrc = Dip2ChrSrc[BedSrc]

HapBedFile = BedFile[HapBedSrc]
DipBedFile = BedFile[DipBedSrc]
Dip1BedFile = BedFile[Dip1BedSrc]
Dip2BedFile = BedFile[Dip2BedSrc]

# TODO this is for more than just "bed files" (right now it basically means "a
# file that is either not zipped or gzipped but not bgzipped")

AnySrcT = TypeVar("AnySrcT", Haploid[BedSrc], Diploid[BedSrc])

# TODO clean this up with real polymorphism when mypy catches up with Haskell
# 98, see https://github.com/python/typing/issues/548
# TODO this union can't be aliased in python 3.10
AnyBedT = TypeVar("AnyBedT", HapBedSrc, Dip1BedSrc | Dip2BedSrc)
AnyBedT_ = TypeVar("AnyBedT_", HapBedSrc, Dip1BedSrc, Dip2BedSrc)
AnyBedFileT = BedFile[AnyBedT]


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
        """Read a bed file at 'path' on disk and return dataframe"""
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
        c = pattern.to_chr_name(i)
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


class GCParams(BaseModel):
    """The params by which to generate GC stratifications.

    Members:
    low: the lower boundaries to use; for instance, a list like [X, Y, Z] will
         correspond to bed files with GC content <X, X-Y, and Y-Z
    high: reverse of 'low'

    The second part of the bound corresponds to whether the boundary should be
    used to create combined boundary (True means yes). The number of True's in
    each list must equal, and will be matched in inverse order in each list.
    Thus something like low = [(X1, True) (X2, True)] and
    high = [(Y1, True), (Y2, True)] will correspond to two bed files with GC
    content <X1 & >Y2 and <X2 & >Y1.
    """

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
    # Defaulting to a finite set means that the user never needs to specify it
    # if they want this for diploid (which they probably do) and they don't need
    # to care in haploid cases. The only issue would be if the user specified
    # this in the haploid case; it technically should be a validation error
    # since it makes no sense in the case of haploid, but here it is setup to
    # not hurt anything.
    hets: set[int] = {10, 20, 50, 100, 500}


class OtherBedFile(BedFile[AnyBedT], Generic[AnyBedT]):
    """A bed file that is imported with minimal processing and included as-is
    in a given stratification package. Useful for one-off bed files made in a
    somewhat hacky (but documented) manner that I don't feel like enshrining
    via code here.

    If 'remove_gaps' is True, subtract that gaps bed if present. This is the
    only processing done to these files. However, they are still checked for
    correctness during the validation stage.
    """

    remove_gaps: bool = False


class Bench(GenericModel, Generic[AnyBedT, AnyBedT_]):
    """Configuration for benchmark to use when validating stratifications.

    Note: the two vcf files need to have haploid/diploid layouts that correspond
    to the target reference (ie if the reference is dip1, these two must also be
    dip1). This is a limitation of happy/vcfeval, which won't know what to do
    if we given them two files with different haplotypes.

    This restriction doesn't apply to the bed file since we can split/combine
    these are necessary to match the reference.
    """

    # TODO I could probably split/compine the VCFs as well...but that sounds
    # like too much work
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


class Malloc(BaseModel):
    """Manual memory allocations for rules in Mb

    Note, this can obviously be done with snakemake itself, but I want to set
    these per-build if needed.

    These are only needed for "high" memory steps, which really only means
    things that involve sorting, hap.py, minimap, or GEM.
    """

    # mappability steps
    gemIndex: int = 16000  # GEM
    gemMappability: int = 12000  # GEM
    gemToWig: int = 8000  # GEM
    mergeNonunique: int = 4000  # sort

    # normalization steps (all of which involve a sort)
    normalizeRmsk: int = 4000
    normalizeSimreps: int = 4000
    normalizeCensat: int = 4000
    normalizeSuperdups: int = 4000
    normalizeCds: int = 4000
    normalizeOther: int = 4000

    # diploid steps
    crossAlignBreaks: int = 16000  # minimap2
    crossAlignVariants: int = 16000  # minimap2
    filterSortVariantCrossAlignment: int = 16000  # samtools sort

    # happy
    runHappy: int = 24000


class Build(GenericModel, Generic[AnyBedT, AnyBedT_]):
    chr_filter: set[ChrIndex]
    comparison: BuildCompare | None = None
    bench: Bench[AnyBedT, AnyBedT_] | None = None
    other_strats: dict[OtherLevelKey, dict[OtherStratKey, OtherBedFile[AnyBedT]]] = {}
    # TODO if I really want I could validate this such that the user would be
    # politely alerted in case they specify any diploid params for a haploid
    # config.
    include: Include = Include()
    malloc: Malloc | None = Malloc()

    @property
    def compare_key(self) -> CompareKey | None:
        return fmap_maybe(lambda x: x.other, self.comparison)


class Functional(GenericModel, Generic[X]):
    """Configuration for Functional stratifications."""

    ftbl_src: X
    gff_src: X


class StratInputs(GenericModel, Generic[AnyBedT, AnySrcT]):
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


HapStratInputs = StratInputs[HapChrSrc[BedSrc], Haploid[BedSrc]]
DipStratInputs = StratInputs[Dip1ChrSrc[BedSrc] | Dip2ChrSrc[BedSrc], Diploid[BedSrc]]

StratInputT = TypeVar("StratInputT", HapStratInputs, DipStratInputs)

RefSourceT = TypeVar(
    "RefSourceT",
    HapChrSrc[RefSrc],
    Dip1ChrSrc[RefSrc],
    Dip2ChrSrc[RefSrc],
)


@dataclass(frozen=True)
class RefData_(Generic[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]):
    """A helper class corresponding a given reference and its builds.

    This is primarily meant to provide a glue layer b/t the configuration
    structure and the functions that consume data from it. Because the logic of
    looking up a refkey and determining if it is hap/dip1/dip2 is tedious and
    annoying, this type will represent the results of such a lookup and provide
    an interface for downstream processing. It also is typed generically such
    that mypy can make inferences regarding its membership in hap/dip1/dip2.

    """

    refkey: RefKey
    ref: RefSourceT
    strat_inputs: StratInputs[AnyBedT, AnySrcT]
    builds: dict[BuildKey, Build[AnyBedT, AnyBedT_]]

    @property
    def ref_refkeys(self) -> list[RefKeyFull]:
        "The list of full refkeys for the reference (either one or two)"
        return self.ref.src.to_refkeys(self.refkey)

    @property
    def ref_str_refkeys(self) -> list[RefKeyFullS]:
        "Like 'ref_refkeys' but returns strings."
        return self.ref.src.to_str_refkeys(self.refkey)

    @property
    def mappability_patterns(self) -> list[str]:
        """List of mappability patterns for use in filtering extra contigs.

        Return an empty list if mappability is not given."""
        return fmap_maybe_def(
            [],
            lambda m: m.unplaced_chr_patterns,
            self.strat_inputs.mappability,
        )

    def to_build_data_unsafe(
        self,
        bk: BuildKey,
    ) -> "BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]":
        "Lookup a given build with a build key (and throw DesignError on fail)"
        bd = self.to_build_data(bk)
        if bd is None:
            raise DesignError(f"Could not create build data from key '{bk}'")
        return bd

    def to_build_data(
        self,
        bk: BuildKey,
    ) -> "BuildData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT] | None":
        "Lookup a given build with a build key"
        try:
            return BuildData_(self, bk, self.builds[bk])
        except KeyError:
            return None

    def get_refkeys_unsafe_(self, f: RefDataToSrc) -> list[RefKeyFullS]:
        """
        Get the list of refkeys (either one or two) given a function
        that retrieves an input file
        """
        return not_none_unsafe(
            f(self),
            lambda s: s.to_str_refkeys(self.refkey),
        )

    def get_refkeys_unsafe(self, f: StratInputToSrc) -> list[RefKeyFullS]:
        """Like 'get_refkeys_unsafe_' but the input function is restricted to
        the 'strat_inputs' member of this object.

        """
        return not_none_unsafe(
            f(self.strat_inputs),
            lambda s: s.to_str_refkeys(self.refkey),
        )

    @property
    def has_low_complexity_rmsk(self) -> bool:
        """Return True if this reference has repeat masker specified."""
        return self.strat_inputs.low_complexity.rmsk is not None

    @property
    def has_low_complexity_simreps(self) -> bool:
        """Return True if this reference has simple repeats specified."""
        return self.strat_inputs.low_complexity.simreps is not None

    @property
    def has_low_complexity_censat(self) -> bool:
        """Return True if this reference has satellites specified."""
        return self.strat_inputs.low_complexity.satellites is not None


HapRefData = RefData_[HapChrSrc[RefSrc], HapBedSrc, HapBedSrc, Haploid[BedSrc]]
Dip1RefData = RefData_[Dip1ChrSrc[RefSrc], DipBedSrc, Dip1BedSrc, Diploid[BedSrc]]
Dip2RefData = RefData_[Dip2ChrSrc[RefSrc], DipBedSrc, Dip2BedSrc, Diploid[BedSrc]]

AnyRefData = HapRefData | Dip1RefData | Dip2RefData


@dataclass(frozen=True)
class BuildData_(Generic[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]):
    """A helper class corresponding a given build.

    This follows a similar motivation as 'RefData_' above.
    """

    refdata: RefData_[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]
    buildkey: BuildKey
    build: Build[AnyBedT, AnyBedT_]

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
    def want_hets(self) -> bool:
        r = self.refdata.ref
        if isinstance(r, HapChrSrc):
            return False
        elif isinstance(r, Dip1ChrSrc) or isinstance(r, Dip2ChrSrc):
            return len(self.build.include.hets) > 0
        else:
            assert_never(r)

    @property
    def mappability_params(
        self,
    ) -> tuple[list[int], list[int], list[int]]:
        ms = self.build.include.mappability
        return unzip3([(m.length, m.mismatches, m.indels) for m in ms])


class Stratification(GenericModel, Generic[RefSourceT, AnyBedT, AnyBedT_, AnySrcT]):
    """Configuration for stratifications for a given reference."""

    ref: RefSourceT
    strat_inputs: StratInputs[AnyBedT, AnySrcT]
    builds: dict[BuildKey, Build[AnyBedT, AnyBedT_]]


HapBuildData = BuildData_[
    HapChrSrc[RefSrc],
    HapBedSrc,
    HapBedSrc,
    Haploid[BedSrc],
]


Dip1BuildData = BuildData_[
    Dip1ChrSrc[RefSrc],
    DipBedSrc,
    Dip1BedSrc,
    Diploid[BedSrc],
]

Dip2BuildData = BuildData_[
    Dip2ChrSrc[RefSrc],
    DipBedSrc,
    Dip2BedSrc,
    Diploid[BedSrc],
]


AnyBuildData = HapBuildData | Dip1BuildData | Dip2BuildData


HapStrat = Stratification[
    HapChrSrc[RefSrc],
    HapBedSrc,
    HapBedSrc,
    Haploid[BedSrc],
]
Dip1Strat = Stratification[
    Dip1ChrSrc[RefSrc],
    DipBedSrc,
    Dip1BedSrc,
    Diploid[BedSrc],
]
Dip2Strat = Stratification[
    Dip2ChrSrc[RefSrc],
    DipBedSrc,
    Dip2BedSrc,
    Diploid[BedSrc],
]


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
    comparison_strats: dict[CompareKey, HttpUrl] = {}
    haploid_stratifications: dict[RefKey, HapStrat] = {}
    diploid1_stratifications: dict[RefKey, Dip1Strat] = {}
    diploid2_stratifications: dict[RefKey, Dip2Strat] = {}
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
    malloc: Malloc = Malloc()

    # TODO validate comparison keys
    @validator(
        "haploid_stratifications",
        "diploid1_stratifications",
        "diploid2_stratifications",
        each_item=True,
    )
    def builds_have_valid_existing(
        cls,
        v: HapStrat,
        values: dict[str, Any],
    ) -> HapStrat:
        try:
            levels = cast(list[OtherLevelKey], values["other_levels"])
            bad = [
                f"level='{lk}'; build='{bk}'"
                for bk, b in v.builds.items()
                for lk in b.other_strats
                if lk not in levels
            ]
            if len(bad) > 0:
                assert (
                    False
                ), f"builds referencing invalid strat categories: {', '.join(bad)}"
        except KeyError:
            pass
        return v

    @validator(
        "haploid_stratifications",
        "diploid1_stratifications",
        "diploid2_stratifications",
        each_item=True,
    )
    def builds_have_valid_old_version(
        cls,
        v: HapStrat,
        values: dict[str, Any],
    ) -> HapStrat:
        try:
            prev: dict[CompareKey, HttpUrl] = values["comparison_strats"]
            bad = [
                f"version='{pk}'; build='{bk}'"
                for bk, b in v.builds.items()
                if b.comparison is not None
                if (pk := b.comparison.other) not in prev
            ]
            assert (
                len(bad) == 0
            ), f"builds referencing invalid previous version keys: {', '.join(bad)}"
        except KeyError:
            pass
        return v

    @validator("diploid2_stratifications")
    def no_overlapping_refkeys(
        cls,
        v: dict[RefKey, Dip2Strat],
        values: dict[str, Any],
    ) -> dict[RefKey, Dip2Strat]:
        try:
            hap: list[RefKey] = list(values["haploid_stratifications"])
            dip1: list[RefKey] = list(values["diploid2_stratifications"])
            dip2 = list(v)
            ds = list(duplicates_everseen(hap + dip1 + dip2))
            assert len(ds) == 0, f"duplicate refkeys: {', '.join(ds)}"
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

    def buildkey_to_malloc(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        f: Callable[[Malloc], int],
    ) -> int:
        bd = self.to_build_data(strip_full_refkey(rk), bk)
        return max(
            fmap_maybe_def(f(self.malloc), lambda m: f(m), bd.build.malloc), 1000
        )

    def buildkey_to_ref_mappers(
        self, rk: RefKeyFullS, bk: BuildKey
    ) -> tuple[bed.InitMapper, bed.FinalMapper]:
        """Lookup a given build and return the init/final mappers
        corresponding to the reference chromosome names.

        This is useful for cases where the reference itself is used to
        generate a bed-like file which then needs to be sorted.
        """
        m = self.with_build_data_full(
            rk,
            bk,
            lambda bd: bd.refdata.ref.noop_conversion(bd.chr_indices),
            lambda bd: bd.refdata.ref.noop_conversion(bd.chr_indices),
            lambda hap, bd: hap.from_either(
                *bd.refdata.ref.noop_conversion(bd.chr_indices)
            ),
        )
        return (m.init_mapper, m.final_mapper)

    def buildkey_to_other_keys(
        self, rk: RefKeyFullS, bk: BuildKey
    ) -> list[tuple[OtherLevelKey, OtherStratKey]]:
        """Lookup a given build and return a list of keys corresponding to the
        external bed files we wish to include in the final package.
        """
        bd = self.to_build_data(strip_full_refkey(rk), bk)
        return [(lk, sk) for lk, s in bd.build.other_strats.items() for sk in s]

    def refsrckey_to_ref_src(self, rsk: RefKeyFullS) -> RefSrc:
        """Lookup a given reference and return its source object (haplotype
        specific)."""
        rk, hap = parse_full_refkey(rsk)
        src = self.to_ref_data(rk).ref.src
        return from_hap_or_dip(src, hap)

    def refkey_to_bed_refsrckeys(
        self, f: StratInputToBed, rk: RefKey
    ) -> list[RefKeyFullS]:
        """Lookup a given reference and return the full refkeys for the
        bed file obtained with the given function.

        This is useful for bed file normalization rules which are all in terms
        of the bare refkey (ie not haplotype specific) but need to somehow
        get a list of inputs which are downloaded. Since there might be one or
        two inputs which may or may not have a haplotype associated with them,
        this function provides the full refkeys for obtaining said inputs.
        """
        # TODO this seems like a useful glue function (the labmda that is)
        return self.to_ref_data(rk).get_refkeys_unsafe(
            lambda si: fmap_maybe(lambda x: x.data.src, f(si))
        )

    def _refkey_to_src(self, f: RefDataToSrc, rk: RefKeyFullS) -> BedSrc:
        rk_, hap = parse_full_refkey(rk)
        src = with_ref_data(
            self.to_ref_data(rk_), lambda rd: f(rd), lambda rd: f(rd), lambda rd: f(rd)
        )
        # TODO mypy doens't like me using my 'maybe' functional functions
        if src is None:
            raise DesignError()
        return from_hap_or_dip(src, hap)

    def refsrckey_to_bed_src(self, f: StratInputToBed, rk: RefKeyFullS) -> BedSrc:
        """Lookup a haplotype-specific bed file source with the given function."""
        return self._refkey_to_src(
            lambda rd: fmap_maybe(lambda x: x.data.src, f(rd.strat_inputs)), rk
        )

    def _refsrckey_to_xy_feature_src(self, rsk: RefKeyFullS, i: ChrIndex) -> BedSrc:
        return (
            self.to_ref_data(strip_full_refkey(rsk))
            .strat_inputs.xy_feature_bed_unsafe(i)
            .data.src.hap
        )

    def refsrckey_to_x_features_src(self, rsk: RefKeyFullS) -> BedSrc:
        """Return the X features source file for a given reference."""
        return self._refsrckey_to_xy_feature_src(rsk, ChrIndex.CHRX)

    def refsrckey_to_y_features_src(self, rsk: RefKeyFullS) -> BedSrc:
        """Return the Y features source file for a given reference."""
        return self._refsrckey_to_xy_feature_src(rsk, ChrIndex.CHRY)

    def refkey_to_xy_ref_chr_pattern(
        self, rk: RefKeyFullS, i: ChrIndex
    ) -> HapChrPattern:
        """Return the XY chr pattern for a given reference and haplotype."""
        return self.with_ref_data_full(
            rk,
            lambda rd: rd.ref.chr_pattern,
            lambda rd: rd.ref.chr_pattern.to_hap_pattern(i.xy_to_hap_unsafe),
            lambda hap, rd: rd.ref.chr_pattern.from_either(hap),
        )

    def buildkey_to_bed_refsrckeys(
        self, f: BuildDataToBed, rk: RefKey, bk: BuildKey
    ) -> list[RefKeyFullS]:
        """Like 'refkey_to_bed_refsrckeys' but build-specific.

        Used for looking up benchmark files for each build.
        """
        # TODO this "update" function is not DRY
        # return self.refkey_to_bed_refsrckeys(lambda rd: f(rd.to_build_data(bk)), rk)
        return self.to_ref_data(rk).get_refkeys_unsafe_(
            lambda rd: fmap_maybe(lambda x: x.data.src, f(rd.to_build_data(bk)))
        )

    def buildkey_to_bed_src(
        self, f: BuildDataToBed, rk: RefKeyFullS, bk: BuildKey
    ) -> BedSrc:
        """Like 'refsrckey_to_bed_src' but build-specific.

        Used for looking up benchmark sources for each build.
        """
        return self._refkey_to_src(
            lambda rd: fmap_maybe(lambda x: x.data.src, f(rd.to_build_data(bk))),
            rk,
        )

    def buildkey_to_vcf_src(
        self, f: BuildDataToVCF, rk: RefKeyFullS, bk: BuildKey
    ) -> BedSrc:
        """Like 'buildkey_to_bed_src' but for benchmark VCF sources."""
        # TODO not DRY
        rk_, hap = parse_full_refkey(rk)
        bd = self.to_build_data(rk_, bk)
        src = with_build_data(bd, lambda bd: f(bd), lambda bd: f(bd), lambda bd: f(bd))
        if src is None:
            raise DesignError()
        return from_hap_or_dip(src.data.src, hap)

    def refkey_to_functional_refsrckeys(
        self, f: StratInputToSrc, rk: RefKey
    ) -> list[RefKeyFullS]:
        """Like 'refkey_to_bed_refsrckeys' but for source files in the
        "Functional" stratification level."""
        return self.to_ref_data(rk).get_refkeys_unsafe(f)

    def refsrckey_to_functional_src(
        self, f: StratInputToSrc, rk: RefKeyFullS
    ) -> BedSrc:
        """Like 'refsrckey_to_bed_src' but for source files in the
        "Functional" stratification level."""
        return self._refkey_to_src(lambda rd: f(rd.strat_inputs), rk)

    def to_ref_data(self, rk: RefKey) -> AnyRefData:
        """Lookup refdata object for a given refkey."""
        if rk in self.haploid_stratifications:
            return to_ref_data_unsafe(self.haploid_stratifications, rk)
        elif rk in self.diploid1_stratifications:
            return to_ref_data_unsafe(self.diploid1_stratifications, rk)
        elif rk in self.diploid2_stratifications:
            return to_ref_data_unsafe(self.diploid2_stratifications, rk)
        else:
            raise DesignError(f"invalid ref key: '{rk}'")

    def to_build_data(self, rk: RefKey, bk: BuildKey) -> AnyBuildData:
        """Lookup builddata object for a given refkey and build key."""

        def hap(rd: HapRefData) -> AnyBuildData:
            return rd.to_build_data_unsafe(bk)

        def dip1(rd: Dip1RefData) -> AnyBuildData:
            return rd.to_build_data_unsafe(bk)

        def dip2(rd: Dip2RefData) -> AnyBuildData:
            return rd.to_build_data_unsafe(bk)

        return with_ref_data(self.to_ref_data(rk), hap, dip1, dip2)

    def with_ref_data(
        self,
        rk: RefKey,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Dip1RefData], X],
        dip2_f: Callable[[Dip2RefData], X],
    ) -> X:
        """Lookup refdata and apply function depending on if it is hap, dip1,
        or dip2.
        """
        return with_ref_data(self.to_ref_data(rk), hap_f, dip1_f, dip2_f)

    def with_ref_data_full(
        self,
        rk: RefKeyFullS,
        hap_f: Callable[[HapRefData], X],
        dip1_f: Callable[[Dip1RefData], X],
        dip2_f: Callable[[Haplotype, Dip2RefData], X],
    ) -> X:
        """Like 'with_ref_data_full' but takes a full refkey and supplies the
        haplotype in the dip2 case.
        """
        rk_, hap = parse_full_refkey(rk)
        return self.with_ref_data(
            rk_,
            lambda rd: none_unsafe(hap, hap_f(rd)),
            lambda rd: none_unsafe(hap, dip1_f(rd)),
            lambda rd: not_none_unsafe(hap, lambda hap: dip2_f(hap, rd)),
        )

    def with_build_data(
        self,
        rk: RefKey,
        bk: BuildKey,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data' but for build data"""
        return with_build_data(self.to_build_data(rk, bk), hap_f, dip1_f, dip2_f)

    def with_build_data_full(
        self,
        rfk: RefKeyFullS,
        bk: BuildKey,
        hap_f: Callable[[HapBuildData], X],
        dip1_f: Callable[[Dip1BuildData], X],
        dip2_f: Callable[[Haplotype, Dip2BuildData], X],
    ) -> X:
        """Like 'with_ref_data_full' but for build data"""
        rk_, hap = parse_full_refkey(rfk)
        return self.with_build_data(
            rk_,
            bk,
            lambda bd: none_unsafe(hap, hap_f(bd)),
            lambda bd: none_unsafe(hap, dip1_f(bd)),
            lambda bd: not_none_unsafe(hap, lambda hap: dip2_f(hap, bd)),
        )

    def with_ref_data_and_bed(
        self,
        rk: RefKey,
        get_bed_f: RefDataToBed,
        hap_f: Callable[[HapRefData, HapBedFile], X],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFile], X],
        dip_1to2_f: Callable[[Dip2RefData, Dip1BedFile], X],
        dip_2to1_f: Callable[[Dip1RefData, Dip2BedFile], X],
        dip_2to2_f: Callable[[Dip2RefData, Dip2BedFile], X],
    ) -> X:
        """Lookup refdata and bed file according to the given function.
        Then apply function depending on if the bed is haploid or diploid and
        if the reference is hap/dip1/dip2.
        """
        return self.with_ref_data(
            rk,
            lambda rd: not_none_unsafe(get_bed_f(rd), lambda bd: hap_f(rd, bd)),
            lambda rd: with_dip_bedfile(
                not_none_unsafe(get_bed_f(rd), noop),
                lambda bf: dip_1to1_f(rd, bf),
                lambda bf: dip_2to1_f(rd, bf),
            ),
            lambda rd: with_dip_bedfile(
                not_none_unsafe(get_bed_f(rd), noop),
                lambda bf: dip_1to2_f(rd, bf),
                lambda bf: dip_2to2_f(rd, bf),
            ),
        )

    def with_ref_data_and_bed_full(
        self,
        rk: RefKeyFullS,
        get_bed_f: RefDataToBed,
        hap_f: Callable[[HapRefData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1RefData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2RefData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1RefData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2RefData, Dip2BedFile], Z],
    ) -> Z:
        """Like 'with_ref_data_and_bed' but also takes a full refkey and
        supplies the haplotype to the dip2-reference case."""
        rk_, hap = parse_full_refkey(rk)
        return self.with_ref_data_and_bed(
            rk_,
            get_bed_f,
            hap_f,
            dip_1to1_f,
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_1to2_f(h, rd, bf)),
            dip_2to1_f,
            lambda rd, bf: not_none_unsafe(hap, lambda h: dip_2to2_f(h, rd, bf)),
        )

    def with_build_data_and_bed(
        self,
        rk: RefKey,
        bk: BuildKey,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        """Like 'with_ref_data_and_bed' but for build data."""
        return self.with_ref_data_and_bed(
            rk,
            lambda rd: get_bed_f(rd.to_build_data_unsafe(bk)),
            lambda rd, bf: hap_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_1to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_1to2_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_2to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_2to2_f(rd.to_build_data_unsafe(bk), bf),
        )

    def with_build_data_and_bed_full(
        self,
        rk: RefKeyFullS,
        bk: BuildKey,
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[Haplotype, Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[Haplotype, Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        """Like 'with_ref_data_and_bed_full' but for build data."""
        return self.with_ref_data_and_bed_full(
            rk,
            lambda rd: get_bed_f(rd.to_build_data_unsafe(bk)),
            lambda rd, bf: hap_f(rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_1to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda hap, rd, bf: dip_1to2_f(hap, rd.to_build_data_unsafe(bk), bf),
            lambda rd, bf: dip_2to1_f(rd.to_build_data_unsafe(bk), bf),
            lambda hap, rd, bf: dip_2to2_f(hap, rd.to_build_data_unsafe(bk), bf),
        )

    def _with_build_data_and_bed_i(
        self,
        rk: RefKey,
        bk: BuildKey,
        inputs: list[X],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[tuple[X, X], Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[tuple[X, X], Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        """Like 'with_build_data_and_bed' but also take a list of input files
        and supply it to one of the supplied higher-order functions.

        Throw DesignError if the list of inputs does not correspond to the
        expected number of inputs (either 1 or 2).

        """
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

    def _with_build_data_and_bed_io(
        self,
        rk: RefKey,
        bk: BuildKey,
        inputs: list[X],
        output_f: Callable[[RefKeyFull], Y],
        write_outputs: Callable[[list[Y]], None],
        get_bed_f: BuildDataToBed,
        hap_f: Callable[[X, Y, HapBuildData, HapBedFile], Z],
        dip_1to1_f: Callable[[X, Y, Dip1BuildData, Dip1BedFile], Z],
        dip_1to2_f: Callable[[X, tuple[Y, Y], Dip2BuildData, Dip1BedFile], Z],
        dip_2to1_f: Callable[[tuple[X, X], Y, Dip1BuildData, Dip2BedFile], Z],
        dip_2to2_f: Callable[[tuple[X, X], tuple[Y, Y], Dip2BuildData, Dip2BedFile], Z],
    ) -> Z:
        """Like '_with_build_data_and_bed_i' but also take a function that
        generates an output file path and function that writes the outputs paths
        (NOT the data itself) to disk. The five higher order functions
        corresponding to each of the haplotype configurations must take the
        correct number of output paths which they are assumed to use for
        writing.
        """

        def out1(src: HapRefSrc) -> Y:
            return with_first(output_f(src.key(rk)), lambda o: write_outputs([o]))

        def out2(src: DipRefSrc) -> tuple[Y, Y]:
            return with_first(
                both(output_f, src.keys(rk)), lambda o: write_outputs([*o])
            )

        return self._with_build_data_and_bed_i(
            rk,
            bk,
            inputs,
            get_bed_f,
            lambda i, bd, bf: hap_f(i, out1(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_1to1_f(i, out1(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_1to2_f(i, out2(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_2to1_f(i, out1(bd.refdata.ref.src), bd, bf),
            lambda i, bd, bf: dip_2to2_f(i, out2(bd.refdata.ref.src), bd, bf),
        )

    def with_build_data_and_bed_io(
        self,
        rk: RefKey,
        bk: BuildKey,
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
        """Like 'with_build_data_and_bed_io' with the following differences.

        * This function takes an input pattern to be used to generate the list
        of output files. This pattern must have a '%s' for where the full ref
        key will be subbed, and must not contain any unexpanded snakemake
        wildcards.

        * This functions takes a single path to which the output list will be
        written instead of a function. The output list will be written as a json
        dump.

        * here, the function corresponding to dip2->dip2 here only takes a
        single haplotype, which is useful since the haplotypes can be processed
        independently and thus the supplied function can be made less
        redundant.

        """

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

        self._with_build_data_and_bed_io(
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

    # final refkey/buildkey lists (for the "all" target and related)

    @property
    def all_build_keys(self) -> tuple[list[RefKey], list[BuildKey]]:
        return unzip2(
            all_build_keys(self.haploid_stratifications)
            + all_build_keys(self.diploid1_stratifications)
            + all_build_keys(self.diploid2_stratifications)
        )

    @property
    def all_full_build_keys(self) -> tuple[list[RefKeyFullS], list[BuildKey]]:
        return unzip2(
            all_ref_build_keys(self.haploid_stratifications)
            + all_ref_build_keys(self.diploid1_stratifications)
            + all_ref_build_keys(self.diploid2_stratifications)
        )

    # source refkey/buildkey lists (for the "all resources" rule)

    @property
    def all_ref_refsrckeys(self) -> list[RefKeyFullS]:
        return (
            all_ref_refsrckeys(self.haploid_stratifications)
            + all_ref_refsrckeys(self.diploid1_stratifications)
            + all_ref_refsrckeys(self.diploid2_stratifications)
        )

    def _all_bed_build_and_refsrckeys(
        self, f: BuildDataToSrc
    ) -> list[tuple[RefKeyFullS, BuildKey]]:
        return (
            all_bed_build_and_refsrckeys(self.haploid_stratifications, f)
            + all_bed_build_and_refsrckeys(self.diploid1_stratifications, f)
            + all_bed_build_and_refsrckeys(self.diploid2_stratifications, f)
        )

    def _all_bed_refsrckeys(self, f: BuildDataToSrc) -> list[RefKeyFullS]:
        return [rk for rk, _ in self._all_bed_build_and_refsrckeys(f)]

    @property
    def all_refkey_gap(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_gaps(bd))
        )

    @property
    def all_refkey_rmsk(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_rmsk(bd))
        )

    @property
    def all_refkey_simreps(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_simreps(bd))
        )

    @property
    def all_refkey_censat(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_satellites(bd))
        )

    @property
    def all_refkey_segdups(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_superdups(bd))
        )

    @property
    def all_buildkey_bench(self) -> list[tuple[RefKeyFullS, BuildKey]]:
        return self._all_bed_build_and_refsrckeys(
            lambda bd: fmap_maybe(lambda x: x.data.src, bd_to_bench_vcf(bd))
        )

    @property
    def all_refkey_ftbl(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(bd_to_ftbl)

    @property
    def all_refkey_gff(self) -> list[RefKeyFullS]:
        return self._all_bed_refsrckeys(bd_to_gff)

    # other nice functions

    def refkey_is_dip1(self, rk: RefKeyFullS) -> bool:
        """Test if refkey is dip1 or dip2.

        Return True if dip1, false if dip2, and error otherwise.

        Incoming key must have a haplotype appended to it. Otherwise error.

        This function is useful for cases where dip1 rules need to be "split"
        into component haplotypes. Normally each rule for a dip1 configuration
        is denoted by only the refkey (since both haplotypes are in each file).
        However, if split in the manner above we need to disambiguate by
        appending the haplotype, in which case the key become indistinguishable
        from a dip2 key unless we run the logic of checking the het1/dip1/dip2
        refkey maps (which is what this function does).
        """
        rk_, hap = parse_full_refkey(rk)
        if hap is None:
            raise DesignError(f"Refkey must have a haplotype {rk_}")
        if rk_ in self.haploid_stratifications:
            raise DesignError(f"Refkey must be either dip1 or dip2: {rk_}")
        elif rk_ in self.diploid1_stratifications:
            return True
        elif rk_ in self.diploid2_stratifications:
            return False
        else:
            raise DesignError(f"Unknown refkey: {rk_}")

    def dip1_either(self, left: X, right: X, rk: RefKeyFullS) -> X:
        """Return left if dip1, right if dip2, and error otherwise."""
        return left if self.refkey_is_dip1(rk) else right

    def thread_per_chromosome(self, rk: RefKeyFullS, bk: BuildKey, n: int) -> int:
        cis = self.to_build_data(strip_full_refkey(rk), bk).chr_indices
        return min(n, len(cis))


################################################################################
# protocols

# hacky rankN type mimicry


class RefDataToBed(Protocol):
    A = TypeVar("A", HapChrSrc[BedSrc], Dip1ChrSrc[BedSrc] | Dip2ChrSrc[BedSrc])

    def __call__(
        self,
        __x: RefData_[RefSourceT, A, AnyBedT_, AnySrcT],
    ) -> BedFile[A] | None:
        pass


class RefDataToSrc(Protocol):
    A = TypeVar("A", Haploid[BedSrc], Diploid[BedSrc])

    def __call__(
        self,
        __x: RefData_[RefSourceT, AnyBedT, AnyBedT_, A],
    ) -> A | None:
        pass


class StratInputToBed(Protocol):
    A = TypeVar("A", HapChrSrc[BedSrc], Dip1ChrSrc[BedSrc] | Dip2ChrSrc[BedSrc])

    def __call__(self, __x: StratInputs[A, AnySrcT]) -> BedFile[A] | None:
        pass


class StratInputToSrc(Protocol):
    A = TypeVar("A", Haploid[BedSrc], Diploid[BedSrc])

    def __call__(self, __x: StratInputs[AnyBedT, A]) -> A | None:
        pass


class BuildDataToBed(Protocol):
    A = TypeVar("A", HapBedSrc, DipBedSrc)

    def __call__(
        self,
        __x: BuildData_[RefSourceT, A, AnyBedT_, AnySrcT],
    ) -> BedFile[A] | None:
        pass


class BuildDataToVCF(Protocol):
    A = TypeVar("A", HapBedSrc, Dip1BedSrc, Dip2BedSrc)

    def __call__(
        self,
        __x: BuildData_[RefSourceT, AnyBedT, A, AnySrcT],
    ) -> BedFile[A] | None:
        pass


class BuildDataToSrc(Protocol):
    A = TypeVar("A", Haploid[BedSrc], Diploid[BedSrc])

    def __call__(
        self,
        __x: BuildData_[RefSourceT, AnyBedT, AnyBedT_, A],
    ) -> A | None:
        pass
