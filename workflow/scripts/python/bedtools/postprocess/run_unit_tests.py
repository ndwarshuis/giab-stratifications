import pandas as pd
from typing import Any, NamedTuple
from pathlib import Path
from pybedtools import BedTool as bt  # type: ignore
from os.path import dirname, basename
from os import scandir
import common.config as cfg
from common.io import setup_logging, is_bgzip
from common.bed import InternalChrIndex
import subprocess as sp

log = setup_logging(snakemake.log[0])  # type: ignore

RevMapper = dict[str, InternalChrIndex]


class GaplessBT(NamedTuple):
    auto: bt
    parY: bt


def test_bgzip(strat_file: Path) -> list[str]:
    """Each bed file should be bgzip'ed (not just gzip'ed)."""
    return [] if is_bgzip(strat_file) else ["is not bgzip file"]


def test_bed(
    strat_file: Path,
    reverse_map: RevMapper,
    gapless: GaplessBT,
) -> list[str]:
    """Each stratification should be a valid bed file."""

    # test bed files in multiple phases, where each phase depends on the
    # previous being valid
    #
    # Phase 1: try to import the dataframe (with almost no assumptions except
    # that it has no header)
    try:
        # no low mem since we don't want to assume what dtypes each column has,
        # and by default pandas will chunk the input which could lead to mixed
        # types
        df = pd.read_table(strat_file, header=0, low_memory=False)
    except pd.errors.EmptyDataError:
        # if df is empty, nothing to do
        return []
    except Exception as e:
        # catch any other weird read errors
        return [str(e)]

    # Phase 2: ensure "bed file" has 3 columns
    if len(df.columns) != 3:
        return ["bed file has wrong number of columns"]

    # Phase 3: ensure "bed file" has column of the correct type
    cols = {"chrom": str, "start": int, "end": int}
    df.columns = pd.Index([*cols])

    try:
        df = df.astype(cols)
    except ValueError as e:
        return [str(e)]

    # Phase 4: we now know we have a real bed file, ensure that we didn't
    # invent any new chromosomes (try to convert chrom column to ints, which
    # will produce NaNs on error)
    chr_ints = df["chrom"].map(reverse_map).astype("Int64")
    invalid_chrs = df["chrom"][chr_ints.isnull()].unique().tolist()
    if len(invalid_chrs) > 0:
        return [f"invalid chrs: {invalid_chrs}"]

    # Phase 5: test remaining assumptions:
    # - chromosome column is sorted
    # - each region is disjoint and ascending (eg the start of a region is at
    #   least 1bp away from the end of the previous)
    # - does not intersect with gap regions (if any) or extend beyond chr bounds
    same_chrom = df["chrom"] == df["chrom"].shift(1)
    distance = (df["start"] - df["end"].shift(1))[same_chrom]
    overlapping_lines = distance[distance < 1].index.tolist()

    # choose gap file depending on if this is an XY strat or not; note that the
    # built-in gaps file (if applicable) is obviously whitelisted from this
    # check
    if "gaps_slop15kb" in basename(strat_file):
        no_gaps = True
    else:
        isXY = basename(dirname(strat_file)) == "XY"
        gapless_bt = gapless.parY if isXY else gapless.auto
        gaps = bt.from_dataframe(df).subtract(gapless_bt).to_dataframe()
        no_gaps = len(gaps) == 0

    return [
        msg
        for result, msg in [
            (
                pd.Index(chr_ints).is_monotonic_increasing,
                "chromosomes not sorted",
            ),
            (
                len(overlapping_lines) == 0,
                f"non-disjoint regions at lines: {overlapping_lines}",
            ),
            (no_gaps, "bed contains gap regions"),
        ]
        if result is False
    ]


def test_checksums(checksums: Path) -> list[str]:
    """The md5sum utility should pass when invoked."""
    res = sp.run(
        ["md5sum", "-c", "--strict", "--quiet", checksums.name],
        cwd=checksums.parent,
        capture_output=True,
        text=True,
    )
    errors = [] if (s := res.stdout) == "" else s.strip().split("\n")
    return [f"checksum error: {e}" for e in errors]


def test_tsv_list(tsv_list: Path) -> list[str]:
    """The final stratifications list should match final beds exactly.

    We are running our hotel on very tight margins; no extra or missing beds
    allowed.
    """

    def strat_set(root: Path, sub: Path) -> set[str]:
        if root.is_dir():
            deeper = (Path(p.path) for p in scandir(root))
            return {str(q) for d in deeper for q in strat_set(d, sub / d.name)}
        elif root.is_file() and root.name.endswith(".bed.gz"):
            return {str(sub)}
        else:
            return set()

    current = strat_set(tsv_list.parent, Path("./"))
    with open(tsv_list, "r") as f:
        listed = {line.strip().split("\t")[1] for line in f}
        missing = [f"not in final directory: {p}" for p in listed - current]
        extra = [f"not in final list: {p}" for p in current - listed]
        return missing + extra


def run_all_tests(
    strat_file: Path,
    reverse_map: RevMapper,
    gapless: GaplessBT,
) -> list[tuple[Path, str]]:
    return [
        (strat_file, msg)
        for msg in test_bgzip(strat_file) + test_bed(strat_file, reverse_map, gapless)
    ]


def strat_files(path: str) -> list[Path]:
    with open(path, "r") as f:
        return [Path(s.strip()) for s in f]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    fm = sconf.with_build_data_final(
        ws["ref_key"],
        ws["build_key"],
        lambda bd: bd.ref_chr_conversion.final_mapper,
        lambda bd: bd.ref_chr_conversion.final_mapper,
        lambda hap, bd: hap.from_either(*bd.ref_chr_conversion).final_mapper,
    )
    reverse_map = {v: k for k, v in fm.items()}

    # check global stuff first (since this is faster)

    global_failures: list[str] = test_checksums(
        Path(smk.input["checksums"])
    ) + test_tsv_list(Path(smk.input["strat_list"]))

    for j in global_failures:
        log.error(j)

    if len(global_failures) > 0:
        exit(1)

    # check individual stratification files

    gapless = GaplessBT(
        auto=bt(str(smk.input["gapless_auto"])),
        parY=bt(str(smk.input["gapless_parY"])),
    )

    strat_failures = [
        res
        for p in strat_files(str(smk.input["strats"]))
        for res in run_all_tests(p, reverse_map, gapless)
    ]

    for i in strat_failures:
        log.error("%s: %s" % i)

    if len(strat_failures) > 0:
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
