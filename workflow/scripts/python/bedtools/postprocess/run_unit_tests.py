from typing import Any, NamedTuple
from pathlib import Path
import pandas as pd
from pybedtools import BedTool as bt  # type: ignore
from os.path import dirname, basename
import common.config as cfg
from common.io import setup_logging, is_bgzip

log = setup_logging(snakemake.log[0])  # type: ignore


class GaplessBT(NamedTuple):
    auto: bt
    parY: bt


def test_bgzip(strat_file: Path) -> list[str]:
    return [] if is_bgzip(strat_file) else ["is not bgzip file"]


def test_bed(
    strat_file: Path,
    reverse_map: dict[str, int],
    gapless: GaplessBT,
) -> list[str]:
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

    # choose gap file depending on if this is an XY strat or not
    isXY = basename(dirname(strat_file)) == "XY"
    gapless_bt = gapless.parY if isXY else gapless.auto
    gaps = bt.from_dataframe(df).subtract(gapless_bt).to_dataframe()

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
            (len(gaps) == 0, "bed contains gap regions"),
        ]
        if result is False
    ]


def run_all_tests(
    strat_file: Path,
    reverse_map: dict[str, int],
    gapless: GaplessBT,
) -> list[tuple[Path, str]]:
    return [
        (strat_file, msg)
        for msg in test_bgzip(strat_file) + test_bed(strat_file, reverse_map, gapless)
    ]


# TODO not DRY...but also this is literally 3 lines so whatever
def strat_files(path: str) -> list[Path]:
    with open(path, "r") as f:
        return [Path(s.strip()) for s in f]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    reverse_map = {v: k for k, v in sconf.buildkey_to_final_chr_mapping(rk, bk).items()}

    gapless = GaplessBT(
        auto=bt(str(smk.input["gapless_auto"])),
        parY=bt(str(smk.input["gapless_parY"])),
    )

    all_failures = [
        res
        for p in strat_files(str(smk.input["strats"]))
        for res in run_all_tests(p, reverse_map, gapless)
    ]

    for res in all_failures:
        log.error("%s: %s" % res)

    if len(all_failures) > 0:
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
