from common.config import is_bgzip
from typing import Any
from os import scandir
from pathlib import Path
import pandas as pd
import common.config as cfg


def test_bgzip(strat_file: Path) -> list[str]:
    return [] if is_bgzip(strat_file) else ["is not bgzip file"]


def test_bed_format(strat_file: Path, reverse_map: dict[str, int]) -> list[str]:
    # TODO this could be more informative
    try:
        df = pd.read_table(
            strat_file,
            names=["chrom", "start", "end"],
            dtype={0: str, 1: int, 2: int},
            header=0,
        )
    except:
        # if we can't read it, there's something wrong with it (likely an
        # incorrect number of columns)
        return ["error when assessing bed format"]

    # if we can read the file without breaking physics, test the following:
    # - no invalid chromosomes (test by converting known chromosomes to
    #   integers, which will produce NaNs for non-matches)
    # - chromosome column is sorted
    # - each region is disjoint and ascending (eg the start of a region is at
    #   least 1bp away from the end of the previous)

    chr_ints = df["chrom"].map(reverse_map).astype("Int64")
    invalid_chrs = df["chrom"][chr_ints.isnull()].unique().tolist()
    df["chrom"] = chr_ints

    same_chrom = (df["chrom"] - df["chrom"].shift(1)) == 0
    gaps = (df["start"] - df["end"].shift(1))[same_chrom]
    gap_lines = gaps[gaps < 1].index.tolist()

    return [
        msg
        for result, msg in [
            (len(invalid_chrs) == 0, f"invalid chrs: {invalid_chrs}"),
            (
                pd.Index(df["chrom"]).is_monotonic_increasing,
                "chromosomes not sorted",
            ),
            (len(gap_lines) == 0, f"disjoint regions at lines: {gap_lines}"),
        ]
        if result is False
    ]


def run_all_tests(
    strat_file: Path, reverse_map: dict[str, int]
) -> list[tuple[Path, str]]:
    return [
        (strat_file, msg)
        for msg in test_bgzip(strat_file) + test_bed_format(strat_file, reverse_map)
    ]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    reverse_map = {v: k for k, v in sconf.buildkey_to_final_chr_mapping(rk, bk).items()}

    strat_files = [
        Path(p.path)
        for i in smk.params["strat_dirs"]
        for p in scandir(i)
        if p.is_file() and p.path.endswith(".gz")
    ]

    all_failures = [res for p in strat_files for res in run_all_tests(p, reverse_map)]

    with open(smk.log[0], "w") as f:
        for res in all_failures:
            f.write("ERROR - %s: %s\n" % res)

    assert len(all_failures) == 0, "failed at least one unit test"


main(snakemake, snakemake.config)  # type: ignore
