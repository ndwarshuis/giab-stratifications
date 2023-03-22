from common.config import is_bgzip
from typing import Any
from os import scandir
from io import TextIOWrapper
from pathlib import Path
import pandas as pd
import common.config as cfg


def log_fail(strat_file: Path, log: TextIOWrapper, msg: str) -> None:
    log.write("ERROR - %s: %s" % (strat_file, msg))


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    strat_files = [
        Path(p.path)
        for i in smk.input
        for p in scandir(i)
        if p.is_file() and p.path.endswith(".gz")
    ]
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    # TODO this is a mess and can be cleaned up :(
    def test_bed_format(strat_file: Path, log: TextIOWrapper) -> bool:
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
            log_fail(strat_file, log, "error when assessing bed format")
            return False

        # check that chromosomes are valid
        reverse_map = {
            v: k for k, v in sconf.buildkey_to_final_chr_mapping(rk, bk).items()
        }

        chr_ints = df["chrom"].map(reverse_map).astype("Int64")
        invalid_chrs = df["chrom"][chr_ints.isnull()].unique().tolist()
        df["chrom"] = chr_ints

        if len(invalid_chrs > 0):
            log_fail(strat_file, log, f"invalid chrs: {invalid_chrs}")
            return False

        # test chromosomes sorted (start/end tested later)
        if not pd.Index(df["chrom"]).is_monotonic_increasing:
            log_fail(strat_file, log, "chromosomes not sorted")
            return False

        # test that start - prev_end >= 1 (which implies sorting)
        same_chrom = (df["chrom"] - df["chrom"].shift(-1)) == 0
        gaps = (df["end"] - df["start"].shift(-1))[same_chrom]
        n_overlapping = len(gaps[gaps < 1])

        if n_overlapping > 0:
            log_fail(strat_file, log, f"{n_overlapping} regions")
            return False

        return True

    def run_tests(strat_file: Path, log: TextIOWrapper) -> bool:
        if not is_bgzip(strat_file):
            log_fail(strat_file, log, "is not bgzip file")
            return False

        return test_bed_format(strat_file, log)

    with open(smk.log[0], "w") as f:
        assert all(run_tests(p, f) for p in strat_files) is True


main(snakemake, snakemake.config)  # type: ignore
