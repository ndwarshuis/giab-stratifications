from typing import Any
import common.config as cfg
from common.bed import read_filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    def get_bed(x: cfg.Stratification) -> cfg.BedFile:
        b = x.segdups.superdups
        assert b is not None
        return b

    read_filter_sort_bed(
        sconf,
        smk.input[0],
        smk.output[0],
        get_bed,
        rk,
        bk,
        [],
    )


main(snakemake, snakemake.config)  # type: ignore
