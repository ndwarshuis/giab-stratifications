from typing import Any
import common.config as cfg
from common.bed import read_filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    ss = sconf.refkey_to_strat(rk).low_complexity.simreps
    assert ss is not None, "this should not happen"
    read_filter_sort_bed(sconf, smk.input[0], smk.output[0], ss.params, rk, bk)


main(snakemake, snakemake.config)  # type: ignore
