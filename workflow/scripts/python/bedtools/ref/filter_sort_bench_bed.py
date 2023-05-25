from typing import Any
import common.config as cfg
from common.bed import read_filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    bench = sconf.stratifications[rk].builds[bk].bench
    assert bench is not None, "this should not happen"
    ps = bench.bench_bed.params
    read_filter_sort_bed(sconf, smk.input[0], smk.output[0], ps, rk, bk)


main(snakemake, snakemake.config)  # type: ignore
