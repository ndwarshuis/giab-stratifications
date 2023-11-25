from typing import Any
import common.config as cfg
from common.bed import read_filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    match smk.input:
        case [hap]:
            p = cfg.HaploidRefPair(
                smk.wildcards["ref_key"],
                smk.wildcards["build_key"],
            )
            rm = sconf.refkey_to_strat(p).low_complexity.rmsk
            assert rm is not None, "this should not happen"
            read_filter_sort_bed(
                sconf, hap, smk.output[0], rm.params, p, [rm.class_col]
            )
        case [dip0, dip1]:
            pass


main(snakemake, snakemake.config)  # type: ignore
