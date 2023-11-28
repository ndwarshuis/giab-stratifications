from pathlib import Path
from typing import Any
import common.config as cfg

# from common.bed import read_filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ons: list[Path] = smk.output[0]
    ins: list[Path] = smk.input
    bd = sconf.to_build_data(ws["ref_key"], ws["build_key"])
    hap = cfg.to_haplotype(ws["hap"])
    match (ins, ons):
        case ([i], [o]):
            # one haplotype for both bed and ref (no combine)
            if isinstance(bd, cfg.HaploidBuildData) and hap is None:
                bd.read_filter_sort_hap_bed(lambda si: si.low_complexity.simreps, i, o)
            # one bed with both haps in it; one reference with both haps (no combine)
            elif isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                bd.read_filter_sort_dip_bed(lambda si: si.low_complexity.simreps, i, o)
            # one bed with both haps in it; two references for both haps (split)
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.simreps, i, o, hap
                )
        case ([i], [o0, o1]):
            # one bed and one ref for a single haplotype in a diploid reference
            if isinstance(bd, cfg.Diploid2BuildData) and hap is None:
                bd.read_filter_sort_dip_bed(
                    lambda si: si.low_complexity.simreps,
                    i,
                    (o0, o1),
                )
            else:
                assert False, "this should not happen"
        case ([i0, i1], [o]):
            # two beds for both haps; one reference with both haps (combine)
            if isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.simreps, (i0, i1), o
                )
            else:
                assert False, "this should not happen"
        case _:
            assert False, "this should not happen"


# def main(smk: Any, sconf: cfg.GiabStrats) -> None:
#     rk = cfg.RefKey(smk.wildcards["ref_key"])
#     bk = cfg.BuildKey(smk.wildcards["build_key"])
#     ss = sconf.refkey_to_strat(rk).low_complexity.simreps
#     assert ss is not None, "this should not happen"
#     read_filter_sort_bed(sconf, smk.input[0], smk.output[0], ss.params, rk, bk)


main(snakemake, snakemake.config)  # type: ignore
