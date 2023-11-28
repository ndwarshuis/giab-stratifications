from typing import Any
from pathlib import Path
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    o = smk.output[0]
    ws: dict[str, str] = smk.wildcards
    ins: list[Path] = smk.input
    bd = sconf.to_build_data(ws["ref_key"], ws["build_key"])
    hap = cfg.to_haplotype(ws["hap"])
    match ins:
        case [i]:
            # one haplotype for both bed and ref (no combine)
            if isinstance(bd, cfg.HaploidBuildData) and hap is None:
                bd.read_filter_sort_hap_bed(lambda si: si.low_complexity.rmsk, i, o)
            # one bed with both haps in it; one reference with both haps (no combine)
            elif isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                bd.read_filter_sort_dip_bed(lambda si: si.low_complexity.rmsk, i, o)
            # two beds for both haps; one reference with both haps (combine)
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.rmsk, i, o, hap
                )
            else:
                assert False, "this should not happen"
        case [i0, i1]:
            # one bed with both haps in it; two references for both haps (split)
            if isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.rmsk, (i0, i1), o
                )
            # one bed and one ref for a single haplotype in a diploid reference
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                bd.read_filter_sort_dip_bed(
                    lambda si: si.low_complexity.rmsk, hap.from_either(i0, i1), o
                )
            else:
                assert False, "this should not happen"
        case _:
            assert False, "this should not happen"


main(snakemake, snakemake.config)  # type: ignore
