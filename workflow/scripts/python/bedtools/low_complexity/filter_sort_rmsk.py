from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import (
    read_filter_sort_hap_bed,
    read_filter_sort_half_hap_bed,
    read_filter_sort_full_hap_bed,
    read_filter_sort_dip1_bed,
    read_filter_sort_dip2_bed,
)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    o = smk.output[0]
    ws: dict[str, str] = smk.wildcards
    ins: list[Path] = smk.input
    rk, bk, sd = sconf.strat_dict(ws["ref_key"], ws["build_key"])
    hap = cfg.to_haplotype(ws["hap"])

    if isinstance(rk, cfg.HaploidRefKey):
        # one haplotype for both bed and ref (no combine)
        match ins:
            case [i]:
                rm = sconf.haploid_stratifications[rk].strat_inputs.low_complexity.rmsk
                assert rm is not None, "this shouldn't happen"
                read_filter_sort_hap_bed(
                    sconf, i, o, rm.params, rk, bk, rm.data.chr_pattern, [rm.class_col]
                )
            case _:
                assert False, "this should not happen"

    elif isinstance(rk, cfg.Diploid1RefKey):
        # one bed with both haps in it; one reference with both haps (no combine)
        rm0 = sconf.diploid1_stratifications[rk].strat_inputs.low_complexity.rmsk
        assert rm0 is not None, "this shouldn't happen"
        pat = rm0.data.chr_pattern
        if len(ins) == 2 and isinstance(pat, cfg.Diploid_):
            read_filter_sort_dip2_bed(
                sconf,
                (ins[0], ins[1]),
                o,
                rm0.params,
                rk,
                bk,
                pat,
                [rm0.class_col],
            )
        # two beds for both haps; one reference with both haps (combine)
        elif len(ins) == 1 and isinstance(pat, cfg.DipChrPattern):
            read_filter_sort_dip1_bed(
                sconf,
                ins[0],
                o,
                rm0.params,
                rk,
                bk,
                pat,
                [rm0.class_col],
            )
        else:
            assert False, "this should not happen"

    elif isinstance(rk, cfg.Diploid2RefKey) and hap is not None:
        # one bed with both haps in it; two references for both haps (split)
        rm1 = sconf.diploid2_stratifications[rk].strat_inputs.low_complexity.rmsk
        assert rm1 is not None, "this should not happen"
        pat = rm1.data.chr_pattern
        if len(ins) == 1 and isinstance(pat, cfg.DipChrPattern):
            read_filter_sort_half_hap_bed(
                sconf, ins[0], o, rm1.params, rk, bk, pat, [rm1.class_col]
            )
        # one bed and one ref for a single haplotype in a diploid reference
        elif len(ins) == 1 and isinstance(pat, cfg.Diploid_):
            # TODO
            read_filter_sort_full_hap_bed(
                sconf, ins[0], o, rm1.params, rk, bk, pat, [rm1.class_col]
            )
        else:
            assert False, "this should not happen"

    else:
        # assert_never(p.ref)
        assert False, "this should not happen"


main(snakemake, snakemake.config)  # type: ignore
