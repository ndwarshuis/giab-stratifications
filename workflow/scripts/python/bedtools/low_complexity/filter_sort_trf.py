from pathlib import Path
from typing import Any, TypeVar, Callable
import common.config as cfg

X = TypeVar("X")
Y = TypeVar("Y")

# from common.bed import read_filter_sort_bed


# TODO ewwwwwww
def giant_switch_thingy(
    inputs: list[X],
    outputs: list[Y],
    hap: cfg.Haplotype | None,
    bd: cfg.HaploidBuildData | cfg.Diploid1BuildData | cfg.Diploid2BuildData,
    hap_f: Callable[[X, Y, cfg.HaploidBuildData], None],
    dip_to_dip_f: Callable[[X, Y, cfg.Diploid1BuildData], None],
    hap_to_hap_f: Callable[[X, Y, cfg.Diploid2BuildData], None],
    dip_to_hap_f: Callable[[X, Y, Y, cfg.Diploid2BuildData], None],
    hap_to_dip_f: Callable[[X, X, Y, cfg.Diploid1BuildData], None],
) -> None:
    match (inputs, outputs):
        case ([i], [o]):
            # one haplotype for both bed and ref (no combine)
            if isinstance(bd, cfg.HaploidBuildData) and hap is None:
                hap_f(i, o, bd)
            # one bed with both haps in it; one reference with both haps (no combine)
            elif isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                dip_to_dip_f(i, o, bd)
            # one bed and one ref for a single haplotype in a diploid reference
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                hap_to_hap_f(i, o, bd)
        case ([i], [o0, o1]):
            # one bed with both haps in it; two references for both haps (split)
            if isinstance(bd, cfg.Diploid2BuildData) and hap is None:
                dip_to_hap_f(i, o0, o1, bd)
            else:
                assert False, "this should not happen"
        case ([i0, i1], [o]):
            # two beds for both haps; one reference with both haps (combine)
            if isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                hap_to_dip_f(i0, i1, o, bd)
            else:
                assert False, "this should not happen"
        case _:
            assert False, "this should not happen"


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
