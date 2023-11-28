import pandas as pd
from pathlib import Path
from typing import Any
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed


def filter_ct(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df[3].str.startswith("ct_")]


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
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.satellites, i, o, filter_ct
                )
            # one bed with both haps in it; one reference with both haps (no combine)
            elif isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                bd.read_filter_sort_dip_bed(
                    lambda si: si.low_complexity.satellites, i, o, filter_ct
                )
            # two beds for both haps; one reference with both haps (combine)
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.satellites, i, o, hap, filter_ct
                )
            else:
                assert False, "this should not happen"
        case [i0, i1]:
            # one bed with both haps in it; two references for both haps (split)
            if isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                bd.read_filter_sort_hap_bed(
                    lambda si: si.low_complexity.satellites, (i0, i1), o, filter_ct
                )
            # one bed and one ref for a single haplotype in a diploid reference
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                bd.read_filter_sort_dip_bed(
                    lambda si: si.low_complexity.satellites,
                    hap.from_either(i0, i1),
                    o,
                    filter_ct,
                )
            else:
                assert False, "this should not happen"
        case _:
            assert False, "this should not happen"


# def main(smk: Any, sconf: cfg.GiabStrats) -> None:
#     rk = cfg.RefKey(smk.wildcards["ref_key"])
#     bk = cfg.BuildKey(smk.wildcards["build_key"])

#     bedfile = sconf.refkey_to_strat(rk).low_complexity.satellites
#     assert bedfile is not None, "this should not happen"

#     ps = bedfile.params
#     conv = sconf.buildkey_to_chr_conversion(rk, bk, ps.chr_pattern)

#     df = read_bed(smk.input[0], ps, [bedfile.sat_col])
#     df = filter_sort_bed(conv, df)
#     df = df[~df[3].str.startswith("ct_")]
#     write_bed(smk.output[0], df)


main(snakemake, snakemake.config)  # type: ignore
