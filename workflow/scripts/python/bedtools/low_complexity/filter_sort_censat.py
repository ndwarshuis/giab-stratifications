from typing import Any
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    bedfile = sconf.refkey_to_strat(rk).low_complexity.satellites
    assert bedfile is not None, "this should not happen"

    ps = bedfile.params
    conv = sconf.buildkey_to_chr_conversion(rk, bk, ps.chr_pattern)

    df = read_bed(smk.input[0], ps, [bedfile.sat_col])
    df = filter_sort_bed(conv, df)
    df = df[~df[3].str.startswith("ct_")]
    write_bed(smk.output[0], df)


main(snakemake, snakemake.config)  # type: ignore
