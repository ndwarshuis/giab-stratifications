from typing import Any
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    prefix = sconf.refkey_to_final_chr_prefix(rk)
    conv = cfg.fullset_conv(prefix)
    df = read_bed(smk.input[0])
    df_ = filter_sort_bed(conv, df)
    write_bed(smk.output[0], df_)


main(snakemake, snakemake.config)  # type: ignore
