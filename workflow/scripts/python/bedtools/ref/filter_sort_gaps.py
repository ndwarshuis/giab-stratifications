from typing import Any
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.bed import filter_sort_bed, read_bed, write_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    bedfile = sconf.stratifications[rk].gap
    assert bedfile is not None, "this should not happen"

    df = read_bed(smk.input[0], bedfile)

    conv = sconf.buildkey_to_chr_conversion(rk, bk, bedfile.chr_prefix)

    filtered = filter_sort_bed(conv, df)
    merged = bt.from_dataframe(filtered).merge(d=100).to_dataframe()
    write_bed(smk.output[0], merged)


main(snakemake, snakemake.config)  # type: ignore
