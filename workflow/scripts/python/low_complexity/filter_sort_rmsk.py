import pandas as pd
from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed


# TODO not dry
def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    bedfile = smk.input[0]
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    rmsk = sconf.stratifications[rk].low_complexity.rmsk
    bedcols = [*rmsk.bed_cols.columns, rmsk.class_col]
    df = pd.read_table(bedfile, header=None, usecols=bedcols)
    df.columns = pd.Index(range(len(bedcols)))
    df_ = filter_sort_bed(sconf, lambda x: x.low_complexity.rmsk, rk, bk, df)
    df_.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
