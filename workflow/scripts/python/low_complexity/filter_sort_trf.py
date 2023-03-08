import pandas as pd
from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    bedfile = smk.input[0]
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    bedcols = [2, 3, 4]  # TODO don't hardcode this
    df = pd.read_table(bedfile, header=None, usecols=bedcols)
    df.columns = pd.Index(range(len(bedcols)))
    df_ = filter_sort_bed(sconf, lambda x: x.low_complexity.simreps, rk, bk, df)
    df_.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
