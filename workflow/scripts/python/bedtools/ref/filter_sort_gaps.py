from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed
import pandas as pd
from pybedtools import BedTool as bt  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    def get_bed(x: cfg.Stratification) -> cfg.BedFile:
        b = x.gap
        assert b is not None
        return b

    # TODO not dry
    bedfile = get_bed(sconf.stratifications[rk])

    df = pd.read_table(
        smk.input[0],
        header=None,
        usecols=bedfile.bed_cols.columns,
        comment="#",
        skiprows=bedfile.skip_lines,
    )

    filtered = filter_sort_bed(sconf, get_bed, rk, bk, df)
    merged = bt.from_dataframe(filtered).merge(d=100).to_dataframe()
    merged.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
