import pandas as pd
import common.config as cfg
from typing import Any, Callable


def sort_bed_numerically(df: pd.DataFrame) -> pd.DataFrame:
    cols = df.columns.tolist()
    return df.sort_values(
        by=[cols[0], cols[1], cols[2]],
        axis=0,
        ignore_index=True,
    )


def filter_sort_bed(
    sconf: cfg.GiabStrats,
    f: Callable[[cfg.Stratification], cfg.BedFile],
    rk: cfg.RefKey,
    bk: cfg.BuildKey,
    df: pd.DataFrame,
) -> pd.DataFrame:
    from_map = sconf.buildkey_to_init_chr_mapping(f, rk, bk)
    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)
    df[0] = df[0].map(from_map)
    df = sort_bed_numerically(df.dropna(subset=[0]))
    df[0] = df[0].map(to_map)
    return df


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    bedfile = smk.input["bed"][0]
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    rmsk = sconf.stratifications[rk].low_complexity.rmsk
    bedcols = [rmsk.chr_col, rmsk.start_col, rmsk.end_col, rmsk.class_col]
    df = pd.read_table(bedfile, header=None, usecols=bedcols)
    df.columns = pd.Index(range(len(bedcols)))
    df_ = filter_sort_bed(sconf, lambda x: x.low_complexity.rmsk, rk, bk, df)
    print(df_)
    df_.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
