import pandas as pd
import common.config as cfg
from typing import Callable
from Bio import bgzf  # type: ignore
import csv


def sort_bed_numerically(df: pd.DataFrame) -> pd.DataFrame:
    cols = df.columns.tolist()
    return df.sort_values(
        by=[cols[0], cols[1], cols[2]],
        axis=0,
        ignore_index=True,
    )


def filter_sort_bed_inner(
    from_map: dict[str, int],
    to_map: dict[int, str],
    df: pd.DataFrame,
) -> pd.DataFrame:
    chr_col = df.columns.tolist()[0]
    df[chr_col] = df[chr_col].map(from_map)
    df = sort_bed_numerically(df.dropna(subset=[chr_col]))
    df[chr_col] = df[chr_col].map(to_map)
    return df


def filter_sort_bed(
    sconf: cfg.GiabStrats,
    f: Callable[[cfg.Stratification], cfg.BedFile],
    rk: cfg.RefKey,
    bk: cfg.BuildKey,
    df: pd.DataFrame,
) -> pd.DataFrame:
    from_map = sconf.buildkey_to_init_chr_mapping(f, rk, bk)
    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)
    return filter_sort_bed_inner(from_map, to_map, df)


def read_bed(path: str, b: cfg.BedFile, more: list[int]) -> pd.DataFrame:
    bedcols = [*b.bed_cols.columns, *more]
    df = pd.read_table(
        path,
        header=None,
        usecols=bedcols,
        sep=b.sep,
        comment="#",
        skiprows=b.skip_lines,
        # satisfy type checker :/
        dtype={k: v for k, v in b.bed_cols.columns.items()},
    )
    df.columns = pd.Index(range(len(bedcols)))
    return df


def write_bed(path: str, df: pd.DataFrame) -> None:
    with bgzf.open(path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        for r in df.itertuples(index=False):
            w.writerow(r)


def read_filter_sort_bed(
    sconf: cfg.GiabStrats,
    ipath: str,
    opath: str,
    f: Callable[[cfg.Stratification], cfg.BedFile],
    rk: cfg.RefKey,
    bk: cfg.BuildKey,
    more: list[int],
) -> None:
    bedfile = f(sconf.stratifications[rk])
    df = read_bed(ipath, bedfile, more)
    df_ = filter_sort_bed(sconf, f, rk, bk, df)
    write_bed(opath, df_)
