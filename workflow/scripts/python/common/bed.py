import pandas as pd
from pathlib import Path
import common.config as cfg
from Bio import bgzf  # type: ignore
import csv


def read_bed(
    path: Path,
    b: cfg.BedFileParams = cfg.BedFileParams(),
    more: list[int] = [],
) -> pd.DataFrame:
    """Read a bed file as a pandas dataframe.

    Return a dataframe where the first three columns are numbered 0, 1, 2 and
    typed str, int, int (first is str regardless of how the chr names are
    formated). Columns from 'more' are appended to the end of the dataframe
    in the order given starting from 3.
    """
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


def write_bed(path: Path, df: pd.DataFrame) -> None:
    """Write a bed file in bgzip format from a dataframe.

    Dataframe is not checked to make sure it is a "real" bed file.
    """
    with bgzf.open(path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        for r in df.itertuples(index=False):
            w.writerow(r)


def sort_bed_numerically(df: pd.DataFrame) -> pd.DataFrame:
    """Sort a bed file encoded by a dataframe.

    Assumes the first three columns correspond to coordinates, and that
    all are integer typed.
    """
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
    """Filter and sort a bed file.

    Arguments:
    from_map - dict containing chr name -> int mappings (int = order)
    to_map - dict containing int -> chr name mappings
    df - dataframe to sort

    Assumes the first three columns correspond to the coordinates of a bed
    file.

    Any chr name not specified in 'from_map' will be removed (hence the filter).
    Furthermore, 'to_map' should contain at least all corresponding entries
    from 'from_map', otherwise the final df will have NaNs.
    """
    chr_col = df.columns.tolist()[0]
    df[chr_col] = df[chr_col].map(from_map)
    df = sort_bed_numerically(df.dropna(subset=[chr_col]))
    df[chr_col] = df[chr_col].map(to_map)
    return df


def filter_sort_bed(
    conv: cfg.ChrConversion,
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Filter and sort a bed file from a dataframe."""
    from_map = cfg.conversion_to_init_mapper(conv)
    to_map = cfg.conversion_to_final_mapper(conv)
    return filter_sort_bed_inner(from_map, to_map, df)


def read_filter_sort_bed(
    sconf: cfg.GiabStrats,
    ipath: Path,
    opath: Path,
    bp: cfg.BedFileParams,
    rk: cfg.RefKey,
    bk: cfg.BuildKey,
    more: list[int] = [],
) -> None:
    """Read a bed file, sort it, and write it in bgzip format."""
    conv = sconf.buildkey_to_chr_conversion(rk, bk, bp.chr_prefix)
    df = read_bed(ipath, bp, more)
    df_ = filter_sort_bed(conv, df)
    write_bed(opath, df_)
