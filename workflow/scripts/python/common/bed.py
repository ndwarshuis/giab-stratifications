import gzip
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
    # NOTE: the 'comment="#"' parameter in pandas will strip everything after
    # the '#' in the line, which if at the beginning will include the entire
    # line and it will be skipped, and if not will only obliterate the
    # remainder. Either way this is a problem, since I only care about initial
    # lines starting with '#'. Bed files shouldn't have comments in the middle,
    # and some 'bed' files use '#' as a delimiter within a field.
    #
    # This hacky bit will count the number of lines starting with a '#' and add
    # to the original "skip_lines" parameter, thereby skipping all starting
    # comments as well as the number of lines we wish to skip with 'skip_lines'.
    total_skip = b.skip_lines
    with gzip.open(path, "rt") as f:
        while line := next(f, None):
            if line.startswith("#"):
                total_skip += 1
            else:
                break
    df = pd.read_table(
        path,
        header=None,
        usecols=bedcols,
        sep=b.sep,
        skiprows=total_skip,
        # satisfy type checker :/
        dtype={
            **{k: v for k, v in b.bed_cols.columns.items()},
            **{m: str for m in more},
        },
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


def sort_bed_numerically(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Sort a bed file encoded by a dataframe.

    Assumes the first three columns correspond to coordinates, and that all are
    integer typed. Use 'n = 2' to sort only by chr/start, and 'n=1' to sort only
    by chr.

    """
    cols = df.columns.tolist()
    bycols = [cols[i] for i in range(0, n)]
    return df.sort_values(
        by=bycols,
        axis=0,
        ignore_index=True,
    )


def filter_sort_bed_inner(
    from_map: dict[str, int],
    to_map: dict[int, str],
    df: pd.DataFrame,
    n: int = 3,
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
    df = sort_bed_numerically(df.dropna(subset=[chr_col]), n)
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
    conv = sconf.buildkey_to_chr_conversion(rk, bk, bp.chr_pattern)
    df = read_bed(ipath, bp, more)
    df_ = filter_sort_bed(conv, df)
    write_bed(opath, df_)
