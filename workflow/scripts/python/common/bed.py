import gzip
import pandas as pd
from typing import Type, NewType
from pathlib import Path
from Bio import bgzf  # type: ignore
import csv

InternalChrIndex = NewType("InternalChrIndex", int)
InitMapper = dict[str, InternalChrIndex]
FinalMapper = dict[InternalChrIndex, str]
SplitMapper = dict[str, bool]


def make_split_mapper(im: InitMapper, fm: FinalMapper) -> SplitMapper:
    return {n: i in fm for n, i in im.items()}


def read_bed(
    path: Path,
    columns: dict[int, Type[int | str]],
    skip_lines: int,
    sep: str,
    more: list[int],
) -> pd.DataFrame:
    """Read a bed file as a pandas dataframe.

    Return a dataframe where the first three columns are numbered 0, 1, 2 and
    typed str, int, int (first is str regardless of how the chr names are
    formated). Columns from 'more' are appended to the end of the dataframe
    in the order given starting from 3.
    """
    bedcols = [*columns, *more]
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
    total_skip = skip_lines
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
        sep=sep,
        skiprows=total_skip,
        # satisfy type checker :/
        dtype={
            **{k: v for k, v in columns.items()},
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


def filter_sort_bed(
    from_map: InitMapper,
    to_map: FinalMapper,
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


def split_bed(
    split_map: SplitMapper,
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # TODO doc string
    chr_col = df.columns.tolist()[0]
    sp = df[chr_col].map(split_map)
    return df[sp], df[~sp]


# def filter_sort_bed(
#     conv: cfg.ChrConversion_,
#     df: pd.DataFrame,
# ) -> pd.DataFrame:
#     """Filter and sort a bed file from a dataframe."""
#     from_map = conv.init_mapper
#     to_map = conv.final_mapper
#     # from_map = cfg.conversion_to_init_mapper(conv)
#     # to_map = cfg.conversion_to_final_mapper(conv)
#     return filter_sort_bed_inner(from_map, to_map, df)


# def read_filter_sort_hap_bed(
#     sconf: cfg.GiabStrats,
#     ipath: Path,
#     opath: Path,
#     bp: cfg.BedFileParams,
#     rk: cfg.HaploidRefKey,
#     bk: cfg.HaploidBuildKey,
#     pat: cfg.HapChrPattern,
#     more: list[int] = [],
# ) -> None:
#     """Read a haploid bed file, sort it, and write it in bgzip format."""
#     conv = sconf.haploid_stratifications.to_build_data_unsafe(rk, bk).chr_conversion(
#         pat
#     )
#     df = read_bed(ipath, bp, more)
#     df_ = filter_sort_bed(conv, df)
#     write_bed(opath, df_)


# # one diploid bed to be split into two haplotypes
# def read_filter_sort_half_hap_bed(
#     sconf: cfg.GiabStrats,
#     ipath: Path,
#     opath: tuple[Path, Path],
#     bp: cfg.BedFileParams,
#     rk: cfg.Diploid2RefKey,
#     bk: cfg.Diploid2BuildKey,
#     pat: cfg.DipChrPattern,
#     more: list[int] = [],
# ) -> None:
#     conv = sconf.diploid2_stratifications.to_build_data_unsafe(
#         rk, bk
#     ).dip_chr_conversion(pat)
#     imap, splitter = conv.init_mapper
#     fmap0, fmap1 = conv.final_mapper

#     def go(
#         o: Path,
#         df: pd.DataFrame,
#         fmap: dict[int, str],
#     ) -> None:
#         df_ = filter_sort_bed_inner(imap, fmap, df)
#         write_bed(o, df_)

#     df = read_bed(ipath, bp, more)
#     df0, df1 = split_bed(splitter, df)
#     go(opath[0], df0, fmap0)
#     go(opath[1], df1, fmap1)


# # def read_filter_sort_full_hap_bed(
# #     sconf: cfg.GiabStrats,
# #     ipath: Path,
# #     opath: Path,
# #     bp: cfg.BedFileParams,
# #     p: cfg.Diploid2BuildPair,
# #     pat: cfg.Diploid_[cfg.HapChrPattern],
# #     more: list[int] = [],
# # ) -> None:
# #     conv = sconf.diploid2_stratifications.buildkey_to_hap_chr_conversion(
# #         p.ref, p.build, pat
# #     )
# #     df = read_bed(ipath, bp, more)
# #     df_ = filter_sort_bed(conv, df)
# #     write_bed(opath, df_)
# #     # conv = sconf.buildkey_to_dip_chr_conversion(p, pat)
# #     # df = read_bed(ipath, bp, more)
# #     # df_ = filter_sort_bed(conv, df)
# #     # write_bed(opath, df_)
# #     pass


# def read_filter_sort_dip1_bed(
#     sconf: cfg.GiabStrats,
#     ipath: Path,
#     opath: Path,
#     bp: cfg.BedFileParams,
#     rk: cfg.Diploid1RefKey,
#     bk: cfg.Diploid1BuildKey,
#     pat: cfg.DipChrPattern,
#     more: list[int] = [],
# ) -> None:
#     """Read a diploid bed file, sort it, and write it in bgzip format."""
#     conv = sconf.diploid1_stratifications.to_build_data_unsafe(
#         rk, bk
#     ).dip_chr_conversion(pat)
#     df = read_bed(ipath, bp, more)
#     df_ = filter_sort_bed(conv, df)
#     write_bed(opath, df_)


# def read_filter_sort_dip2_bed(
#     sconf: cfg.GiabStrats,
#     ipath: tuple[Path, Path],
#     opath: Path,
#     bp: cfg.BedFileParams,
#     rk: cfg.Diploid1RefKey,
#     bk: cfg.Diploid1BuildKey,
#     pat: cfg.Diploid_[cfg.HapChrPattern],
#     more: list[int] = [],
# ) -> None:
#     """Read two haploid bed files, combine and sort them as diploid, and write
#     it in bgzip format.
#     """
#     conv = sconf.diploid1_stratifications.to_build_data_unsafe(
#         rk, bk
#     ).hap_chr_conversion(pat)
#     imap1, imap2 = conv.init_mapper
#     fmap = conv.final_mapper

#     def go(i: Path, imap: dict[str, int]) -> pd.DataFrame:
#         df = read_bed(i, bp, more)
#         return filter_sort_bed_inner(imap, fmap, df)

#     df = pd.concat(
#         [
#             go(*x)
#             for x in [
#                 (ipath[0], imap1),
#                 (ipath[1], imap2),
#             ]
#         ]
#     )
#     write_bed(opath, df)
