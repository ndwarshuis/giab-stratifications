import pandas as pd
from pathlib import Path
import json
from typing import Any, NamedTuple
import common.config as cfg
from common.bed import (
    filter_sort_bed,
    read_bed,
    InitMapper,
    FinalMapper,
    make_split_mapper,
    split_bed,
)


class CDSInput(NamedTuple):
    bed: Path
    mapper: Path


def load_init_mapper(p: Path) -> InitMapper:
    with open(p, "r") as f:
        m: InitMapper = json.load(f)
    assert all([isinstance(k, str) and isinstance(v, int) for k, v in m.items()])
    return m


def _read_df(i: Path) -> pd.DataFrame:
    df = read_bed(i, {0: str, 1: int, 2: int}, 0, "\t", more=[3, 4])
    return df[df[3].str.contains("RefSeq") & (df[4] == "CDS")][[0, 1, 2]].copy()


def read_df(i: Path, m: Path, fmap: FinalMapper) -> pd.DataFrame:
    from_map = load_init_mapper(m)
    df = _read_df(i)
    return filter_sort_bed(from_map, fmap, df)


def write_df(path: Path, df: pd.DataFrame) -> None:
    df.to_csv(path, header=False, index=False, sep="\t")


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ins = [CDSInput(b, m) for b, m in zip(smk.input["bed"], smk.input["mapper"])]
    ons: list[Path] = smk.output

    def hap(i: CDSInput, o: Path, bd: cfg.HaploidBuildData) -> None:
        write_df(o, read_df(i.bed, i.mapper, bd.final_mapper))

    def dip_to_dip(i: CDSInput, o: Path, bd: cfg.Diploid1BuildData) -> None:
        write_df(o, read_df(i.bed, i.mapper, bd.final_mapper))

    def hap_to_hap(
        i: CDSInput,
        o: Path,
        hap: cfg.Haplotype,
        bd: cfg.Diploid2BuildData,
    ) -> None:
        fmap = bd.final_mapper
        write_df(o, read_df(i.bed, i.mapper, hap.from_either(fmap[0], fmap[1])))

    def dip_to_hap(
        i: CDSInput,
        os: tuple[Path, Path],
        bd: cfg.Diploid2BuildData,
    ) -> None:
        df = _read_df(i.bed)
        im = load_init_mapper(i.mapper)
        fmap0, fmap1 = bd.final_mapper
        sm = make_split_mapper(im, fmap0)
        df0, df1 = split_bed(sm, df)
        for o, df_, fmap_ in zip(os, (df0, df1), (fmap0, fmap1)):
            write_df(o, filter_sort_bed(im, fmap_, df_))

    def hap_to_dip(
        i: tuple[CDSInput, CDSInput],
        o: Path,
        bd: cfg.Diploid1BuildData,
    ) -> None:
        df = pd.concat([read_df(bed, mapper, bd.final_mapper) for bed, mapper in i])
        write_df(o, df)

    sconf.with_build_src_data_unsafe(
        ws["ref_key"],
        ws["build_key"],
        ins,
        ons,
        cfg.to_haplotype(ws["hap"]),
        hap,
        dip_to_dip,
        hap_to_hap,
        dip_to_hap,
        hap_to_dip,
    )


main(snakemake, snakemake.config)  # type: ignore
