import pandas as pd
from pathlib import Path
import json
from typing import Any
import common.config as cfg
from common.bed import (
    filter_sort_bed,
    read_bed,
    InitMapper,
    FinalMapper,
    make_split_mapper,
    split_bed,
)


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
    ons: list[Path] = smk.output[0]

    ws: dict[str, str] = smk.wildcards
    bed_ins: list[Path] = smk.input["bed"]
    map_ins: list[Path] = smk.input["mapper"]
    bd = sconf.to_build_data(ws["ref_key"], ws["build_key"])
    hap = cfg.to_haplotype(ws["hap"])
    match (bed_ins, map_ins, ons):
        case ([bi], [mi], [o]):
            # one haplotype for both bed and ref (no combine)
            if isinstance(bd, cfg.HaploidBuildData) and hap is None:
                df = read_df(bi, mi, bd.final_mapper)
            # one bed with both haps in it; one reference with both haps (no combine)
            elif isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                df = read_df(bi, mi, bd.final_mapper)
            # one bed with both haps in it; two references for both haps (split)
            elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                fmap = bd.final_mapper
                df = read_df(bi, mi, hap.from_either(fmap[0], fmap[1]))
            else:
                assert False, "this should not happen"
            write_df(o, df)
        case ([bi], [mi], [o0, o1]):
            # one bed and one ref for a single haplotype in a diploid reference
            if isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
                # TODO bleh....
                im = load_init_mapper(mi)
                df = _read_df(bi)
                fmap0, fmap1 = bd.final_mapper
                sm = make_split_mapper(im, fmap0)
                df0, df1 = split_bed(sm, df)
                df0_ = filter_sort_bed(im, fmap0, df0)
                df1_ = filter_sort_bed(im, fmap1, df1)
                write_df(o0, df0_)
                write_df(o1, df1_)
            else:
                assert False, "this should not happen"

        case ([bi0, bi1], [mi0, mi1], [o]):
            # two beds for both haps; one reference with both haps (combine)
            if isinstance(bd, cfg.Diploid1BuildData) and hap is None:
                df0 = read_df(bi0, mi0, bd.final_mapper)
                df1 = read_df(bi1, mi1, bd.final_mapper)
                write_df(o, pd.concat([df0, df1]))
            else:
                assert False, "this should not happen"
        case _:
            assert False, "this should not happen"


main(snakemake, snakemake.config)  # type: ignore
