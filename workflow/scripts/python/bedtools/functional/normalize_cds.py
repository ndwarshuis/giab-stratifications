import pandas as pd
from pathlib import Path
import json
from typing import Any, Callable, TypeVar
import common.config as cfg
from common.functional import DesignError, match1_unsafe, match12_unsafe, both
from common.bed import filter_sort_bed, InitMapper, split_bed, write_bed

X = TypeVar("X")

# This seems a bit wonky since it is unlike many of the other chromosome name
# mapping operations. The FTBL file has a mapping b/t bare chromosome names (ie
# 1-22, X, Y) and accession numbers (ie "NC_gibberish.420"), and the accession
# numbers are what is reported in GFF files. Thus we need to translate these
# accession numbers to the final chromosomal names on the target reference.
#
# Thus the mapper created here connects accession numbers to each chromosomal
# index; the main difference b/t this and the other "initial mappers" is that
# the "accession numbers" are the "chromosome names" found in a normal input bed
# file.
FTBLMapper = InitMapper

VDJ_PAT = "^ID=gene-(IGH|IGK|IGL|TRA|TRB|TRG);"


# def _read_df(i: Path) -> pd.DataFrame:
#     df = read_bed(i, {0: str, 1: int, 2: int}, 0, "\t", more=[3, 4])
#     return df[df[3].str.contains("RefSeq") & (df[4] == "CDS")][[0, 1, 2]].copy()


def read_gff(i: Path) -> pd.DataFrame:
    """Read a gff file and return a condensed bed-like dataframe."""
    # Pull the following columns and rearrange like so:
    # 0 -> 0: accession number
    # 3 -> 1: start pos
    # 4 -> 2: end pos
    # 1 -> 3: source
    # 2 -> 4: type (eg gene vs exon)
    # 8 -> 5: attributes
    #
    # NOTE that source and type are used for creating the CDS bed, and
    # attributes are used for creating the VDJ bed (hence why all are included)
    columns = {0: str, 3: int, 4: int, 1: str, 2: str, 8: str}
    df = pd.read_table(
        i,
        header=None,
        comment="#",
        dtype={k: v for k, v in columns.items()},
        usecols=list(columns),
    )[list(columns)]
    df = df.set_axis(range(len(df.columns)), axis=1)
    # for some reason there are lots of rows where the start/end are the same;
    # these are useless to us so remove them
    return df[df[1] != df[2]].copy()


def write_gff(
    o: Path,
    mask_fun: Callable[[pd.DataFrame], "pd.Series[bool]"],
    df: pd.DataFrame,
) -> None:
    write_bed(o, df[mask_fun(df)][[0, 1, 2]])


def cds_mask(df: pd.DataFrame) -> "pd.Series[bool]":
    refseq_mask = df[3].str.contains("RefSeq")
    cds_mask = df[4] == "CDS"
    return refseq_mask & cds_mask


def vdj_mask(df: pd.DataFrame) -> "pd.Series[bool]":
    return df[5].str.match(VDJ_PAT)


def read_ftbl(path: Path, cis: set[cfg.ChrIndex], hap: cfg.Haplotype) -> FTBLMapper:
    filter_cols = ["assembly_unit", "seq_type"]
    map_cols = ["chromosome", "genomic_accession"]
    df = pd.read_table(path, header=0, usecols=filter_cols + map_cols, dtype=str)
    chr_mask = df["seq_type"] == "chromosome"
    asm_mask = df["assembly_unit"] == "Primary Assembly"
    ser = (
        df[chr_mask & asm_mask][map_cols]
        .drop_duplicates()
        .copy()
        .set_index("chromosome")["genomic_accession"]
    )
    try:
        return {ser[i.chr_name]: i.to_internal_index(hap) for i in cis}
    except KeyError:
        raise DesignError("Feature table has wonky chromosome names, fixmeplz")


def iamnotlivingimasleep(_: Any) -> Any:
    """Don't go down the rabbit whole"""
    raise DesignError("NOT IMPLEMENTED")


def write_outputs(p: Path, xs: list[Path]) -> None:
    with open(p, "w") as f:
        json.dump([str(x) for x in xs], f)


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ftbl_inputs: list[Path] = [Path(i) for i in smk.input["ftbl"]]
    gff_inputs: list[Path] = [Path(i) for i in smk.input["gff"]]
    cds_pattern: str = ps["cds_output"]
    vdj_pattern: str = ps["vdj_output"]

    def write_vdj_maybe1(
        switch: bool,
        rfk: cfg.RefKeyFull,
        df: pd.DataFrame,
    ) -> list[Path]:
        if switch:
            v = cfg.sub_output_path(vdj_pattern, rfk)
            write_gff(v, vdj_mask, df)
            return [v]
        else:
            return []

    def write_vdj_maybe2(
        switch: bool,
        rfks: tuple[cfg.RefKeyFull, cfg.RefKeyFull],
        df: tuple[pd.DataFrame, pd.DataFrame],
    ) -> list[Path]:
        if switch:
            vs: tuple[Path, Path] = both(
                lambda r: cfg.sub_output_path(vdj_pattern, r), rfks
            )
            write_gff(vs[0], vdj_mask, df[0])
            write_gff(vs[1], vdj_mask, df[1])
            return [*vs]
        else:
            return []

    def hap(bd: cfg.HapBuildData) -> tuple[list[Path], list[Path]]:
        rd = bd.refdata
        rk = rd.ref.src.key(rd.refkey)
        c = cfg.sub_output_path(cds_pattern, rk)

        def go(f: Path, g: Path) -> pd.DataFrame:
            im = read_ftbl(f, bd.chr_indices, cfg.Haplotype.HAP1)
            fm = bd.refdata.ref.chr_pattern.final_mapper(
                bd.chr_indices, cfg.Haplotype.HAP1
            )
            gff = read_gff(g)
            gff_ = filter_sort_bed(im, fm, gff)
            write_gff(c, cds_mask, gff_)
            return gff_

        gff = match1_unsafe(
            list(zip(ftbl_inputs, gff_inputs)),
            lambda x: go(*x),
        )

        v = write_vdj_maybe1(bd.build.include.vdj, rk, gff)
        return [c], v

    def dip1(bd: cfg.Dip1BuildData) -> tuple[list[Path], list[Path]]:
        fm = bd.refdata.ref.chr_pattern.final_mapper(bd.chr_indices)
        rd = bd.refdata
        rk = rd.ref.src.key(rd.refkey)

        im = match12_unsafe(
            ftbl_inputs,
            iamnotlivingimasleep,
            lambda f0, f1: {
                **read_ftbl(f0, bd.chr_indices, cfg.Haplotype.HAP1),
                **read_ftbl(f1, bd.chr_indices, cfg.Haplotype.HAP2),
            },
        )

        gff = match12_unsafe(
            gff_inputs,
            lambda g: filter_sort_bed(im, fm, read_gff(g)),
            # NOTE this should be ok to do since the accession numbers should be
            # unique
            lambda g0, g1: filter_sort_bed(
                im, fm, pd.concat([read_gff(g) for g in [g0, g1]])
            ),
        )

        c = cfg.sub_output_path(cds_pattern, rk)
        write_gff(c, cds_mask, gff)

        v = write_vdj_maybe1(bd.build.include.vdj, rk, gff)
        return [c], v

    def dip2(bd: cfg.Dip2BuildData) -> tuple[list[Path], list[Path]]:
        fm0, fm1 = bd.refdata.ref.chr_pattern.both(
            lambda x, hap: x.final_mapper(bd.chr_indices, hap)
        )
        im0, im1 = match12_unsafe(
            ftbl_inputs,
            iamnotlivingimasleep,
            lambda f0, f1: (
                read_ftbl(f0, bd.chr_indices, cfg.Haplotype.HAP1),
                read_ftbl(f1, bd.chr_indices, cfg.Haplotype.HAP2),
            ),
        )
        gffs = match12_unsafe(
            gff_inputs,
            # TODO set a new PR for number of characters in one lambda :)
            lambda g: (
                filter_sort_bed(
                    im0,
                    fm0,
                    (split := split_bed({k: True for k in im0}, read_gff(g)))[0],
                ),
                filter_sort_bed(im1, fm1, split[1]),
            ),
            lambda g0, g1: (
                filter_sort_bed(im0, fm0, read_gff(g0)),
                filter_sort_bed(im1, fm1, read_gff(g1)),
            ),
        )

        rd = bd.refdata
        rks = rd.ref.src.keys(rd.refkey)

        cs = both(lambda r: cfg.sub_output_path(cds_pattern, r), rks)

        write_gff(cs[0], cds_mask, gffs[0])
        write_gff(cs[1], cds_mask, gffs[1])

        vs = write_vdj_maybe2(bd.build.include.vdj, rks, gffs)
        return [*cs], vs

    cs, vs = sconf.with_build_data(
        cfg.wc_to_refkey(ws),
        cfg.wc_to_buildkey(ws),
        hap,
        dip1,
        dip2,
    )

    write_outputs(smk.output["cds"], cs)
    write_outputs(smk.output["vdj"], vs)


main(snakemake, snakemake.config)  # type: ignore
