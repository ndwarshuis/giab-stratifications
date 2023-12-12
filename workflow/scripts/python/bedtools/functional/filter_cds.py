import pandas as pd
from pathlib import Path
import json
from typing import Any, Callable
from typing_extensions import assert_never
import common.config as cfg
from common.io import DesignError, match1_unsafe, match12_unsafe, match2_unsafe
from common.bed import filter_sort_bed, read_bed, InitMapper, split_bed, write_bed

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


def _read_df(i: Path) -> pd.DataFrame:
    df = read_bed(i, {0: str, 1: int, 2: int}, 0, "\t", more=[3, 4])
    return df[df[3].str.contains("RefSeq") & (df[4] == "CDS")][[0, 1, 2]].copy()


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
    return read_bed(i, {0: str, 3: int, 4: int}, 0, "\t", more=[1, 2, 8])


def write_gff(
    o: Path,
    mask_fun: Callable[[pd.DataFrame], "pd.Series[bool]"],
    df: pd.DataFrame,
) -> None:
    write_bed(o, df[mask_fun(df)])


def cds_mask(df: pd.DataFrame) -> "pd.Series[bool]":
    refseq_mask = df[3].str.contains("RefSeq")
    cds_mask = df[4] == "CDS"
    return refseq_mask & cds_mask


def vdj_mask(df: pd.DataFrame) -> "pd.Series[bool]":
    return df[3].str.match(VDJ_PAT)


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


def write_vdj1_maybe(os: list[Path], want_vdj: bool, gff: pd.DataFrame) -> None:
    match (os, want_vdj):
        case ([v], True):
            write_gff(v, vdj_mask, gff)
        case ([], False):
            pass
        case _:
            raise DesignError(f"Invalid VDJ combination: {os}, {want_vdj}")


def iamnotlivingimasleep(_: Any) -> Any:
    """Don't go down the rabbit whole"""
    raise DesignError("NOT IMPLEMENTED")


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ftbl_inputs: list[Path] = [Path(i) for i in smk.input["ftbl"]]
    gff_inputs: list[Path] = [Path(i) for i in smk.input["gff"]]
    cds_out: Path = smk.output["cds"]
    vdj_out: Path = smk.output["cds"]
    cds_outputs = [Path(i) for i in ps["cds_outputs"]]
    vdj_outputs = [Path(i) for i in ps["vdj_outputs"]]

    rk_ = cfg.strip_refkey(ws["ref_final_key"])
    bd = sconf.to_build_data(rk_, ws["build_key"])
    want_vdj = bd.build.include.vdj

    def hap(bd: cfg.HaploidBuildData) -> None:
        def go(f: Path, g: Path, c: Path) -> pd.DataFrame:
            im = read_ftbl(f, bd.chr_indices, cfg.Haplotype.HAP1)
            gff = read_gff(g)
            gff_ = filter_sort_bed(im, bd.final_mapper, gff)
            write_gff(c, cds_mask, gff_)
            return gff_

        gff = match1_unsafe(
            list(zip(ftbl_inputs, gff_inputs, cds_outputs)),
            lambda x: go(*x),
        )

        write_vdj1_maybe(vdj_outputs, want_vdj, gff)

    def dip1(bd: cfg.Diploid1BuildData) -> None:
        fm = bd.final_mapper

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

        match1_unsafe(cds_outputs, lambda c: write_gff(c, cds_mask, gff))

        write_vdj1_maybe(vdj_outputs, want_vdj, gff)

    def dip2(bd: cfg.Diploid2BuildData) -> None:
        fm0, fm1 = bd.final_mapper
        im0, im1 = match12_unsafe(
            ftbl_inputs,
            iamnotlivingimasleep,
            lambda f0, f1: (
                read_ftbl(f0, bd.chr_indices, cfg.Haplotype.HAP1),
                read_ftbl(f1, bd.chr_indices, cfg.Haplotype.HAP2),
            ),
        )
        gff0, gff1 = match12_unsafe(
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

        # TODO ...because lambdas can't have two statements (and python doesn't
        # have the >> operator)
        def go(c0: Path, c1: Path) -> None:
            write_gff(c0, cds_mask, gff0)
            write_gff(c1, cds_mask, gff1)

        match2_unsafe(cds_outputs, go)

        match (vdj_outputs, want_vdj):
            case ([v0, v1], True):
                write_gff(v0, vdj_mask, gff0)
                write_gff(v1, vdj_mask, gff1)
            case ([], False):
                pass
            case _:
                raise DesignError(f"Invalid VDJ combination: {vdj_outputs}, {want_vdj}")

    if isinstance(bd, cfg.HaploidBuildData):
        return hap(bd)
    elif isinstance(bd, cfg.Diploid1BuildData):
        return dip1(bd)
    elif isinstance(bd, cfg.Diploid2BuildData):
        return dip2(bd)
    else:
        assert_never(bd)

    with open(cds_out, "w") as f:
        json.dump(cds_outputs, f)

    with open(vdj_out, "w") as f:
        json.dump(vdj_outputs, f)


main(snakemake, snakemake.config)  # type: ignore
