import pandas as pd
from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed_inner

JOIN_KEY = "genomic_accession"

# ASSUME chromosomes in the ftbl file are like "1, 2, 3...X, Y"
FROM_MAP = {x.chr_name: x.value for x in cfg.ChrIndex}


def read_ftbl(path: str) -> pd.DataFrame:
    cols = ["assembly_unit", "seq_type", "chromosome", JOIN_KEY]
    df = pd.read_table(
        path,
        header=0,
        usecols=cols,
        dtype=str,
    ).drop_duplicates()
    chr_mask = df["seq_type"] == "chromosome"
    asm_mask = df["assembly_unit"] == "Primary Assembly"
    return df[chr_mask & asm_mask].set_index(JOIN_KEY)[["chromosome"]].copy()


def read_gff(path: str) -> pd.DataFrame:
    cols = [JOIN_KEY, "source", "type", "start", "end"]
    df = pd.read_table(
        path,
        header=None,
        comment="#",
        names=cols,
        usecols=cols,
    )
    return (
        df[(df["type"] == "CDS") & df["source"].str.contains("RefSeq")]
        .copy()
        .set_index(JOIN_KEY)
    )


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ftbl_df = read_ftbl(smk.input["ftbl"][0])
    gff_df = read_gff(smk.input["gff"][0])

    df = (
        ftbl_df.join(gff_df, how="left")[["chromosome", "start", "end"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)
    # TODO this seems hacky
    from_map = {k: v for k, v in FROM_MAP.items() if v in to_map}

    filter_sort_bed_inner(from_map, to_map, df).to_csv(
        smk.output[0],
        header=False,
        index=False,
        sep="\t",
    )


main(snakemake, snakemake.config)  # type: ignore
