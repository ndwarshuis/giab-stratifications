import json
import pandas as pd
from typing import Any
import common.config as cfg


def read_ftbl(path: str) -> "pd.Series[str]":
    filter_cols = ["assembly_unit", "seq_type"]
    map_cols = ["chromosome", "genomic_accession"]
    df = pd.read_table(path, header=0, usecols=filter_cols + map_cols, dtype=str)
    chr_mask = df["seq_type"] == "chromosome"
    asm_mask = df["assembly_unit"] == "Primary Assembly"
    return (
        df[chr_mask & asm_mask][map_cols]
        .drop_duplicates()
        .copy()
        .set_index("chromosome")["genomic_accession"]
    )


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    cis = sconf.buildkey_to_chr_indices(rk, bk)

    # if this fails then something is wrong with the ftbl (eg it doesn't have
    # a complete set of chromosomes)
    ser = read_ftbl(smk.input[0])
    mapper = {ser[i.chr_name]: i.value for i in cis}

    with open(smk.output[0], "w") as f:
        json.dump(mapper, f)


main(snakemake, snakemake.config)  # type: ignore
