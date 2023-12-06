from pathlib import Path
import json
import pandas as pd
from typing import Any
import common.config as cfg

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


def read_ftbl(path: Path) -> "pd.Series[str]":
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


# TODO make this work for dip/hap
def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ins: list[Path] = smk.input
    ons: list[Path] = smk.output

    # if this fails then something is wrong with the ftbl (eg it doesn't have
    # a complete set of chromosomes)
    ser = read_ftbl(ins)
    try:
        mapper = {ser[i.chr_name]: i.value for i in cis}
    except KeyError:
        assert False, "FTBL file does not have a complete set of chromosomes"

    with open(smk.output[0], "w") as f:
        json.dump(mapper, f)

    sconf.with_build_src_data_unsafe(
        ws["ref_key"],
        ws["build_key"],
        cfg.to_haplotype(ws["hap"]),
        ins,
        ons,
    )


main(snakemake, snakemake.config)  # type: ignore
