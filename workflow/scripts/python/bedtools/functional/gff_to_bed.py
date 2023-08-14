import pandas as pd
from typing import Any
import common.config as cfg
from common.bed import write_bed

JOIN_KEY = "genomic_accession"


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
    df = read_gff(smk.input[0])
    write_bed(smk.output[0], df)


main(snakemake, snakemake.config)  # type: ignore
