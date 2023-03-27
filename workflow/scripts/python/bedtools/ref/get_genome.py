from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed_inner
import pandas as pd


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)
    # since this is generated from the reference itself, the chr -> int map
    # is just the reverse of the final int -> chr map
    from_map = {v: k for k, v in to_map.items()}

    # ASSUME the input for this is a .fa.fai file (columns = chr, length)
    df = pd.read_table(smk.input[0], header=None, dtype={0: str, 1: int})

    filtered = filter_sort_bed_inner(from_map, to_map, df)
    filtered.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
