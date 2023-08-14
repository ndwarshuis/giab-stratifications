import json
from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed_inner, read_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)

    with open(smk.input["mapper"], "r") as f:
        from_map = json.load(f)

    df = read_bed(smk.input["bed"], more=[3, 4])
    df = df[df[3].str.contains("RefSeq") & (df[4] == "CDS")][[0, 1, 2]].copy()
    df = filter_sort_bed_inner(from_map, to_map, df)

    df.to_csv(
        smk.output[0],
        header=False,
        index=False,
        sep="\t",
    )


main(snakemake, snakemake.config)  # type: ignore
