import json
from typing import Any
import common.config as cfg
from common.bed import filter_sort_bed_inner, read_bed, write_bed

VDJ_PAT = "^ID=gene-(IGH|IGK|IGL|TRA|TRB|TRG);"


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    to_map = sconf.buildkey_to_final_chr_mapping(rk, bk)

    with open(smk.input["mapper"], "r") as f:
        from_map = json.load(f)

    df = read_bed(smk.input["bed"], more=[5])
    df = df[df[3].str.match(VDJ_PAT)][[0, 1, 2]].copy()
    df = filter_sort_bed_inner(from_map, to_map, df)

    write_bed(smk.output[0], df)


main(snakemake, snakemake.config)  # type: ignore
