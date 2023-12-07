import pandas as pd
import json
from pathlib import Path
from typing import Any
import common.config as cfg


def filter_ct(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df[3].str.startswith("ct_")]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins: list[Path] = smk.input
    output_list: Path = smk.output
    bed_outputs = list(map(Path, ps["bed_outputs"]))

    def go(x: cfg.StratInputs_[cfg.AnyBedT]) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.low_complexity.simreps

    sconf.with_build_data_and_bed_io(
        ws["ref_final_key"],
        ws["build_key"],
        ins,
        bed_outputs,
        go,
        lambda i, o, bd, b: bd.read_filter_sort_hap_bed(b, i, o, filter_ct),
        lambda i, o, bd, b: bd.read_filter_sort_dip_bed(b, i, o, filter_ct),
        lambda i, o, bd, b: bd.read_filter_sort_dip_bed(b, i, o, filter_ct),
        lambda i, o, bd, b: bd.read_filter_sort_hap_bed(b, i, o, filter_ct),
        lambda i, o, hap, bd, b: bd.read_filter_sort_hap_bed(b, i, o, hap, filter_ct),
    )

    with open(output_list, "w") as f:
        json.dump(bed_outputs, f)


main(snakemake, snakemake.config)  # type: ignore
