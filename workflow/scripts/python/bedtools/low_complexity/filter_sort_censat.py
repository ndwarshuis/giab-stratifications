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

    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.IncludeT]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.strat_inputs.low_complexity.rmsk
        return x.strat_inputs.low_complexity.simreps

    sconf.with_build_data_and_bed_io_(
        ws["ref_final_key"],
        ws["build_key"],
        ins,
        bed_outputs,
        go,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, filter_ct),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o, filter_ct),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o, filter_ct),
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, filter_ct),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_hap_bed(
            b, i, o, hap, filter_ct
        ),
    )

    with open(output_list, "w") as f:
        json.dump(bed_outputs, f)


main(snakemake, snakemake.config)  # type: ignore
