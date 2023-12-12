from pathlib import Path
import json
from typing import Any
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins: list[Path] = smk.input
    output_list = smk.output[0]
    output_pattern: str = ps["output_pattern"]

    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.IncludeT]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.strat_inputs.low_complexity.simreps

    sconf.with_build_data_and_bed_io_(
        ws["ref_final_key"],
        ws["build_key"],
        ins,
        bed_outputs,
        go,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, hap),
    )

    with open(output_list, "w") as f:
        json.dump(bed_outputs, f)


main(snakemake, snakemake.config)  # type: ignore
