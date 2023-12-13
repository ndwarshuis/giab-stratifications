from pathlib import Path
from typing import Any
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins: list[Path] = smk.input
    output_pattern: str = ps["output_pattern"]

    def go(
        x: cfg.BuildData_[
            cfg.RefKeyT,
            cfg.BuildKeyT,
            cfg.RefSourceT,
            cfg.AnyBedT,
            cfg.AnyBedT_,
            cfg.IncludeT,
        ]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.refdata.strat_inputs.low_complexity.simreps

    sconf.with_build_data_and_bed_io_(
        ws["ref_key"],
        ws["build_key"],
        ins,
        smk.output[0],
        output_pattern,
        go,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, hap),
    )


main(snakemake, snakemake.config)  # type: ignore
