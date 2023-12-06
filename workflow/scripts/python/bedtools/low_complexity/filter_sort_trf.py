from pathlib import Path
from typing import Any
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins: list[Path] = smk.input
    o = smk.output[0]

    def go(x: cfg.StratInputs_[cfg.AnyBedT]) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.low_complexity.simreps

    sconf.with_build_data_and_bed_io(
        ws["ref_key"],
        ws["build_key"],
        ins,
        [],
        go,
        lambda i, o, bd, b: bd.read_filter_sort_hap_bed(b, i, o),
        lambda i, o, bd, b: bd.read_filter_sort_dip_bed(b, i, o),
        lambda i, o, hap, bd, b: bd.read_filter_sort_hap_bed(b, i, o, hap),
        lambda i, o, bd, b: bd.read_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_filter_sort_hap_bed(b, i, o),
    )


main(snakemake, snakemake.config)  # type: ignore
