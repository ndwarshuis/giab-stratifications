from typing import Any
import common.config as cfg


def main(smk: Any) -> None:
    cfg.filter_sort_bed_main(cfg.bd_to_bench_bed, smk)


main(snakemake)  # type: ignore
