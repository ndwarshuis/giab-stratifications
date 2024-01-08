from typing import Any
import common.config as cfg


def main(smk: Any) -> None:
    cfg.filter_sort_bed_main(lambda bd: cfg.bd_to_si(cfg.si_to_superdups, bd), smk)


main(snakemake)  # type: ignore
