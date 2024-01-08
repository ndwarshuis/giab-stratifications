import pandas as pd
from typing import Any
import common.config as cfg


def filter_ct(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df[3].str.startswith("ct_")]


def main(smk: Any) -> None:
    cfg.filter_sort_bed_main(
        lambda bd: cfg.bd_to_si(cfg.si_to_satellites, bd),
        smk,
        filter_ct,
    )


main(snakemake)  # type: ignore
