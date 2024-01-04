from typing import Any
import pandas as pd
import common.config as cfg
from common.bed import filter_sort_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, Any] = smk.wildcards

    # TODO not DRY
    im, fm = sconf.with_build_data_final(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
        lambda bd: ((c := bd.ref_chr_conversion).init_mapper, c.final_mapper),
        lambda bd: ((c := bd.ref_chr_conversion).init_mapper, c.final_mapper),
        lambda hap, bd: (
            (c := hap.from_either(*bd.ref_chr_conversion)).init_mapper,
            c.final_mapper,
        ),
    )

    # ASSUME the input for this is a .fa.fai file (columns = chr, length)
    df = pd.read_table(
        smk.input[0],
        header=None,
        dtype={0: str, 1: int},
        usecols=[0, 1],
    )

    filtered = filter_sort_bed(im, fm, df, 2)
    filtered.to_csv(smk.output[0], sep="\t", header=False, index=False)


main(snakemake, snakemake.config)  # type: ignore
