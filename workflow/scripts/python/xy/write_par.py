from typing import Any
import common.config as cfg
from Bio import bgzf  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    i = cfg.ChrIndex.from_name(smk.wildcards["chr"])
    k = cfg.RefKey(smk.wildcards["ref_key"])
    prefix = sconf.refkey_to_final_chr_prefix(k)
    if i == cfg.ChrIndex.CHRX:
        par = sconf.stratifications[k].xy.fmt_x_par(prefix)
    elif i == cfg.ChrIndex.CHRY:
        par = sconf.stratifications[k].xy.fmt_y_par(prefix)
    with bgzf.BgzfWriter(smk.output[0], "w") as f:
        f.write(par)


main(snakemake, snakemake.config)  # type: ignore
