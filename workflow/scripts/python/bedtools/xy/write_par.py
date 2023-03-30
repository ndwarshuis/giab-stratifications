from typing import Any
from Bio import bgzf  # type: ignore
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    i = cfg.ChrIndex.from_name(smk.wildcards["chr"])
    k = cfg.RefKey(smk.wildcards["ref_key"])
    prefix = sconf.refkey_to_final_chr_prefix(k)
    cxy = sconf.stratifications[k].xy

    assert i in [cfg.ChrIndex.CHRX, cfg.ChrIndex.CHRY], "invalid sex chr"

    par_fun = cxy.fmt_x_par if i == cfg.ChrIndex.CHRX else cxy.fmt_y_par

    with bgzf.BgzfWriter(smk.output[0], "w") as f:
        f.write(par_fun(prefix))


main(snakemake, snakemake.config)  # type: ignore
