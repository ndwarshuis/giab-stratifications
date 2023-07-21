from typing import Any
from Bio import bgzf  # type: ignore
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    i = cfg.ChrIndex.from_name(smk.wildcards["sex_chr"])
    k = cfg.RefKey(smk.wildcards["ref_key"])
    pattern = sconf.refkey_to_final_chr_pattern(k)
    cxy = sconf.stratifications[k].xy

    assert i in [cfg.ChrIndex.CHRX, cfg.ChrIndex.CHRY], "invalid sex chr"

    par_fun = cxy.fmt_x_par if i == cfg.ChrIndex.CHRX else cxy.fmt_y_par
    out = par_fun(pattern)
    assert out is not None, "this should not happen"

    with bgzf.BgzfWriter(smk.output[0], "w") as f:
        f.write(out)


main(snakemake, snakemake.config)  # type: ignore
