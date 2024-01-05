from typing import Any
from Bio import bgzf  # type: ignore
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, Any] = smk.wildcards
    i = cfg.ChrIndex.from_name(cfg.wc_lookup(ws, "sex_chr"))

    # TODO this pattern is DRY?
    rfk = cfg.wc_to_reffinalkey(ws)
    rk, hap = cfg.parse_full_refkey(rfk)
    pat = sconf.refkey_to_xy_ref_chr_pattern(rfk, i)
    cxy = sconf.to_ref_data(rk).strat_inputs.xy

    par_fun = cfg.choose_xy_unsafe(i, cxy.fmt_x_par_unsafe, cxy.fmt_y_par_unsafe)

    # TODO not dry?
    with bgzf.BgzfWriter(smk.output[0], "w") as f:
        f.write(par_fun(pat))


main(snakemake, snakemake.config)  # type: ignore
