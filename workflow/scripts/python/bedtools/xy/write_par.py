from typing import Any
from Bio import bgzf  # type: ignore
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    i = cfg.ChrIndex.from_name(smk.wildcards["sex_chr"])

    cxy, pat = sconf.with_ref_data_unsafe(
        smk.wildcards["ref_key"],
        lambda rd: (
            rd.strat_inputs.xy,
            rd.ref.chr_pattern,
        ),
        lambda rd: (
            rd.strat_inputs.xy,
            rd.ref.chr_pattern.to_hap_pattern(i.xy_to_hap_unsafe),
        ),
        lambda hap, rd: (
            rd.strat_inputs.xy,
            rd.ref.chr_pattern.from_either(hap),
        ),
    )

    par_fun = i.choose_xy_unsafe(cxy.fmt_x_par_unsafe, cxy.fmt_y_par_unsafe)

    with bgzf.BgzfWriter(smk.output[0], "w") as f:
        f.write(par_fun(pat))


main(snakemake, snakemake.config)  # type: ignore
