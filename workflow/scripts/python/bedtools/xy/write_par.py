from typing import Any
from Bio import bgzf  # type: ignore
import common.config as cfg
from common.functional import not_none_unsafe, none_unsafe


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    i = cfg.ChrIndex.from_name(smk.wildcards["sex_chr"])

    # TODO this pattern is DRY?
    rk, hap = cfg.parse_final_refkey(
        smk.wildcards["ref_final_key"],
    )
    # NOTE tuple thing is to appease mypy...be nice to mypy...respect mypy
    cxy, pat = sconf.with_ref_data(
        rk,
        lambda rd: (
            rd.strat_inputs.xy,
            none_unsafe(hap, rd.ref.chr_pattern),
        ),
        lambda rd: (
            rd.strat_inputs.xy,
            none_unsafe(hap, rd.ref.chr_pattern.to_hap_pattern(i.xy_to_hap_unsafe)),
        ),
        lambda rd: (
            rd.strat_inputs.xy,
            not_none_unsafe(hap, lambda h: rd.ref.chr_pattern.from_either(h)),
        ),
    )

    par_fun = i.choose_xy_unsafe(cxy.fmt_x_par_unsafe, cxy.fmt_y_par_unsafe)

    # TODO not dry?
    with bgzf.BgzfWriter(smk.output[0], "w") as f:
        f.write(par_fun(pat))


main(snakemake, snakemake.config)  # type: ignore
