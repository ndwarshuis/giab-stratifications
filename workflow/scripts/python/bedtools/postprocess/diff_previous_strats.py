from pathlib import Path
from typing import Any
import common.config as cfg
from common.stratdiff.lib.diff import compare_all
from common.io import setup_logging
from common.functional import DesignError


log = setup_logging(snakemake.log[0])  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    rkf = cfg.RefKeyFullS(ws["ref_final_key"])
    rk = cfg.parse_final_refkey(rkf)[0]
    bk = cfg.BuildKey(rkf)
    bd = sconf.to_build_data(rk, bk)
    comparison = bd.build.comparison

    fm = sconf.with_build_data_final(
        rkf,
        bk,
        lambda bd: bd.ref_chr_conversion.final_mapper,
        lambda bd: bd.ref_chr_conversion.final_mapper,
        lambda hap, bd: hap.from_either(*bd.ref_chr_conversion).final_mapper,
    )
    chr_names = list(fm.values())

    if comparison is None:
        raise DesignError("comparison should not be none")

    outdir = Path(smk.output[0]).parent

    logged = compare_all(
        Path(smk.input["new_list"]).parent,
        smk.input["old"],
        outdir,
        comparison.path_mapper,
        comparison.replacements,
        chr_names,
        comparison.ignore_generated,
        comparison.ignore_other,
    )
    for x in logged:
        log.info(x)


main(snakemake, snakemake.config)  # type: ignore
