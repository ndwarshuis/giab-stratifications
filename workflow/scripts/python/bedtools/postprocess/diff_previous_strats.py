from pathlib import Path
from typing import Any
import common.config as cfg
from common.stratdiff.lib.diff import compare_all
from common.io import setup_logging
from common.functional import DesignError


log = setup_logging(snakemake.log[0])  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, Any] = smk.wildcards
    rkf = cfg.wc_to_reffinalkey(ws)
    bk = cfg.wc_to_buildkey(ws)
    bd = sconf.to_build_data(cfg.strip_full_refkey(rkf), bk)
    comparison = bd.build.comparison

    fm = sconf.buildkey_to_ref_mappers(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
    )[1]
    chr_names = list(fm.values())

    if comparison is None:
        raise DesignError("comparison should not be None")

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
