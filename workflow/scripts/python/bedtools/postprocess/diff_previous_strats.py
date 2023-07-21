from pathlib import Path
from typing import Any
import common.config as cfg
from common.stratdiff.lib.diff import compare_all
from common.io import setup_logging


log = setup_logging(snakemake.log[0])  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards.ref_key)
    bk = cfg.BuildKey(smk.wildcards.build_key)
    comparison = sconf.buildkey_to_build(rk, bk).comparison
    pattern = sconf.refkey_to_final_chr_pattern(rk)
    ixs = sconf.buildkey_to_chr_indices(rk, bk)

    assert comparison is not None, "this should not happen"

    outdir = Path(smk.output[0]).parent

    logged = compare_all(
        Path(smk.input["new_list"]).parent,
        smk.input["old"],
        outdir,
        comparison.path_mapper,
        comparison.replacements,
        [i.chr_name_full(pattern) for i in ixs],
        comparison.ignore_generated,
        comparison.ignore_other,
    )
    for x in logged:
        log.info(x)


main(snakemake, snakemake.config)  # type: ignore
