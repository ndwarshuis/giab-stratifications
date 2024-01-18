from pathlib import Path
from typing import Any
import common.config as cfg
from common.functional import DesignError
import subprocess as sp


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: cfg.SmkWildcards = smk.wildcards
    i = Path(smk.input[0])
    o = Path(smk.output[0])

    # This looks a bit weird since this script is only supposed to work on dip1
    # references, but the output looks like that for a dip2 reference since it
    # has the haplotype
    bk = cfg.wc_to_buildkey(ws)
    rk, hap = cfg.parse_full_refkey(cfg.wc_to_reffinalkey(ws))
    bd = sconf.to_build_data(rk, bk)

    if not isinstance(bd, cfg.Dip1BuildData):
        raise DesignError("This code was a masterpiece 'til you refactored it")

    if hap is None:
        raise DesignError("Just between us did the typo maim you too?")

    ns = bd.refdata.ref.chr_pattern.to_hap_pattern(hap).to_names(bd.chr_indices)
    sp.run(["samtools", "faidx", str(i), *ns, "-o", str(o)])


main(snakemake, snakemake.config)  # type: ignore
