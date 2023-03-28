from typing import Any
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed
from pybedtools import BedTool as bt  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    level = smk.params["level"]
    i = cfg.ChrIndex.from_name(smk.wildcards["chr"])

    assert i in [cfg.ChrIndex.CHRX, cfg.ChrIndex.CHRY], "invalid sex chr"

    fs = sconf.refkey_to_strat(rk).xy.features
    bedfile = fs.x_bed if i is cfg.ChrIndex.CHRX else fs.y_bed

    assert bedfile is not None, "this should not happen"

    level_col = bedfile.level_col
    conv = sconf.buildkey_to_chr_conversion(rk, bk, bedfile.chr_prefix)

    df = read_bed(smk.input["bed"], bedfile, [level_col])
    filtsort = filter_sort_bed(conv, df)
    filtsort[level_col] = filtsort[level_col].str.contains(level)
    gapless = bt(smk.input["gapless"])
    nogaps = (
        bt.from_dataframe(filtsort.drop(columns=[level_col]))
        .intersect(b=gapless, sorted=True)
        .to_dataframe()
    )
    write_bed(smk.output[0], nogaps)


main(snakemake, snakemake.config)  # type: ignore
