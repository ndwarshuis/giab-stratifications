from typing import Any
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.bed import read_bed, filter_sort_bed, write_bed


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    level = smk.params["level"]
    i = cfg.ChrIndex.from_name(smk.wildcards["sex_chr"])

    assert i in [cfg.ChrIndex.CHRX, cfg.ChrIndex.CHRY], "invalid sex chr"

    fs = sconf.refkey_to_strat(rk).xy.features
    assert fs is not None, "this should not happen"
    bedfile = fs.x_bed if i is cfg.ChrIndex.CHRX else fs.y_bed
    ps = bedfile.params

    level_col = bedfile.level_col
    conv = sconf.buildkey_to_chr_conversion(rk, bk, ps.chr_pattern)

    df = read_bed(smk.input["bed"], ps, [level_col])
    filtsort = filter_sort_bed(conv, df)
    filtsort = filtsort[filtsort[level_col].str.contains(level)].drop(
        columns=[level_col]
    )
    gapless = bt(smk.input["gapless"])
    nogaps = (
        bt.from_dataframe(filtsort)
        .intersect(b=gapless, sorted=True, g=str(smk.input["genome"]))
        .to_dataframe()
    )
    write_bed(smk.output[0], nogaps)


main(snakemake, snakemake.config)  # type: ignore
