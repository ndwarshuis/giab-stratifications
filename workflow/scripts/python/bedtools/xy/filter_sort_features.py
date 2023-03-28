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

    def get_bed(x: cfg.Stratification) -> cfg.BedFile:
        f = x.xy.features
        b = f.x_bed if i is cfg.ChrIndex.CHRX else f.y_bed
        assert b is not None
        return b

    # TODO don't hardcode class column
    df = read_bed(smk.input["bed"], get_bed(sconf.stratifications[rk]), [3])
    filtsort = filter_sort_bed(sconf, get_bed, rk, bk, df)
    filtsort[3] = filtsort[3].str.contains(level)
    gapless = bt(smk.input["gapless"])
    nogaps = (
        bt.from_dataframe(filtsort.drop(columns=[3]))
        .intersect(b=gapless, sorted=True)
        .to_dataframe()
    )
    write_bed(smk.output[0], nogaps)


main(snakemake, snakemake.config)  # type: ignore
