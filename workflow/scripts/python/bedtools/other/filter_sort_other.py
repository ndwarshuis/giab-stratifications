from typing import Any
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed
from pybedtools import BedTool as bt  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    lk = cfg.OtherLevelKey(smk.wildcards["other_level_key"])
    sk = cfg.OtherStratKey(smk.wildcards["other_strat_key"])
    bed = sconf.otherkey_to_bed(rk, bk, lk, sk)
    ps = bed.params
    conv = sconf.buildkey_to_chr_conversion(rk, bk, ps.chr_pattern)
    df = filter_sort_bed(conv, read_bed(smk.input["bed"][0], ps))
    if bed.remove_gaps:
        df = (
            bt()
            .from_dataframe(df)
            .intersect(
                smk.input["gapless"],
                sorted=True,
                g=smk.input["genome"][0],
            )
            .to_dataframe()
        )
    write_bed(smk.output[0], df)


main(snakemake, snakemake.config)  # type: ignore
