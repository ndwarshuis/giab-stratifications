from typing import Any
from pathlib import Path
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.bed import filter_sort_bed, read_bed, write_bed
import pandas as pd


def read_genome_bed(p: Path) -> pd.DataFrame:
    df = pd.read_table(
        p,
        header=None,
        names=["chrom", "end"],
        dtype={"chrom": str, "end": int},
    )
    df["start"] = 0
    return df[["chrom", "start", "end"]].copy()


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    inputs = smk.input
    auto_out = Path(smk.output["auto"])
    parY_out = Path(smk.output["parY"])

    # convert genome to bed file (where each region is just the length of
    # one chromosome)
    genome_bed = read_genome_bed(Path(inputs["genome"]))

    # If we have gap input, make the gapless file, otherwise just symlink to the
    # genome bed file (which just means the entire genome is gapless)
    if hasattr(inputs, "gaps"):
        rk = cfg.RefKey(smk.wildcards["ref_key"])
        bk = cfg.BuildKey(smk.wildcards["build_key"])

        bedfile = sconf.stratifications[rk].gap
        assert bedfile is not None, "this should not happen"
        ps = bedfile.params

        conv = sconf.buildkey_to_chr_conversion(rk, bk, ps.chr_pattern)

        gaps_src = Path(inputs["gaps"])

        gaps = read_bed(gaps_src, ps)
        gaps = filter_sort_bed(conv, gaps)
        gaps = bt.from_dataframe(gaps).merge(d=100).to_dataframe()
        gaps_bed = bt().from_dataframe(gaps)
        gaps_with_parY = bt().from_dataframe(genome_bed).subtract(gaps_bed)

        # If we have a parY bed, subtract parY from the gaps bed, otherwise
        # just link them since we have nothing to subtract off
        if hasattr(inputs, "parY"):
            parY_src = Path(inputs["parY"])
            gaps_no_parY = gaps_with_parY.subtract(bt(parY_src))

            write_bed(parY_out, gaps_with_parY.to_dataframe())
            write_bed(auto_out, gaps_no_parY.to_dataframe())
        else:
            write_bed(auto_out, gaps)
            parY_out.symlink_to(auto_out.resolve())
    else:
        write_bed(auto_out, genome_bed)
        parY_out.symlink_to(auto_out.resolve())


main(snakemake, snakemake.config)  # type: ignore
