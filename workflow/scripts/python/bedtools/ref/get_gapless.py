from typing import Any
from pathlib import Path
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.bed import write_bed
from common.functional import match1_unsafe, match2_unsafe
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
    ws: dict[str, str] = smk.wildcards

    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.AnySrcT]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.refdata.strat_inputs.gap

    # convert genome to bed file (where each region is just the length of
    # one chromosome)
    genome_bed = read_genome_bed(Path(inputs["genome"]))

    # If we have gap input, make the gapless file, otherwise just symlink to the
    # genome bed file (which just means the entire genome is gapless)
    if not hasattr(inputs, "gaps"):
        write_bed(auto_out, genome_bed)
        parY_out.symlink_to(auto_out.resolve())
    else:
        gap_inputs: list[Path] = inputs["gaps"]

        gaps_df = sconf.with_build_data_and_bed_hap(
            cfg.RefKeyFullS(ws["ref_final_key"]),
            cfg.BuildKey(ws["build_key"]),
            go,
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: bd.read_filter_sort_hap_bed(bf, i)
            ),
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: bd.read_filter_sort_dip1_bed(bf, i)
            ),
            lambda hap, bd, bf: match1_unsafe(
                gap_inputs,
                lambda i: hap.from_either(*bd.read_filter_sort_dip1_bed(bf, i)),
            ),
            lambda bd, bf: match2_unsafe(
                gap_inputs,
                lambda i0, i1: bd.read_filter_sort_dip2_bed(bf, (i0, i1)),
            ),
            lambda hap, bd, bf: match2_unsafe(
                gap_inputs,
                lambda i0, i1: bd.read_filter_sort_dip2_bed(
                    bf,
                    *hap.from_either((i0, cfg.Haplotype.HAP1), (i1, cfg.Haplotype.HAP2))
                ),
            ),
        )

        gaps: pd.DataFrame = bt.from_dataframe(gaps_df).merge(d=100).to_dataframe()
        gaps_bed = bt().from_dataframe(gaps)
        gaps_with_parY = bt().from_dataframe(genome_bed).subtract(gaps_bed)

        # If we have a parY bed, subtract parY from the gaps bed, otherwise
        # just link them since we have nothing to subtract off
        # TODO add logic to ignore parY when on maternal hap (not necessary but
        # makes things cleaner)
        if hasattr(inputs, "parY"):
            parY_src = Path(inputs["parY"])
            gaps_no_parY = gaps_with_parY.subtract(bt(parY_src))

            write_bed(parY_out, gaps_with_parY.to_dataframe())
            write_bed(auto_out, gaps_no_parY.to_dataframe())
        else:
            write_bed(auto_out, gaps)
            parY_out.symlink_to(auto_out.resolve())


main(snakemake, snakemake.config)  # type: ignore
