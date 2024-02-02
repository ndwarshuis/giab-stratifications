from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import (
    write_bed,
    bed_to_stream,
    read_bed_default,
    complementBed,
    mergeBed,
    subtractBed,
)
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
    ws: dict[str, Any] = smk.wildcards

    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.AnySrcT]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.refdata.strat_inputs.gap

    # convert genome to bed file (where each region is just the length of
    # one chromosome)
    genome_path = Path(inputs["genome"])
    genome_bed = read_genome_bed(genome_path)

    # If we have gap input, make the gapless file, otherwise just symlink to the
    # genome bed file (which just means the entire genome is gapless)
    if not hasattr(inputs, "gaps"):
        write_bed(auto_out, genome_bed)
        parY_out.symlink_to(auto_out.resolve())
    else:
        gap_inputs: list[Path] = inputs["gaps"]

        gaps_df = sconf.with_build_data_and_bed_full(
            cfg.wc_to_reffinalkey(ws),
            cfg.wc_to_buildkey(ws),
            go,
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: cfg.read_filter_sort_hap_bed(bd, bf, i)
            ),
            lambda bd, bf: match1_unsafe(
                gap_inputs, lambda i: cfg.read_filter_sort_dip1to1_bed(bd, bf, i)
            ),
            lambda hap, bd, bf: match1_unsafe(
                gap_inputs,
                lambda i: hap.from_either(*cfg.read_filter_sort_dip1to2_bed(bd, bf, i)),
            ),
            lambda bd, bf: match2_unsafe(
                gap_inputs,
                lambda i0, i1: cfg.read_filter_sort_dip2to1_bed(bd, bf, (i0, i1)),
            ),
            lambda hap, bd, bf: match2_unsafe(
                gap_inputs,
                lambda i0, i1: cfg.read_filter_sort_dip2to2_bed(
                    bd,
                    bf,
                    *hap.from_either((i0, cfg.Haplotype.HAP1), (i1, cfg.Haplotype.HAP2))
                ),
            ),
        )

        with bed_to_stream(gaps_df) as s:
            _, o = mergeBed(s, ["-d", "100"])
            gaps_bed = read_bed_default(o)

        with bed_to_stream(gaps_bed) as s:
            _, o = complementBed(s, genome_path)
            gaps_with_parY = read_bed_default(o)

        # gaps: pd.DataFrame = bt.from_dataframe(gaps_df).merge(d=100).to_dataframe()
        # gaps_bed = bt().from_dataframe(gaps)
        # gaps_with_parY = bt().from_dataframe(genome_bed).subtract(gaps_bed)

        # If we have a parY bed, subtract parY from the gaps bed, otherwise
        # just link them since we have nothing to subtract off
        # TODO add logic to ignore parY when on maternal hap (not necessary but
        # makes things cleaner)
        if hasattr(inputs, "parY"):
            parY_src = Path(inputs["parY"])
            with bed_to_stream(gaps_with_parY) as s:
                _, o = subtractBed(s, parY_src)
                gaps_no_parY = read_bed_default(o)
            write_bed(parY_out, gaps_with_parY)
            write_bed(auto_out, gaps_no_parY)

        else:
            write_bed(auto_out, gaps_bed)
            parY_out.symlink_to(auto_out.resolve())


main(snakemake, snakemake.config)  # type: ignore
