from pathlib import Path
from typing import Any, NamedTuple
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore
from common.bed import filter_sort_bed, write_bed
from common.functional import match1_unsafe, match2_unsafe


class SortInputs(NamedTuple):
    gapless: Path
    genome: Path


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    bed_input: Path = smk.input["bed"]
    gapless_inputs: list[Path] = smk.input["gapless"]
    genome_inputs: list[Path] = smk.input["genome"]
    sort_inputs = [SortInputs(x, y) for x, y in zip(gapless_inputs, genome_inputs)]
    bed_output = smk.output[0]

    level = smk.params["level"]
    i = cfg.ChrIndex.from_name_unsafe(ws["sex_chr"])

    # TODO kinda not DRY (this gapless/genome thing shows up lots of places)
    def read_write_xy(
        s: SortInputs,
        bd: cfg.BuildData_[
            cfg.RefKeyT,
            cfg.BuildKeyT,
            cfg.RefSourceT,
            cfg.AnyBedT,
            cfg.AnyBedT_,
            cfg.AnySrcT,
            cfg.IncludeT,
        ],
        refPat: cfg.HapChrPattern,
    ) -> None:
        bf = bd.refdata.strat_inputs.xy_feature_bed_unsafe(i)
        conv = cfg.HapToHapChrConversion(bf.data.chr_pattern, refPat, bd.chr_indices)
        df = bf.read(bed_input)
        df_sorted = filter_sort_bed(conv.init_mapper, conv.final_mapper, df)
        level_mask = df_sorted[bf.level_col].str.contains(level)
        df_filtered = df_sorted[level_mask].drop(columns=[bf.level_col])
        # TODO put this in its own rule to simplify script?
        df_gapless = (
            bt()
            .from_dataframe(df_filtered)
            .intersect(
                s.gapless,
                sorted=True,
                g=s.genome,
            )
            .to_dataframe()
        )
        write_bed(bed_output, df_gapless)

    def hap_f(bd: cfg.HapBuildData) -> None:
        match1_unsafe(
            sort_inputs,
            lambda s: read_write_xy(s, bd, bd.refdata.ref.chr_pattern),
        )

    def dip1_f(bd: cfg.Dip1BuildData) -> None:
        pat = bd.refdata.ref.chr_pattern.to_hap_pattern(i.xy_to_hap_unsafe)
        match1_unsafe(sort_inputs, lambda s: read_write_xy(s, bd, pat))

    def dip2_f(h: cfg.Haplotype, bd: cfg.Dip2BuildData) -> None:
        pat = bd.refdata.ref.chr_pattern.from_either(h)
        match2_unsafe(
            sort_inputs,
            lambda s1, s2: read_write_xy(h.from_either(s1, s2), bd, pat),
        )

    sconf.with_build_data_final(
        ws["ref_final_key"],
        ws["build_key"],
        hap_f,
        dip1_f,
        dip2_f,
    )


main(snakemake, snakemake.config)  # type: ignore
