import json
import pandas as pd
from pathlib import Path
from typing import Any, NamedTuple, cast
import common.config as cfg
from pybedtools import BedTool as bt  # type: ignore


class PostPaths(NamedTuple):
    gapless: Path
    genome: Path
    output: Path


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins_bed: list[Path] = list(map(Path, smk.input["bed"]))
    ins_gapless: list[Path] = list(map(Path, smk.input["gapless"]))
    ins_genome: list[Path] = list(map(Path, smk.input["genome"]))
    output_list = smk.output[0]
    bed_outputs = list(map(Path, ps["bed_outputs"]))

    lk = cfg.OtherLevelKey(ws["other_level_key"])
    sk = cfg.OtherStratKey(ws["other_strat_key"])

    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.IncludeT]
    ) -> cfg.OtherBedFile[cfg.AnyBedT] | None:
        return x.build.other_strats[lk][sk]

    def postprocess(
        gapless_path: Path,
        genome_path: Path,
        bf: cfg.OtherBedFile[cfg.AnyBedT],
        df: pd.DataFrame,
    ) -> pd.DataFrame:
        if bf.remove_gaps:
            return cast(
                pd.DataFrame,
                bt()
                .from_dataframe(df)
                .intersect(
                    gapless_path,
                    sorted=True,
                    g=genome_path,
                )
                .to_dataframe(),
            )
        else:
            return df

    # TODO make types here better once mypy gets higher order type var support
    sconf.with_build_data_and_bed_io_(
        ws["ref_final_key"],
        ws["build_key"],
        ins_bed,
        # ASSUME the cardinality of these is the the same; if not this function
        # will create a black hole and consume the entire galaxy...kinda bad but
        # at least you will know ;)
        [PostPaths(x, y, z) for x, y, z in zip(ins_gapless, ins_genome, bed_outputs)],
        go,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(
            b,
            i,
            o.output,
            lambda df: postprocess(
                o.gapless, o.genome, cast(cfg.OtherBedFile[cfg.HapBedSrc], b), df
            ),
        ),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(
            b,
            i,
            o.output,
            lambda df: postprocess(
                o.gapless, o.genome, cast(cfg.OtherBedFile[cfg.DipBedSrc], b), df
            ),
        ),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(
            b,
            i,
            (o[0].output, o[1].output),
            lambda df: postprocess(
                o[0].gapless, o[0].genome, cast(cfg.OtherBedFile[cfg.DipBedSrc], b), df
            ),
            lambda df: postprocess(
                o[1].gapless, o[1].genome, cast(cfg.OtherBedFile[cfg.DipBedSrc], b), df
            ),
        ),
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(
            b,
            i,
            o.output,
            lambda df: postprocess(
                o.gapless, o.genome, cast(cfg.OtherBedFile[cfg.DipBedSrc], b), df
            ),
        ),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_hap_bed(
            b,
            i,
            o.output,
            hap,
            lambda df: postprocess(
                o.gapless, o.genome, cast(cfg.OtherBedFile[cfg.DipBedSrc], b), df
            ),
        ),
    )

    with open(output_list, "w") as f:
        json.dump(bed_outputs, f)


main(snakemake, snakemake.config)  # type: ignore
