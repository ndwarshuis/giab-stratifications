from pathlib import Path
from typing import Any, NamedTuple
import common.config as cfg


class PostPaths(NamedTuple):
    gapless: Path
    genome: Path
    output: Path


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins_bed: list[Path] = list(map(Path, smk.input["bed"]))
    output_pattern: str = ps["output_pattern"]

    lk = cfg.OtherLevelKey(ws["other_level_key"])
    sk = cfg.OtherStratKey(ws["other_strat_key"])

    def go(
        x: cfg.BuildData_[
            cfg.RefKeyT,
            cfg.BuildKeyT,
            cfg.RefSourceT,
            cfg.AnyBedT,
            cfg.AnyBedT_,
            cfg.IncludeT,
        ]
    ) -> cfg.OtherBedFile[cfg.AnyBedT] | None:
        return x.build.other_strats[lk][sk]

    # def postprocess(
    #     gapless_path: Path,
    #     genome_path: Path,
    #     bf: cfg.OtherBedFile[cfg.AnyBedT],
    #     df: pd.DataFrame,
    # ) -> pd.DataFrame:
    #     if bf.remove_gaps:
    #         return cast(
    #             pd.DataFrame,
    #             bt()
    #             .from_dataframe(df)
    #             .intersect(
    #                 gapless_path,
    #                 sorted=True,
    #                 g=genome_path,
    #             )
    #             .to_dataframe(),
    #         )
    #     else:
    #         return df

    sconf.with_build_data_and_bed_io_(
        ws["ref_key"],
        ws["build_key"],
        ins_bed,
        smk.output[0],
        output_pattern,
        go,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, hap),
    )


main(snakemake, snakemake.config)  # type: ignore
