from typing import Any
import common.config as cfg


def main(smk: Any) -> None:
    ws: dict[str, str] = smk.wildcards
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

    cfg.filter_sort_bed_main(go, smk)


main(snakemake)  # type: ignore
