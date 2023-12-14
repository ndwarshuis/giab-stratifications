from typing import Any
import common.config as cfg


def main(smk: Any) -> None:
    def go(
        x: cfg.BuildData_[
            cfg.RefKeyT,
            cfg.BuildKeyT,
            cfg.RefSourceT,
            cfg.AnyBedT,
            cfg.AnyBedT_,
            cfg.IncludeT,
        ]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return x.refdata.strat_inputs.low_complexity.rmsk

    cfg.filter_sort_bed_main(go, smk)


main(snakemake)  # type: ignore
