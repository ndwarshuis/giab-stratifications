from typing import Any
import common.config as cfg
from common.functional import fmap_maybe


def main(smk: Any) -> None:
    def go(
        x: cfg.BuildData_[
            cfg.RefKeyT,
            cfg.BuildKeyT,
            cfg.RefSourceT,
            cfg.AnyBedT,
            cfg.AnyBedT_,
            cfg.AnySrcT,
            cfg.IncludeT,
        ]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return fmap_maybe(lambda y: y.bench_bed, x.build.bench)

    cfg.filter_sort_bed_main(go, smk)


main(snakemake)  # type: ignore
