from typing import Any
import common.config as cfg


def main(smk: Any) -> None:
    ws: dict[str, Any] = smk.wildcards
    lk = cfg.OtherLevelKey(cfg.wc_lookup(ws, "other_level_key"))
    sk = cfg.OtherStratKey(cfg.wc_lookup(ws, "other_strat_key"))
    cfg.filter_sort_bed_main(lambda bd: cfg.bd_to_other(lk, sk, bd), smk)


main(snakemake)  # type: ignore
