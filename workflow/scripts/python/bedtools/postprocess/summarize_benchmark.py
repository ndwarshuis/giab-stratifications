from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import read_filter_sort_bed
import pandas as pd

COLS = {
    "Type": str,
    "Subtype": str,
    "Subset": str,
    "Filter": str,
    "METRIC.Recall": float,
    "METRIC.Precision": float,
    "METRIC.Frac_NA": float,
    "Subset.Size": float,
    "Subset.IS_CONF.Size": float,
    "TRUTH.TOTAL": float,
    "TRUTH.TP": float,
    "TRUTH.TP.het": float,
    "TRUTH.TP.homalt": float,
    "TRUTH.FN": float,
    "TRUTH.FN.het": float,
    "TRUTH.FN.homalt": float,
    "TRUTH.TOTAL.het_hom_ratio": float,
    "QUERY.TOTAL": float,
    "QUERY.FP": float,
    "QUERY.FP.het": float,
    "QUERY.FP.homalt": float,
    "QUERY.TOTAL.het_hom_ratio": float,
    "FP.gt": float,
    "FP.al": float,
}


def read_happy_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=[*COLS], dtype=COLS)
    df = df[
        (df["Filter"] == "PASS") & df["Subtype"].isin(["*", "I16_PLUS", "D16_PLUS"])
    ].drop(columns=["Filter"])
    df["Query"] = path.parent.parent.parent.name
    df["Truth/kb"] = df["TRUTH.TOTAL"] / df["Subset.IS_CONF.Size"] * 1000
    df["Query/kb"] = df["QUERY.TOTAL"] / df["Subset.Size"] * 1000
    return df.copy()


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    bench = sconf.stratifications[rk].builds[bk].bench
    assert bench is not None, "this should not happen"
    ps = bench.bench_bed.params
    read_filter_sort_bed(sconf, smk.input[0], smk.output[0], ps, rk, bk)


main(snakemake, snakemake.config)  # type: ignore
