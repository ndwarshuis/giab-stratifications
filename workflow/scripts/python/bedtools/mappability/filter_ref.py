import re
from typing import Any
import pandas as pd
import common.config as cfg
import subprocess as sp
from pathlib import Path

# This is the one exception to the pattern wherein we filter all input data to
# the chromosomes we actually care about. In the case of mappability, we want to
# run GEM against all random/unplaced contigs, since these are hard to map to
# begin with (otherwise they wouldn't be unplaced). Here, we need to filter
# the fasta to include all the primary contigs we care about + all the other
# unplaced/random contigs. Later we will filter out all the non-primary contigs.


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    idx = pd.read_table(
        Path(smk.input["idx"]),
        header=None,
        dtype={0: str},
        usecols=[0],
    )

    def run_samtools(
        bd: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.AnySrcT],
        pat: cfg.ChrPattern,
    ) -> None:
        main_chrs = pat.to_names(bd.chr_indices)
        chrs: list[str] = [
            c
            for c in idx[0].tolist()
            if c in main_chrs
            or any(re.match(p, c) is not None for p in bd.refdata.mappability_patterns)
        ]
        with open(Path(smk.output[0]), "w") as f:
            # ASSUME sort order doesn't matter and ref is bgzip'd or unzip'd
            # NOTE the output is not bgzip'd because GEM doesn't like it
            p = sp.Popen(["samtools", "faidx", smk.input["fa"], *chrs], stdout=f)
            p.wait()

    sconf.with_build_data_full(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
        lambda bd: run_samtools(bd, bd.refdata.ref.chr_pattern),
        lambda bd: run_samtools(bd, bd.refdata.ref.chr_pattern),
        lambda hap, bd: run_samtools(bd, bd.refdata.ref.chr_pattern.from_either(hap)),
    )


main(snakemake, snakemake.config)  # type: ignore
