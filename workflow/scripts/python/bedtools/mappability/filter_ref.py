from typing import Any
import pandas as pd
import common.config as cfg
import subprocess as sp
from pathlib import Path

# This is the one exception to the pattern wherein we filter all input data to
# the chromosomes we actually care about. In the case of mappability, we want to
# run GEM against all random/unplaced contigs, since these are hard to map to
# begin with (otherwise they wouldn't be unplaced). Here, we need to filter
# the fasta to include all the primary contigs we care about, all random contigs
# pertinent to the chromosomes we want, and all unknown contigs. Later, we will
# filter out all these extra contigs.


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])
    prefix = sconf.refkey_to_final_chr_prefix(rk)
    cis = sconf.buildkey_to_chr_indices(rk, bk)
    idx = pd.read_table(
        Path(smk.input["idx"]),
        header=None,
        dtype={0: str},
        usecols=[0],
    )
    chrs = [
        c
        for c in idx[0].tolist()
        if any([i.matches(prefix, c) for i in cis])
        or cfg.chr_is_unknown(prefix, c)
        or cfg.chr_is_mito(prefix, c)
    ]
    with open(Path(smk.output[0]), "w") as f:
        # ASSUME sort order doesn't matter
        # ASSUME ref is bgziped or unzip'd
        # NOTE the output is not bgzip'd because GEM doesn't like it
        p = sp.Popen(["samtools", "faidx", Path(smk.input["fa"]), *chrs], stdout=f)
        p.wait()


main(snakemake, snakemake.config)  # type: ignore
