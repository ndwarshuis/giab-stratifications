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
    hap = cfg.to_haplotype(ws["hap"])
    bd = sconf.to_build_data(ws["ref_key"], ws["build_key"])

    # TODO this may be two files in the diploid2 case
    idx = pd.read_table(
        Path(smk.input["idx"]),
        header=None,
        dtype={0: str},
        usecols=[0],
    )

    def run_samtools(main_chrs: list[str]) -> None:
        chrs = [
            *filter(
                lambda c: c in main_chrs
                or any(re.match(p, c) for p in bd.mappability_patterns),
                idx[0].tolist(),
            )
        ]
        with open(Path(smk.output[0]), "w") as f:
            # ASSUME sort order doesn't matter and ref is bgzip'd or unzip'd
            # NOTE the output is not bgzip'd because GEM doesn't like it
            p = sp.Popen(["samtools", "faidx", smk.input["fa"], *chrs], stdout=f)
            p.wait()

    if isinstance(bd, cfg.HaploidBuildData) and hap is None:
        main_chrs = bd.ref.chr_pattern.to_names(bd.chr_indices)
        run_samtools(main_chrs)
    elif isinstance(bd, cfg.Diploid1BuildData) and hap is None:
        main_chrs = bd.ref.chr_pattern.to_names(bd.chr_indices)
        run_samtools(main_chrs)
    elif isinstance(bd, cfg.Diploid2BuildData) and hap is not None:
        pat = bd.ref.chr_pattern
        main_chrs = hap.from_either(pat.hap1, pat.hap2).to_names(bd.chr_indices)
        run_samtools(main_chrs)
    else:
        assert False


main(snakemake, snakemake.config)  # type: ignore
