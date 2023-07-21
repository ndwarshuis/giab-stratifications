import re
from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed
from pybedtools import BedTool as bt  # type: ignore
import json


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])

    pattern = sconf.refkey_to_final_chr_pattern(rk)
    conv = cfg.fullset_conv(pattern)

    inputs = smk.input["bed"]
    genome = Path(smk.input["genome"][0])
    gapless = bt(smk.input["gapless"])

    def final_path(name: str) -> Path:
        p = Path(str(smk.params.path_pattern).format(name))
        p.parent.mkdir(exist_ok=True, parents=True)
        return p

    def to_single_output(basename: str) -> Path:
        m = re.match(".*(l\\d+_m\\d+_e\\d+).*", basename)
        assert m is not None, f"malformed mappability file name: {basename}"
        return final_path(f"nonunique_{m[1]}")

    def read_sort_bed(p: Path) -> bt:
        df = read_bed(p)
        # Sort here because we can't assume wig2bed sorts its output. Also,
        # filtering is necessary because the output should have unplaced contigs
        # in it that we don't want.
        return bt().from_dataframe(filter_sort_bed(conv, df))

    def merge_bed(bed: bt, out: Path) -> None:
        df = bed.merge(d=100).intersect(b=gapless, sorted=True, g=genome).to_dataframe()
        write_bed(out, df)

    def merge_single(bed: bt, out: Path) -> None:
        comp = bed.complement(g=str(genome))
        merge_bed(comp, out)

    all_lowmap = final_path("lowmappabilityall")
    single_lowmap = []

    # If there is only one input, merge this to make one "all_nonunique" bed.
    # Otherwise, merge each individual input and then combine these with
    # multi-intersect to make the "all_nonunique" bed.
    if len(inputs) == 1:
        bed = read_sort_bed(Path(inputs[0]))
        merge_single(bed, all_lowmap)
    else:
        single = [
            (
                to_single_output((p := Path(i)).name),
                read_sort_bed(p),
            )
            for i in inputs
        ]
        for sout, bed in single:
            merge_single(bed, sout)
        mint = bt().multi_intersect(i=[s[0] for s in single])
        merge_bed(mint, all_lowmap)
        single_lowmap = [str(s[0]) for s in single]

    with open(smk.output[0], "w") as f:
        obj = {"all_lowmap": str(all_lowmap), "single_lowmap": single_lowmap}
        json.dump(obj, f)


main(snakemake, snakemake.config)  # type: ignore
