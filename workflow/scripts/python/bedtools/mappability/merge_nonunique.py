import re
from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import read_bed, filter_sort_bed, write_bed
from snakemake.io import expand  # type: ignore
from pybedtools import BedTool as bt  # type: ignore


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards["ref_key"])
    bk = cfg.BuildKey(smk.wildcards["build_key"])

    prefix = sconf.refkey_to_final_chr_prefix(rk)
    conv = cfg.fullset_conv(prefix)

    inputs = smk.input["bed"]
    output = Path(smk.output[0])
    genome = Path(smk.input["genome"][0])
    gapless = bt(smk.input["gapless"])
    mapdir = smk.params["map_dir"]

    def read_sort_bed(p: Path) -> bt:
        df = read_bed(p)
        # sort here because we can't assume wig2bed sorts its output
        return bt().from_dataframe(filter_sort_bed(conv, df))

    def merge_bed(bed: bt, out: Path) -> None:
        df = bed.merge(d=100).intersect(b=gapless, sorted=True, g=genome).to_dataframe()
        write_bed(out, df)

    def merge_single(bed: bt, out: Path) -> None:
        comp = bed.complement(g=str(genome))
        merge_bed(comp, out)

    # dirty hack since we can't use functions in snakemake output files
    def to_single_output(basename: str) -> Path:
        m = re.match(".*(l\\d+_m\\d+_e\\d+).*", basename)
        assert m is not None, f"malformed mappability file name: {basename}"
        return Path(
            expand(
                sconf.build_strat_path(mapdir, f"nonunique_{m[1]}"),
                ref_key=rk,
                build_key=bk,
            )[0]
        )

    # If there is only one input, merge this to make one "all_nonunique" bed.
    # Otherwise, merge each individual input and then combine these with
    # multi-intersect to make the "all_nonunique" bed.
    if len(inputs) == 1:
        bed = read_sort_bed(Path(inputs[0]))
        merge_single(bed, output)
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
        merge_bed(mint, output)


main(snakemake, snakemake.config)  # type: ignore
