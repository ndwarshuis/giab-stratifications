import re
import pandas as pd
from typing import Any, IO
from pathlib import Path
import common.config as cfg
from common.bed import (
    read_bed,
    filter_sort_bed,
    bed_to_stream,
    mergeBed,
    intersectBed,
    complementBed,
    multiIntersectBed,
    bgzip_file,
)
import json


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards

    im, fm = sconf.buildkey_to_ref_mappers(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
    )

    inputs = smk.input["bed"]
    genome = Path(smk.input["genome"][0])
    gapless = Path(smk.input["gapless"])

    def final_path(name: str) -> Path:
        p = Path(str(smk.params.path_pattern).format(name))
        p.parent.mkdir(exist_ok=True, parents=True)
        return p

    def to_single_output(basename: str) -> Path:
        m = re.match(".*(l\\d+_m\\d+_e\\d+).*", basename)
        assert m is not None, f"malformed mappability file name: {basename}"
        return final_path(f"nonunique_{m[1]}")

    def read_sort_bed(p: Path) -> pd.DataFrame:
        # Sort here because we can't assume wig2bed sorts its output. Also,
        # filtering is necessary because the output should have unplaced contigs
        # in it that we don't want.
        return filter_sort_bed(
            im, fm, read_bed(p, {0: str, 1: int, 2: int}, 0, "\t", [])
        )

    def merge_bed(bed: IO[bytes], out: Path) -> None:
        _, o1 = mergeBed(bed, ["-d", "100"])
        _, o2 = intersectBed(o1, gapless, genome)
        bgzip_file(o2, out)

    def merge_single(i: Path, o: Path) -> None:
        bed = read_sort_bed(i)
        with bed_to_stream(bed) as s:
            _, o1 = complementBed(s, genome)
            merge_bed(o1, o)

    all_lowmap = final_path("lowmappabilityall")
    single_lowmap = []

    # If there is only one input, merge this to make one "all_nonunique" bed.
    # Otherwise, merge each individual input and then combine these with
    # multi-intersect to make the "all_nonunique" bed.
    if len(inputs) == 1:
        merge_single(Path(inputs[0]), all_lowmap)
    else:
        # first read/sort/write all single bed files
        single_lowmap = [to_single_output(Path(i).name) for i in inputs]
        for i, o in zip(inputs, single_lowmap):
            merge_single(i, o)
        # once all the single files are on disk, stream them together; this
        # allows us to avoid keeping multiple dataframes in memory at once
        _, mi_out = multiIntersectBed(single_lowmap)
        merge_bed(mi_out, all_lowmap)

    with open(smk.output[0], "w") as f:
        obj = {"all_lowmap": str(all_lowmap), "single_lowmap": single_lowmap}
        json.dump(obj, f)


main(snakemake, snakemake.config)  # type: ignore
