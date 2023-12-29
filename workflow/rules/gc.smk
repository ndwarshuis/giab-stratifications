from functools import partial
from common.config import CoreLevel
from more_itertools import unzip
import json

gc = config.to_bed_dirs(CoreLevel.GC)


def seqtk_args(wildcards):
    _frac = int(wildcards["frac"])
    rk, _ = parse_final_refkey(wildcards.ref_final_key)
    bk = wildcards.build_key
    gps = config.to_build_data(rk, bk).build.include.gc

    if _frac in gps.low_fractions:
        switch, frac = ("w", 100 - _frac)
    elif _frac in gps.high_fractions:
        switch, frac = ("", _frac)
    else:
        assert False, "this should not happen"

    return f"-{switch}f 0.{frac}"


# TODO this can take a gzipped fa file
rule find_gc_content:
    input:
        ref=rules.filter_sort_ref.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        gc.inter.postsort.data / "gc{frac}.bed.gz",
    params:
        args=seqtk_args,
    conda:
        "../envs/seqtk.yml"
    wildcard_constraints:
        frac="\d+",
    shell:
        """
        seqtk gc {params.args} -l 100 {input.ref} | \
        cut -f1-3 | \
        slopBed -i stdin -g {input.genome} -b 50 | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c \
        > {output}
        """


use rule find_gc_content as find_gc_content_final with:
    output:
        gc.final("gc{frac}_slop50"),


def range_inputs(wildcards):
    def _expand(p, frac):
        return expand(p, allow_missing=True, frac=frac)

    def expand_inter(frac):
        return _expand(rules.find_gc_content.output, frac)

    def expand_final(frac):
        return _expand(rules.find_gc_content_final.output, frac)

    # TODO hmmmm...sketchy
    rk, _ = parse_final_refkey(wildcards.ref_final_key)
    bk = wildcards.build_key
    gps = config.to_build_data(rk, bk).build.include.gc
    lowest, lower = gps.low_bounds
    highest, higher = gps.high_bounds

    return {
        "low": expand_final(lowest) + expand_inter(lower),
        "high": expand_inter(higher) + expand_final(highest),
    }


checkpoint intersect_gc_ranges:
    input:
        unpack(range_inputs),
        genome=rules.get_genome.output[0],
        gapless=rules.get_gapless.output.auto,
    output:
        gc.inter.postsort.data / "intersect_output.json",
    conda:
        "../envs/bedtools.yml"
    # hack together a format pattern that will be used for output
    params:
        path_pattern=lambda w: expand(
            gc.final("{{}}"),
            ref_final_key=w.ref_final_key,
            build_key=w.build_key,
        )[0],
    script:
        "../scripts/python/bedtools/gc/intersect_ranges.py"


def gc_inputs(ref_final_key, build_key):
    c = checkpoints.intersect_gc_ranges.get(
        ref_final_key=ref_final_key,
        build_key=build_key,
    )
    with c.output[0].open() as f:
        return json.load(f)


def gc_inputs_flat(ref_final_key, build_key):
    print(ref_final_key)
    res = gc_inputs(ref_final_key, build_key)
    return [*res["gc_ranges"], res["widest_extreme"], *res["other_extremes"]]
