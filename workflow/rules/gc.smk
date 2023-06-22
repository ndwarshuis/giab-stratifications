from functools import partial
from common.config import CoreLevel
from more_itertools import unzip
import json

gc_dir = CoreLevel.GC
gc_inter_dir = config.intermediate_build_dir / gc_dir.value


def gc_final_path(name):
    return config.build_strat_path(gc_dir, name)


def seqtk_args(wildcards):
    _frac = int(wildcards["frac"])
    rk = wildcards.ref_key
    bk = wildcards.build_key
    gps = config.buildkey_to_include(rk, bk).gc

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
        gc_inter_dir / "gc{frac}.bed.gz",
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
        gc_final_path("gc{frac}_slop50"),


def range_inputs(wildcards):
    def _expand(p, frac):
        return expand(p, allow_missing=True, frac=frac)

    def expand_inter(frac):
        return _expand(rules.find_gc_content.output, frac)

    def expand_final(frac):
        return _expand(rules.find_gc_content_final.output, frac)

    rk = wildcards.ref_key
    bk = wildcards.build_key
    gps = config.buildkey_to_include(rk, bk).gc
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
        gc_inter_dir / "intersect_output.json",
    conda:
        "../envs/bedtools.yml"
    # hack together a format pattern that will be used for output
    params:
        path_pattern=lambda w: expand(
            gc_final_path("{{}}"),
            ref_key=w.ref_key,
            build_key=w.build_key,
        )[0],
    script:
        "../scripts/python/bedtools/gc/intersect_ranges.py"


def gc_inputs(wildcards):
    with checkpoints.intersect_gc_ranges.get(**wildcards).output[0].open() as f:
        return json.load(f)


def gc_inputs_flat(wildcards):
    res = gc_inputs(wildcards)
    return [*res["gc_ranges"], res["widest_extreme"], *res["other_extremes"]]
