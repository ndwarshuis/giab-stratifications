from functools import partial
from common.config import CoreLevel
from more_itertools import unzip

gc_dir = CoreLevel.GC
gc_inter_dir = config.intermediate_build_dir / gc_dir.value


def gc_final_path(name):
    return config.build_strat_path(gc_dir, name)


GC_LOW = [
    (15, None),
    (20, None),
    (25, True),
    (30, True),
]
GC_HIGH = [
    (55, True),
    (60, None),
    (65, True),
    (70, None),
    (75, None),
    (80, None),
    (85, None),
]

GC_BEDS = [15, 20, 25, 30, 55, 60, 65, 70, 75, 80, 85]
GC_LIMIT = 40


def subset_fractions(x):
    x_ = int(x)
    ys = (
        [i for i in GC_BEDS if i <= x_]
        if x_ < GC_LIMIT
        else [i for i in GC_BEDS if i >= x_]
    )
    return list(zip(ys, ys[1:]))


def expand_frac(f):
    # final = f == GC_BEDS[0] or f == GC_BEDS[-1]
    final = f == GC_LOW[0][0] or f == GC_HIGH[-1][0]
    p = rules.find_gc_content_final.output if final else rules.find_gc_content.output
    return expand(p, allow_missing=True, frac=f)


def subtract_inputs(wildcards):
    lf = int(wildcards["lower_frac"])
    uf = int(wildcards["upper_frac"])

    # hacky sanity check
    # assert abs(lf - uf) == 5, "invalid GC fraction combo"

    bed_a, bed_b = (
        (expand_frac(uf), expand_frac(lf))
        if lf < GC_LIMIT
        else (expand_frac(lf), expand_frac(uf))
    )
    return {"bed_a": bed_a, "bed_b": bed_b}


def seqtk_args(wildcards):
    _frac = int(wildcards["frac"])
    switch, frac = ("w", 100 - _frac) if _frac < GC_LIMIT else ("", _frac)
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


# rule subtract_gc_content:
#     input:
#         unpack(subtract_inputs),
#     output:
#         gc_final_path("gc{lower_frac}to{upper_frac}_slop50"),
#     conda:
#         "../envs/bedtools.yml"
#     wildcard_constraints:
#         lower_frac="\d+",
#         upper_frac="\d+",
#     shell:
#         """
#         subtractBed -a {input.bed_a} -b {input.bed_b} | \
#         bgzip -c > {output}
#         """


# def range_inputs(wildcards):
#     lower, upper = map(
#         list,
#         unzip(
#             subset_fractions(wildcards["lower"]) + subset_fractions(wildcards["upper"])
#         ),
#     )
#     return (
#         expand(
#             rules.subtract_gc_content.output,
#             zip,
#             allow_missing=True,
#             lower_frac=lower,
#             upper_frac=upper,
#         )
#         + expand_frac(lower)
#         + expand_frac(upper)
#     )


# rule intersect_gc_ranges:
#     input:
#         range_inputs,
#     output:
#         gc_final_path("gclt{lower,[0-9]+}orgt{upper,[0-9]+}_slop50"),
#     conda:
#         "../envs/bedtools.yml"
#     shell:
#         """
#         multiIntersectBed -i {input} | \
#         mergeBed -i stdin | \
#         bgzip -c > {output}
#         """


# rule subtract_middle:
#     input:
#         bed_high=lambda w: expand_frac(w.upper),
#         bed_low=lambda w: expand_frac(w.lower),
#         genome=rules.get_genome.output,
#     output:
#         gc_final_path("gclt{lower,[0-9]+}orgt{upper,[0-9]+}_slop50"),
#     shell:
#         """
#         complementBed -i {input.bed_high} -g {input.genome}  | \
#         subtractBed -a stdin -b {input.bed_low} \
#         > {output}
#         """


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


# NOTE: this makes a dummy file to give snakemake a build target; the logic of
# combining these files is complex enough that I would rather have it in a
# python script
rule intersect_gc_ranges:
    input:
        unpack(range_inputs),
        genome=rules.get_genome.output[0],
        gapless=rules.get_gapless.output.auto,
    # randomly need this output to use in union strats
    output:
        gc_inter_dir / "gc_wider_range.bed.gz",
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


# ASSUME this alone will pull in all the individual range beds
# rule all_gc:
#     input:
# **{
#     key: expand(
#         rules.intersect_gc_ranges.output,
#         allow_missing=True,
#         lower=lwr,
#         upper=upr,
#     )
#     for key, lwr, upr in [("wide", 25, 65), ("narrow", 30, 55)]
# },
# # TODO not DRY
# middle=expand(
#     rules.intersect_gc_ranges.output,
#     lower=30,
#     upper=55,
#     allow_missing=True,
# ),
