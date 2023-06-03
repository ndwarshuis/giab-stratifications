from functools import partial
from common.config import CoreLevel
from more_itertools import unzip

gc_dir = CoreLevel.GC
gc_inter_dir = config.intermediate_build_dir / gc_dir.value


def gc_final_path(name):
    return config.build_strat_path(gc_dir, name)


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
    final = f == GC_BEDS[0] or f == GC_BEDS[-1]
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


rule subtract_gc_content:
    input:
        unpack(subtract_inputs),
    output:
        gc_final_path("gc{lower_frac}to{upper_frac}_slop50"),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        lower_frac="\d+",
        upper_frac="\d+",
    shell:
        """
        subtractBed -a {input.bed_a} -b {input.bed_b} | \
        bgzip -c > {output}
        """


def range_inputs(wildcards):
    lower, upper = map(
        list,
        unzip(
            subset_fractions(wildcards["lower"]) + subset_fractions(wildcards["upper"])
        ),
    )
    return (
        expand(
            rules.subtract_gc_content.output,
            zip,
            allow_missing=True,
            lower_frac=lower,
            upper_frac=upper,
        )
        + expand_frac(lower)
        + expand_frac(upper)
    )


rule intersect_gc_ranges:
    input:
        range_inputs,
    output:
        gc_final_path("gclt{lower,[0-9]+}orgt{upper,[0-9]+}_slop50"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


# ASSUME this alone will pull in all the individual range beds
rule all_gc:
    input:
        **{
            key: expand(
                rules.intersect_gc_ranges.output,
                allow_missing=True,
                lower=lwr,
                upper=upr,
            )
            for key, lwr, upr in [("wide", 25, 65), ("narrow", 30, 55)]
        },
        # TODO not DRY
        middle=expand(
            rules.subtract_gc_content.output,
            lower_frac=30,
            upper_frac=55,
            allow_missing=True,
        ),
