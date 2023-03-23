from functools import partial
from more_itertools import unzip

gc_inter_dir = intermediate_dir / "GCcontent"
gc_final_dir = final_dir / "GCcontent"
gc_log_dir = log_dir / "GCcontent"


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
    assert abs(lf - uf) == 5, "invalid GC fraction combo"

    bed_a, bed_b = (
        (expand_frac(uf), expand_frac(lf))
        if lf < GC_LIMIT
        else (expand_frac(lf), expand_frac(uf))
    )
    return {"bed_a": bed_a, "bed_b": bed_b}


rule find_gc_content:
    input:
        ref=rules.filter_sort_ref.output,
        genome=rules.get_genome.output,
    output:
        gc_inter_dir / "l100_gc{frac,[0-9]+}.bed.gz",
    params:
        switch=lambda w: "-wf" if int(w["frac"]) < GC_LIMIT else "-f",
    conda:
        envs_path("seqtk.yml")
    shell:
        """
        seqtk gc {params.switch} 0.{wildcards.frac} -l 100 {input.ref} | \
        cut -f1-3 | \
        slopBed -i stdin -g {input.genome} -b 50 | \
        mergeBed -i stdin | \
        bgzip -c \
        > {output}
        """


use rule find_gc_content as find_gc_content_final with:
    output:
        gc_final_dir / "GRCh38_l100_gc{frac,[0-9]+}.bed.gz",


rule subtract_gc_content:
    input:
        unpack(subtract_inputs),
    output:
        gc_final_dir
        / "GRCh38_l100_gc{lower_frac,[0-9]+}to{upper_frac,[0-9]+}_slop50.bed.gz",
    conda:
        envs_path("bedtools.yml")
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
    return expand(
        rules.subtract_gc_content.output,
        zip,
        allow_missing=True,
        lower_frac=lower,
        upper_frac=upper,
    )


rule intersect_gc_ranges:
    input:
        range_inputs,
    output:
        gc_final_dir / "GRCh38_l100_gclt{lower,[0-9]+}orgt{upper,[0-9]+}_slop50.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "multiIntersectBed -i {input} | mergeBed -i stdin | bgzip -c > {output}"


# ASSUME this alone will pull in all the individual range beds
rule all_gc:
    input:
        expand(
            rules.intersect_gc_ranges.output,
            zip,
            allow_missing=True,
            lower=[25, 30],
            upper=[65, 50],
        ),
