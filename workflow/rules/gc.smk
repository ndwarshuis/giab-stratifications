from functools import partial
from more_itertools import unzip

# gc_src_dir = ref_src_dir / "GCcontent"
gc_inter_dir = intermediate_dir / "GCcontent"
gc_final_dir = final_dir / "GCcontent"
gc_log_dir = log_dir / "GCcontent"


GC_BEDS = [15, 20, 25, 30, 55, 60, 65, 70, 75, 80, 85]
GC_LIMIT = 40


def subset_fractions(x):
    ys = (
        [i for i in GC_BEDS if i <= x]
        if x < GC_LIMIT
        else [i for i in GC_BEDS if i >= x]
    )
    return zip(ys, ys[1:])


def expand_frac(f):
    final = f == GC_BEDS[0] or f == GC_BEDS[-1]
    p = rules.find_gc_content_final.output if final else rules.find_gc_content.output
    return expand(p, allow_missing=True, frac=f)


def range_inputs(wildcards):
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
        gc_inter_dir / "l100_gc{frac}.bed.gz",
    params:
        switch=lambda w: "-wf" if int(wildcards["frac"]) < GC_LIMIT else "-f",
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
        gc_final_dir / "l100_gc{frac}.bed.gz",


rule subtract_gc_content:
    input:
        range_inputs,
    output:
        gc_final_dir / "GRCh38_l100_gc{lower_frac}to{upper_frac}_slop50.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        subtractBed -a {input.bigger_bed} -b {input.smaller_bed} | \
        bgzip -c > {output}
        """


rule intersect_gc_ranges:
    input:
        lambda w: (
            fracs := unzip(subset_fraction(w["lower"]) + subset_fractions(w["upper"])),
            expand(
                rules.subtract_gc_content.output,
                zip,
                lower_frac=fracs[0],
                upper_frac=fracs[1],
            ),
        ),
    output:
        gc_final_dir / "GRCh38_l100_gclt{lower}orgt{upper}_slop50.bed",
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
            lower=[25, 30],
            upper=[65, 50],
        ),
