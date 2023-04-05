segdup_src_dir = config.ref_src_dir / "SegmentalDuplications"
segdup_inter_dir = config.intermediate_build_dir / "SegmentalDuplications"


def segdup_final_path(name):
    return config.build_strat_path("SegmentalDuplications", name)


use rule download_ref as download_superdups with:
    output:
        segdup_src_dir / "selfChain.txt.gz",
    params:
        src=lambda w: config.refkey_to_superdups_src(w.ref_key),


rule filter_sort_superdups:
    input:
        rules.download_superdups.output,
    output:
        segdup_inter_dir / "filter_sorted.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/segdups/filter_sort_superdups.py"


rule merge_superdups:
    input:
        bed=rules.filter_sort_superdups.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        segdup_final_path("segdups"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} -d 100 | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_long_superdups:
    input:
        rules.merge_superdups.output,
    output:
        segdup_final_path("segdups_gt10kb"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input} | \
        awk '($3-$2 > 10000)' | \
        bgzip -c > {output}
        """


rule notin_superdups:
    input:
        bed=rules.merge_superdups.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        segdup_final_path("notinsegdups"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


use rule notin_superdups as notin_long_superdups with:
    input:
        bed=rules.filter_long_superdups.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        segdup_final_path("notinsegdups_gt10kb"),


rule all_segdups:
    input:
        rules.notin_superdups.output,
        rules.notin_long_superdups.output,
