segdup_src_dir = ref_src_dir / "SegmentalDuplications"
segdup_inter_dir = intermediate_dir / "SegmentalDuplications"
segdup_final_dir = final_dir / "SegmentalDuplications"


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
        envs_path("bedtools.yml")
    script:
        scripts_path("python/bedtools/segdups/filter_sort_superdups.py")


rule merge_superdups:
    input:
        bed=rules.filter_sort_superdups.output,
        gapless=rules.get_gapless.output.auto,
    output:
        segdup_final_dir / "GRCh38_segdups.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        mergeBed -i {input.bed} -d 100 | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


rule filter_long_superdups:
    input:
        rules.merge_superdups.output,
    output:
        segdup_final_dir / "GRCh38_segdups_gt10kb.bed.gz",
    conda:
        envs_path("bedtools.yml")
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
        segdup_final_dir / "GRCh38_notinsegdups.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


use rule notin_superdups as notin_long_superdups with:
    input:
        bed=rules.filter_long_superdups.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        segdup_final_dir / "GRCh38_notinsegdups_gt10kb.bed.gz",


rule all_segdups:
    input:
        rules.notin_superdups.output,
        rules.notin_long_superdups.output,
