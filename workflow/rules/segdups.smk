from common.config import CoreLevel, si_to_superdups

segdup = config.to_bed_dirs(CoreLevel.SEGDUPS)


use rule download_ref as download_superdups with:
    output:
        segdup.src.data / "superdups.txt.gz",
    log:
        segdup.src.log / "superdups.log",
    params:
        src=lambda w: config.refsrckey_to_bed_src(si_to_superdups, w.ref_src_key),
    localrule: True


checkpoint normalize_superdups:
    input:
        lambda w: bed_src_inputs(rules.download_superdups.output, si_to_superdups, w),
    output:
        segdup.inter.filtersort.data / "segdups.json",
    params:
        output_pattern=lambda w: to_output_pattern(segdup, "segdups", w),
    conda:
        "../envs/bedtools.yml"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeSuperdups
        ),
    benchmark:
        segdup.inter.filtersort.bench / "normalize_segdups.txt"
    script:
        "../scripts/python/bedtools/segdups/normalize_superdups.py"


rule merge_superdups:
    input:
        bed=lambda w: read_checkpoint("normalize_superdups", w),
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        segdup.final("segdups"),
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
        segdup.final("segdups_gt10kb"),
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
        segdup.final("notinsegdups"),
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
        segdup.final("notinsegdups_gt10kb"),


rule all_segdups:
    input:
        rules.merge_superdups.output,
        rules.filter_long_superdups.output,
        rules.notin_superdups.output,
        rules.notin_long_superdups.output,
    localrule: True
