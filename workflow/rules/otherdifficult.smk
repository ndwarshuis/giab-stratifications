from common.config import CoreLevel

odiff = config.to_bed_dirs(CoreLevel.OTHER_DIFFICULT)


rule get_gaps:
    input:
        gapless=rules.get_gapless.output.auto,
        genome=rules.get_genome.output,
    output:
        odiff.final("gaps_slop15kb"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.gapless} -g {input.genome} | \
        slopBed -i stdin -b 15000 -g {input.genome} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule remove_vdj_gaps:
    input:
        bed=lambda w: read_named_checkpoint("normalize_cds", "vdj", w),
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        odiff.final("VDJ"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """
