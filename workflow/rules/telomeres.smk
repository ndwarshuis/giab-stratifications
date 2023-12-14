from common.config import CoreLevel

telo = config.to_bed_dirs(CoreLevel.TELOMERES)


rule find_telomeres:
    input:
        ref=rules.filter_sort_ref.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        telo.final("telomeres"),
    conda:
        "../envs/seqtk.yml"
    log:
        telo.inter.postsort.log / "seqtk.log",
    shell:
        """
        seqtk telo {input.ref} 2> {log} | \
        cut -f1-3 | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c \
        > {output}
        """
