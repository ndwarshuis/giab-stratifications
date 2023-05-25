from common.config import CoreLevel

telo_dir = CoreLevel.TELOMERES
telo_inter_dir = config.intermediate_build_dir / telo_dir.value
telo_log_dir = config.log_build_dir / telo_dir.value


def telo_final_path(name):
    return config.build_strat_path(telo_dir, name)


rule find_telomeres:
    input:
        ref=rules.filter_sort_ref.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        telo_final_path("telomeres"),
    conda:
        "../envs/seqtk.yml"
    log:
        telo_log_dir / "seqtk.log",
    shell:
        """
        seqtk telo {input.ref} 2> {log} | \
        cut -f1-3 | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c \
        > {output}
        """
