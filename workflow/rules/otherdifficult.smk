from common.config import CoreLevel


def other_difficult_final_path(name):
    return config.build_strat_path(CoreLevel.OtherDifficult, name)


rule get_gaps:
    input:
        gapless=rules.get_gapless.output.auto,
        genome=rules.get_genome.output,
    output:
        other_difficult_final_path("gaps_slop15kb"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.gapless} -g {input.genome} | \
        slopBed -i stdin -b 15000 -g {input.genome} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """
