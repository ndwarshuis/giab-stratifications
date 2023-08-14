from common.config import CoreLevel

otherdiff_dir = CoreLevel.OtherDifficult
otherdiff_inter_dir = config.intermediate_build_dir / otherdiff_dir.value


def other_difficult_final_path(name):
    return config.build_strat_path(otherdiff_dir, name)


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
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_vdj:
    input:
        mapper=rules.ftbl_to_mapper.output[0],
        bed=rules.gff_to_bed.output[0],
    output:
        otherdiff_inter_dir / "vdj.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/otherdifficult/filter_vdj.py"


rule remove_vdj_gaps:
    input:
        bed=rules.filter_vdj.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        other_difficult_final_path("VDJ"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """
