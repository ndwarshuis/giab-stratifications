from common.config import CoreLevel

func_dir = CoreLevel.FUNCTIONAL
func_src_dir = config.ref_src_dir / func_dir.value
func_inter_dir = config.intermediate_build_dir / func_dir.value
func_log_build_dir = config.log_src_dir / func_dir.value


def func_final_path(name):
    return config.build_strat_path(func_dir, name)


use rule download_ref as download_ftbl with:
    output:
        func_src_dir / "ftbl.txt.gz",
    params:
        src=lambda w: config.refkey_to_functional_ftbl_src(w.ref_key),
    localrule: True
    log:
        func_log_build_dir / "ftbl.log",


use rule download_ref as download_gff with:
    output:
        func_src_dir / "gff.txt.gz",
    params:
        src=lambda w: config.refkey_to_functional_gff_src(w.ref_key),
    localrule: True
    log:
        func_log_build_dir / "gff.log",


rule combine_ftbl_and_gff:
    input:
        ftbl=rules.download_ftbl.output,
        gff=rules.download_gff.output,
    output:
        func_inter_dir / "refseq_cds.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/functional/combine_ftbl_and_gff.py"


rule merge_functional:
    input:
        bed=rules.combine_ftbl_and_gff.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        func_final_path("refseq_cds"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule invert_functional:
    input:
        bed=rules.merge_functional.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        func_final_path("notinrefseq_cds"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


rule all_functional:
    input:
        rules.merge_functional.output,
        rules.invert_functional.output,
    localrule: True
