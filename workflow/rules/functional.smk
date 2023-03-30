func_src_dir = config.ref_src_dir / "FunctionalRegions"
func_inter_dir = config.intermediate_build_dir / "FunctionalRegions"
# func_log_dir = config.log_dir / "FunctionalRegions"


def func_final_path(name):
    return config.build_strat_path("FunctionalRegions", name)


use rule download_ref as download_ftbl with:
    output:
        func_src_dir / "ftbl.txt.gz",
    params:
        src=lambda w: config.refkey_to_functional_ftbl_src(w.ref_key),


use rule download_ref as download_gff with:
    output:
        func_src_dir / "func.gff.gz",
    params:
        src=lambda w: config.refkey_to_functional_gff_src(w.ref_key),


rule combine_ftbl_and_gff:
    input:
        ftbl=rules.download_ftbl.output,
        gff=rules.download_gff.output,
    output:
        func_inter_dir / "refseq_cds.bed.gz",
    conda:
        config.env_path("bedtools")
    script:
        config.python_script("bedtools/functional/combine_ftbl_and_gff.py")


rule merge_functional:
    input:
        bed=rules.combine_ftbl_and_gff.output,
        gapless=rules.get_gapless.output.auto,
    output:
        func_final_path("refseq_cds"),
    conda:
        config.env_path("bedtools")
    shell:
        """
        mergeBed -i {input.bed} | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
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
        config.env_path("bedtools")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


rule all_functional:
    input:
        rules.invert_functional.output,
