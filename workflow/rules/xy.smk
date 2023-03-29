xy_src_dir = config.ref_src_dir / "XY"
xy_inter_dir = config.build_intermediate_dir / "XY"


def xy_final_path(name):
    return config.build_strat_path("XY", name)


use rule download_ref as download_genome_features_bed with:
    output:
        xy_src_dir / "genome_features_{chr}.bed.gz",
    params:
        src=lambda w: (
            config.refkey_to_y_features_src
            if w.chr == "Y"
            else config.refkey_to_x_features_src
        )(w.ref_key),


use rule write_PAR_intermediate as write_PAR_final with:
    output:
        xy_inter_dir / "GRCh38_chr{chr}_PAR.bed.gz",


rule filter_XTR_features:
    input:
        bed=rules.download_genome_features_bed.output[0],
        gapless=rules.get_gapless.output.auto,
    output:
        xy_final_path("chr{chr}_XTR"),
    conda:
        envs_path("bedtools.yml")
    params:
        level="XTR",
    script:
        scripts_path("python/bedtools/xy/filter_sort_features.py")


use rule filter_XTR_features as filter_ampliconic_features with:
    input:
        bed=rules.download_genome_features_bed.output[0],
        gapless=rules.get_gapless.output.auto,
    output:
        xy_final_path("chr{chr}_ampliconic"),
    params:
        level="Ampliconic",


def par_input(wildcards):
    test_fun = config.want_xy_y if wildcards.chr == "Y" else config.want_xy_x
    return (
        rules.write_PAR_final.output
        if test_fun(wildcards.ref_key, wildcards.build_key)
        else rules.write_PAR_intermediate.output
    )


rule invert_PAR:
    input:
        bed=par_input,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        xy_final_path("chr{chr}_nonPAR"),
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        grep {wildcards.chr} | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


rule filter_autosomes:
    input:
        bed=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        xy_final_path("AllAutosomes"),
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        awk -v OFS='\t' {{'print $1,\"0\",$2'}} {input.bed} | \
        grep -v \"X\|Y\" | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


# NOTE the XTR X region is different relative to old strats because it was
# updated a few months ago
rule all_xy_sex:
    input:
        rules.filter_XTR_features.output,
        rules.filter_ampliconic_features.output,
        rules.invert_PAR.output,


# # TODO the post processing scripts for these have merge steps; I probably don't
# # need them (at least for hg38) but they were likely added for a reason


rule all_xy_auto:
    input:
        rules.filter_autosomes.output,


# def xy_inputs(wildcards):
#     rk = wildcards.ref_key
#     bk = wildcards.build_key
#     auto = rules.filter_autosomes.output if config.want_xy_auto(rk, bk) else []
#     sex = expand(
#         [
#             *rules.filter_XTR_features.output,
#             *rules.filter_ampliconic_features.output,
#             *rules.invert_PAR.output,
#         ],
#         allow_missing=True,
#         chr=config.wanted_xy_chr_names(rk, bk),
#     )
#     return auto + sex
# rule all_xy:
#     input:
#         xy_inputs,
#     output:
#         directory(xy_final_dir),
