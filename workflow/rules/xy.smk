xy_src_dir = ref_src_dir / "XY"
xy_inter_dir = intermediate_dir / "XY"
xy_final_dir = final_dir / "XY"


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
        rules.download_genome_features_bed.output,
    output:
        xy_final_dir / "GRCh38_chr{chr}_XTR.bed.gz",
    params:
        filt="XTR",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input} | \
        grep {params.filt} | \
        cut -f 1-3 | \
        sort -k2,2n -k3,3n | \
        bgzip -c > {output}
        """


use rule filter_XTR_features as filter_ampliconic_features with:
    input:
        rules.download_genome_features_bed.output,
    output:
        xy_final_dir / "GRCh38_chr{chr}_ampliconic.bed.gz",
    params:
        filt="Ampliconic",


# rule filter_xy_features:
#     input:
#         rules.download_genome_features_bed.output,
#     output:
#         xy_final_dir / "GRCh38_chr{chr}_{region,XTR|ampliconic}.bed.gz",
#     params:
#         filt=lambda wildcards: "Ampliconic"
#         if wildcards.region == "ampliconic"
#         else wildcards.region,
#     conda:
#         envs_path("bedtools.yml")
#     shell:
#         """
#         grep {params.filt} {input} | \
#         cut -f 1-3 | \
#         sort -k2,2n -k3,3n | \
#         bgzip -c > {output}
#         """


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
    output:
        xy_final_dir / "GRCh38_chr{chr}_nonPAR.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        grep {wildcards.chr} | \
        bgzip -c > {output}
        """


rule filter_autosomes:
    input:
        rules.get_genome.output,
    output:
        xy_final_dir / "GRCh38_AllAutosomes.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        awk -v OFS='\t' {{'print $1,\"0\",$2'}} {input} | \
        grep -v \"X\|Y\" | \
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
