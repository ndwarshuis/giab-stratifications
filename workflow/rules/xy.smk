xy_src_dir = ref_src_dir / "XY"
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


rule write_PAR:
    output:
        xy_final_dir / "GRCh38_chr{chr}_PAR.bed.gz",
    conda:
        envs_path("bedtools.yml")
    script:
        scripts_path("python/bedtools/xy/write_par.py")


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


rule invert_PAR:
    input:
        bed=rules.write_PAR.output,
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
rule all_xy:
    input:
        expand(
            [
                *rules.write_PAR.output,
                *rules.filter_XTR_features.output,
                *rules.filter_ampliconic_features.output,
                *rules.invert_PAR.output,
            ],
            allow_missing=True,
            chr=["X", "Y"],
        ),


# TODO the post processing scripts for these have merge steps; I probably don't
# need them (at least for hg38) but they were likely added for a reason


rule all_auto:
    input:
        rules.filter_autosomes.output,
