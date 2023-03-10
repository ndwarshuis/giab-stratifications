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


use rule download_genome_features_bed as download_X_PAR with:
    output:
        xy_src_dir / "GRCh38_chrX_PAR.bed.gz",
    params:
        src=lambda w: config.refkey_to_x_par_src(w.ref_key),


rule compress_X_PAR:
    input:
        rules.download_X_PAR.output,
    output:
        xy_final_dir / "GRCh38_chrX_PAR.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "cat {input} | bgzip -c > {output}"


# TODO not sure where the actual PAR is, but this will do for now
# TODO the "chr" in front isn't constant across all refs
rule write_Y_PAR:
    output:
        xy_final_dir / "GRCh38_chrY_PAR.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        echo "chrY\t10000\t2781479\nchrY\t56887902\t57217415" | \
        bgzip -c > {output}
        """


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
        bed=lambda wildcards: rules.download_X_PAR.output
        if wildcards.chr == "X"
        else rules.write_Y_PAR.output,
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
                *rules.filter_XTR_features.output,
                *rules.filter_ampliconic_features.output,
                *rules.invert_PAR.output,
            ],
            allow_missing=True,
            chr=["X", "Y"],
        ),
        rules.compress_X_PAR.output,
        rules.write_Y_PAR.output,


# TODO the post processing scripts for these have merge steps; I probably don't
# need them (at least for hg38) but they were likely added for a reason


rule all_auto:
    input:
        rules.filter_autosomes.output,
