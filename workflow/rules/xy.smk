xy_src_dir = ref_src_dir / "XY"
xy_final_dir = final_dir / "XY"


rule download_genome_features_bed:
    output:
        xy_src_dir / "genome_features_{chr}.bed",
    params:
        url=lambda wildcards: lookup_ref_wc(
            ["XY", f"{wildcards.chr}_features"], wildcards
        ),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule download_genome_features_bed as download_X_PAR with:
    output:
        xy_final_dir / "GRCh38_chrX_PAR.bed",
    params:
        url=partial(lookup_ref_wc, ["XY", "X_PAR"]),


# TODO not sure where the actual PAR is, but this will do for now
# TODO the "chr" in front isn't constant across all refs
rule write_Y_PAR:
    output:
        xy_final_dir / "GRCh38_chrY_PAR.bed",
    shell:
        """
        echo "chrY\t10000\t2781479" >> {output}
        echo "chrY\t56887902\t57217415" >> {output}
        """


rule filter_xy_features:
    input:
        rules.download_genome_features_bed.output,
    output:
        xy_final_dir / "GRCh38_chr{chr}_{region,XTR|ampliconic}.bed",
    params:
        filt=lambda wildcards: "Ampliconic"
        if wildcards.region == "ampliconic"
        else wildcards.region,
    shell:
        """
        grep {params.filt} {input} | \
        cut -f 1-3 | \
        sort -k2,2n -k3,3n > \
        {output}
        """


rule invert_PAR:
    input:
        bed=lambda wildcards: rules.download_X_PAR.output
        if wildcards.chr == "X"
        else rules.write_Y_PAR.output,
        genome=rules.get_genome.output,
    output:
        xy_final_dir / "GRCh38_chr{chr}_nonPAR.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        grep {wildcards.chr} > \
        {output}
        """


rule filter_autosomes:
    input:
        rules.get_genome.output,
    output:
        xy_final_dir / "GRCh38_AllAutosomes.bed",
    shell:
        """awk -v OFS='\t' {{'print $1,\"0\",$2'}} {input} | \
        grep -v \"X\|Y\" > {output}
        """


# NOTE the XTR X region is different relative to old strats because it was
# updated a few months ago
rule all_xy:
    input:
        expand(
            rules.filter_xy_features.output,
            allow_missing=True,
            chr=["X", "Y"],
            region=["XTR", "ampliconic"],
        ),
        rules.download_X_PAR.output,
        rules.write_Y_PAR.output,
        expand(rules.invert_PAR.output, allow_missing=True, chr=["X", "Y"]),


# TODO the post processing scripts for these have merge steps; I probably don't
# need them (at least for hg38) but they were likely added for a reason


rule all_auto:
    input:
        rules.filter_autosomes.output,
