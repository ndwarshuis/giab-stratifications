ref_src_dir = resources_dir / "{ref_key}"
ref_master_dir = results_dir / "ref" / "{ref_key}"
ref_inter_dir = intermediate_dir / "ref"


rule download_ref:
    output:
        ref_src_dir / "ref.fna.gz",
    params:
        src=lambda w: config.refkey_to_ref_src(w.ref_key),
    conda:
        envs_path("bedtools.yml")
    script:
        scripts_path("python/bedtools/get_file.py")


rule unzip_ref:
    input:
        rules.download_ref.output,
    output:
        ref_master_dir / "ref.fna",
    conda:
        envs_path("utils.yml")
    shell:
        "gunzip -c {input} > {output}"


# TODO don't hardcode chr here
rule index_ref:
    input:
        rules.download_ref.output,
    output:
        ref_master_dir / "ref.fna.fai",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        samtools faidx - -o - | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        sed 's/^chr//' | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n -k2,2n -k3,3n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' > \
        {output}
        """


# TODO don't use filt here
rule get_genome:
    input:
        rules.index_ref.output,
    output:
        ref_inter_dir / "genome.txt",
    params:
        filt=lambda wildcards: config.buildkey_to_chr_pattern(
            wildcards.ref_key, wildcards.build_key
        ),
    shell:
        """
        cut -f 1,2 {input} | \
        sed -n '/^\(#\|{params.filt}\)\t/p' \
        > {output}
        """


rule genome_to_bed:
    input:
        rules.get_genome.output,
    output:
        ref_inter_dir / "genome.bed",
    shell:
        "awk 'BEGIN {{ FS = OFS = \"\t\"}} {{ print $1, 0, $2 }}' {input} > {output}"


rule filter_sort_ref:
    input:
        fa=rules.unzip_ref.output,
        genome=rules.get_genome.output,
    output:
        ref_inter_dir / "ref_filtered.fa",
    conda:
        envs_path("utils.yml")
    shell:
        """
        samtools faidx {input.fa} $(cut -f1 {input.genome} | tr '\n' ' ') > \
        {output}
        """


use rule download_ref as download_gaps with:
    output:
        ref_src_dir / "gap.bed.gz",
    params:
        src=lambda w: config.refkey_to_gap_src(w.ref_key),


# TODO don't use filt here
rule merge_gaps:
    input:
        gaps=rules.download_gaps.output,
        genome=rules.get_genome.output,
    output:
        ref_inter_dir / "gaps_merged.bed",
    conda:
        envs_path("bedtools.yml")
    params:
        filt=lambda wildcards: config.buildkey_to_chr_pattern(
            wildcards.ref_key, wildcards.build_key
        ),
    shell:
        """
        gunzip -c {input.gaps} | \
        cut -f2-4 | \
        sed -n '/^\(#\|{params.filt}\)\t/p' | \
        bedtools sort -i stdin -g {input.genome} | \
        bedtools merge -i stdin -d 100 \
        > {output}
        """
