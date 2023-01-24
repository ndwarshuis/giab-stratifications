ref_src_dir = resources_dir / "{ref_key}"
ref_master_dir = results_dir / "ref" / "{ref_key}"
ref_inter_dir = intermediate_dir / "ref"


rule download_ref:
    output:
        ref_src_dir / "ref.fna.gz",
    params:
        url=partial(lookup_ref_wc, ["ref_url"]),
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unzip_ref:
    input:
        rules.download_ref.output,
    output:
        ref_master_dir / "ref.fna",
    conda:
        envs_path("utils.yml")
    shell:
        "gunzip -c {input} > {output}"


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


rule get_genome:
    input:
        rules.index_ref.output,
    output:
        ref_inter_dir / "genome.txt",
    params:
        filt=lookup_filter_wc,
    shell:
        """
        cut -f 1,2 {input} | \
        sed -n '/^\(#\|{params.filt}\)\t/p' \
        > {output}
        """


# TODO filter a subset of this for testing
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
        url=partial(lookup_ref_wc, ["gap_url"]),
    conda:
        envs_path("utils.yml")


rule merge_gaps:
    input:
        gaps=rules.download_gaps.output,
        genome=rules.get_genome.output,
    output:
        ref_inter_dir / "gaps_merged.bed",
    conda:
        envs_path("bedtools.yml")
    params:
        filt=lookup_filter_wc,
    shell:
        """
        gunzip -c {input.gaps} | \
        cut -f2-4 | \
        sed -n '/^\(#\|{params.filt}\)\t/p' | \
        bedtools sort -t stdin -g {input.genome} | \
        bedtools merge -i stdin -d 100 \
        > {output}
        """
