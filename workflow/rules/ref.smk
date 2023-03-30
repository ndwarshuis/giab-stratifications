ref_master_dir = config.intermediate_root_dir / "ref" / "{ref_key}"
ref_master_log_dir = config.log_root_dir / "ref" / "{ref_key}"

ref_inter_dir = config.intermediate_build_dir / "ref"
ref_log_dir = config.log_build_dir / "ref"


# lots of things depend on PAR so move this out of the XY ruleset
rule write_PAR_intermediate:
    output:
        ref_inter_dir / "chr{chr}_PAR.bed.gz",
    conda:
        config.env_path("bedtools")
    script:
        config.python_script("bedtools/xy/write_par.py")


rule download_ref:
    output:
        config.ref_src_dir / "ref.fna.gz",
    params:
        src=lambda w: config.refkey_to_ref_src(w.ref_key),
    conda:
        config.env_path("bedtools")
    script:
        config.python_script("bedtools/misc/get_file.py")


rule unzip_ref:
    input:
        rules.download_ref.output,
    output:
        ref_master_dir / "ref.fna",
    conda:
        config.env_path("utils")
    shell:
        "gunzip -c {input} > {output}"


# TODO is this only used for getting the genome? if so combine with rule 'get_genome'
rule index_ref:
    input:
        rules.download_ref.output,
    output:
        ref_master_dir / "ref.fna.fai",
    conda:
        config.env_path("utils")
    log:
        ref_master_log_dir / "index_ref.log",
    shell:
        """
        gunzip -c {input} | \
        samtools faidx - -o - \
        2> {log} > {output}
        """


rule get_genome:
    input:
        rules.index_ref.output,
    output:
        ref_inter_dir / "genome.txt",
    conda:
        config.env_path("bedtools")
    log:
        ref_log_dir / "get_genome.log",
    script:
        config.python_script("bedtools/ref/get_genome.py")


rule genome_to_bed:
    input:
        rules.get_genome.output,
    output:
        ref_inter_dir / "genome.bed",
    shell:
        "awk 'BEGIN {{ FS = OFS = \"\t\"}} {{ print $1, 0, $2 }}' {input} > {output}"


rule filter_sort_ref:
    input:
        fa=rules.download_ref.output,
        genome=rules.get_genome.output,
    output:
        ref_inter_dir / "ref_filtered.fa",
    conda:
        config.env_path("utils")
    log:
        ref_log_dir / "filter_sort_ref.log",
    shell:
        """
        samtools faidx {input.fa} $(cut -f1 {input.genome} | tr '\n' ' ') 2> {log} > \
        {output}
        """


use rule download_ref as download_gaps with:
    output:
        config.ref_src_dir / "gap.bed.gz",
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
        config.env_path("bedtools")
    script:
        scripts_path("python/bedtools/ref/filter_sort_gaps.py")


# This is super convoluted since many rules use the gap file, but I only want
# to do all this logic once. Shell script it is...I guess this works
rule get_gapless:
    input:
        unpack(
            lambda w: {
                "gaps": rules.merge_gaps.output[0],
                "parY": expand(
                    rules.write_PAR_intermediate.output,
                    allow_missing=True,
                    chr="Y",
                ),
            }
            if config.refkey_to_gap_src(w.ref_key)
            else {"gaps": [], "parY": []}
        ),
        genome=rules.genome_to_bed.output,
    output:
        auto=ref_inter_dir / "genome_gapless.bed",
        parY=ref_inter_dir / "genome_gapless_parY.bed",
    conda:
        config.env_path("bedtools")
    params:
        hasgaps=lambda _, input: 1 if "gaps" in input else 0,
    shell:
        """
        gapfile="{input.gaps}"
        if [[ -z "$gapfile" ]]; then
            ln -sfr {input.genome} {output.auto}
            ln -sfr {input.genome} {output.parY}
        else
            bedtools subtract -a {input.genome} -b {input.gaps} > {output.parY}
            bedtools subtract -a {output.parY} -b {input.parY} > {output.auto}
        fi
        """
