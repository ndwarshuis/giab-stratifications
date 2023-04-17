ref_master_dir = config.intermediate_root_dir / "ref" / "{ref_key}"
ref_master_log_dir = config.log_root_dir / "ref" / "{ref_key}"

ref_inter_dir = config.intermediate_build_dir / "ref"
ref_log_dir = config.log_build_dir / "ref"


# lots of things depend on PAR so move this out of the XY ruleset
rule write_PAR_intermediate:
    output:
        ref_inter_dir / "chr{sex_chr}_PAR.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/xy/write_par.py"


rule download_ref:
    output:
        config.ref_src_dir / "ref.fna.gz",
    params:
        src=lambda w: config.refkey_to_ref_src(w.ref_key),
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/misc/get_file.py"


# TODO is this only used for getting the genome? if so combine with rule 'get_genome'
rule index_ref:
    input:
        rules.download_ref.output,
    output:
        ref_master_dir / "ref.fna.fai",
    conda:
        "../envs/utils.yml"
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
        "../envs/bedtools.yml"
    log:
        ref_log_dir / "get_genome.log",
    script:
        "../scripts/python/bedtools/ref/get_genome.py"


rule filter_sort_ref:
    input:
        fa=rules.download_ref.output,
        genome=rules.get_genome.output,
    output:
        ref_inter_dir / "ref_filtered.fa",
    conda:
        "../envs/utils.yml"
    log:
        ref_log_dir / "filter_sort_ref.log",
    shell:
        """
        samtools faidx {input.fa} $(cut -f1 {input.genome} | tr '\n' ' ') 2> {log} | \
        bgzip -c > {output}
        """


use rule download_ref as download_gaps with:
    output:
        config.ref_src_dir / "gap.bed.gz",
    params:
        src=lambda w: config.refkey_to_gap_src(w.ref_key),
    localrule: True


rule get_gapless:
    input:
        unpack(
            lambda w: {
                "gaps": rules.download_gaps.output[0],
                "parY": expand(
                    rules.write_PAR_intermediate.output[0],
                    allow_missing=True,
                    sex_chr="Y",
                )[0],
            }
            if config.refkey_to_gap_src(w.ref_key)
            else {}
        ),
        genome=rules.get_genome.output[0],
    output:
        auto=ref_inter_dir / "genome_gapless.bed",
        parY=ref_inter_dir / "genome_gapless_parY.bed",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/ref/get_gapless.py"
