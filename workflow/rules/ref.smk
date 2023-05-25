ref_dir = Path("ref")
ref_inter_dir = config.intermediate_build_dir / ref_dir
ref_log_src_dir = config.log_src_dir / ref_dir
ref_log_build_dir = config.log_build_dir / ref_dir

_ref_master_dir = ref_dir / "{ref_key}"
ref_master_dir = config.intermediate_root_dir / _ref_master_dir
ref_master_log_dir = config.log_results_dir / _ref_master_dir

bench_src_dir = config.ref_src_dir / "bench"
bench_src_log_dir = ref_log_src_dir / "bench"


# lots of things depend on PAR which is why this isn't part of the XY module
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
    log:
        ref_log_src_dir / "download_ref.log",
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/misc/get_file.py"


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
        ref_log_build_dir / "get_genome.log",
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
        ref_log_build_dir / "filter_sort_ref.log",
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
    log:
        ref_log_src_dir / "download_gaps.log",


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
        auto=ref_inter_dir / "genome_gapless.bed.gz",
        parY=ref_inter_dir / "genome_gapless_parY.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/ref/get_gapless.py"


use rule download_ref as download_bench_vcf with:
    output:
        bench_src_dir / "{build_key}_bench.vcf.gz",
    params:
        src=lambda w: config.refkey_to_bench_vcf_src(w.ref_key, w.build_key),
    localrule: True
    log:
        bench_src_log_dir / "{build_key}_download_bench_vcf.log",


use rule download_ref as download_bench_bed with:
    output:
        bench_src_dir / "{build_key}_bench.bed.gz",
    params:
        src=lambda w: config.refkey_to_bench_bed_src(w.ref_key, w.build_key),
    localrule: True
    log:
        bench_src_log_dir / "{build_key}_download_bench_bed.log",


use rule download_ref as download_query_vcf with:
    output:
        bench_src_dir / "{build_key}_query.vcf.gz",
    params:
        src=lambda w: config.refkey_to_query_vcf_src(w.ref_key, w.build_key),
    localrule: True
    log:
        bench_src_log_dir / "{build_key}_download_query_vcf.log",


rule filter_sort_bench_bed:
    input:
        rules.download_bench_bed.output,
    output:
        ref_inter_dir / "bench_filtered.bed.gz",
    log:
        ref_log_build_dir / "bench_filtered.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/ref/filter_sort_bench_bed.py"
