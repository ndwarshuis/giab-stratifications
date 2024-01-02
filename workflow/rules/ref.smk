from common.config import si_to_gaps, bd_to_bench_bed, bd_to_bench_vcf, bd_to_query_vcf

ref = config.ref_dirs


# lots of things depend on PAR which is why this isn't part of the XY module
rule write_PAR_intermediate:
    output:
        ref.inter.build.data / "chr{sex_chr}_PAR.bed.gz",
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/xy/write_par.py"


rule download_ref:
    output:
        ref.src.reference.data / "ref.fna.gz",
    params:
        src=lambda w: config.refsrckey_to_ref_src(w.ref_src_key),
    log:
        ref.src.reference.log / "download_ref.log",
    conda:
        "../envs/bedtools.yml"
    localrule: True
    script:
        "../scripts/python/bedtools/misc/get_file.py"


rule index_ref:
    input:
        lambda w: expand_final_to_src(rules.download_ref.output, w, bla=1),
    output:
        ref.inter.prebuild.data / "ref.fna.fai",
    conda:
        "../envs/utils.yml"
    log:
        ref.inter.prebuild.log / "index_ref.log",
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
        ref.inter.build.data / "genome.txt",
    conda:
        "../envs/bedtools.yml"
    log:
        ref.inter.build.data / "get_genome.log",
    script:
        "../scripts/python/bedtools/ref/get_genome.py"


rule filter_sort_ref:
    input:
        fa=lambda w: expand_final_to_src(rules.download_ref.output, w, bla=2),
        genome=rules.get_genome.output,
    output:
        ref.inter.build.data / "ref_filtered.fa",
    conda:
        "../envs/utils.yml"
    log:
        ref.inter.build.log / "filter_sort_ref.log",
    shell:
        """
        samtools faidx {input.fa} $(cut -f1 {input.genome} | tr '\n' ' ') 2> {log} | \
        bgzip -c > {output}
        """


use rule download_ref as download_gaps with:
    output:
        ref.src.reference.data / "gap.bed.gz",
    params:
        src=lambda w: config.refsrckey_to_bed_src(si_to_gaps, w.ref_src_key),
    localrule: True
    log:
        ref.src.reference.log / "download_gaps.log",


def gapless_input(wildcards):
    rk, _ = parse_final_refkey(wildcards.ref_final_key)
    si = config.to_ref_data(rk).strat_inputs
    if si.gap is not None:
        gaps = {
            "gaps": expand(
                rules.download_gaps.output,
                ref_src_key=config.refkey_to_bed_refsrckeys(si_to_gaps, rk),
            )
        }
        if si.xy.y_par is None:
            return gaps
        else:
            return {
                **gaps,
                "parY": expand(
                    rules.write_PAR_intermediate.output[0],
                    allow_missing=True,
                    sex_chr="Y",
                )[0],
            }
    else:
        return {}


rule get_gapless:
    input:
        unpack(gapless_input),
        genome=rules.get_genome.output[0],
    output:
        auto=ref.inter.build.data / "genome_gapless.bed.gz",
        parY=ref.inter.build.data / "genome_gapless_parY.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/ref/get_gapless.py"


use rule download_ref as download_bench_vcf with:
    output:
        ref.src.benchmark.data / "bench.vcf.gz",
    params:
        src=lambda w: config.buildkey_to_bed_src(
            bd_to_bench_bed,
            w.ref_src_key,
            w.build_key,
        ),
    localrule: True
    log:
        ref.src.benchmark.log / "download_bench_vcf.log",


use rule download_ref as download_bench_bed with:
    output:
        ref.src.benchmark.data / "bench.bed.gz",
    params:
        src=lambda w: config.buildkey_to_vcf_src(
            bd_to_bench_vcf,
            w.ref_src_key,
            w.build_key,
        ),
    localrule: True
    log:
        ref.src.benchmark.log / "download_bench_bed.log",


use rule download_ref as download_query_vcf with:
    output:
        ref.src.benchmark.data / "query.vcf.gz",
    params:
        src=lambda w: config.buildkey_to_vcf_src(
            bd_to_query_vcf,
            w.ref_src_key,
            w.build_key,
        ),
    localrule: True
    log:
        ref.src.benchmark.log / "download_query_vcf.log",


rule filter_sort_bench_bed:
    input:
        lambda w: expand_final_to_src(rules.download_bench_bed.output, w),
    output:
        ref.inter.build.data / "bench_filtered.bed.gz",
    log:
        ref.inter.build.log / "bench_filtered.bed.gz",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/ref/filter_sort_bench_bed.py"
