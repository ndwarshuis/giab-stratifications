from os.path import dirname, basename
from more_itertools import unique_everseen, unzip
from os import scandir
from common.config import CoreLevel

post_inter_dir = config.intermediate_build_dir / "postprocess"
post_log_dir = config.log_build_dir / "postprocess"


# NOTE This will only give a limited set of "targets" from each of the strat
# level types. Because of the way snakemake works, it is not practical or
# maintainable to specify every single stratification file; therefore we only
# use the "toplevel" targets which will pull in all others. In downstream rules
# from this, it is much easier to condense this to a list of parent directories
# then manually iterate all strat files in these directories.
def expand_strat_targets(wildcards):
    rk = wildcards.ref_key
    bk = wildcards.build_key

    # all targets except sex (which needs an additional wildcard)
    targets = [
        (rules.all_low_complexity.input, config.want_low_complexity),
        (rules.all_xy_auto.input, config.want_xy_auto),
        (rules.all_map.input, config.want_mappability),
        (
            [
                rules.all_gc.input.wide[0],
                rules.all_gc.input.narrow[0],
                rules.all_gc.input.middle[0],
            ],
            config.want_gc,
        ),
        (rules.all_functional.input, config.want_functional),
        (rules.all_segdups.input, config.want_segdups),
        (rules.find_telomeres.output, config.want_telomeres),
        (rules.invert_segdup_and_map.output, config.want_segdup_and_map),
        (rules.invert_alldifficult.output, config.want_alldifficult),
    ]
    auto = [t for tgt, test in targets if test(rk, bk) for t in tgt]
    other = all_other(rk, bk)

    # xy (expand target depending on which chromosomes we have selected)
    sex = all_xy_features(wildcards) + all_xy_PAR(wildcards)

    # combine and ensure that all "targets" refer to final bed files
    all_targets = auto + sex + other
    invalid = [f for f in all_targets if not f.startswith(str(config.final_root_dir))]
    assert len(invalid) == 0, f"invalid targets: {invalid}"
    return all_targets


rule list_all_strats:
    input:
        expand_strat_targets,
    output:
        post_inter_dir / "all_strats.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/list_strats.py"


rule generate_md5sums:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "{ref_key}-genome-stratifications-md5s.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/list_md5.py"


rule generate_tsv_list:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "{ref_key}-all-stratifications.tsv",
    script:
        "../scripts/python/bedtools/postprocess/generate_tsv.py"


rule unit_test_strats:
    input:
        strats=rules.list_all_strats.output,
        gapless_auto=rules.get_gapless.output.auto,
        gapless_parY=rules.get_gapless.output.parY,
    output:
        touch(post_inter_dir / "unit_tests.done"),
    log:
        post_log_dir / "unit_tests.log",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/run_unit_tests.py"


rks, bks = map(
    list,
    unzip([(rk, bk) for rk, r in config.stratifications.items() for bk in r.builds]),
)


rule download_comparison_strat_tarball:
    output:
        directory(config.resources_dir / "comparisons" / "{compare_key}"),
    conda:
        "../envs/utils.yml"
    # ASSUME this is not None (we control for this when building the "all" rule)
    params:
        url=lambda w: config.comparison_strats[w.compare_key],
        tar=lambda w: f"/tmp/comparison_{w.compare_key}.tar.gz",
    shell:
        """
        curl -fsSq -o {params.tar} {params.url} && \
        mkdir {output} && \
        tar xzf {params.tar} --directory {output} --strip-components=1
        """


rule compare_strats:
    input:
        old=lambda w: expand(
            rules.download_comparison_strat_tarball.output,
            compare_key=config.buildkey_to_comparekey(w.ref_key, w.build_key),
        )[0],
        # use this to target a specific rule to satisfy the snakemake scheduler,
        # the thing I actually need here is the parent directory
        new_list=rules.generate_tsv_list.output[0],
    output:
        anti=post_inter_dir / "comparison_anti.tsv.gz",
        diagnostics=post_inter_dir / "comparison_diagnostics.tsv",
    log:
        post_log_dir / "comparison.log",
    conda:
        "../envs/bedtools.yml"
    threads: 8
    script:
        "../scripts/python/bedtools/postprocess/diff_previous_strats.py"


rule all_comparisons:
    input:
        [
            expand(rules.compare_strats.output, ref_key=rk, build_key=bk)[0]
            for rk, r in config.stratifications.items()
            for bk in r.builds
            if config.buildkey_to_comparekey(rk, bk) is not None
        ],


# Silly rule to make a dataframe which maps chromosome names to a common
# nomenclature (for the validation markdown). The only reason this exists is
# because snakemake apparently mangles this elegant list of tuples when I pass
# it as a param to the rmd script...dataframe in IO monad it is :/
rule write_chr_name_mapper:
    output:
        config.final_root_dir / ".validation" / "chr_mapper.tsv",
    run:
        with open(output[0], "w") as f:
            for ref_key, build_key in zip(rks, bks):
                for i in config.buildkey_to_chr_indices(ref_key, build_key):
                    prefix = config.refkey_to_final_chr_prefix(ref_key)
                    line = [
                        i.chr_name_full(prefix),
                        i.chr_name,
                        f"{ref_key}@{build_key}",
                    ]
                    f.write("\t".join(line) + "\n")


rule validate_strats:
    input:
        # this first input isn't actually used, but ensures the unit tests pass
        # before running the rmd script
        **{
            "_test": expand(
                rules.unit_test_strats.output,
                zip,
                ref_key=rks,
                build_key=bks,
            ),
            "strats": expand(
                rules.list_all_strats.output,
                zip,
                ref_key=rks,
                build_key=bks,
            ),
            "nonN": expand(
                rules.get_gapless.output.auto,
                zip,
                ref_key=rks,
                build_key=bks,
            ),
        },
        chr_mapper=rules.write_chr_name_mapper.output,
    output:
        config.final_root_dir / ".validation" / "coverage_plots.html",
    conda:
        "../envs/rmarkdown.yml"
    params:
        core_levels=[c.value for c in CoreLevel],
        other_levels=config.other_levels,
    script:
        "../scripts/rmarkdown/rmarkdown/validate.Rmd"


rule unzip_ref:
    input:
        rules.filter_sort_ref.output,
    output:
        post_log_dir / "happy" / "ref.fa",
    shell:
        "gunzip -c {input} > {output}"


rule index_unzipped_ref:
    input:
        rules.unzip_ref.output,
    output:
        rules.unzip_ref.output[0] + ".fai",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        samtools faidx {input}
        """


rule run_happy:
    input:
        refi=rules.index_unzipped_ref.output,
        ref=rules.unzip_ref.output,
        bench_vcf=rules.download_bench_vcf.output,
        bench_bed=rules.filter_sort_bench_bed.output,
        query_vcf=rules.download_query_vcf.output,
        # query_bed=rules.run_dipcall.output.bed,
        strats=rules.generate_tsv_list.output,
    output:
        post_inter_dir / "happy" / "happy.extended.csv",
    params:
        prefix=lambda _, output: str(output[0]).replace(".extended.csv", ""),
    conda:
        "../envs/happy.yml"
    log:
        post_log_dir / "happy" / "happy.log",
    threads: 8
    script:
        "../scripts/python/bedtools/postprocess/run_happy.py"


rule summarize_happy:
    input:
        [
            expand(
                rules.run_happy.output,
                ref_key=rk,
                build_key=bk,
            )
            for rk, bk in zip(rks, bks)
            if config.want_benchmark(rk, bk)
        ],
    output:
        config.final_root_dir / ".validation" / "benchmark_summary.html",
    conda:
        "../envs/rmarkdown.yml"
    script:
        "../scripts/rmarkdown/rmarkdown/benchmark.Rmd"
