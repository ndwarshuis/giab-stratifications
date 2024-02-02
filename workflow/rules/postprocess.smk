from os.path import dirname, basename
from more_itertools import unique_everseen, unzip
from os import scandir
from common.config import CoreLevel, strip_full_refkey
from common.functional import DesignError

post_inter_dir = config.intermediate_build_dir / "postprocess"
post_log_dir = config.log_build_dir / "postprocess"
post_bench_dir = config.bench_build_dir / "postprocess"
validation_dir = config.final_root_dir / "validation"


def expand_strat_targets_inner(ref_final_key, build_key):
    # TODO this is crude
    # NOTE we need to do this because despite the final key potentially having
    # one or two haps, the configuration applies to both haps simultaneously
    bd = config.to_build_data(strip_full_refkey(ref_final_key), build_key)
    function_targets = [
        (all_low_complexity, bd.want_low_complexity),
        (gc_inputs_flat, bd.want_gc),
        (mappabilty_inputs, bd.want_mappability),
        (het_hom_inputs, bd.want_hets),
        (all_xy_sex, True),
        (all_other, True),
    ]
    rule_targets = [
        (rules.filter_autosomes.output, bd.want_xy_auto),
        (rules.all_functional.input, bd.want_functional),
        (rules.all_segdups.input, bd.want_segdups),
        (rules.find_telomeres.output, bd.want_telomeres),
        (rules.all_segdup_and_map.input, bd.want_segdup_and_map),
        (rules.all_alldifficult.input, bd.want_alldifficult),
        (rules.get_gaps.output, lambda r, _: bd.want_gaps(r)),
        (rules.remove_vdj_gaps.output, bd.want_vdj),
    ]
    all_function = [f(ref_final_key, build_key) for f, test in function_targets if test]
    all_rule = [tgt for tgt, test in rule_targets if test]

    # combine and ensure that all "targets" refer to final bed files
    all_targets = [y for xs in all_function + all_rule for y in xs]
    invalid = [f for f in all_targets if not f.startswith(str(config.final_root_dir))]
    if len(invalid) > 0:
        raise DesignError(f"invalid targets: {invalid}")
    return all_targets


def expand_strat_targets(wildcards):
    return expand_strat_targets_inner(wildcards.ref_final_key, wildcards.build_key)


rule list_all_strats:
    input:
        expand_strat_targets,
    output:
        post_inter_dir / "all_strats.txt",
    localrule: True
    script:
        "../scripts/python/bedtools/postprocess/list_strats.py"


rule generate_md5sums:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "{ref_final_key}-genome-stratifications-md5s.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/list_md5.py"


rule generate_tsv_list:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "{ref_final_key}-all-stratifications.tsv",
    localrule: True
    script:
        "../scripts/python/bedtools/postprocess/generate_tsv.py"


rule unit_test_strats:
    input:
        strats=rules.list_all_strats.output,
        gapless_auto=rules.get_gapless.output.auto,
        gapless_parY=rules.get_gapless.output.parY,
        strat_list=rules.generate_tsv_list.output[0],
        checksums=rules.generate_md5sums.output[0],
    output:
        touch(post_inter_dir / "unit_tests.done"),
    log:
        post_log_dir / "unit_tests.log",
    benchmark:
        post_bench_dir / "unit_test_strats.txt"
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/run_unit_tests.py"


rule get_coverage_table:
    input:
        _test=rules.unit_test_strats.output,
        bedlist=rules.list_all_strats.output[0],
        gapless=rules.get_gapless.output.auto,
    output:
        touch(post_inter_dir / "coverage.tsv.gz"),
    benchmark:
        post_bench_dir / "get_coverage_table.txt"
    script:
        "../scripts/python/bedtools/postprocess/get_coverage_table.py"


rule download_comparison_strat_tarball:
    output:
        directory(config.resources_dir / "comparisons" / "{compare_key}"),
    conda:
        "../envs/utils.yml"
    # ASSUME this is not None (we control for this when building the "all" rule)
    params:
        url=lambda w: config.comparison_strats[w.compare_key],
        tar=lambda w: f"/tmp/comparison_{w.compare_key}.tar.gz",
    localrule: True
    shell:
        """
        curl -fsSq -o {params.tar} {params.url} && \
        mkdir {output} && \
        tar xzf {params.tar} --directory {output} --strip-components=1
        """


rule compare_strats:
    input:
        # ensure tests are run before this rule
        _test=rules.unit_test_strats.output,
        # NOTE: in the case id dip2 configurations, each comparison key will
        # correspond to one half of the diploid dataset (each tarball is one
        # haplotype)
        old=lambda w: expand(
            rules.download_comparison_strat_tarball.output,
            compare_key=config.to_build_data(
                w.ref_final_key, w.build_key
            ).build.compare_key,
        )[0],
        # use this to target a specific rule to satisfy the snakemake scheduler,
        # the thing I actually need here is the parent directory
        new_list=rules.generate_tsv_list.output[0],
    output:
        validation_dir / "{ref_final_key}@{build_key}" / "diagnostics.tsv",
    log:
        post_log_dir / "comparison.log",
    conda:
        "../envs/diff.yml"
    threads: 8
    script:
        "../scripts/python/bedtools/postprocess/diff_previous_strats.py"


rule all_comparisons:
    input:
        [
            expand(rules.compare_strats.output, ref_final_key=rk, build_key=bk)[0]
            for rk, bk in zip(*config.all_build_keys)
            if config.to_build_data(rk, bk).build.compare_key is not None
        ],
    localrule: True


rule make_coverage_plots:
    input:
        # this first input isn't actually used, but ensures the unit tests pass
        # before running the rmd script
        **{
            "_test": expand(
                rules.unit_test_strats.output,
                zip,
                ref_final_key=(t := config.all_full_build_keys)[0],
                build_key=t[1],
            ),
            "coverage": expand(
                rules.get_coverage_table.output,
                zip,
                ref_final_key=t[0],
                build_key=t[1],
            ),
        },
    output:
        validation_dir / "coverage_plots.html",
    conda:
        "../envs/rmarkdown.yml"
    params:
        core_levels=[c.value for c in CoreLevel],
        other_levels=config.other_levels,
    script:
        "../scripts/rmarkdown/rmarkdown/coverage_plots.Rmd"


rule run_happy:
    input:
        _test=rules.unit_test_strats.output,
        refi=rules.index_unzipped_ref.output,
        ref=rules.unzip_ref.output,
        bench_vcf=lambda w: expand_final_to_src(rules.download_bench_vcf.output, w),
        bench_bed=lambda w: read_checkpoint("normalize_bench_bed", w),
        query_vcf=lambda w: expand_final_to_src(rules.download_query_vcf.output, w),
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
    benchmark:
        post_bench_dir / "run_happy.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.runHappy
        ),
    script:
        "../scripts/python/bedtools/postprocess/run_happy.py"


rule summarize_happy:
    input:
        [
            expand(
                rules.run_happy.output,
                ref_final_key=rk,
                build_key=bk,
            )
            for rk, bk in zip(*config.all_build_keys)
            if config.to_build_data(rk, bk).want_benchmark
        ],
    output:
        validation_dir / "benchmark_summary.html",
    params:
        subsets=config.benchmark_subsets,
    conda:
        "../envs/rmarkdown.yml"
    script:
        "../scripts/rmarkdown/rmarkdown/benchmark.Rmd"


rule copy_READMEs:
    input:
        main="workflow/files/README_main.md",
        validation="workflow/files/README_validation.md",
    output:
        validation=validation_dir / "README.md",
        main=config.final_root_dir / "README.md",
    shell:
        """
        cp {input.main} {output.main}
        cp {input.validation} {output.validation}
        """


rule generate_tarballs:
    input:
        all_strats=rules.generate_tsv_list.output,
        _checksums=rules.generate_md5sums.output,
    output:
        config.final_root_dir
        / "genome-stratifications-{ref_final_key}@{build_key}.tar.gz",
    params:
        parent=lambda _, input: Path(input.all_strats[0]).parent.parent,
        target=lambda _, input: Path(input.all_strats[0]).parent.name,
    shell:
        """
        tar czf {output} -C {params.parent} {params.target}
        """


rule checksum_everything:
    input:
        rules.copy_READMEs.output,
        rules.make_coverage_plots.output,
        rules.summarize_happy.output,
        rules.all_comparisons.input,
        [
            expand(rules.generate_tarballs.output, ref_final_key=rk, build_key=bk)
            for rk, bk in zip(*config.all_full_build_keys)
        ],
    output:
        config.final_root_dir / "genome-stratifications-md5s.txt",
    params:
        root=config.final_root_dir,
    shell:
        """
        find {params.root} -type f -exec md5sum {{}} + | \
        sed 's|{params.root}|\.|' | \
        grep -v "\./genome-stratifications-md5s\.txt" | \
        sort -k2,2 > {output} && \
        cd {params.root} && \
        md5sum -c --quiet genome-stratifications-md5s.txt
        """
