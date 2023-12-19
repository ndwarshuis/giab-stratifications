from os.path import dirname, basename
from more_itertools import unique_everseen, unzip
from os import scandir
from common.config import CoreLevel

post_inter_dir = config.intermediate_build_dir / "postprocess"
post_log_dir = config.log_build_dir / "postprocess"
validation_dir = config.final_root_dir / "validation"


def expand_strat_targets_inner(ref_final_key, build_key):
    bd = config.to_build_data(ref_final_key, build_key)
    function_targets = [
        (all_low_complexity, bd.want_low_complexity),
        (gc_inputs_flat, bd.want_gc),
        (mappabilty_inputs, bd.want_mappability),
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
    assert len(invalid) == 0, f"invalid targets: {invalid}"
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
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/run_unit_tests.py"


# TODO this is totally incorrect
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
    localrule: True
    shell:
        """
        curl -fsSq -o {params.tar} {params.url} && \
        mkdir {output} && \
        tar xzf {params.tar} --directory {output} --strip-components=1
        """


rule compare_strats:
    input:
        _test=expand(
            rules.unit_test_strats.output,
            zip,
            ref_final_key=rks,
            build_key=bks,
        ),
        old=lambda w: expand(
            rules.download_comparison_strat_tarball.output,
            compare_key=config.buildkey_to_comparekey(w.ref_final_key, w.build_key),
        )[0],
        # use this to target a specific rule to satisfy the snakemake scheduler,
        # the thing I actually need here is the parent directory
        new_list=rules.generate_tsv_list.output[0],
    output:
        validation_dir / "{ref_final_key}@{build_key}" / "diagnostics.tsv",
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
            expand(rules.compare_strats.output, ref_final_key=rk, build_key=bk)[0]
            for rk, r in config.stratifications.items()
            for bk in r.builds
            if config.buildkey_to_comparekey(rk, bk) is not None
        ],
    localrule: True


# Rule to make a dataframe which maps chromosome names to a common nomenclature
# (for the validation markdown). The only reason this exists is because
# snakemake flattens this list of tuples when I pass it as a param to the rmd
# script...dataframe in IO monad it is :/
rule write_chr_name_mapper:
    output:
        config.intermediate_root_dir / ".validation" / "chr_mapper.tsv",
    localrule: True
    run:
        with open(output[0], "w") as f:
            for ref_final_key, build_key in zip(rks, bks):
                for i in config.buildkey_to_chr_indices(ref_final_key, build_key):
                    pattern = config.refkey_to_final_chr_pattern(ref_final_key)
                    line = [
                        str(i.value),
                        f"{ref_final_key}@{build_key}",
                        i.chr_name_full(pattern),
                    ]
                    f.write("\t".join(line) + "\n")


rule validate_strats:
    input:
        # this first input isn't actually used, but ensures the unit tests pass
        # before running the rmd script
        _test=expand(
            rules.unit_test_strats.output,
            zip,
            ref_final_key=rks,
            build_key=bks,
        ),
        strats=expand(
            rules.list_all_strats.output,
            zip,
            ref_final_key=rks,
            build_key=bks,
        ),
        nonN=expand(
            rules.get_gapless.output.auto,
            zip,
            ref_final_key=rks,
            build_key=bks,
        ),
        chr_mapper=rules.write_chr_name_mapper.output,
    output:
        validation_dir / "coverage_plots.html",
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
        _test=expand(
            rules.unit_test_strats.output,
            zip,
            ref_final_key=rks,
            build_key=bks,
        ),
        refi=rules.index_unzipped_ref.output,
        ref=rules.unzip_ref.output,
        bench_vcf=rules.download_bench_vcf.output,
        bench_bed=rules.filter_sort_bench_bed.output,
        query_vcf=rules.download_query_vcf.output,
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
                ref_final_key=rk,
                build_key=bk,
            )
            for rk, bk in zip(rks, bks)
            if config.want_benchmark(rk, bk)
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


# TODO expand rule is wrong
rule checksum_everything:
    input:
        rules.copy_READMEs.output,
        rules.validate_strats.output,
        rules.summarize_happy.output,
        rules.all_comparisons.input,
        [
            expand(rules.generate_tarballs.output, ref_final_key=rk, build_key=bk)
            for rk, r in config.stratifications.items()
            for bk in r.builds
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
