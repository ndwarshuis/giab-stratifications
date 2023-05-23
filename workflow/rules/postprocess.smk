from os.path import dirname, basename
from more_itertools import unique_everseen, unzip
from os import scandir
from common.config import CoreLevel

post_inter_dir = config.intermediate_build_dir / "postprocess"


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
        ([rules.all_gc.input.wide[0], rules.all_gc.input.narrow[0]], config.want_gc),
        (rules.all_functional.input, config.want_functional),
        (rules.all_segdups.input, config.want_segdups),
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
        post_inter_dir / "unit_tests.log",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/run_unit_tests.py"


rks, bks = map(
    list,
    unzip([(rk, bk) for rk, r in config.stratifications.items() for bk in r.builds]),
)


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
    output:
        config.final_root_dir / "coverage_plots.html",
    conda:
        "../envs/rmarkdown.yml"
    params:
        core_levels=[c.value for c in CoreLevel],
        other_levels=config.other_levels,
        chr_mapper=[
            [
                i.chr_name_full(config.refkey_to_final_chr_prefix(ref_key)),
                i.chr_name,
                f"{ref_key}@{build_key}",
            ]
            for ref_key, build_key in zip(rks, bks)
            for i in config.buildkey_to_chr_indices(ref_key, build_key)
        ],
    script:
        "../scripts/rmarkdown/rmarkdown/validate.Rmd"
