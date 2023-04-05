from os.path import dirname, basename
from more_itertools import unique_everseen
from os import scandir

# TODO add a nice header to the top informing user that "these are strats"?

post_inter_dir = config.intermediate_build_dir / "postprocess"


# NOTE This will only give a limited set of "targets" from each of the strat
# level types. Because of the way snakemake works, it is not practical or
# maintainable to specify every single stratification file; therefore we only
# use the "toplevel" targets which will pull in all others. In downstream rules
# from this, it is much easier to use condense this to a list of parent
# directories then manually iterate all strat files in these directories.
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

    # xy (expand target depending on which chromosomes we have selected)
    sex = all_xy_features(wildcards) + all_xy_PAR(wildcards)

    # combine and ensure that all "targets" refer to final bed files
    all_targets = auto + sex
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


# TODO don't hardcode version
rule generate_md5sums:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "v3.1-genome-stratifications-{ref_key}-md5s.txt",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/postprocess/list_md5.py"


rule generate_tsv_list:
    input:
        rules.list_all_strats.output,
    output:
        config.final_build_dir / "v3.1-{ref_key}-all-stratifications.tsv",
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


rule validate_strats:
    input:
        # this input isn't actually used, but ensures the unit tests pass
        # before running the rmd script
        _test=rules.unit_test_strats.output,
        strats=rules.list_all_strats.output,
        nonN=rules.get_gapless.output.auto,
    output:
        config.final_build_dir / "validation.html",
    conda:
        "../envs/rmarkdown.yml"
    script:
        "../scripts/rmarkdown/rmarkdown/validate.Rmd"
