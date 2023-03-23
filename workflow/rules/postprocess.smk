from os.path import dirname, basename
from more_itertools import unique_everseen
from os import scandir

# TODO merge outputs? (this is something done in postprocessing but seems like
# it should be part of the validation)

# TODO need to get a bed file that has all Ns from the reference (so they can be subtracted)
# TODO need a PSA_Y_GRCh38.bed file (assuming to subtract off pseudo-autosomal regions)
# TODO also somehow need to generate the hap.py tables (the tsvs in the root)
# TODO add a nice header to the top informing user that "these are strats"?

post_inter_dir = intermediate_dir / "postprocess"


rule remove_gaps:
    input:
        genome=rules.genome_to_bed.output,
        gaps=rules.merge_gaps.output,
        parY=expand(rules.write_PAR_intermediate.output, allow_missing=True, chr="Y"),
    output:
        ref_inter_dir / "genome_gapless.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        cat {input.genome} | \
        bedtools subtract -a stdin -b {input.gaps} | \
        bedtools subtract -a stdin -b {input.parY} \
        > {output}
        """


# NOTE This will only give a limited set of "targets" from each of the strat
# level types. Because of the way snakemake works, it is not practical or
# maintainable to specify every single stratification file; therefore we only
# use the "toplevel" targets which will pull in all others. In downstream rules
# from this, it is much easier to use condense this to a list of parent
# directories then manually iterate all strat files in these directories.
def expand_strat_targets(wildcards):
    rk = wildcards.ref_key
    bk = wildcards.build_key
    targets = [
        (rules.all_low_complexity.input, config.want_low_complexity),
        (rules.all_xy_auto.input, config.want_xy_auto),
        (rules.all_map.input, config.want_map),
        (rules.all_gc.input, config.want_gc),
        (rules.all_functional.input, config.want_functional),
        (rules.all_segdups.input, config.want_segdups),
    ]
    auto = [t for tgt, test in targets if test(rk, bk) for t in tgt]
    sex = expand(
        rules.all_xy_sex.input,
        allow_missing=True,
        chr=config.wanted_xy_chr_names(rk, bk),
    )
    all_targets = auto + sex
    invalid = [f for f in all_targets if not f.startswith("results/builds/final")]
    assert len(invalid) == 0, f"invalid targets: {invalid}"
    return all_targets


def common_dirs(files):
    return list(unique_everseen([dirname(str(f)) for f in files]))


# TODO don't hardcode version
rule generate_md5sums:
    input:
        expand_strat_targets,
    output:
        final_dir / "v3.1-genome-stratifications-{ref_key}-md5s.txt",
    params:
        root=lambda wildcards, output: dirname(str(output[0])),
        all_strats=lambda _, input: [
            p.path
            for i in common_dirs(input)
            for p in scandir(i)
            if p.path.endswith(".bed.gz")
        ],
    shell:
        """
        md5sum {params.all_strats} | \
        sed 's|{params.root}|{wildcards.ref_key}|' \
        > {output}
        """


rule unit_test_strats:
    input:
        expand_strat_targets,
    output:
        touch(post_inter_dir / "verify.done"),
    log:
        post_inter_dir / "verify_strats.log",
    conda:
        envs_path("bedtools.yml")
    params:
        strat_dirs=lambda _, input: common_dirs(input),
    script:
        scripts_path("python/bedtools/verify.py")


rule validate_strats:
    input:
        # this input isn't actually used, but ensures the unit tests pass
        # before running the rmd script
        _test=rules.unit_test_strats.output,
        targets=expand_strat_targets,
        nonN=lambda w: rules.genome_to_bed.output
        if config.refkey_to_gap_src(w.ref_key) is None
        else rules.remove_gaps.output,
    output:
        final_dir / "validation.html",
    conda:
        envs_path("rmarkdown.yml")
    params:
        expected_chroms=lambda wildcards: config.buildkey_to_chr_names(
            wildcards.ref_key, wildcards.build_key
        ),
        strat_dirs=lambda _, input: common_dirs(input.targets),
    script:
        scripts_path("rmarkdown/validate.Rmd")
