from os.path import dirname
from scripts.python.config import StratOutputs


rule remove_gaps:
    input:
        genome=rules.genome_to_bed.output,
        gaps=rules.merge_gaps.output,
        parY=rules.write_Y_PAR.output,
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


# TODO need to get a bed file that has all Ns from the reference (so they can be subtracte

# TODO merge outputs? (this is something done in postprocessing but seems like
# it should be part of the validation)

# TODO need a PSA_Y_GRCh38.bed file (assuming to subtract off pseudo-autosomal regions)

# lets worry about this last, I might end up replace lots of the merge commands
# with python, which will make this easy
# rule add_strat_header:
#     input:
#         "",
#     output:
#         "",

# TODO also somehow need to generate the hap.py tables (the tsvs in the root)


# def expand_strat_targets(ref_key, build_key):
#     inc = lookup_build(["include"], ref_key, build_key)
#     chr_indices = lookup_chr_indices(ref_key, build_key)

#     # XY is special since we should not include the XY strats if we don't also
#     # have X and Y in the filter. Likewise we should not include an autosomes
#     # if they are not in the filter
#     want_xy = "X" in chr_indices and "Y" in chr_indices and inc["xy"]
#     want_autosomes = len(set(chr_indices) - set(["X", "Y"])) > 0

#     all_targets = [
#         (rules.all_low_complexity.input, inc["low_complexity"]),
#         (rules.all_xy.input, want_xy),
#         (rules.all_auto.input, want_autosomes),
#         (rules.all_map.input, inc["map"]),
#     ]

#     return chain(
#         *[
#             expand(target, allow_missing=True, ref_key=ref_key, build_key=build_key)
#             for target, wants in all_targets
#             if wants
#         ]
#     )


def expand_strat_targets_wc(wildcards):
    s = StratOutputs(
        low_complexity=rules.all_low_complexity.input,
        xy_sex=rules.all_xy.input,
        xy_auto=rules.all_auto.input,
        map=rules.all_map.input,
    )
    return config.strat_targets(wildcards.ref_key, wildcards.build_key, s)


# TODO don't hardcode version
rule generate_md5sums:
    input:
        expand_strat_targets_wc,
    output:
        final_dir / "v3.1-genome-stratifications-{ref_key}-md5s.txt",
    params:
        root=lambda wildcards, output: dirname(str(output[0])),
    shell:
        """md5sum {input} | \
        sed 's|{params.root}|{wildcards.ref_key}|' \
        > {output}
        """


rule validate_strats:
    input:
        strats=expand_strat_targets_wc,
        nonN=lambda wildcards: rules.genome_to_bed.output
        if config.refkey_to_gap_url(wildcards.ref_key) is None
        else rules.remove_gaps.output,
    output:
        final_dir / "validation.html",
    conda:
        envs_path("rmarkdown.yml")
    params:
        expected_chroms=lambda wildcards: config.buildkey_to_chr_names(
            wildcards.ref_key, wildcards.build_key
        ),
    script:
        scripts_path("rmarkdown/validate.Rmd")
