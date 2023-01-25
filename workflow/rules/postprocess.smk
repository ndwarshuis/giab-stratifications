from os.path import dirname
from itertools import chain


rule remove_genome_gaps:
    input:
        gaps=rules.merge_gaps.output,
        genome=rules.get_genome.output,
    output:
        ref_inter_dir / "genome_noN.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        bedtools complement -i {input.gaps} -g {input.genome} > {output}
        """


rule remove_gaps_PAR_Y:
    input:
        gaps=rules.remove_genome_gaps.output,
        par=rules.write_Y_PAR.output,
    output:
        ref_inter_dir / "genome_noN_noPARY.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        bedtools subtract -a {input.gaps} -b {input.par} > {output}
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


def expand_strat_targets(ref_key, build_key):
    inc = lookup_build(["include"], ref_key, build_key)
    chr_indices = lookup_chr_indices(ref_key, build_key)

    # XY is special since we should not include the XY strats if we don't also
    # have X and Y in the filter. Likewise we should not include an autosomes
    # if they are not in the filter
    want_xy = "X" in chr_indices and "Y" in chr_indices and inc["xy"]
    want_autosomes = len(set(chr_indices) - set(["X", "Y"])) > 0

    all_targets = [
        (rules.all_low_complexity.input, inc["low_complexity"]),
        (rules.all_xy.input, want_xy),
        (rules.all_auto.input, want_autosomes),
        (rules.all_map.input, inc["map"]),
    ]

    return chain(
        *[
            expand(target, allow_missing=True, ref_key=ref_key, build_key=build_key)
            for target, wants in all_targets
            if wants
        ]
    )


def expand_strat_targets_wc(wildcards):
    return expand_strat_targets(wildcards.ref_key, wildcards.build_key)


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
        gaps=rules.remove_gaps_PAR_Y.output,
    output:
        final_dir / "validation.html",
    conda:
        envs_path("rmarkdown.yml")
    params:
        expected_chroms=lookup_chrs_wc,
    script:
        scripts_path("rmarkdown/validate.Rmd")
