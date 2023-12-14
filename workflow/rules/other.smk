from more_itertools import unzip

# make sure the wildcards here can match everything except the "built-in" other-
# difficult beds (listed here)
builtin_other_targets = ["gaps_slop15kb", "VDJ"]

other_constraints = {
    "other_level_key": f"({'|'.join(config.other_levels)})",
    "other_strat_key": f"(?!({'|'.join(builtin_other_targets)}))[A-Za-z0-9-._]+",
}


use rule download_ref as download_other with:
    output:
        config.intermediate_build_dir
        / "other"
        / "{other_level_key}_{other_strat_key}.bed.gz",
    log:
        config.log_build_dir / "existing" / "{other_level_key}_{other_strat_key}.log",
    params:
        src=lambda w: config.buildkey_to_bed_src(
            lambda bd: bd_to_other(w.other_level_key, w.other_strat_key, bd),
            w.ref_src_key,
            w.build_key,
        ),
    wildcard_constraints:
        **other_constraints,
    localrule: True


rule filter_sort_other:
    input:
        lambda w: expand(
            rules.download_other.output,
            ref_src_key=config.buildkey_to_bed_refsrckeys(
                lambda bd: bd_to_other(w.other_level_key, w.other_strat_key, bd),
                w.ref_key,
                w.build_key,
            ),
        ),
    output:
        self.intermediate_build_hapless_dir
        / "{other_level_key}"
        / "{other_strat_key}.json",
    params:
        output_pattern=lambda w: expand(
            self.intermediate_build_dir
            / "{other_level_key}"
            / "{other_strat_key}.bed.gz",
            ref_final_key="%s",
            build_key=w.build_key,
            other_level_key=w.other_level_key,
            other_strat_key=w.other_strat_key,
        ),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        **other_constraints,
    script:
        "../scripts/python/bedtools/other/filter_sort_other.py"


rule remove_gaps_other:
    input:
        bed=lambda w: read_checkpoint("filter_sort_other", w),
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        config.build_final_strat_path("{other_level_key}", "{other_strat_key}"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -a {input.bed} -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


def all_other(ref_final_key, build_key):
    other = config.buildkey_to_build(ref_final_key, build_key).other_strats
    return [
        expand(
            rules.remove_gaps_other.output,
            ref_final_key=ref_final_key,
            build_key=build_key,
            other_level_key=lk,
            other_strat_key=sk,
        )[0]
        for lk, s in other.items()
        for sk in s
    ]
