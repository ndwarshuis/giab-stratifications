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
        src=lambda w: config.otherkey_to_src(
            w.ref_key, w.build_key, w.other_level_key, w.other_strat_key
        ),
    wildcard_constraints:
        **other_constraints,
    localrule: True


rule filter_sort_other:
    input:
        bed=rules.download_other.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        config.build_final_strat_path("{other_level_key}", "{other_strat_key}"),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        **other_constraints,
    script:
        "../scripts/python/bedtools/other/filter_sort_other.py"


def all_other(ref_key, build_key):
    other = config.buildkey_to_build(ref_key, build_key).other_strats
    return [
        expand(
            rules.filter_sort_other.output,
            ref_key=ref_key,
            build_key=build_key,
            other_level_key=lk,
            other_strat_key=sk,
        )[0]
        for lk, s in other.items()
        for sk in s
    ]
