from more_itertools import unzip, flatten
from collections import namedtuple
from common.config import (
    CoreLevel,
    si_to_rmsk,
    si_to_simreps,
    si_to_satellites,
    strip_full_refkey,
)
from functools import partial

lc = config.to_bed_dirs(CoreLevel.LOWCOMPLEXITY)

################################################################################
## uniform repeats

URep = namedtuple("URep", ["unit_len", "range_indices", "total_lens"])

uniform_repeats = {
    "homopolymer": URep(1, 3, [4, 7, 12, 21]),
    "diTR": URep(2, 3, [10, 50, 150]),
    "triTR": URep(3, 3, [14, 50, 150]),
    "quadTR": URep(4, 3, [19, 50, 150]),
}

unit_name_constraint = f"({'|'.join(uniform_repeats)})"

bases_constraint = "[ATGC]{2}"

COMPLEMENTS = ["AT", "GC"]
IMPERFECT_LENS = [11, 21]


# NOTE: weird sed command ensures all bases are capitalized, otherwise for some
# refs we might miss a repeat like "AAaa"
rule find_perfect_uniform_repeats:
    input:
        ref=rules.filter_sort_ref.output,
        bin=rules.build_repseq.output,
    output:
        lc.inter.postsort.data / "uniform_repeats_R{unit_len}_T{total_len}.bed",
    log:
        lc.inter.postsort.log / "uniform_repeats_R{unit_len}_T{total_len}.log",
    wildcard_constraints:
        unit_len="\d+",
        total_len="\d+",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input.ref} | \
        sed 's/^[A-Za-z]\\+$/\\U&/' | \
        {input.bin} {wildcards.unit_len} {wildcards.total_len} - 2> {log} \
        > {output} 
        """


rule filter_perfect_uniform_repeats:
    input:
        rules.find_perfect_uniform_repeats.output,
    output:
        lc.inter.postsort.data / "uniform_repeats_{bases}_R{unit_len}_T{total_len}.bed",
    log:
        lc.inter.postsort.log / "uniform_repeats_{bases}_R{unit_len}_T{total_len}.log",
    wildcard_constraints:
        unit_len="\d+",
        total_len="\d+",
        bases=bases_constraint,
    shell:
        "sed -n '/unit=[{wildcards.bases}]\+/p' {input} > {output}"


# bundle all these together in a dummy rule so I can access them later
rule all_perfect_uniform_repeats:
    input:
        **{
            f"R{u.unit_len}_T{t}": expand(
                rules.find_perfect_uniform_repeats.output,
                allow_missing=True,
                unit_len=u.unit_len,
                total_len=t,
            )
            for u in uniform_repeats.values()
            for t in u.total_lens
        },
        **{
            f"R{u.unit_len}_T{t}_{bs}": expand(
                rules.filter_perfect_uniform_repeats.output,
                allow_missing=True,
                unit_len=u.unit_len,
                total_len=t,
                bases=bs,
            )
            for u in uniform_repeats.values()
            for t in u.total_lens
            for bs in COMPLEMENTS
        },
    localrule: True


def lookup_perfect_uniform_repeat(unit_name, total_len):
    ul = uniform_repeats[unit_name].unit_len
    return rules.all_perfect_uniform_repeats.input[f"R{ul}_T{total_len}"]


def lookup_perfect_uniform_repeat_complement(unit_name, total_len, bases):
    ul = uniform_repeats[unit_name].unit_len
    return rules.all_perfect_uniform_repeats.input[f"R{ul}_T{total_len}_{bases}"]


def repeat_range_inputs(wildcards):
    return {
        k: lookup_perfect_uniform_repeat(wildcards.unit_name, t)
        for k, t in zip(
            "ab",
            [wildcards.total_lenA, wildcards.total_lenB],
        )
    }


rule subtract_uniform_repeats:
    input:
        unpack(repeat_range_inputs),
    output:
        lc.inter.postsort.data
        / "uniform_repeat_range_{unit_name}_{total_lenA}to{total_lenB}.bed",
    conda:
        "../envs/bedtools.yml"
    benchmark:
        lc.inter.postsort.bench / "subtract_uniform_repeats_{unit_name}_{total_lenA}to{total_lenB}.txt"
    wildcard_constraints:
        unit_name=unit_name_constraint,
        total_lenA="\d+",
        total_lenb="\d+",
    shell:
        "subtractBed -a {input.a} -b {input.b} -sorted > {output}"


def repeat_range_complement_inputs(wildcards):
    return {
        k: lookup_perfect_uniform_repeat_complement(
            wildcards.unit_name, t, wildcards.bases
        )
        for k, t in zip(
            "ab",
            [wildcards.total_lenA, wildcards.total_lenB],
        )
    }


rule subtract_uniform_repeat_complement:
    input:
        unpack(repeat_range_complement_inputs),
    output:
        lc.inter.postsort.data
        / "uniform_repeat_range_{bases}_{unit_name}_{total_lenA}to{total_lenB}.bed",
    conda:
        "../envs/bedtools.yml"
    benchmark:
        lc.inter.postsort.bench / "subtract_uniform_repeats_{bases}_{unit_name}_{total_lenA}to{total_lenB}.txt"
    wildcard_constraints:
        unit_name=unit_name_constraint,
        total_lenA="\d+",
        total_lenb="\d+",
        bases=bases_constraint,
    shell:
        "subtractBed -a {input.a} -b {input.b} -sorted > {output}"


rule slop_uniform_repeats:
    input:
        lambda w: lookup_perfect_uniform_repeat(w.unit_name, int(w.total_len)),
    output:
        lc.final("SimpleRepeat_{unit_name}_ge{total_len}_slop5"),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        unit_name=unit_name_constraint,
        total_len="\d+",
    params:
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    shell:
        """
        slopBed -i {input} -b 5 -g {params.genome} | \
        cut -f1-3 | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {params.gapless} -sorted -g {params.genome} | \
        bgzip -c \
        > {output}
        """


# +1 to lenB since the filename is [X, Y] and not [X, Y)
use rule slop_uniform_repeats as slop_uniform_repeat_ranges with:
    input:
        lambda w: expand(
            rules.subtract_uniform_repeats.output,
            allow_missing=True,
            unit_name=w.unit_name,
            total_lenA=int(w.total_lenA),
            total_lenB=int(w.total_lenB) + 1,
        ),
    output:
        lc.final("SimpleRepeat_{unit_name}_{total_lenA}to{total_lenB}_slop5"),
    wildcard_constraints:
        unit_name=unit_name_constraint,
        total_lenA="\d+",
        total_lenB="\d+",


use rule slop_uniform_repeats as slop_uniform_repeats_complement with:
    input:
        lambda w: lookup_perfect_uniform_repeat_complement(
            w.unit_name,
            int(w.total_len),
            w.bases,
        ),
    output:
        lc.final("SimpleRepeat_{unit_name}_ge{total_len}_{bases}_slop5"),
    wildcard_constraints:
        unit_name=unit_name_constraint,
        total_len="\d+",
        bases=bases_constraint,


use rule slop_uniform_repeats as slop_uniform_repeat_ranges_complement with:
    input:
        lambda w: expand(
            rules.subtract_uniform_repeat_complement.output,
            allow_missing=True,
            unit_name=w.unit_name,
            total_lenA=int(w.total_lenA),
            total_lenB=int(w.total_lenB) + 1,
            bases=w.bases,
        ),
    output:
        lc.final("SimpleRepeat_{unit_name}_{total_lenA}to{total_lenB}_{bases}_slop5"),
    wildcard_constraints:
        unit_name=unit_name_constraint,
        total_lenA="\d+",
        total_lenB="\d+",
        bases=bases_constraint,


rule merge_perfect_uniform_repeats:
    input:
        rules.all_perfect_uniform_repeats.input.R1_T4,
    output:
        lc.inter.postsort.data / "repeats_imp_hp_T{merged_len}_B{base}.bed",
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        merged_len="\d+",
    shell:
        """
        grep 'unit={wildcards.base}' {input} | \
        mergeBed -i stdin -d 1 | \
        awk '$3-$2>={wildcards.merged_len}' > \
        {output}
        """


rule merge_imperfect_uniform_repeats:
    input:
        expand(
            rules.merge_perfect_uniform_repeats.output,
            allow_missing=True,
            base=["A", "C", "G", "T"],
        ),
    output:
        lc.final("SimpleRepeat_imperfecthomopolge{merged_len}_slop5"),
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        merged_len="\d+",
    params:
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    shell:
        """
        multiIntersectBed -i {input} | \
        slopBed -i stdin -b 5 -g {params.genome} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {params.gapless} -sorted | \
        bgzip -c > {output}
        """


use rule merge_imperfect_uniform_repeats as merge_imperfect_uniform_repeats_complement with:
    input:
        lambda w: expand(
            rules.merge_perfect_uniform_repeats.output,
            allow_missing=True,
            base=list(w.bases),
        ),
    output:
        lc.final("SimpleRepeat_imperfecthomopolge{merged_len}_{bases}_slop5"),
    wildcard_constraints:
        merged_len="\d+",
        bases=bases_constraint,


rule all_uniform_repeats:
    input:
        # Perfect (greater than X)
        *[
            p
            for k, v in uniform_repeats.items()
            for x in v.total_lens[v.range_indices - 1 :]
            for p in expand(
                rules.slop_uniform_repeats.output,
                allow_missing=True,
                unit_name=k,
                total_len=x,
            )
            + (
                expand(
                    rules.slop_uniform_repeats_complement.output,
                    allow_missing=True,
                    unit_name=k,
                    total_len=x,
                    bases=COMPLEMENTS,
                )
                if k == "homopolymer"
                else []
            )
        ],
        # Perfect (between X and Y)
        *[
            p
            for k, v in uniform_repeats.items()
            for a, b in zip(
                v.total_lens[0 : v.range_indices - 1],
                v.total_lens[1 : v.range_indices],
            )
            for p in expand(
                rules.slop_uniform_repeat_ranges.output,
                allow_missing=True,
                unit_name=k,
                total_lenA=a,
                total_lenB=b - 1,
            )
            + (
                expand(
                    rules.slop_uniform_repeat_ranges_complement.output,
                    allow_missing=True,
                    unit_name=k,
                    total_lenA=a,
                    total_lenB=b - 1,
                    bases=COMPLEMENTS,
                )
                if k == "homopolymer"
                else []
            )
        ],
        # Imperfect (greater than X)
        expand(
            rules.merge_imperfect_uniform_repeats_complement.output,
            allow_missing=True,
            merged_len=IMPERFECT_LENS,
            bases=COMPLEMENTS,
        ),
        **{
            f"imperfect_ge{x}": expand(
                rules.merge_imperfect_uniform_repeats.output,
                allow_missing=True,
                merged_len=x,
            )
            for x in IMPERFECT_LENS
        },
    localrule: True


################################################################################
## simple repeats


use rule download_ref as download_simreps with:
    output:
        lc.src.data / "simreps.txt.gz",
    log:
        lc.src.log / "simreps.log",
    params:
        src=lambda w: to_bed_src(si_to_simreps, w),
    localrule: True


checkpoint normalize_simreps:
    input:
        lambda w: bed_src_inputs(rules.download_simreps.output, si_to_simreps, w),
    output:
        lc.inter.filtersort.data / "simreps.json",
    params:
        output_pattern=lambda w: to_output_pattern(lc, "simreps", w),
    conda:
        "../envs/bedtools.yml"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeSimreps
        ),
    benchmark:
        lc.inter.filtersort.bench / "normalize_simreps.txt"
    script:
        "../scripts/python/bedtools/low_complexity/normalize_simreps.py"


rule merge_simreps:
    input:
        partial(read_checkpoint, "normalize_simreps"),
    output:
        lc.inter.postsort.data / "simreps_merged.bed.gz",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


################################################################################
## rmsk


use rule download_ref as download_rmsk with:
    output:
        lc.src.data / "rmsk.txt.gz",
    log:
        lc.src.log / "rmsk.log",
    params:
        src=lambda w: to_bed_src(si_to_rmsk, w),
    localrule: True


checkpoint normalize_rmsk:
    input:
        lambda w: bed_src_inputs(rules.download_rmsk.output, si_to_rmsk, w),
    output:
        lc.inter.filtersort.data / "rmsk.txt.gz",
    params:
        output_pattern=lambda w: to_output_pattern(lc, "rmsk", w),
    conda:
        "../envs/bedtools.yml"
    benchmark:
        lc.inter.filtersort.bench / "normalize_rmsk.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeRmsk
        ),
    script:
        "../scripts/python/bedtools/low_complexity/normalize_rmsk.py"


rule merge_rmsk_class:
    input:
        lambda w: read_checkpoint("normalize_rmsk", w),
    output:
        lc.inter.postsort.data / "rmsk_class_{rmsk_class}.bed.gz",
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        rmsk_class="\w+",
    shell:
        """
        gunzip -c {input} | \
        grep {wildcards.rmsk_class} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


all_rmsk_classes = {
    c: expand(
        rules.merge_rmsk_class.output,
        allow_missing=True,
        rmsk_class=c,
    )
    for c in ["Low_complexity", "Simple_repeat", "Satellite"]
}


################################################################################
## Satellites (censat, alternative to RMSK as in above)


use rule download_ref as download_censat with:
    output:
        lc.src.data / "censat.txt.gz",
    log:
        lc.src.log / "censat.log",
    params:
        src=lambda w: to_bed_src(si_to_satellites, w),
    localrule: True


# TODO add config directives to override the mem_mb in each of these since the
# input sizes can vary wildly depending on how many columns/rows/datatypes are
# actually in them. 2-4k is probably a safe default
checkpoint normalize_censat:
    input:
        lambda w: bed_src_inputs(rules.download_censat.output, si_to_satellites, w),
    output:
        lc.inter.filtersort.data / "censat.json",
    params:
        output_pattern=lambda w: to_output_pattern(lc, "censat", w),
    benchmark:
        lc.inter.filtersort.bench / "normalize_censat.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_key, w.build_key, lambda m: m.normalizeCensat
        ),
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/low_complexity/normalize_censat.py"


# split this from the final rule since other rules depend on the satellites
# bed file but add slop of their own, so this avoids adding slop twice
rule merge_satellites_intermediate:
    input:
        lambda w: read_checkpoint("normalize_censat", w)
        if config.to_build_data(
            strip_full_refkey(w.ref_final_key), w.build_key
        ).refdata.has_low_complexity_censat
        else all_rmsk_classes["Satellite"],
    output:
        lc.inter.postsort.data / "merged_satellites.bed.gz",
    conda:
        "../envs/bedtools.yml"
    shell:
        "mergeBed -i {input} | bgzip -c > {output}"


rule merge_satellites:
    input:
        bed=rules.merge_satellites_intermediate.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        lc.final("satellites_slop5"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        slopBed -i {input.bed} -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule invert_satellites:
    input:
        bed=rules.merge_satellites.output,
    output:
        lc.final("notinsatellites_slop5"),
    conda:
        "../envs/bedtools.yml"
    # this is a nice trick to avoid specifying input files for rule overrides
    # when they never change
    params:
        gapless=rules.get_gapless.output.auto,
        genome=rules.get_genome.output,
    shell:
        """
        complementBed -i {input.bed} -g {params.genome} |
        intersectBed -a stdin -b {params.gapless} -sorted | \
        bgzip -c > {output}
        """


################################################################################
## Tandem Repeats


rule merge_all_uniform_repeats:
    input:
        imperfect=rules.all_uniform_repeats.input.imperfect_ge11,
        perfect=rules.all_perfect_uniform_repeats.input.R1_T7,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        lc.final("AllHomopolymers_ge7bp_imperfectge11bp_slop5"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        slopBed -i {input.perfect} -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        multiIntersectBed -i stdin {input.imperfect} | \
        mergeBed -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c \
        > {output}
        """


use rule invert_satellites as invert_all_uniform_repeats with:
    input:
        bed=rules.merge_all_uniform_repeats.output,
    output:
        lc.final("notinAllHomopolymers_ge7bp_imperfectge11bp_slop5"),


rule merge_repeats:
    input:
        beds=all_rmsk_classes["Low_complexity"]
        + all_rmsk_classes["Simple_repeat"]
        + rules.merge_simreps.output
        + rules.all_perfect_uniform_repeats.input.R2_T10
        + rules.all_perfect_uniform_repeats.input.R3_T14
        + rules.all_perfect_uniform_repeats.input.R4_T19
        + rules.merge_satellites_intermediate.output,
        genome=rules.get_genome.output,
    output:
        lc.inter.postsort.data / "AllTandemRepeats_intermediate.bed",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input.beds} | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        mergeBed -i stdin > {output}
        """


# NOTE: this is pre-slop
tr_bounds = {
    "le50": {"lower": 0, "upper": 51},
    "51to200": {"lower": 51, "upper": 201},
    "201to10000": {"lower": 201, "upper": 10001},
    "ge10001": {"lower": 10001, "upper": 1e10},  # NOTE 1e10 ~ Inf
    "ge101": {"lower": 101, "upper": 1e10},
}


rule filter_TRs:
    input:
        tr=rules.merge_repeats.output,
        hp=rules.merge_all_uniform_repeats.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        lc.final("AllTandemRepeats_{tr_bound}bp_slop5"),
    params:
        lower=lambda w: tr_bounds[w.tr_bound]["lower"],
        upper=lambda w: tr_bounds[w.tr_bound]["upper"],
    conda:
        "../envs/bedtools.yml"
    wildcard_constraints:
        tr_bound=f"({'|'.join(tr_bounds)})",
    # NOTE +10 since this is processing a bed file that had 5bp slop added
    shell:
        """
        awk '({params.lower}+10)<=($3-$2) && ($3-$2)<({params.upper}+10)' {input.tr} | \
        subtractBed -a stdin -b {input.hp} -sorted | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule all_TRs:
    # funny prefix since snakemake doesn't like keys that start with digits
    input:
        **{
            f"_{k}": expand(
                rules.filter_TRs.output,
                allow_missing=True,
                tr_bound=k,
            )
            for k in tr_bounds
        },
    localrule: True


rule merge_filtered_TRs:
    input:
        rules.all_TRs.input._le50,
        rules.all_TRs.input._51to200,
        rules.all_TRs.input._201to10000,
        rules.all_TRs.input._ge10001,
    output:
        lc.final("AllTandemRepeats"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


use rule invert_satellites as invert_TRs with:
    input:
        bed=rules.merge_filtered_TRs.output,
    output:
        lc.final("notinallTandemRepeats"),


################################################################################
## Combine all the beds to make a Pink Floyd album cover


use rule merge_filtered_TRs as merge_HPs_and_TRs with:
    input:
        rules.merge_filtered_TRs.input,
        rules.merge_all_uniform_repeats.output,
    output:
        lc.final("AllTandemRepeatsandHomopolymers_slop5"),


use rule invert_satellites as invert_HPs_and_TRs with:
    input:
        bed=rules.merge_HPs_and_TRs.output,
    output:
        lc.final("notinAllTandemRepeatsandHomopolymers_slop5"),


# rule all_low_complexity:
#     input:
#         # Uniform repeats
#         rules.all_uniform_repeats.input,
#         rules.merge_all_uniform_repeats.output,
#         rules.invert_all_uniform_repeats.output,
#         # Satellites
#         rules.merge_satellites.output,
#         rules.invert_satellites.output,
#         # Tandem Repeats
#         rules.all_TRs.input,
#         rules.merge_filtered_TRs.output,
#         rules.invert_TRs.output,
#         # "Everything" (in theory)
#         rules.merge_HPs_and_TRs.output,
#         rules.invert_HPs_and_TRs.output,
#     localrule: True


def all_low_complexity(ref_final_key, _):
    rd = config.to_ref_data(strip_full_refkey(ref_final_key))
    rmsk = rd.has_low_complexity_rmsk
    simreps = rd.has_low_complexity_simreps
    censat = rd.has_low_complexity_censat
    has_sats = rmsk or censat

    # include uniform repeats no matter what
    urs = (
        rules.all_uniform_repeats.input
        + rules.merge_all_uniform_repeats.output
        + rules.invert_all_uniform_repeats.output
    )

    # include satellites only if we have rmsk or censat
    sats = (
        rules.merge_satellites.output + rules.invert_satellites.output
        if has_sats
        else []
    )

    # include tandem repeats and merged output if we have rmsk/censat and simreps
    trs = (
        rules.all_TRs.input + rules.merge_filtered_TRs.output + rules.invert_TRs.output
    )
    merged = rules.merge_HPs_and_TRs.output + rules.invert_HPs_and_TRs.output
    all_trs_and_hps = trs + merged if has_sats and simreps else []

    return all_trs_and_hps + sats + urs
