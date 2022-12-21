from more_itertools import unzip
from collections import namedtuple
from functools import partial


lc_inter_dir = intermediate_dir / "LowComplexity"
lc_final_dir = final_dir / "LowComplexity"
lc_log_dir = log_dir / "LowComplexity"


rule download_ref:
    output:
        resources_dir / "ref.fna.gz",
    params:
        url=config["ref_url"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unzip_ref:
    input:
        rules.download_ref.output,
    output:
        resources_dir / "ref.fna",
    conda:
        envs_path("utils.yml")
    shell:
        "gunzip -c {input} > {output}"


rule index_ref:
    input:
        rules.download_ref.output,
    output:
        ref_dir / "ref.fna.fai",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        samtools faidx - -o - | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        sed 's/^chr//' | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n -k2,2n -k3,3n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' > \
        {output}
        """


rule get_genome:
    input:
        rules.index_ref.output,
    output:
        ref_dir / "genome.txt",
    shell:
        "cut -f 1,2 {input} > {output}"


rule filter_sort_ref:
    input:
        fa=rules.unzip_ref.output,
        genome=rules.get_genome.output,
    output:
        ref_dir / "ref_filtered.fa",
    conda:
        envs_path("utils.yml")
    shell:
        """
        samtools faidx {input.fa} $(cut -f1 {input.genome} | tr '\n' ' ') > \
        {output}
        """


################################################################################
## uniform repeats

URep = namedtuple("URep", ["unit_len", "range_indices", "total_lens"])

uniform_repeats = {
    "homopolymer": URep(1, 3, [4, 7, 12, 21]),
    "diTR": URep(2, 3, [10, 50, 200]),
    "triTR": URep(3, 3, [14, 50, 200]),
    "quadTR": URep(4, 3, [19, 50, 200]),
}


rule find_perfect_uniform_repeats:
    input:
        ref=rules.filter_sort_ref.output,
        bin=rules.build_repseq.output,
    output:
        lc_inter_dir / "uniform_repeats_R{unit_len}_L{total_len}.bed",
    log:
        lc_log_dir / "uniform_repeats_R{unit_len}_L{total_len}.log",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        {input.bin} {wildcards.unit_len} {wildcards.total_len} {input.ref} \
        > {output} 2> {log}
        """


# bundle all these together in a dummy rule so I can access them later
rule all_perfect_uniform_repeats:
    input:
        **{
            f"R{u.unit_len}_T{t}": expand(
                rules.find_perfect_uniform_repeats.output,
                unit_len=u.unit_len,
                total_len=t,
            )
            for u in uniform_repeats.values()
            for t in u.total_lens
        },


def lookup_perfect_uniform_repeat(unit_name, total_len):
    ul = uniform_repeats[unit_name].unit_len
    return rules.all_perfect_uniform_repeats.input[f"R{ul}_T{total_len}"]


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
        lc_inter_dir / "uniform_repeat_range_{unit_name}_{total_lenA}to{total_lenB}.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        "subtractBed -a {input.a} -b {input.b} > {output}"


def get_offset(unit_name):
    return 0 if unit_name == "homopolymer" else 1


rule slop_uniform_repeats:
    input:
        bed=lambda wildcards: lookup_perfect_uniform_repeat(
            wildcards.unit_name,
            int(wildcards.total_len) + 1 - get_offset(wildcards.unit_name),
        ),
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_SimpleRepeat_{unit_name}_gt{total_len}_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        slopBed -i {input.bed} -b 5 -g {input.genome} | \
        cut -f1-3 | \
        bgzip -c \
        > {output}
        """


use rule slop_uniform_repeats as slop_uniform_repeat_ranges with:
    input:
        bed=lambda wildcards: expand(
            rules.subtract_uniform_repeats.output,
            unit_name=wildcards.unit_name,
            total_lenA=int(wildcards.total_lenA)
            - (offset := get_offset(wildcards.unit_name)),
            total_lenB=int(wildcards.total_lenB) - offset + 1,
        ),
        genome=rules.get_genome.output,
    output:
        lc_final_dir
        / "GRCh38_SimpleRepeat_{unit_name}_{total_lenA}to{total_lenB}_slop5.bed.gz",


rule merge_perfect_uniform_repeats:
    input:
        rules.all_perfect_uniform_repeats.input.R1_T4,
    output:
        lc_inter_dir / "repeats_imp_hp_L{merged_len}_B{base}.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        grep 'unit={wildcards.base}' {input} | \
        mergeBed -i stdin -d 1 | \
        awk '$3-$2>{wildcards.merged_len}' > \
        {output}
        """


rule merge_imperfect_uniform_repeats:
    input:
        beds=expand(
            rules.merge_perfect_uniform_repeats.output,
            allow_missing=True,
            base=["A", "C", "G", "T"],
        ),
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_SimpleRepeat_imperfecthomopolgt{merged_len}_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        multiIntersectBed -i {input.beds} | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


rule all_uniform_repeats:
    input:
        # Perfect (greater than X)
        **{
            f"perfect_{k}_gt{x - 1}": expand(
                rules.slop_uniform_repeats.output,
                unit_name=k,
                total_len=x - 1 + (offset := get_offset(k)),
            )
            for k, v in uniform_repeats.items()
            for x in v.total_lens[v.range_indices - 1 :]
        },
        # Perfect (between X and Y)
        **{
            f"perfect_{k}_{a}to{b}": expand(
                rules.slop_uniform_repeat_ranges.output,
                unit_name=k,
                total_lenA=a + (offset := get_offset(k)),
                total_lenB=b + offset,
            )
            for k, v in uniform_repeats.items()
            for a, b in zip(
                v.total_lens[0 : v.range_indices - 1],
                map(lambda x: x - 1, v.total_lens[1 : v.range_indices]),
            )
        },
        # Imperfect (greater than X)
        **{
            f"imperfect_gt{x}": expand(
                rules.merge_imperfect_uniform_repeats.output, merged_len=x
            )
            for x in [10, 20]
        },


################################################################################
## trf


use rule download_ref as download_trf with:
    output:
        resources_dir / "trf_simreps.txt.gz",
    params:
        url=config["simreps_url"],


rule filter_sort_trf:
    input:
        bed=rules.download_trf.output,
        genome=rules.get_genome.output,
    output:
        ref_dir / "trf.txt.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input.bed} | \
        cut -f 2,3,4 | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        bedtools sort -i stdin -g {input.genome} | \
        bgzip -c > {output}
        """


rule merge_trf:
    input:
        rules.filter_sort_trf.output,
    output:
        lc_inter_dir / "trf_simpleRepeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input} | \
        mergeBed -i stdin |  \
        bgzip -c > {output}
        """


################################################################################
## rmsk


use rule download_ref as download_rmsk with:
    output:
        resources_dir / "rmsk.txt.gz",
    params:
        url=config["rmsk_url"],


rule filter_sort_rmsk:
    input:
        bed=rules.download_rmsk.output,
        genome=rules.get_genome.output,
    output:
        ref_dir / "rmsk.txt.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input.bed} | \
        cut -f 6,7,8,12 | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        bedtools sort -i stdin -g {input.genome} | \
        bgzip -c > {output}
        """


rule merge_rmsk_class:
    input:
        rmsk=rules.filter_sort_rmsk.output,
        idx=rules.index_ref.output,
    output:
        lc_inter_dir / "rmsk_{rmsk_class}.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        zgrep {wildcards.rmsk_class} {input.rmsk} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


rule all_rmsk_classes:
    input:
        **{
            c: expand(
                rules.merge_rmsk_class.output,
                rmsk_class=c,
            )
            for c in ["Low_complexity", "Simple_repeat", "Satellite"]
        },


rule merge_satellites:
    input:
        bed=rules.all_rmsk_classes.input.Satellite,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_satellites_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        slopBed -i {input.bed} -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


# TODO the previous version of this had chrM in it? (seems like a mistake)
rule invert_satellites:
    input:
        bed=rules.merge_satellites.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_notinsatellites_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "complementBed -i {input.bed} -g {input.genome} | bgzip -c > {output}"


################################################################################
## Tandem Repeats


rule merge_all_uniform_repeats:
    input:
        rules.all_uniform_repeats.input.imperfect_gt10,
        perfect=rules.all_perfect_uniform_repeats.input.R1_T7,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    threads: 4
    shell:
        """
        slopBed -i {input.perfect} -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        multiIntersectBed -i stdin {input.imperfect} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


rule invert_all_uniform_repeats:
    input:
        bed=rules.merge_all_uniform_repeats.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "complementBed -i {input.bed} -g {input.genome} | bgzip -c > {output}"


rule merge_repeats:
    input:
        beds=rules.all_rmsk_classes.input
        + rules.merge_trf.output
        + [
            rules.all_perfect_uniform_repeats.input.R2_T10,
            rules.all_perfect_uniform_repeats.input.R3_T14,
            rules.all_perfect_uniform_repeats.input.R4_T19,
        ],
        genome=rules.get_genome.output,
    output:
        lc_inter_dir / "AllTandemRepeats_intermediate.bed",
    conda:
        envs_path("bedtools.yml")
    threads: 3
    shell:
        """
        multiIntersectBed -i {input.beds} | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        mergeBed -i stdin > {output}
        """


# NOTE these look funny because they have slop added
tr_bounds = {
    "lt51": {"lower": 0, "upper": 61},
    "51to200": {"lower": 60, "upper": 211},
    "201to10000": {"lower": 210, "upper": 10011},
    "gt10000": {"lower": 10010, "upper": 1e10},
    "gt100": {"lower": 110, "upper": 1e10},
}


rule filter_TRs:
    input:
        tr=rules.merge_repeats.output,
        hp=rules.merge_all_uniform_repeats.output,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeats_{tr_bound}bp_slop5.bed.gz",
    params:
        lower=lambda wildcards: tr_bounds[wildcards.tr_bound]["lower"],
        upper=lambda wildcards: tr_bounds[wildcards.tr_bound]["upper"],
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        awk '$3-$2>{params.lower} && $3-$2<{params.upper}' {input.tr} | \
        subtractBed -a stdin -b {input.hp} | \
        bgzip -c > {output}
        """


rule all_TRs:
    input:
        **{f"_{k}": expand(rules.filter_TRs.output, tr_bound=k) for k in tr_bounds},


rule merge_filtered_TRs:
    input:
        rules.all_TRs.input._lt51,
        rules.all_TRs.input._51to200,
        rules.all_TRs.input._201to10000,
        rules.all_TRs.input._gt10000,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        multiIntersectBed -i {input.beds} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


rule invert_TRs:
    input:
        beds=rules.merge_filtered_TRs.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_notinallTandemRepeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "complementBed -i {input.beds} -g {input.genome} | bgzip -c > {output}"


################################################################################
## Combine all the beds to make a Pink Floyd album cover


rule merge_HPs_and_TRs:
    input:
        rules.merge_filtered_TRs.input,
        rules.merge_all_uniform_repeats.output,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    threads: 2
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


# NOTE chrM in the original
rule invert_HPs_and_TRs:
    input:
        beds=rules.merge_HPs_and_TRs.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "complementBed -i {input.beds} -g {input.genome} | bgzip -c > {output}"


rule all_low_complexity:
    input:
        # Uniform repeats
        rules.all_uniform_repeats.input,
        # Satellites
        rules.merge_satellites.output,
        rules.invert_satellites.output,
        # Tandem Repeats
        rules.all_TRs.input,
        rules.merge_filtered_TRs.output,
        rules.invert_TRs.output,
        # "Everything" (in theory)
        rules.merge_HPs_and_TRs.output,
        rules.invert_HPs_and_TRs.output,
        rules.invert_all_uniform_repeats.output,
