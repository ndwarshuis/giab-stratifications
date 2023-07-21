from common.config import CoreLevel

xy_dir = CoreLevel.XY
xy_src_dir = config.ref_src_dir / xy_dir.value
xy_inter_dir = config.intermediate_build_dir / xy_dir.value
xy_log_src_dir = config.log_src_dir / xy_dir.value


def xy_final_path(name):
    return config.build_strat_path(xy_dir, name)


use rule download_ref as download_genome_features_bed with:
    output:
        xy_src_dir / "genome_features_{sex_chr}.bed.gz",
    log:
        xy_log_src_dir / "genome_features_{sex_chr}.log",
    params:
        src=lambda w: (
            config.refkey_to_y_features_src
            if w.sex_chr == "Y"
            else config.refkey_to_x_features_src
        )(w.ref_key),
    localrule: True


rule write_PAR_final:
    input:
        bed=rules.write_PAR_intermediate.output,
        gapless=rules.get_gapless.output.parY,
        genome=rules.get_genome.output,
    output:
        xy_final_path("chr{sex_chr}_PAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -a {input.bed} -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_XTR_features:
    input:
        bed=rules.download_genome_features_bed.output[0],
        gapless=rules.get_gapless.output.parY,
        genome=rules.get_genome.output,
    output:
        xy_final_path("chr{sex_chr}_XTR"),
    conda:
        "../envs/bedtools.yml"
    params:
        level="XTR",
    script:
        "../scripts/python/bedtools/xy/filter_sort_features.py"


use rule filter_XTR_features as filter_ampliconic_features with:
    input:
        bed=rules.download_genome_features_bed.output[0],
        gapless=rules.get_gapless.output.parY,
        genome=rules.get_genome.output,
    output:
        xy_final_path("chr{sex_chr}_ampliconic"),
    params:
        level="Ampliconic",


# def par_input(wildcards):
#     test_fun = config.want_xy_y if wildcards.sex_chr == "Y" else config.want_xy_x
#     return (
#         rules.write_PAR_final.output
#         if test_fun(wildcards.ref_key, wildcards.build_key)
#         else rules.write_PAR_intermediate.output
#     )


rule invert_PAR:
    input:
        bed=rules.write_PAR_final.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.parY,
    output:
        xy_final_path("chr{sex_chr}_nonPAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        grep {wildcards.sex_chr} | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_autosomes:
    input:
        bed=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
        genome=rules.get_genome.output,
    output:
        xy_final_path("AllAutosomes"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk -v OFS='\t' {{'print $1,\"0\",$2'}} {input.bed} | \
        grep -v \"X\|Y\" | \
        intersectBed -a stdin -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


# helper functions for build targets
def all_xy_features(ref_key, build_key):
    targets = [
        (rules.filter_XTR_features.output[0], config.want_xy_XTR),
        (rules.filter_ampliconic_features.output[0], config.want_xy_ampliconic),
    ]
    cns = config.wanted_xy_chr_names(ref_key, build_key)
    return [
        t
        for p, f in targets
        if f(ref_key)
        for t in expand(
            p,
            allow_missing=True,
            sex_chr=cns,
            ref_key=ref_key,
            build_key=build_key,
        )
    ]


def all_xy_PAR(ref_key, build_key):
    wanted_chrs = [
        c
        for (c, fun) in [("X", config.want_x_PAR), ("Y", config.want_y_PAR)]
        if fun(ref_key, build_key)
    ]
    return expand(
        rules.invert_PAR.output + rules.write_PAR_final.output,
        allow_missing=True,
        sex_chr=wanted_chrs,
        ref_key=ref_key,
        build_key=build_key,
    )


def all_xy_sex(ref_key, build_key):
    return all_xy_PAR(ref_key, build_key) + all_xy_features(ref_key, build_key)
