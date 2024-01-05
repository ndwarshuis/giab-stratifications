from common.config import CoreLevel, strip_full_refkey

xy = config.to_bed_dirs(CoreLevel.XY)


use rule download_ref as download_genome_features_bed with:
    output:
        xy.src.data / "genome_features_{sex_chr}.bed.gz",
    log:
        xy.src.log / "genome_features_{sex_chr}.log",
    params:
        src=lambda w: (
            config.refsrckey_to_y_features_src
            if w.sex_chr == "Y"
            else config.refsrckey_to_x_features_src
        )(w.ref_src_key),
    localrule: True


rule write_PAR_final:
    input:
        bed=rules.write_PAR_intermediate.output,
        gapless=rules.get_gapless.output.parY,
        genome=rules.get_genome.output,
    output:
        xy.final("chr{sex_chr}_PAR"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -a {input.bed} -b {input.gapless} -sorted -g {input.genome} | \
        bgzip -c > {output}
        """


rule filter_XTR_features:
    input:
        bed=lambda w: expand_final_to_src(
            rules.download_genome_features_bed.output, w
        )[0],
        gapless=rules.get_gapless.output.parY,
        genome=rules.get_genome.output,
    output:
        xy.final("chr{sex_chr}_XTR"),
    conda:
        "../envs/bedtools.yml"
    params:
        level="XTR",
    script:
        "../scripts/python/bedtools/xy/filter_sort_features.py"


use rule filter_XTR_features as filter_ampliconic_features with:
    output:
        xy.final("chr{sex_chr}_ampliconic"),
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
        xy.final("chr{sex_chr}_nonPAR"),
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
        xy.final("AllAutosomes"),
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
def all_xy_features(ref_final_key, build_key):
    bd = config.to_build_data(strip_full_refkey(ref_final_key), build_key)
    targets = [
        (rules.filter_XTR_features.output[0], bd.want_xy_XTR),
        (rules.filter_ampliconic_features.output[0], bd.want_xy_ampliconic),
    ]
    return [
        t
        for p, test in targets
        if test
        for t in expand(
            p,
            allow_missing=True,
            sex_chr=bd.wanted_xy_chr_names,
            ref_final_key=ref_final_key,
            build_key=build_key,
        )
    ]


def all_xy_PAR(ref_final_key, build_key):
    bd = config.to_build_data(strip_full_refkey(ref_final_key), build_key)
    wanted_chrs = [c for (c, t) in [("X", bd.want_x_PAR), ("Y", bd.want_y_PAR)] if t]
    return expand(
        rules.invert_PAR.output + rules.write_PAR_final.output,
        allow_missing=True,
        sex_chr=wanted_chrs,
        ref_final_key=ref_final_key,
        build_key=build_key,
    )


def all_xy_sex(ref_final_key, build_key):
    return all_xy_PAR(ref_final_key, build_key) + all_xy_features(
        ref_final_key, build_key
    )
