from common.config import CoreLevel


def uni_final_path(name):
    return config.build_strat_path(CoreLevel.UNION, name)


rule intersect_segdup_and_map:
    input:
        rules.merge_superdups.output,
        lambda w: nonunique_inputs(w.ref_key, w.build_key)["all_lowmap"],
    output:
        uni_final_path("alllowmapandsegdupregions"),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        bgzip -c \
        > {output}
        """


rule invert_segdup_and_map:
    input:
        bed=rules.intersect_segdup_and_map.output,
    output:
        uni_final_path("notinalllowmapandsegdupregions"),
    conda:
        "../envs/bedtools.yml"
    params:
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    shell:
        """
        complementBed -i {input.bed} -g {params.genome} |
        intersectBed -a stdin -b {params.gapless} -sorted -g {params.genome} | \
        bgzip -c \
        > {output}
        """


use rule intersect_segdup_and_map as intersect_alldifficult with:
    input:
        rules.intersect_segdup_and_map.output,
        rules.merge_HPs_and_TRs.output,
        lambda w: all_xy_features(w.ref_key, w.build_key),
        lambda w: gc_inputs(w.ref_key, w.build_key)["widest_extreme"],
    output:
        uni_final_path("alldifficultregions"),


use rule invert_segdup_and_map as invert_alldifficult with:
    input:
        bed=rules.intersect_alldifficult.output,
    output:
        uni_final_path("notinalldifficultregions"),


rule all_segdup_and_map:
    input:
        rules.intersect_segdup_and_map.output,
        rules.invert_segdup_and_map.output,
    localrule: True


rule all_alldifficult:
    input:
        rules.intersect_alldifficult.output,
        rules.invert_alldifficult.output,
    localrule: True
