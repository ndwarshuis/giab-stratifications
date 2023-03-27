uni_final_dir = final_dir / "Union"


rule intersect_segdup_and_map:
    input:
        rules.merge_superdups.output,
        rules.merge_nonunique.output,
    output:
        uni_final_dir / "GRCh38_alllowmapandsegdupregions.bed.gz",
    conda:
        envs_path("bedtools.yml")
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
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        uni_final_dir / "GRCh38_notinalllowmapandsegdupregions.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        complementBed -i {input.bed} -g {params.genome} |
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c \
        > {output}
        """


# TODO still need to add xtr + ampliconic (as applicable) to this
use rule intersect_segdup_and_map as intersect_alldifficult with:
    input:
        rules.intersect_segdup_and_map.output,
        rules.all_gc.input.wide,
        rules.merge_HPs_and_TRs.output,
    output:
        uni_final_dir / "GRCh38_alldifficultregions.bed.gz",


use rule invert_segdup_and_map as invert_alldifficult with:
    input:
        bed=rules.intersect_alldifficult.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        uni_final_dir / "GRCh38_notinalldifficultregions.bed.gz",
