from common.config import (
    CoreLevel,
    parse_full_refkey_class,
    strip_full_refkey,
    flip_full_refkey,
)

dip = config.to_bed_dirs(CoreLevel.DIPLOID)

# TODO don't hardcode minimap2 params (which might be changed if we move to a
# different asm)


def minimap_inputs(wildcards):
    rk = wildcards["ref_final_key"]
    other_rk = flip_full_refkey(rk)
    fa, idx = config.dip1_either(
        (rules.split_ref.output, rules.index_split_ref.output),
        (rules.unzip_ref.output, rules.index_unzipped_ref.output),
        rk,
    )

    def expand_rk(path, rk):
        return expand(path, allow_missing=True, ref_final_key=rk)

    return {
        "this_hap": expand_rk(fa, rk),
        "_this_idx": expand_rk(idx, rk),
        "other_hap": expand_rk(fa, other_rk),
        "_other_idx": expand_rk(idx, other_rk),
    }


# Dipcall normally outputs a bed file that roughly corresponds to "regions with
# low divergence". Since we are interested in large structural variation in
# addition to small variants, obtain this bed file using the same process and
# the flip it.
rule cross_align_large:
    input:
        unpack(minimap_inputs),
    output:
        dip.inter.postsort.data / "large_cross_align.paf.gz",
    # TODO I can cheat a bit here and cap the thread count to the number of
    # chromosomes desired (since mm2 runs with one thread/chromosome)
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 8)
    log:
        dip.inter.postsort.log / "cross_align_large.log",
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        minimap2 -c --paf-no-hit -t{threads} --cs -z200000,10000,200 -xasm5 \
          {input.this_hap} \
          {input.other_hap} \
          2> {log} | \
        bgzip -c > {output}
        """


rule large_cross_alignment_to_bed:
    input:
        paf=rules.cross_align_large.output,
        paftools_bin=rules.download_paftools.output,
        genome=lambda w: config.dip1_either(
            rules.get_split_genome.output,
            rules.get_genome.output,
            w["ref_final_key"],
        ),
    output:
        dip.inter.postsort.data / "large_cross_align.bed.gz",
    log:
        dip.inter.postsort.log / "large_cross_alignment_to_bed.log",
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        gunzip -c {input.paf} | \
        sort -k6,6 -k8,8n | \
        k8 {input.paftools_bin} call - 2> {log} | \
        grep ^R | \
        cut -f2- | \
        bedtools complement -i stdin -g {input.genome} | \
        gzip -c > {output}
        """


rule cross_align_small:
    input:
        unpack(minimap_inputs),
    output:
        dip.inter.postsort.data / "small_cross_align.sam.gz",
    # TODO I can cheat a bit here and cap the thread count to the number of
    # chromosomes desired (since mm2 runs with one thread/chromosome)
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 8)
    log:
        dip.inter.postsort.log / "cross_align_small.log",
    # this is what dipcall does to produce the pair file in one step
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        minimap2 -a -t{threads} --cs -z200000,10000,200 -xasm5 \
          {input.this_hap} \
          {input.other_hap} \
          2> {log} | \
        bgzip -c > {output}
        """


rule filter_sort_small_cross_alignment:
    input:
        aux_bin=rules.download_dipcall_aux.output,
        sam=rules.cross_align_small.output,
    output:
        dip.inter.postsort.data / "sorted_small_cross_alignments.bam",
    conda:
        "../envs/quasi-dipcall.yml"
    resources:
        # TODO don't hardcode
        mem_mb=1000,
    shell:
        """
        k8 {input.aux_bin} samflt {input.sam} | \
        samtools sort -m{resources.mem_mb}M --threads {threads} \
        > {output}
        """


# TODO not sure if I also need to filter out alignments where one chr is
# clearly not mapped to the other, something like this: awk '{if ((gensub("_PATERNAL", "", 1, $1) == gensub("_MATERNAL", "", 1, $3)) || ("@" == substr($0,0,1))) { print $0 }}'


rule small_cross_alignment_to_bed:
    input:
        bam=rules.filter_sort_small_cross_alignment.output,
        hap=lambda w: config.dip1_either(
            rules.split_ref.output,
            rules.unzip_ref.output,
            w["ref_final_key"],
        ),
    output:
        dip.inter.postsort.data / "small_cross_align.bed.gz",
    conda:
        "../envs/quasi-dipcall.yml"
    shell:
        """
        htsbox pileup -q5 -evcf {input.hap} {input.bam} | \
        grep -v '^#' | \
        awk 'OFS="\t" {{print $1, $2-1, $2+length($4)-1}}' | \
        gzip -c > {output}
        """


rule merge_large_and_small_hets:
    input:
        small=rules.small_cross_alignment_to_bed.output,
        large=rules.large_cross_alignment_to_bed.output,
    output:
        dip.inter.postsort.data / "het_regions.bed.gz",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        multiIntersectBed -i {input.small} {input.large} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


rule combine_dip1_hets:
    # NOTE: the wildcard value of ref_final_key will not have a haplotype (as
    # per logic of dip1 references)
    input:
        lambda w: {
            h.name: expand(
                rules.merge_large_and_small.output,
                ref_final_key=cfg.RefKeyFull(w.ref_final_key, h).name,
            )
            for h in cfg.Haplotype
        },
    output:
        dip.inter.postsort.data / "combined_het_regions.bed.gz",
    shell:
        """
        cat {input.hap1} {input.hap2} > {output}
        """


def het_region_inputs(wildcards):
    return (
        rules.merge_large_and_small_hets.output
        if parse_full_refkey_class(wildcards.ref_final_key).has_hap
        else rules.combine_dip1_hets.output
    )


rule filter_snp_hets:
    input:
        het_region_inputs,
    output:
        dip.inter.postsort.data / "snp_het_regions.bed.gz",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        gunzip -c {input} | \
        awk '{{if($3-$2==1) {{print $0}}}}' | \
        bgzip -c > {output}
        """


rule merge_het_regions:
    input:
        het_region_inputs,
    output:
        dip.final("het_regions_{merge_len}k"),
    conda:
        "../envs/bedtools.yml"
    params:
        gapless=rules.get_gapless.output.auto,
    wildcard_constraints:
        merge_len=f"\d+",
    shell:
        """
        mergeBed -i {input} -d $(({wildcards.merge_len}*1000)) | \
        intersectBed -a stdin -b {params.gapless} -sorted | \
        bgzip -c > {output}
        """


rule invert_het_regions:
    input:
        rules.merge_het_regions.output,
    params:
        gapless=rules.get_gapless.output.auto,
        genome=rules.get_genome.output,
    conda:
        "../envs/bedtools.yml"
    output:
        dip.final("hom_regions_{merge_len}k"),
    wildcard_constraints:
        merge_len=f"\d+",
    shell:
        """
        complementBed -i {input} -g {params.genome} | \
        intersectBed -a stdin -b {params.gapless} -sorted | \
        bgzip -c > {output}
        """


use rule merge_het_regions as merge_het_snp_regions with:
    input:
        rules.filter_snp_hets.output,
    output:
        dip.final("het_snp_regions_{merge_len}k"),
    wildcard_constraints:
        merge_len=f"\d+",


use rule invert_het_regions as invert_het_snp_regions with:
    input:
        rules.merge_het_snp_regions.output,
    output:
        dip.final("hom_snp_regions_{merge_len}k"),
    wildcard_constraints:
        merge_len=f"\d+",


def het_hom_inputs(ref_final_key, build_key):
    bd = config.to_build_data(strip_full_refkey(ref_final_key), build_key)
    return expand(
        rules.invert_het_regions.output
        + rules.invert_het_snp_regions.output
        + rules.merge_het_regions.output
        + rules.merge_het_snp_regions.output,
        merge_len=bd.build.include.hets,
        ref_final_key=ref_final_key,
        build_key=build_key,
    )
