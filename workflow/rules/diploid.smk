# TODO add rule to split ref if dip1


# htsbox doesn't work if the ref is zipped
rule unzip_ref:
    input:
        rules.filter_sort_ref.output,
    output:
        TODO,
    conda:
        "../envs/utils.yml"
    shell:
        """
        gunzip -c {input} > {output}
        """


rule cross_align_hap:
    input:
        this_hap=TODO,
        other_hap=TODO,
    output:
        TODO,
    # TODO I can cheat a bit here and cap the thread count to the number of
    # chromosomes desired (since mm2 runs with one thread/chromosome)
    threads: 8
    params:
        # hardcoded for now, we probably don't need to change this until we
        # move to another asm that's more heterozygous (ie African) if ever
        mm2_args=lambda _, threads: f"-z200000,10000,200 -xasm5 --cs -t{threads}",
    log:
        TODO,
    # this is what dipcall does to produce the pair file in one step
    shell:
        """
        minimap2 -a {params.mm2_args} {input.this_hap} {input.other_hap} 2> {log} | \
        bgzip -c > {output}
        """


rule filter_sort_cross_alignment:
    input:
        aux_bin=TODO,
        # TODO this needs to be unzipped
        hap=TODO,
        sam=rules.cross_align_hap.output,
    output:
        TODO,
    shell:
        """
        k8 {input.aux_bin} samflt {input.sam} | \
        samtools sort -m{resources.mem_mb}M --threads {threads} | \
        """


# TODO not sure if I also need to filter out alignments where one chr is
# clearly not mapped to the other, something like this: awk '{if ((gensub("_PATERNAL", "", 1, $1) == gensub("_MATERNAL", "", 1, $3)) || ("@" == substr($0,0,1))) { print $0 }}'


rule cross_alignment_to_hets:
    input:
        bam=rule.filter_sort_cross_alignment.output,
        hap=TODO,
    output:
        TODO,
    shell:
        """
        htsbox pileup -q5 -evcf {input.hap} {input.bam} | \
        grep -v '^#' | \
        awk 'OFS="\t" {print $1, $2-1, $2+length($4)-1}' | \
        bgzip -c > {output}
        """


rule zip_hets:
    input:
        rules.cross_alignment_to_bed.output,
    output:
        TODO,
    shell:
        """
        bgzip -c {input} > {output}
        """


rule merge_dip1_hets:
    # NOTE: the wildcard value of ref_final_key will not have a haplotype (as
    # per logic of dip1 references)
    input:
        lambda w: {
            h.name: expand(
                rules.cross_alignment_to_bed.output,
                ref_final_key=cfg.RefKeyFull(w.ref_final_key, h).name,
            )
            for h in cfg.Haplotype
        },
    output:
        TODO,
    shell:
        # TODO inputs are not bgzipped
        """
        cat {input.hap1} {input.hap2} | bgzip -c > {output}
        """


def dip1_or_dip2(dip2, dip1, wildcards):
    return dip2 if cfg.parse_full_refkey(w.ref_key_full).has_hap else dip1


rule filter_snp_hets:
    input:
        lambda w: dip1_or_dip2(rules.zip_hets.output, rules.merge_dip1_hets.output, w),
    output:
        TODO,
    shell:
        """
        gunzip -c xalign_same.bed.gz | \
        awk '{if($3-$2==1) {print $0}}' | \
        bgzip -c > {output}
        """


# TODO add rule to recombine ref if dip1
