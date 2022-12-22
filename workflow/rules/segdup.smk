segdup_res_dir = resources_dir / "SegmentalDuplications"
segdup_final_dir = final_dir / "SegmentalDuplications"


rule download_self_chain:
    output:
        segdup_res_dir / "selfChain.txt.gz",
    params:
        url=config["segdups"]["sefl_chain"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule download_self_chain_link:
    output:
        segdup_res_dir / "selfChain_link.txt.gz",
    params:
        url=config["segdups"]["sefl_chain"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


rule filter_trivial:
    input:
        rules.download_self_chain.output,
    output:
        "GRCh38_chainSelf_notrivial.txt.gz",
    conda:
        envs_path("utils.yml")
    shell:
        "gzcat {input} | awk '$3!=$7 || $5!=$10' | bgzip > {output}"


# TODO what is this "datamash" thing?

# gzcat GRCh38_chainSelfLink.txt.gz | \
# awk '{ print 3 "\t" 5 "\t" 4-$3; }' | \
# sort -k1,1 -k2,2g > \
# ./script_intermediates/38_chainSelfLink_intermediate.bed

# gzcat GRCh38_chainSelfLink.txt.gz | \
# awk '{ print 3 "\t" $4; }' | \
# bedtools merge -d 100 -i stdin | \
# sort -k1,1 -k2,2g | \
# bedtools intersect -sorted -wa -wb -a stdin -b ./script_intermediates/38_chainSelfLink_intermediate.bed | \
# sort -k1,1 -k2,2 -k3,3 | \
# datamash -g1,2,3 min 7 max 8  >\
# ./script_intermediates/38_datamash_intermediate.bed

# gzcat GRCh38_chainSelf_notrivial.txt.gz | \
# awk '{ print 3 "\t" 8 "\t" 13; }' | \
# sort -k1,1 >
# ./script_intermediates/38_chainSelf_intermediate.bed

# join -t $'\t' ./script_intermediates/38_datamash_intermediate.bed ./script_intermediates/38_chainSelf_intermediate.bed > \
# ./script_intermediates/38_joined_intermediate.bed

# awk '{ qs=(4:5; qe=(5:4; printf "%s\t%d\t%d\t%s:%'\''d-%'\''d\t%s\t%d\n", 2, 7, qs+1, qe, 10); }' ./script_intermediates/38_joined_intermediate.bed | \
# sort -k1,1 -k2,2g | \
# bgzip -c > \
# /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_beds/selfchain/new_selfchains/annotation/GRCh38_Aaron_code_chainSelf.bed.gz
# tabix -f ./annotation/GRCh38_Aaron_code_chainSelf.bed.gz


# NOTE new in v3.1 (subtract the PAR region)
rule self_chain_thing:
    input:
        chain="../v3.0-carry-over-stratifications/GRCh38/SegmentalDuplications/GRCh38_chainSelf.bed.gz",
        par=rules.download_x_PAR.output,
        genome=rules.get_genome.output,
    output:
        "SegmentalDuplications/GRCh38_chainSelf.bed.gz",
    shell:
        """
        gzcat {input} | \
        grep -v '#' | \
        subtractBed -a stdin -b {input.par} | \
        sortBed -faidx {input.genome} -i stdin | \
        mergeBed -i stdin -d 100 | \
        bgzip -c > \
        {output}
        """


# TODO do I actually need to sort?
rule notin_self_chain:
    input:
        bed=rules.self_chain_thing.output,
        genome=rules.get_genome.output,
    output:
        "SegmentalDuplications/GRCh38_notinchainSelf.bed.gz",
    shell:
        """
        subtractBed -a {input.genome} -b {input.bed} | \
        sortBed -faidx ref-files/GRCh38.fa.fai -i stdin | \
        mergeBed -i stdin -d 100 | \
        bgzip -c > \
        {output}
        """


# rinse/repeat the above two steps for gt10kb self chains
# TODO add "all rule"
