from os.path import splitext, basename
from pathlib import Path

map_inter_dir = config.intermediate_build_dir / "mappability"
map_log_dir = config.log_build_dir / "mappability"


def map_final_path(name):
    return config.build_strat_path("mappability", name)


################################################################################
# download a bunch of stuff to run GEM
#
# NOTE: this is in bioconda, but the bioconda version does not have gem-2-wig
# for some reason


gemlib_bin = Path("GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin")


rule download_gem:
    output:
        config.tools_src_dir / "gemlib.tbz2",
    params:
        url=config.tools.gemlib,
    conda:
        "../envs/utils.yml"
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unpack_gem:
    input:
        rules.download_gem.output,
    output:
        # called by other binaries
        config.tools_bin_dir / "gem-indexer_fasta2meta+cont",
        config.tools_bin_dir / "gem-indexer_bwt-dna",
        config.tools_bin_dir / "gem-indexer_generate",
        # the things I actually need
        indexer=config.tools_bin_dir / "gem-indexer",
        mappability=config.tools_bin_dir / "gem-mappability",
        gem2wig=config.tools_bin_dir / "gem-2-wig",
    params:
        bins=lambda wildcards, output: " ".join(
            str(gemlib_bin / basename(o)) for o in output
        ),
    shell:
        """
        mkdir -p {config.tools_bin_dir} && \
        tar xjf {input} \
        --directory {config.tools_bin_dir} \
        --strip-components=2 \
        {params.bins}
        """


################################################################################
# index/align


# this seems to be the only place that requires the fa to be unzipped
rule unzip_ref:
    input:
        rules.filter_sort_ref.output,
    output:
        map_inter_dir / "ref.fa",
    shell:
        "gunzip -c {input} > {output}"


rule gem_index:
    input:
        fa=rules.unzip_ref.output,
        bin=rules.unpack_gem.output.indexer,
    output:
        map_inter_dir / "index.gem",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    threads: 8
    log:
        map_log_dir / "index.log",
    shell:
        """
        PATH={config.tools_bin_dir}:$PATH
        {input.bin} \
        --complement emulate \
        -T {threads} \
        -i {input.fa} \
        -o {params.base} > {log} 2>&1
        """


rule gem_mappability:
    input:
        fa=rules.gem_index.output,
        bin=rules.unpack_gem.output.mappability,
    output:
        map_inter_dir / "GRCh38_unique_l{l}_m{m}_e{e}.mappability",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    threads: 8
    log:
        map_log_dir / "mappability_l{l}_m{m}_e{e}.log",
    shell:
        """
        {input.bin} \
        -m {wildcards.m} \
        -e {wildcards.e} \
        -l {wildcards.l} \
        -T {threads} \
        -I {input.fa} \
        -o {params.base} > {log} 2>&1
        """


rule gem_to_wig:
    input:
        idx=rules.gem_index.output,
        map=rules.gem_mappability.output,
        bin=rules.unpack_gem.output.gem2wig,
    output:
        map_inter_dir / "GRCh38_unique_l{l}_m{m}_e{e}.wig",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    log:
        map_log_dir / "gem2wig_l{l}_m{m}_e{e}.log",
    shell:
        """
        {input.bin} \
        -I {input.idx} \
        -i {input.map} \
        -o {params.base} 2> {log}
        """


rule wig_to_bed:
    input:
        rules.gem_to_wig.output,
    output:
        map_inter_dir / "GRCh38_unique_l{l}_m{m}_e{e}.bed.gz",
    resources:
        mem_mb=8000,
    conda:
        "../envs/map.yml"
    shell:
        """
        sed 's/ AC//' {input} | \
        wig2bed -m {resources.mem_mb}M | \
        awk '$5>0.9' | \
        gzip -c > \
        {output}
        """


################################################################################
# create stratifications


rule get_nonunique:
    input:
        rules.wig_to_bed.output,
    output:
        map_final_path("nonunique_l{l}_m{m}_e{e}"),
    conda:
        "../envs/bedtools.yml"
    params:
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    shell:
        """
        complementBed -i {input} -g {params.genome} | \
        mergeBed -d 100 -i stdin | \
        intersectBed -a stdin -b {params.gapless} -sorted | \
        bgzip -c > \
        {output}
        """


# same as above but don't put the file in the final directory
use rule get_nonunique as get_nonunique_single with:
    output:
        map_inter_dir / "nonunique_l{l}_m{m}_e{e}",


# This rule is complicated because we need to deal with the case of single vs
# multiple parameter specifications. If there is only one, we only need the
# "lowmappabilityall" all bed file as it includes all the data we want. If
# multiple, we need multiple bed files for each individual spec (eg
# "nonunique_lx_my_ez") plus the "lowmappabilityall". This hacky script handles
# this logic.
def nonunique_inputs(wildcards):
    rk = wildcards.ref_key
    bk = wildcards.build_key
    l, m, e = config.buildkey_to_mappability(rk, bk)
    n = len(l)
    if n == 0:
        assert False, "this should not happen"
    elif n == 1:
        path = rules.get_nonunique_single.output
    else:
        path = rules.get_nonunique.output
    return expand(path, zip, allow_missing=True, l=l, m=m, e=e)


rule merge_nonunique:
    input:
        bed=nonunique_inputs,
        gapless=rules.get_gapless.output.auto,
    output:
        map_final_path("lowmappabilityall"),
    conda:
        "../envs/bedtools.yml"
    params:
        n=lambda _, input: len(input.bed),
    shell:
        """
        n={params.n}
        if [ $n == 1 ]; then 
            cp {input.bed} {output}
        else
            multiIntersectBed -i {input.bed} | \
            mergeBed -d 100 -i stdin | \
            intersectBed -a stdin -b {input.gapless} -sorted | \
            bgzip -c > \
            {output}
        fi
        """


rule invert_merged_nonunique:
    input:
        rules.merge_nonunique.output,
    output:
        map_final_path("notinlowmappabilityall"),
    conda:
        "../envs/bedtools.yml"
    params:
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    shell:
        """
        complementBed -i {input} -g {params.genome} | \
        intersectBed -a stdin -b {params.gapless} -sorted | \
        bgzip -c > {output}
        """


rule all_map:
    input:
        rules.invert_merged_nonunique.output,
