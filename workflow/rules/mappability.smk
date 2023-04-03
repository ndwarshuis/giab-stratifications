from os.path import splitext, basename
from pathlib import Path

map_inter_dir = config.intermediate_build_dir / "mappability"
map_log_dir = config.log_build_dir / "mappability"


def map_final_path(name):
    return config.build_strat_path("mappability", name)


rule download_gem:
    output:
        config.tools_src_dir / "gemlib.tbz2",
    params:
        url=config.tools.gemlib,
    conda:
        config.env_path("utils")
    shell:
        "curl -sS -L -o {output} {params.url}"


gemlib_bin = Path("GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin")


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


rule gem_index:
    input:
        fa=rules.filter_sort_ref.output,
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
        envs_path("map.yml")
    shell:
        """
        sed 's/ AC//' {input} | \
        wig2bed -m {resources.mem_mb}M | \
        awk '$5>0.9' | \
        gzip -c > \
        {output}
        """


rule invert_unique:
    input:
        bed=rules.wig_to_bed.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        map_final_path("nonunique_l{l}_m{m}_e{e}"),
    conda:
        config.env_path("bedtools")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        mergeBed -d 100 -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > \
        {output}
        """


rule all_nonunique:
    input:
        expand(
            rules.invert_unique.output,
            zip,
            allow_missing=True,
            l=[100, 250],
            m=[2, 0],
            e=[1, 0],
        ),


rule merge_nonunique:
    input:
        bed=rules.all_nonunique.input,
        gapless=rules.get_gapless.output.auto,
    output:
        map_final_path("lowmappabilityall"),
    conda:
        config.env_path("bedtools")
    shell:
        """
        multiIntersectBed -i {input.bed} | \
        mergeBed -d 100 -i stdin | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > \
        {output}
        """


rule invert_merged_nonunique:
    input:
        bed=rules.merge_nonunique.output,
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    output:
        map_final_path("notinlowmappabilityall"),
    conda:
        config.env_path("bedtools")
    shell:
        """
        complementBed -i {input.bed} -g {input.genome} | \
        intersectBed -a stdin -b {input.gapless} -sorted | \
        bgzip -c > {output}
        """


rule all_map:
    input:
        rules.all_nonunique.input,
        rules.merge_nonunique.output,
        rules.invert_merged_nonunique.output,
