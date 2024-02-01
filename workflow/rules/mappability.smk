from os.path import splitext, basename
from pathlib import Path
from common.config import CoreLevel, strip_full_refkey

# TODO this entire thing needs to be split apart when run with dip1 references

mlty = config.to_bed_dirs(CoreLevel.MAPPABILITY)

################################################################################
# download a bunch of stuff to run GEM
#
# NOTE: this is in bioconda, but the bioconda version does not have gem-2-wig
# for some reason


gemlib_bin = Path("GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin")

gem_wc_constraints = {
    "l": "\d+",
    "m": "\d+",
    "e": "\d+",
}


rule download_gem:
    output:
        config.tools_src_dir / "gemlib.tbz2",
    params:
        url=config.tools.gemlib,
    conda:
        "../envs/utils.yml"
    localrule: True
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


rule filter_mappability_ref:
    input:
        fa=lambda w: expand_final_to_src(rules.download_ref.output[0], w)[0],
        idx=rules.index_ref.output[0],
    output:
        mlty.inter.postsort.data / "ref.fa",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/mappability/filter_ref.py"


rule gem_index:
    input:
        fa=rules.filter_mappability_ref.output[0],
        bin=rules.unpack_gem.output.indexer,
    output:
        mlty.inter.postsort.data / "index.gem",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    threads: lambda w: config.thread_per_chromosome(w.ref_final_key, w.build_key, 4) * 1.5
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.gemIndex
        ),
    log:
        mlty.inter.postsort.log / "index.log",
    benchmark:
        mlty.inter.postsort.bench / "index.txt"
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
        idx=rules.gem_index.output[0],
        bin=rules.unpack_gem.output.mappability,
    output:
        mlty.inter.postsort.data / "unique_l{l}_m{m}_e{e}.mappability",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    threads: lambda w: int(config.thread_per_chromosome(w.ref_final_key, w.build_key, 4) * 1.5)
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.gemMappability
        ),
    log:
        mlty.inter.postsort.data / "mappability_l{l}_m{m}_e{e}.log",
    benchmark:
        mlty.inter.postsort.bench / "mappability_l{l}_m{m}_e{e}.txt"
    wildcard_constraints:
        **gem_wc_constraints,
    shell:
        """
        {input.bin} \
        -m {wildcards.m} \
        -e {wildcards.e} \
        -l {wildcards.l} \
        -T {threads} \
        -I {input.idx} \
        -o {params.base} > {log} 2>&1
        """


rule gem_to_wig:
    input:
        idx=rules.gem_index.output,
        map=rules.gem_mappability.output[0],
        bin=rules.unpack_gem.output.gem2wig,
    output:
        mlty.inter.postsort.data / "unique_l{l}_m{m}_e{e}.wig",
    params:
        base=lambda wildcards, output: splitext(output[0])[0],
    log:
        mlty.inter.postsort.log / "gem2wig_l{l}_m{m}_e{e}.log",
    benchmark:
        mlty.inter.postsort.bench / "gem2wig_l{l}_m{m}_e{e}.txt"
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.gemToWig
        ),
    wildcard_constraints:
        **gem_wc_constraints,
    shell:
        """
        {input.bin} \
        -I {input.idx} \
        -i {input.map} \
        -o {params.base} 2> {log}
        """


# NOTE: -d option to not sort and save some speed/memory (we shall sort later,
# because this command does not do it the way I want)
rule wig_to_bed:
    input:
        rules.gem_to_wig.output,
    output:
        mlty.inter.postsort.data / "unique_l{l}_m{m}_e{e}.bed.gz",
    conda:
        "../envs/map.yml"
    wildcard_constraints:
        **gem_wc_constraints,
    shell:
        """
        sed 's/ AC//' {input} | \
        wig2bed -d | \
        awk '$5>0.9' | \
        cut -f1-3 | \
        gzip -c > \
        {output}
        """


################################################################################
# create stratifications


def nonunique_inputs(wildcards):
    rk = strip_full_refkey(wildcards.ref_final_key)
    bk = wildcards.build_key
    l, m, e = config.to_build_data(rk, bk).mappability_params
    return expand(rules.wig_to_bed.output, zip, allow_missing=True, l=l, m=m, e=e)


checkpoint merge_nonunique:
    input:
        bed=nonunique_inputs,
        gapless=rules.get_gapless.output.auto,
        genome=rules.get_genome.output,
    output:
        mlty.inter.postsort.data / "nonunique_output.json",
    params:
        path_pattern=lambda w: expand(
            mlty.final("{{}}"),
            ref_final_key=w.ref_final_key,
            build_key=w.build_key,
        )[0],
    resources:
        mem_mb=lambda w: config.buildkey_to_malloc(
            w.ref_final_key, w.build_key, lambda m: m.mergeNonunique
        ),
    benchmark:
        mlty.inter.postsort.bench / "merge_nonunique.txt"
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/python/bedtools/mappability/merge_nonunique.py"


def nonunique_inputs(ref_final_key, build_key):
    c = checkpoints.merge_nonunique.get(
        ref_final_key=ref_final_key, build_key=build_key
    )
    with c.output[0].open() as f:
        return json.load(f)


rule invert_merged_nonunique:
    input:
        lambda w: nonunique_inputs(w.ref_final_key, w.build_key)["all_lowmap"],
    output:
        mlty.final("notinlowmappabilityall"),
    conda:
        "../envs/bedtools.yml"
    params:
        genome=rules.get_genome.output,
        gapless=rules.get_gapless.output.auto,
    shell:
        """
        complementBed -i {input} -g {params.genome} | \
        intersectBed -a stdin -b {params.gapless} -sorted -g {params.genome} | \
        bgzip -c > {output}
        """


def nonunique_inputs_flat(ref_final_key, build_key):
    res = nonunique_inputs(ref_final_key, build_key)
    return [res["all_lowmap"], *res["single_lowmap"]]


def mappabilty_inputs(ref_final_key, build_key):
    return nonunique_inputs_flat(ref_final_key, build_key) + expand(
        rules.invert_merged_nonunique.output,
        ref_final_key=ref_final_key,
        build_key=build_key,
    )
