rule download_ref:
    output:
        resources_dir / "ref.fna.gz",
    params:
        url=config["ref_url"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule download_ref as download_rmsk with:
    output:
        resources_dir / "rmsk.txt.gz",
    params:
        url=config["rmsk_url"],


use rule download_ref as download_simreps with:
    output:
        resources_dir / "simreps.txt.gz",
    params:
        url=config["simreps_url"],


rule unzip_ref:
    input:
        rules.download_ref.output,
    output:
        resources_dir / "ref.fna",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        grep -Ev '^chr[0-9XYM]_|^chr[0-9][0-9]_|^chrUn_|^chrEBV' \
        > {output}
        """


rule index_ref:
    input:
        rules.download_ref.output,
    output:
        resources_dir / "ref.fna.fai",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        samtools faidx - -o - | \
        grep -Ev '^chr[0-9XYM]_|^chr[0-9][0-9]_|^chrUn_|^chrEBV' > \
        {output}
        """


use rule download_ref as download_repseq with:
    output:
        resources_dir / "tools" / "repseq.tar.gz",
    params:
        url=config["tools"]["repseq"],


rule unpack_repseq:
    input:
        rules.download_repseq.output,
    output:
        directory(results_dir / "build" / "repseq"),
    shell:
        """
        mkdir {output} && \
        tar xzf {input} --directory {output} --strip-components=1
        """


# TODO add logging
rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        results_dir / "bin" / "repseq",
    conda:
        envs_path("build.yml")
    shell:
        "make -C {input} && mv {input}/repseq {output}"


hp_grid = {
    1: [4, 6, 12, 20],
    2: [10, 50, 200],
    3: [14, 50, 200],
    4: [19, 50, 200],
}


rule find_repeats:
    input:
        ref=rules.unzip_ref.output,
        bin=rules.build_repseq.output,
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{len}.bed",
    log:
        log_dir / "low_complexity" / "repeats_R{rep}_L{len}.log",
    shell:
        "{input.bin} {rep} {len} {input.ref} > {output} 2> {log}"


rule subtract_repeats:
    input:
        a=results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}.bed",
        b=results_dir / "low_complexity" / "repeats_R{rep}_L{lenB}.bed",
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}to{lenB}.bed",
    shell:
        "subtractBed -a {input.a} -b {input.b} > {output}"


rule add_repeat_slop:
    input:
        bed=rules.subtract_repeats.output,
        # genome=TODO
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}to{lenB}.bed",


# TODO get multi-inter of this
rule get_imperfect_homopolymers:
    input:
        results_dir / "low_complexity" / "repeats_R1_L4.bed",
    output:
        results_dir / "low_complexity" / "repeats_imp_hp_{len}_{base}.bed",
    shell:
        """
        grep 'unit=C' {input} | \
        mergeBed -i stdin -d 1 | \
        awk '$3-$2>10' > \
        {output}
        """
