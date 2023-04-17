rule download_repseq:
    output:
        config.tools_src_dir / "repseq.tar.gz",
    params:
        url=config.tools.repseq,
    conda:
        "../envs/utils.yml"
    localrule: True
    shell:
        "curl -sS -L -o {output} {params.url}"


rule unpack_repseq:
    input:
        rules.download_repseq.output,
    output:
        directory(config.tools_make_dir / "repseq"),
    shell:
        """
        mkdir {output} && \
        tar xzf {input} --directory {output} --strip-components=1
        """


rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        config.tools_bin_dir / "repseq",
    conda:
        "../envs/build.yml"
    log:
        config.log_root_dir / "tools" / "repseq_build.log",
    shell:
        "make -C {input} 2>&1 > {log} && mv {input}/repseq {output}"
