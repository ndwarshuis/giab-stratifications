rule download_repseq:
    output:
        config.tools_src_dir / "repseq.tar.gz",
    params:
        url=config.tools.repseq,
    conda:
        config.env_path("utils")
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


# TODO add logging
rule build_repseq:
    input:
        rules.unpack_repseq.output,
    output:
        config.tools_bin_dir / "repseq",
    conda:
        config.env_path("build")
    shell:
        "make -C {input} && mv {input}/repseq {output}"
