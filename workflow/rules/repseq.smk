rule download_repseq:
    output:
        resources_dir / "tools" / "repseq.tar.gz",
    params:
        url=config["tools"]["repseq"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


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
