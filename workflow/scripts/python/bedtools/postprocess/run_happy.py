import gzip
import os
import subprocess as sp
from typing import Any
from pathlib import Path


def gzip_empty(path: str) -> bool:
    with gzip.open(path, "rb") as f:
        return len(f.read(1)) == 0


def main(smk: Any) -> None:
    # ASSUME if filtered benchmark bed is empty then we have nothing to
    # benchmark because the desired chromosomes in this build are not included
    # in the benchmark
    with open(smk.log[0], "w") as f:
        if gzip_empty(smk.input["bench_bed"][0]):
            Path(smk.output[0]).touch()
            f.write("no overlapping chromosomes, writing empty dummy file\n")
        else:
            cmd = [
                "hap.py",
                *["--engine", "vcfeval"],
                "--verbose",
                *["--threads", str(smk.threads)],
                *["--stratification", smk.input["strats"][0]],
                *["-f", smk.input["bench_bed"][0]],
                *["-o", smk.params["prefix"]],
                smk.input["bench_vcf"][0],
                smk.input["query_vcf"][0],
            ]
            res = sp.run(
                cmd,
                stderr=sp.STDOUT,
                stdout=f,
                env={"HGREF": smk.input["ref"][0], **os.environ},
            )
            if res.returncode != 0:
                exit(1)


main(snakemake)  # type: ignore
