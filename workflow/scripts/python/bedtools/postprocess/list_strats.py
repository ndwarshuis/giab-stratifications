from typing import Any
from os.path import dirname
from os import scandir
from more_itertools import unique_everseen


def common_dirs(files: list[str]) -> list[str]:
    return list(unique_everseen([dirname(str(f)) for f in files]))


def main(smk: Any) -> None:
    strats = [
        p.path
        for i in common_dirs(smk.input.check)
        for p in scandir(i)
        if p.path.endswith(".bed.gz")
    ]

    with open(smk.output[0], "w") as f:
        for s in sorted(strats):
            f.write(s + "\n")


main(snakemake)  # type: ignore
