from common.config import is_bgzip
from typing import Any
from os import scandir


def main(smk: Any) -> None:
    nonbgzip = [
        p
        for p in scandir(smk.input[0])
        if p.is_file() and p.path.endswith(".gz") and not is_bgzip(p.path)
    ]
    assert len(nonbgzip) == 0, f"non bgzip paths: {nonbgzip}"


main(snakemake)  # type: ignore
