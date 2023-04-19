from typing import Any
from os.path import dirname
from common.io import get_md5


def strat_files(path: str) -> list[str]:
    with open(path, "r") as f:
        return [s.strip() for s in f]


def main(smk: Any) -> None:
    out = str(smk.output)
    root = dirname(out)
    ss = strat_files(str(smk.input))
    with open(out, "w") as op:
        for s in ss:
            h = get_md5(s)
            p = s.replace(root, "")
            op.write(f"{h}  {p}\n")


main(snakemake)  # type: ignore
