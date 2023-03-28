import hashlib
from typing import Any
from os.path import dirname


def get_md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


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
