from pathlib import Path
from typing import Any
import re


def main(smk: Any) -> None:
    with open(smk.input[0], "r") as i, open(smk.output[0], "w") as o:
        for f in i:
            fp = Path(f.strip())
            level = re.sub("^[^_]+_", "", fp.name.replace(".bed.gz", ""))
            o.write(f"{level}\t{Path(fp.parent.name) / fp.name}\n")


main(snakemake)  # type: ignore
