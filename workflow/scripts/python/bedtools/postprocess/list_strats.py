from typing import Any


def main(smk: Any) -> None:
    with open(smk.output[0], "w") as f:
        for s in sorted(smk.input):
            f.write(s + "\n")


main(snakemake)  # type: ignore
