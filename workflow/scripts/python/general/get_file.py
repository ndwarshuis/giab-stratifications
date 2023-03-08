from pathlib import Path
import requests as req  # type: ignore
from typing import Any
import common.config as cfg
from Bio import bgzf  # type: ignore
import zlib


def main(smk: Any, src: cfg.FileSrc | cfg.HttpSrc | None) -> None:
    o = smk.output[0]
    # ASSUME anything here is bgzip compressed
    if isinstance(src, cfg.FileSrc):
        Path(o).symlink_to(src.filepath)
    elif isinstance(src, cfg.HttpSrc):
        res = req.get(src.url, stream=True)
        sc = res.status_code
        assert sc == 200, f"status code was {sc}"
        if src.zipfmt == cfg.ZipFmt.BGZIP:
            with open(o, "wb") as f:
                f.write(res.raw.read())
        elif src.zipfmt == cfg.ZipFmt.GZIP:
            with bgzf.writer(o, "wb") as f:
                # TODO this is probably wrong
                f.write(zlib.decompress(res.raw.read(), 16 + 31))
        elif src.zipfmt == cfg.ZipFmt.NOZIP:
            with bgzf.writer(o, "wb") as f:
                f.write(res.raw.read())
    else:
        assert False, "file src is null"


main(snakemake, snakemake.params.url)  # type: ignore
