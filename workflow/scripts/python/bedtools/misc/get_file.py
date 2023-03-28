from pathlib import Path
import common.config as cfg
import subprocess as sp
from typing_extensions import assert_never

# hacky curl/gzip wrapper; this exists because I got tired of writing
# specialized rules to convert gzip/nozip files to bgzip and back :/
# Solution: force bgzip for references and gzip for bed

GUNZIP = ["gunzip", "-c"]
BGZIP = ["bgzip", "-c"]
GZIP = ["gzip", "-c"]


def main(opath: str, src: cfg.BedSrc | cfg.RefSrc | None) -> None:
    if isinstance(src, cfg.FileSrc_):
        Path(opath).symlink_to(Path(src.filepath).resolve())
    elif isinstance(src, cfg.HttpSrc_):
        curlcmd = ["curl", "-Ss", "-L", "-q", src.url]

        with open(opath, "wb") as f:

            def curl() -> None:
                sp.Popen(curlcmd, stdout=f).wait()

            def curl_gzip(cmd: list[str]) -> None:
                p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                p2 = sp.Popen(cmd, stdin=p1.stdout, stdout=f)
                p2.wait()

            if isinstance(src, cfg.BedHttpSrc):
                if src.fmt is cfg.BedFmt.GZIP:
                    curl()
                elif src.fmt is cfg.BedFmt.NOZIP:
                    curl_gzip(GZIP)
                else:
                    assert_never(src.fmt)
            elif isinstance(src, cfg.RefHttpSrc):
                if src.fmt is cfg.RefFmt.BGZIP:
                    curl()
                elif src.fmt is cfg.RefFmt.NOZIP:
                    curl_gzip(BGZIP)
                elif src.fmt is cfg.RefFmt.GZIP:
                    p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                    p2 = sp.Popen(GUNZIP, stdin=p1.stdout, stdout=sp.PIPE)
                    p3 = sp.Popen(BGZIP, stdin=p2.stdout, stdout=f)
                    p3.wait()
                else:
                    assert_never(src.fmt)
            else:
                assert_never(src)

    elif src is None:
        assert False, "file src is null; this should not happen"
    else:
        assert_never(src)


main(snakemake.output[0], snakemake.params.src)  # type: ignore
