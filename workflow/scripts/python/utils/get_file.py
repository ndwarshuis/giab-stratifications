from pathlib import Path
import common.config as cfg
import subprocess as sp

# hacky curl/gzip wrapper; this exists because I got tired of writing
# specialized rules to convert gzip/nozip files to bgzip and back :/
# Solution: just use bgzip for everything

GUNZIP = ["gunzip", "-c"]
BGZIP = ["bgzip", "-c"]


def main(opath: str, src: cfg.AnySrc | None) -> None:
    # ASSUME anything here is bgzip compressed
    if isinstance(src, cfg.FileSrc):
        Path(opath).symlink_to(src.filepath)
    elif isinstance(src, cfg.HttpSrc):
        curlcmd = ["curl", "-Ss", "-L", "-q", src.url]
        with open(opath, "wb") as f:
            if src.zipfmt == cfg.ZipFmt.BGZIP:
                sp.Popen(curlcmd, stdout=f).wait()
            elif src.zipfmt == cfg.ZipFmt.GZIP:
                p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                p2 = sp.Popen(GUNZIP, stdin=p1.stdout, stdout=sp.PIPE)
                p3 = sp.Popen(BGZIP, stdin=p2.stdout, stdout=f)
                p3.wait()
            elif src.zipfmt == cfg.ZipFmt.NOZIP:
                p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                p2 = sp.Popen(BGZIP, stdin=p1.stdout, stdout=f)
                p2.wait()
    else:
        assert False, "file src is null"


main(snakemake.output[0], snakemake.params.src)  # type: ignore
