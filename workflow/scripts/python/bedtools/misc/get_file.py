from pathlib import Path
import subprocess as sp
from typing import Callable
from typing_extensions import assert_never
from tempfile import NamedTemporaryFile as Tmp
from common.io import is_bgzip, is_gzip
import common.config as cfg
from common.io import get_md5, setup_logging

# hacky curl/gzip wrapper; this exists because I got tired of writing
# specialized rules to convert gzip/nozip files to bgzip and back :/
# Solution: force bgzip for references and gzip for bed


log = setup_logging(snakemake.log[0])  # type: ignore

GUNZIP = ["gunzip", "-c"]
BGZIP = ["bgzip", "-c"]
GZIP = ["gzip", "-c"]
CURL = ["curl", "-Ss", "-L", "-q"]


def main(opath: str, src: cfg.BedSrc | cfg.RefSrc | None) -> None:
    if isinstance(src, cfg.FileSrc_):
        # ASSUME these are already tested via the pydantic class for the
        # proper file format
        Path(opath).symlink_to(Path(src.filepath).resolve())

    elif isinstance(src, cfg.HttpSrc_):
        curlcmd = [*CURL, src.url]

        # to test the format of downloaded files, sample the first 65000 bytes
        # (which should be enough to get one block of a bgzip file, which will
        # allow us to test for it)
        curltestcmd = [*CURL, "-r", "0-65000", src.url]

        with open(opath, "wb") as f, Tmp() as tf:

            def curl() -> None:
                sp.Popen(curlcmd, stdout=f).wait()

            def curl_test(testfun: Callable[[Path], bool]) -> bool:
                sp.Popen(curltestcmd, stdout=tf).wait()
                return testfun(Path(tf.name))

            def curl_gzip(cmd: list[str]) -> None:
                p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                p2 = sp.Popen(cmd, stdin=p1.stdout, stdout=f)
                p2.wait()

            # if we are getting a bed file (or related) ensure it is in gzip
            # format
            if isinstance(src, cfg.BedHttpSrc):
                if curl_test(is_gzip):
                    curl()
                else:
                    curl_gzip(GZIP)

            # if we are getting a fasta file ensure it is in bgzip format
            elif isinstance(src, cfg.RefHttpSrc):
                if curl_test(is_bgzip):
                    curl()
                elif curl_test(is_gzip):
                    p1 = sp.Popen(curlcmd, stdout=sp.PIPE)
                    p2 = sp.Popen(GUNZIP, stdin=p1.stdout, stdout=sp.PIPE)
                    p3 = sp.Popen(BGZIP, stdin=p2.stdout, stdout=f)
                    p3.wait()
                else:
                    curl_gzip(BGZIP)
            else:
                assert_never(src)

    elif src is None:
        assert False, "file src is null; this should not happen"
    else:
        assert_never(src)

    if src.md5 is not None and src.md5 != (actual := get_md5(opath)):
        log.error("md5s don't match; wanted %s, actual %s", src.md5, actual)
        exit(1)


main(snakemake.output[0], snakemake.params.src)  # type: ignore
