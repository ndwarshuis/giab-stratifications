import hashlib
from logging import Logger
from pathlib import Path
from Bio import bgzf  # type: ignore
import gzip


def get_md5(path: str, unzip: bool = False) -> str:
    h = hashlib.md5()
    do_unzip = path.endswith(".gz") and unzip is True
    with gzip.open(path, "rb") if do_unzip else open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


def setup_logging(path: str, console: bool = False) -> Logger:
    import logging

    logging.basicConfig(filename=path, level=logging.INFO)
    logging.captureWarnings(True)
    logger = logging.getLogger()
    if console:
        logger.addHandler(logging.StreamHandler())
    return logger


def is_gzip(p: Path) -> bool:
    # test if gzip by trying to read first byte
    with gzip.open(p, "r") as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False


def is_bgzip(p: Path) -> bool:
    # since bgzip is in blocks (vs gzip), determine if in bgzip by
    # attempting to seek first block
    with open(p, "rb") as f:
        try:
            next(bgzf.BgzfBlocks(f), None)
            return True
        except ValueError:
            return False
