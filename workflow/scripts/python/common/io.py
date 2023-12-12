import hashlib
from typing import TypeVar, Callable
from logging import Logger
from pathlib import Path
from Bio import bgzf  # type: ignore
import gzip

X = TypeVar("X")
Y = TypeVar("Y")


class DesignError(Exception):
    """Exception raised when the code is designed incorrectly (ie the 'this
    should not happen' error)"""

    pass


def none_unsafe(x: X | None, f: Y, msg: None | str = None) -> Y:
    if x is not None:
        raise DesignError(msg if msg is not None else "Should not be None")
    return f


def not_none_unsafe(x: X | None, f: Callable[[X], Y], msg: None | str = None) -> Y:
    if x is None:
        raise DesignError(msg if msg is not None else "Should not be None")
    return f(x)


def match1_unsafe(xs: list[X], f: Callable[[X], Y], msg: None | str = None) -> Y:
    match xs:
        case [x]:
            return f(x)
        case _:
            raise DesignError(
                msg if msg is not None else f"One input expected, got {len(xs)}"
            )


def match2_unsafe(xs: list[X], f: Callable[[X, X], Y], msg: None | str = None) -> Y:
    match xs:
        case [x1, x2]:
            return f(x1, x2)
        case _:
            raise DesignError(
                msg if msg is not None else f"Two inputs expected, got {len(xs)}"
            )


def match12_unsafe(
    xs: list[X],
    f1: Callable[[X], Y],
    f2: Callable[[X, X], Y],
    msg: None | str = None,
) -> Y:
    match xs:
        case [x1]:
            return f1(x1)
        case [x1, x2]:
            return f2(x1, x2)
        case _:
            raise DesignError(
                msg if msg is not None else f"Two inputs expected, got {len(xs)}"
            )


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


def test_not_none(log: Logger, thing: X | None, what: str | None = None) -> X:
    if thing is None:
        log.error(
            "This should not happen: %s is None" % what
            if what is not None
            else "object"
        )
        exit(1)
    return thing


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
