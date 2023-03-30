from typing import TypeVar
from logging import Logger

X = TypeVar("X")


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
