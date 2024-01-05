import common.config as cfg
from common.functional import DesignError


def seqtk_args(
    sconf: cfg.GiabStrats,
    rfk: cfg.RefKeyFullS,
    bk: cfg.BuildKey,
    frac: str,
) -> str:
    try:
        _frac = int(frac)
    except ValueError:
        raise DesignError(f"GC fraction is not an int; got {frac}")

    gps = sconf.to_build_data(cfg.strip_full_refkey(rfk), bk).build.include.gc

    if gps is None:
        raise DesignError("GC include not available")

    if _frac in gps.low_fractions:
        s, f = ("w", 100 - _frac)
    elif _frac in gps.high_fractions:
        s, f = ("", _frac)
    else:
        raise DesignError(
            f"GC fraction is not valid; should be in {gps.low_fractions} or {gps.high_fractions}"
        )

    return f"-{s}f 0.{f}"


def range_bounds(
    sconf: cfg.GiabStrats, rfk: cfg.RefKeyFullS, bk: cfg.BuildKey
) -> tuple[int, list[int], list[int], int]:
    gps = sconf.to_build_data(cfg.strip_full_refkey(rfk), bk).build.include.gc

    if gps is None:
        raise DesignError("GC include not available")

    lowest, lower = gps.low_bounds
    highest, higher = gps.high_bounds
    return (lowest, lower, higher, highest)
