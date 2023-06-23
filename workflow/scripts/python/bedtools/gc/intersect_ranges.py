from pathlib import Path
from typing import Any, NamedTuple, Callable
import subprocess as sp
import common.config as cfg
import json

# use subprocess here because I don't see an easy way to stream a bedtools
# python object to a bgzip file (it only seems to do gzip...oh well)


class GCInput(NamedTuple):
    bed: Path
    fraction: int
    is_range_bound: bool


def write_simple_range_beds(
    final_path: Callable[[str], Path],
    gs: list[GCInput],
    is_low: bool,
) -> list[str]:
    def fmt_out(bigger_frac: int, smaller_frac: int) -> Path:
        lower_frac, upper_frac = (
            (smaller_frac, bigger_frac) if is_low else (bigger_frac, smaller_frac)
        )
        return final_path(f"gc{lower_frac}to{upper_frac}_slop50")

    pairs = zip(gs[1:], gs[:-1]) if is_low else zip(gs[:-1], gs[1:])
    torun = [
        (fmt_out(bigger.fraction, smaller.fraction), bigger, smaller)
        for bigger, smaller in pairs
    ]
    for out, bigger, smaller in torun:
        with open(out, "wb") as f:
            p0 = sp.Popen(
                ["subtractBed", "-a", bigger.bed, "-b", smaller.bed],
                stdout=sp.PIPE,
            )
            p1 = sp.run(["bgzip", "-c"], stdin=p0.stdout, stdout=f)
            p0.wait()
            if not (p0.returncode == p1.returncode == 0):
                exit(1)
    return [str(t[0]) for t in torun]


def write_middle_range_bed(
    final_path: Callable[[str], Path],
    lower: GCInput,
    upper: GCInput,
    genome: Path,
    gapless: Path,
) -> str:
    out = final_path(f"gc{lower.fraction}to{upper.fraction}_slop50")
    with open(out, "wb") as f:
        p0 = sp.Popen(
            ["complementBed", "-i", upper.bed, "-g", genome],
            stdout=sp.PIPE,
        )
        p1 = sp.Popen(
            ["subtractBed", "-a", "stdin", "-b", lower.bed],
            stdin=p0.stdout,
            stdout=sp.PIPE,
        )
        p2 = sp.Popen(
            ["intersectBed", "-a", "stdin", "-b", gapless, "-sorted", "-g", genome],
            stdin=p1.stdout,
            stdout=sp.PIPE,
        )
        p3 = sp.run(["bgzip", "-c"], stdin=p2.stdout, stdout=f)
        p0.wait()
        p1.wait()
        p2.wait()
        if not (p0.returncode == p1.returncode == p2.returncode == p3.returncode == 0):
            exit(1)
    return str(out)


def write_intersected_range_beds(
    final_path: Callable[[str], Path],
    low: list[GCInput],
    high: list[GCInput],
) -> list[str]:
    pairs = zip(
        [x for x in low if x.is_range_bound],
        [x for x in reversed(high) if x.is_range_bound],
    )
    torun = [
        (i, final_path(f"gclt{b1.fraction}orgt{b2.fraction}_slop50"), b1, b2)
        for i, (b1, b2) in enumerate(pairs)
    ]
    for i, bed_out, b1, b2 in torun:
        with open(bed_out, "wb") as f:
            p0 = sp.Popen(
                ["multiIntersectBed", "-i", b1.bed, b2.bed],
                stdout=sp.PIPE,
            )
            p1 = sp.Popen(
                ["mergeBed", "-i", "stdin"],
                stdin=p0.stdout,
                stdout=sp.PIPE,
            )
            p2 = sp.run(["bgzip", "-c"], stdin=p1.stdout, stdout=f)
            p0.wait()
            p1.wait()
            if not (p0.returncode == p1.returncode == p2.returncode == 0):
                exit(1)

    return [str(t[1]) for t in torun]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    rk = cfg.RefKey(smk.wildcards.ref_key)
    bk = cfg.BuildKey(smk.wildcards.build_key)
    gps = sconf.buildkey_to_include(rk, bk).gc
    assert gps is not None, "this should not happen"
    # ASSUME both of these input lists are sorted by GC fraction
    low = [GCInput(p, f, r) for p, (f, r) in zip(smk.input.low, gps.low)]
    high = [GCInput(p, f, r) for p, (f, r) in zip(smk.input.high, gps.high)]
    genome = Path(smk.input.genome)
    gapless = Path(smk.input.gapless)

    def final_path(name: str) -> Path:
        p = Path(str(smk.params.path_pattern).format(name))
        p.parent.mkdir(exist_ok=True, parents=True)
        return p

    low_strats = write_simple_range_beds(final_path, low, True)
    high_strats = write_simple_range_beds(final_path, high, False)
    range_strat = write_middle_range_bed(
        final_path,
        low[-1],
        high[0],
        genome,
        gapless,
    )
    inter_strats = write_intersected_range_beds(
        final_path,
        low,
        high,
    )

    with open(smk.output[0], "w") as f:
        # put the first low and last high input here since these are already
        # in the final directory
        obj = {
            "gc_ranges": [
                low[0][0],
                *low_strats,
                range_strat,
                *high_strats,
                high[-1][0],
            ],
            "widest_extreme": inter_strats[0],
            "other_extremes": inter_strats[1:],
        }
        json.dump(obj, f, indent=4)


main(snakemake, snakemake.config)  # type: ignore
