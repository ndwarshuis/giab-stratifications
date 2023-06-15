from pathlib import Path
from typing import Any, NamedTuple, Callable
import subprocess as sp
import common.config as cfg

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
) -> None:
    pairs = zip(gs[1:], gs[:-1]) if is_low else zip(gs[:-1], gs[1:])
    for bigger, smaller in pairs:
        lower_frac, upper_frac = (
            (smaller.fraction, bigger.fraction)
            if is_low
            else (bigger.fraction, smaller.fraction)
        )
        out = final_path(f"gc{lower_frac}to{upper_frac}_slop50")
        with open(out, "wb") as f:
            p0 = sp.Popen(
                ["subtractBed", "-a", bigger.bed, "-b", smaller.bed],
                stdout=sp.PIPE,
            )
            p1 = sp.run(["bgzip", "-c"], stdin=p0.stdout, stdout=f)
            p0.wait()
            if not (p0.returncode == p1.returncode == 0):
                exit(1)


def write_middle_range_bed(
    final_path: Callable[[str], Path],
    lower: GCInput,
    upper: GCInput,
    genome: Path,
    gapless: Path,
) -> None:
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


def write_intersected_range_beds(
    final_path: Callable[[str], Path],
    low: list[GCInput],
    high: list[GCInput],
    wider_out: Path,
) -> None:
    pairs = zip(
        [x for x in low if x.is_range_bound],
        [x for x in reversed(high) if x.is_range_bound],
    )
    for i, (b1, b2) in enumerate(pairs):
        bed_out = final_path(f"gclt{b1.fraction}orgt{b2.fraction}_slop50")
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

        # This is necessary because unions require the widest of the meta-range
        # GC intersection beds. This is hacky AF but easiest way to do this is
        # to "return" a link that points to the file I want. Need to link b/c
        # snakemake doesn't allow me to dynamically declare output files in a
        # rule, which is basically what this symlink is simulating.
        if i == 0:
            Path(wider_out).symlink_to(Path(bed_out).resolve())


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

    write_simple_range_beds(final_path, low, True)
    write_simple_range_beds(final_path, high, False)
    write_middle_range_bed(final_path, low[-1], high[0], genome, gapless)
    write_intersected_range_beds(final_path, low, high, Path(smk.output[0]))


main(snakemake, snakemake.config)  # type: ignore
