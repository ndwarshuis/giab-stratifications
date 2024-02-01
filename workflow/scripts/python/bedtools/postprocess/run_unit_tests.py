import gzip
from typing import Any, NamedTuple
from pathlib import Path
from os.path import dirname, basename
from os import scandir
import common.config as cfg
from common.io import setup_logging, is_bgzip
from common.bed import InternalChrIndex
import subprocess as sp

log = setup_logging(snakemake.log[0])  # type: ignore

RevMapper = dict[str, InternalChrIndex]


class GaplessBT(NamedTuple):
    auto: Path
    parY: Path


def test_bgzip(strat_file: Path) -> list[str]:
    """Each bed file should be bgzip'ed (not just gzip'ed)."""
    return [] if is_bgzip(strat_file) else ["is not bgzip file"]


def test_bed_format(strat_file: Path, reverse_map: RevMapper) -> str | None:
    with gzip.open(strat_file) as f:
        prevChrom = ""
        prevChromIdx = -1
        prevStart = -1
        prevEnd = -1
        for i in f:
            # should have three columns separated by tabs
            cells = i.split(b"\t")
            if len(cells) != 3:
                return "bed file has wrong number of columns"

            # types of each column should be str/int/int
            chrom = cells[0].decode()
            try:
                start = int(cells[1])
                end = int(cells[2])
            except ValueError as e:
                return str(e)

            # chrom should be valid
            try:
                chromIdx = reverse_map[chrom]
            except KeyError:
                return f"invalid chr: {chrom}"

            # chrom column should be sorted in ascending order
            if not prevChromIdx <= chromIdx:
                return "chrom column not sorted"

            # end should be greater than start
            if not start < end:
                return f"invalid region: {chrom} {start} {end}"

            # regions should be separated by at least one bp
            if not ((prevEnd < start) or (prevChromIdx < chromIdx)):
                return (
                    "non-disjoint regions: "
                    f"{prevChrom} {prevStart} {prevEnd}"
                    f"and {chrom} {start} {end}"
                )

            prevChrom = chrom
            prevChromIdx = chromIdx
            prevStart = start
            prevEnd = end

    # we made it, this file is legit :)
    return None


def test_bed_not_in_gaps(strat_file: Path, gapless: GaplessBT) -> bool:
    # ignore the actual gapless file for obvious reasons
    if "gaps_slop15kb" in basename(strat_file):
        return True
    else:
        isXY = basename(dirname(strat_file)) == "XY"
        gapless_path = gapless.parY if isXY else gapless.auto
        # gunzip needed here because subtract bed will output gibberish if we
        # give it an empty gzip file
        p0 = sp.Popen(["gunzip", "-c", strat_file], stdout=sp.PIPE)
        p1 = sp.run(
            ["subtractBed", "-a", "stdin", "-b", gapless_path, "-sorted"],
            stdin=p0.stdout,
            stdout=sp.PIPE,
        )
        return len(p1.stdout) == 0


def test_bed(
    strat_file: Path,
    reverse_map: RevMapper,
    gapless: GaplessBT,
) -> list[str]:
    if (format_error := test_bed_format(strat_file, reverse_map)) is not None:
        return [format_error]
    else:
        if test_bed_not_in_gaps(strat_file, gapless):
            return []
        else:
            return ["has gaps"]


def test_checksums(checksums: Path) -> list[str]:
    """The md5sum utility should pass when invoked."""
    res = sp.run(
        ["md5sum", "-c", "--strict", "--quiet", checksums.name],
        cwd=checksums.parent,
        capture_output=True,
        text=True,
    )
    errors = [] if (s := res.stdout) == "" else s.strip().split("\n")
    return [f"checksum error: {e}" for e in errors]


def test_tsv_list(tsv_list: Path) -> list[str]:
    """The final stratifications list should match final beds exactly.

    We are running our hotel on very tight margins; no extra or missing beds
    allowed.
    """

    def strat_set(root: Path, sub: Path) -> set[str]:
        if root.is_dir():
            deeper = (Path(p.path) for p in scandir(root))
            return {str(q) for d in deeper for q in strat_set(d, sub / d.name)}
        elif root.is_file() and root.name.endswith(".bed.gz"):
            return {str(sub)}
        else:
            return set()

    current = strat_set(tsv_list.parent, Path("./"))
    with open(tsv_list, "r") as f:
        listed = {line.strip().split("\t")[1] for line in f}
        missing = [f"not in final directory: {p}" for p in listed - current]
        extra = [f"not in final list: {p}" for p in current - listed]
        return missing + extra


def run_all_tests(
    strat_file: Path,
    reverse_map: RevMapper,
    gapless: GaplessBT,
) -> list[tuple[Path, str]]:
    return [
        (strat_file, msg)
        for msg in test_bgzip(strat_file) + test_bed(strat_file, reverse_map, gapless)
    ]


def strat_files(path: str) -> list[Path]:
    with open(path, "r") as f:
        return [Path(s.strip()) for s in f]


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    fm = sconf.buildkey_to_ref_mappers(
        cfg.wc_to_reffinalkey(ws),
        cfg.wc_to_buildkey(ws),
    )[1]
    reverse_map = {v: k for k, v in fm.items()}

    # check global stuff first (since this is faster)

    global_failures: list[str] = test_checksums(
        Path(smk.input["checksums"])
    ) + test_tsv_list(Path(smk.input["strat_list"]))

    for j in global_failures:
        log.error(j)

    if len(global_failures) > 0:
        exit(1)

    # check individual stratification files

    gapless = GaplessBT(
        auto=Path(smk.input["gapless_auto"]),
        parY=Path(smk.input["gapless_parY"]),
    )

    strat_failures = [
        res
        for p in strat_files(str(smk.input["strats"]))
        for res in run_all_tests(p, reverse_map, gapless)
    ]

    for i in strat_failures:
        log.error("%s: %s" % i)

    if len(strat_failures) > 0:
        exit(1)


main(snakemake, snakemake.config)  # type: ignore
