#! /usr/bin/env python3

import re
import sys
from enum import Enum
import argparse
from typing import TypedDict, Callable
from textwrap import indent
import yaml  # type: ignore


BASEURL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/"
    "release/genome-stratifications/v3.1"
)


class Ref(Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


class Url(TypedDict):
    url: str
    md5: str


class ChrPattern(TypedDict):
    template: str


class Src(TypedDict):
    hap: Url
    chr_pattern: ChrPattern


class Bed(TypedDict):
    bed: Src


class Entry(TypedDict):
    data: Bed
    description: str
    remove_gaps: bool


Hashes = dict[str, str]


def to_entry(
    subdir: str,
    key: str,
    desc: str,
    ref: Ref,
    md5s: Hashes,
) -> tuple[str, Entry]:
    url = f"{BASEURL}/{ref.value}/{subdir}/{ref.value}_{key}.bed.gz"
    return (
        key,
        {
            "data": {
                "bed": {
                    "hap": {
                        "url": str(url),
                        "md5": md5s[f"{subdir}/{key}"],
                    },
                    "chr_pattern": {
                        "template": "%i" if ref is Ref.GRCH37 else "chr%i",
                    },
                },
            },
            "description": desc,
            "remove_gaps": True,
        },
    )


def read_md5s(r: Ref) -> Hashes:
    with open(f"config/{r.value}_genome_specific_md5.tsv", "r") as f:
        return {(s := i.strip().split("\t"))[0]: s[1] for i in f}


def to_genome_specific(ref: Ref, md5s: Hashes) -> dict[str, Entry]:

    def go(key: str, desc: str) -> tuple[str, Entry]:
        return to_entry("GenomeSpecific", key, desc, ref, md5s)

    all_genomes = [
        go(f"{g}_v4.2.1_{t}", d)
        for g in [f"HG00{i}" for i in range(1, 8)]
        for (t, d) in [
            (
                "CNV_CCSandONT_elliptical_outlier",
                (
                    f"Potential duplications in {g} relative to the "
                    "reference, detected as higher than normal coverage in "
                    "PacBio CCS and/or ONT."
                ),
            ),
            (
                "CNV_mrcanavarIllumina_CCShighcov_ONThighcov_intersection",
                (
                    f"Potential duplications in {g} relative to the "
                    "reference, detected as higher than normal coverage in "
                    "PacBio CCS and/or ONT and as segmental duplications by "
                    "mrcanavar from Illumina."
                ),
            ),
            (
                "comphetindel10bp_slop50",
                (
                    f"Regions in {g} containing at least one variant on each "
                    "haplotype within 10bp of each other, and at least one of "
                    "the variants is an INDEL, with 50bp slop added on each "
                    "side."
                ),
            ),
            (
                "comphetsnp10bp_slop50",
                (
                    f"Regions in {g} containing at least one variant on each "
                    "haplotype within 10bp of each other, and all variants "
                    "are SNVs, with 50bp slop added on each side."
                ),
            ),
            (
                "complexindel10bp_slop50",
                f"Regions in {g} containing at least two variants on one "
                "haplotype within 10bp of each other, and at least one of the "
                "variants is an INDEL, with 50bp slop added on each side.",
            ),
            (
                "notin_complexandSVs_alldifficultregions",
                f"Complement of the above for {g}",
            ),
            (
                "othercomplexwithin10bp_slop50",
                f"Any other regions in {g} containing at least two variants "
                "within 10bp of each other, with 50bp slop added on each "
                "side.",
            ),
            ("snpswithin10bp_slop50", "FIXME"),
            (
                "CNVsandSVs",
                f"Union of the CNV and SV bed files above for {g}.",
            ),
            (
                "complexandSVs",
                "Union of the above SV, CNV, complex, and compound "
                f"heterozygous variant bed files for {g}.",
            ),
            (
                "complexandSVs_alldifficultregions",
                f"Union of {ref.value}_alldifficultregions.bed.gz from "
                f"`OtherDifficult` and complexandSVs above for {g}.",
            ),
        ]
    ]

    hg002_only = [
        go(f"HG002_{t}", d)
        for (t, d) in [
            (
                "v4.2.1_Tier1plusTier2_v0.6.1",
                "Regions containing HG002 v0.6 Tier1 or Tier2 insertions or "
                "deletions >=50bp, expanded to include overlapping tandem "
                "repeats, with and without expansion by 25% on each side.",
            ),
            (
                "v4.2.1_Tier1plusTier2_v0.6.1_slop25percent",
                "The above with 25% slop added to each side of each region.",
            ),
            (
                "v4.2.1_CNV_gt2assemblycontigs_ONTCanu_ONTFlye_CCSCanu",
                "Potential duplications relative to the reference, detected "
                "as more than 2 contigs aligning in 3 ONT and CCS Trio-binned "
                "assemblies of HG002",
            ),
            (
                "hifiasmv0.11_ComplexVar_in_TRgt100",
                "Complex variants in tandem repeats > 100 bp in HG002",
            ),
        ]
    ]

    parents_only = [
        go(
            f"HG00{g}_v4.2.1_SV_pbsv_slop25percent",
            "FIXME",
        )
        for g in [3, 4, 6, 7]
    ]

    children_only = [
        go(
            f"HG00{g}_v4.2.1_inversions_slop25percent",
            (
                "Putative inversions detected in either haplotype of the "
                "trio-hifiasm assembly using svanalyzer, including regions of "
                "breakpoint homology, expanded by 25% of the region size on "
                "each side"
            ),
        )
        for g in [1, 2, 5]
    ]

    misc = [
        go(
            f"HG00{g}_v4.2.1_SV_pbsv_hifiasm_dipcall_svanalyzer_slop25percent",
            "FIXME",
        )
        for g in [1, 5]
    ]
    return dict(
        [
            *all_genomes,
            *hg002_only,
            *parents_only,
            *children_only,
            *misc,
        ]
    )


def to_functional_tech_diff(ref: Ref, md5s: Hashes) -> dict[str, Entry]:
    genes = (
        "MRC1 and part of CNR2" if ref is Ref.GRCH37 else "CBS, CRYAA, KCNE1, and H19"
    )
    es = [
        to_entry("FunctionalTechnicallyDifficultRegions", *x, ref, md5s)
        for x in [
            (
                "BadPromoters",
                "Identified transcription start sites or first exons that "
                "have systematically low coverage as described in "
                '["Characterizing and measuring bias in sequence data"]'
                "(https://doi.org/10.1186/gb-2013-14-5-r51)",
            ),
            (
                "CMRGv1.00_duplicationinKMT2C",
                "Regions of KMT2C that are duplicated in most or all "
                f"individuals relative to {ref.value}.",
            ),
            (
                "CMRGv1.00_falselyduplicatedgenes",
                f"Genes that are falsely duplicated in {ref.value} ({genes}).",
            ),
        ]
    ]
    return dict(es)


def to_ancestry(ref: Ref, md5s: Hashes) -> dict[str, Entry]:
    # NOTE the subdirectory for v3.1 is lower-case
    es = [
        to_entry(
            "ancestry",
            f"ancestry_{a}",
            f"Regions in {aa} ancestry that closely match GRCh38",
            ref,
            md5s,
        )
        for (a, aa) in [
            ("AFR", "African"),
            ("AMR", "American"),
            ("EAS", "East-Asian"),
            ("EUR", "European"),
            ("Neanderthal", "Neanderthal-introgressed"),
            ("SAS", "South0Asian"),
        ]
    ]
    return dict(es)


def to_otherdifficult(ref: Ref, md5s: Hashes) -> dict[str, Entry]:
    def go(k: str, d: str) -> tuple[str, Entry]:
        return to_entry("OtherDifficult", k, d, ref, md5s)

    common = [
        go(*x)
        for x in [
            (
                "L1H_gt500",
                "L1Hs greater than 500 base pairs.",
            ),
            (
                "contigs_lt500kb",
                f"{ref.value} contigs smaller than 500kb.",
            ),
            (
                "allOtherDifficultregions",
                "Union of all other `OtherDifficult` regions.",
            ),
        ]
    ]

    def only37() -> list[tuple[str, Entry]]:
        return [
            go(*x)
            for x in [
                (
                    "hg38_minimap2_asm20_N10_gt1contig_gt1kb",
                    "GRCh37 regions covered by at least one contig from "
                    "GRCh38 using minimap2.",
                ),
                (
                    "hg38_minimap2_asm20_N10_nocovgt1kb",
                    "GRCh37 regions covered by no contigs from GRCh38 using "
                    "minimap2.",
                ),
                (
                    "hs37d5_decoy_alignments",
                    "Alignments of the hs37d5 decoy sequences to GRCh37, "
                    "potentially duplicated regions.",
                ),
                (
                    "missing_and_multiple_alignments_of_GRCh38",
                    "GRCh37 regions covered by >1 contig or no contigs from "
                    "GRCh38 as defined by GRC as SP or SPonly.",
                ),
            ]
        ]

    def only38() -> list[tuple[str, Entry]]:
        return [
            go(*x)
            for x in [
                (
                    "LD_discordant_haplotypes_slop5bp",
                    "Rare haplotye boundries in GRCh38.",
                ),
                (
                    "collapsed_duplication_FP_regions",
                    "Conservative collapsed errors with clusters of CHM13 "
                    "hets in GRCh38.",
                ),
                (
                    "false_duplications_correct_copy",
                    "Correct copy of falsely duplicated region.",
                ),
                (
                    "false_duplications_incorrect_copy",
                    "Incorrect copy of falsely duplicated region.",
                ),
                (
                    "gnomAD_InbreedingCoeff_slop1bp_merge1000bp",
                    "gnomAD inbreedingcoeff variants.",
                ),
                (
                    "population_CNV_FP_regions",
                    "Collapses in GRCh38 with clusters of CHM13 hets that are "
                    "variable in the population so many not errors.",
                ),
            ]
        ]

    return dict(common + (only37() if ref is Ref.GRCH37 else only38()))


def to_entries(r: Ref) -> str:
    md5s = read_md5s(r)
    return yaml.dump(
        {
            "GenomeSpecific": to_genome_specific(r, md5s),
            "FunctionalTechnicallyDifficult": to_functional_tech_diff(r, md5s),
            "OtherDifficult": to_otherdifficult(r, md5s),
            **({"Ancestry": to_ancestry(r, md5s)} if r is Ref.GRCH38 else {}),
        },
        sort_keys=False,
    )


def sub_pattern(suffix: str, r: Ref, s: str, f: Callable[[Ref], str]) -> str:
    pat = f" +\\${r.value}_{suffix}"

    if (m := re.search(pat, s)) is not None:
        spaces = len(m[0]) - len(m[0].lstrip())
    else:
        print(f"could not find {pat}")
        exit(1)

    return re.sub(pat, indent(f(r), " " * spaces), s)


def main() -> None:
    parser = argparse.ArgumentParser("Build the main config file")
    parser.add_argument("path", metavar="PATH", help="path to template")
    args = parser.parse_args()

    with open(args.path, "r") as f:
        config = f.read()

    for r in Ref:
        config = sub_pattern("external", r, config, to_entries)

    sys.stdout.write(config)


main()
