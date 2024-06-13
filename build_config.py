#! /usr/bin/env python3

import re
from pathlib import Path
from typing import TypedDict
from textwrap import indent
import yaml  # type: ignore


REFS = ["GRCh37", "GRCh38"]

BASEURL = Path(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1"
)


class Url(TypedDict):
    url: str
    md5: str


class Bed(TypedDict):
    hap: Url


class ChrPattern(TypedDict):
    template: str


class Entry(TypedDict):
    bed: Bed
    remove_gaps: bool
    chr_pattern: ChrPattern
    description: str


def to_entries(ref: str) -> str:
    with open(f"config/{ref}_genome_specific_md5.tsv", "rt") as f:
        md5s = {
            Path((s := i.strip().split("\t"))[1])
            .name.replace(f"{ref}_", "")
            .replace(".bed.gz", ""): s[0]
            for i in f
        }

    # dump each of these individually to preserve order
    def to_entry(key: str, desc: str) -> str:
        url = BASEURL / ref / f"{ref}_{key}.bed.gz"
        x: dict[str, Entry] = {
            key: {
                "bed": {
                    "hap": {
                        "url": str(url),
                        "md5": md5s[key],
                    },
                },
                "chr_pattern": {"template": "%i" if ref == "GRCh37" else "chr%i"},
                "remove_gaps": True,
                "description": desc,
            }
        }
        return yaml.dump(x)

    all_genomes = [
        to_entry(f"{g}_v4.2.1_{t}", d)
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
                f"Union of {ref}_alldifficultregions.bed.gz from "
                f"`OtherDifficult` and complexandSVs above for {g}.",
            ),
        ]
    ]

    hg002_only = [
        to_entry(f"HG002_{t}", d)
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
        to_entry(f"HG00{g}_v4.2.1_SV_pbsv_slop25percent", "FIXME") for g in [3, 4, 6, 7]
    ]

    children_only = {
        to_entry(
            f"HG00{g}_v4.2.1_inversions_slop25percent",
            (
                "Putative inversions detected in either haplotype of the "
                "trio-hifiasm assembly using svanalyzer, including regions of "
                "breakpoint homology, expanded by 25% of the region size on "
                "each side"
            ),
        )
        for g in [1, 2, 5]
    }

    misc = {
        to_entry(
            f"HG00{g}_v4.2.1_SV_pbsv_hifiasm_dipcall_svanalyzer_slop25percent", "FIXME"
        )
        for g in [1, 5]
    }
    return "".join(
        [
            *all_genomes,
            *hg002_only,
            *parents_only,
            *children_only,
            *misc,
        ]
    )


with open("config/all.yml", "r") as f:
    config = f.read()

for r in REFS:
    pat = f" +\\${r}_genome_specific"

    if (m := re.search(pat, config)) is not None:
        spaces = len(m[0]) - len(m[0].lstrip())
    else:
        print(f"could not find {pat}")
        exit(1)

    config = re.sub(pat, indent(to_entries("GRCh38"), " " * spaces), config)

print(config)
