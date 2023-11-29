from typing import Any
import common.config as cfg


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    def hap_to_bed(
        si: cfg.HaploidStratInputs | cfg.DiploidStratInputs,
    ) -> (
        cfg.BedFile_[cfg.HapChrSource[cfg.BedSrc]]
        | cfg.BedFile_[cfg.DipChrSource1[cfg.BedSrc]]
        | cfg.BedFile_[cfg.DipChrSource2[cfg.BedSrc]]
        | None
    ):
        return si.low_complexity.simreps

    def hap_to_bed(
        si: cfg.HaploidStratInputs,
    ) -> cfg.BedFile[cfg.HapChrSource[cfg.BedSrc]] | None:
        return si.low_complexity.simreps

    def dip_to_bed(
        si: cfg.DiploidStratInputs,
    ) -> cfg.BedFile[cfg.DipChrSource[cfg.BedSrc]] | None:
        return si.low_complexity.simreps

    ws: dict[str, str] = smk.wildcards
    sconf.with_build_data_unsafe(
        ws["ref_key"],
        ws["build_key"],
        smk.input,
        smk.output,
        cfg.to_haplotype(ws["hap"]),
        lambda i, o, bd: bd.read_filter_sort_hap_bed(hap_to_bed, i, o),
        lambda i, o, bd: bd.read_filter_sort_dip_bed(dip_to_bed, i, o),
        lambda i, o, hap, bd: bd.read_filter_sort_hap_bed(dip_to_bed, i, o, hap),
        lambda i, o, bd: bd.read_filter_sort_dip_bed(dip_to_bed, i, o),
        lambda i, o, bd: bd.read_filter_sort_hap_bed(dip_to_bed, i, o),
    )


main(snakemake, snakemake.config)  # type: ignore
