from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import write_bed, filter_sort_bed, FinalMapper, InitMapper


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    # ws: dict[str, str] = smk.wildcards
    # i: Path = smk.input[0]
    # o: Path = smk.output[0]

    ws: dict[str, str] = smk.wildcards
    ps: dict[str, str] = smk.params
    ins: list[Path] = smk.input
    output_list = smk.output[0]
    bed_outputs = list(map(Path, ps["bed_outputs"]))

    # def read_write(b: cfg.BenchT, im: InitMapper, fm: FinalMapper) -> None:
    #     df = b.bench_bed.read(i)
    #     df_ = filter_sort_bed(im, fm, df)
    #     write_bed(o, df_)

    # def hap_f(bd: cfg.HaploidBuildData) -> None:
    #     b = bd.build.bench
    #     assert b is not None
    #     conv = bd.chr_conversion(b.bench_bed.data.chr_pattern)
    #     read_write(b, conv.init_mapper, conv.final_mapper)

    # def dip1_f(bd: cfg.Diploid1BuildData) -> None:
    #     b = bd.build.bench
    #     assert b is not None
    #     conv = bd.dip_chr_conversion(b.bench_bed.data.chr_pattern)
    #     read_write(b, conv.init_mapper, conv.final_mapper)

    # def dip2_f(hap: cfg.Haplotype, bd: cfg.Diploid2BuildData) -> None:
    #     b = bd.build.bench
    #     assert b is not None
    #     conv = hap.from_either(*bd.hap_chr_conversion(b.bench_bed.data.chr_pattern))
    #     read_write(b, conv.init_mapper, conv.final_mapper)

    # sconf.with_build_data(
    #     ws["ref_key"],
    #     ws["build_key"],
    #     cfg.to_haplotype(ws["hap"]),
    #     hap_f,
    #     dip1_f,
    #     dip2_f,
    # )

    # TODO flatten builddata to make this sort of thing easier (and type the
    # types with less typing)
    def go(
        x: cfg.BuildData_[cfg.RefSourceT, cfg.AnyBedT, cfg.AnyBedT_, cfg.IncludeT]
    ) -> cfg.BedFile[cfg.AnyBedT] | None:
        return cfg.fmap_maybe(lambda y: y.bench_bed, x.build.bench)

    sconf.with_build_data_and_bed_io_(
        ws["ref_final_key"],
        ws["build_key"],
        ins,
        bed_outputs,
        go,
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_dip_bed(b, i, o),
        lambda i, o, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o),
        lambda i, o, hap, bd, b: bd.read_write_filter_sort_hap_bed(b, i, o, hap),
    )


main(snakemake, snakemake.config)  # type: ignore
