from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import write_bed, filter_sort_bed, FinalMapper, InitMapper


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    i: Path = smk.input[0]
    o: Path = smk.output[0]

    def read_write(b: cfg.BenchT, im: InitMapper, fm: FinalMapper) -> None:
        df = b.bench_bed.read(i)
        df_ = filter_sort_bed(im, fm, df)
        write_bed(o, df_)

    def hap_f(bd: cfg.HaploidBuildData) -> None:
        b = bd.build.bench
        assert b is not None
        conv = bd.chr_conversion(b.bench_bed.data.chr_pattern)
        read_write(b, conv.init_mapper, conv.final_mapper)

    def dip1_f(bd: cfg.Diploid1BuildData) -> None:
        b = bd.build.bench
        assert b is not None
        conv = bd.dip_chr_conversion(b.bench_bed.data.chr_pattern)
        read_write(b, conv.init_mapper, conv.final_mapper)

    def dip2_f(hap: cfg.Haplotype, bd: cfg.Diploid2BuildData) -> None:
        b = bd.build.bench
        assert b is not None
        conv = hap.from_either(*bd.hap_chr_conversion(b.bench_bed.data.chr_pattern))
        read_write(b, conv.init_mapper, conv.final_mapper)

    sconf.with_build_data(
        ws["ref_key"],
        ws["build_key"],
        cfg.to_haplotype(ws["hap"]),
        hap_f,
        dip1_f,
        dip2_f,
    )


main(snakemake, snakemake.config)  # type: ignore
