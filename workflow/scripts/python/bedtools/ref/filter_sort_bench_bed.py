from typing import Any
from pathlib import Path
import common.config as cfg
from common.bed import write_bed, filter_sort_bed, FinalMapper, InitMapper


def main(smk: Any, sconf: cfg.GiabStrats) -> None:
    ws: dict[str, str] = smk.wildcards
    # bd = sconf.to_build_data()
    # hap =
    # bench = bd.build.bench
    # assert bench is not None, "this should not happen"
    # b = bench.bench_bed
    i: Path = smk.input[0]
    o: Path = smk.output[0]

    def read_write(
        bd: cfg.BuildData_[cfg.W, cfg.StratInputT, cfg.BuildT],
        im: InitMapper,
        fm: FinalMapper,
    ) -> None:
        b = bd.build.bench
        assert b is not None
        df = b.bench_bed.read(i)
        df_ = filter_sort_bed(im, fm, df)
        write_bed(o, df_)

    def hap_f(bd: cfg.HaploidBuildData) -> None:
        pat0 = bd.build.bench.bench_bed.data.chr_pattern
        conv = bd.chr_conversion(pat0)
        read_write(bd, conv.init_mapper, conv.final_mapper)

    def dip1_f(bd: cfg.Diploid1BuildData) -> None:
        pass

    def dip2_f(hap: cfg.Haplotype, bd: cfg.Diploid2BuildData) -> None:
        pass

    sconf.with_build_data(
        ws["ref_key"],
        ws["build_key"],
        cfg.to_haplotype(ws["hap"]),
        hap_f,
        dip1_f,
        dip2_f,
    )

    # conv = cfg.HapToHapChrConversion(
    #     b.data.chr_pattern, bd.refdata.ref.chr_pattern, bd.chr_indices
    # )


main(snakemake, snakemake.config)  # type: ignore
