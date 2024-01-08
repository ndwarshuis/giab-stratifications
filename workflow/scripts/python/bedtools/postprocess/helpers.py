from pathlib import Path
import common.config as cfg


def write_chr_mapper(sconf: cfg.GiabStrats, o: Path) -> None:
    with open(o, "w") as f:
        for rk, bk in zip(*sconf.all_full_build_keys):
            cis = sconf.to_build_data(cfg.strip_full_refkey(rk), bk).chr_indices
            xs = sconf.with_ref_data_full(
                rk,
                lambda rd: [
                    (i, cfg.Haplotype.HAP1, rd.ref.chr_pattern.to_chr_name(i))
                    for i in cis
                ],
                lambda rd: [
                    (i, h, rd.ref.chr_pattern.to_chr_name(i, h))
                    for i in cis
                    for h in cfg.Haplotype
                ],
                lambda hap, rd: [
                    (i, hap, rd.ref.chr_pattern.from_either(hap).to_chr_name(i))
                    for i in cis
                ],
            )
            for i, h, n in xs:
                if n is not None:
                    # chr number, ref_final_key@build_key, haplotype, chr name
                    line = [str(i.value), f"{rk}@{bk}", str(h.value + 1), n]
                    f.write("\t".join(line) + "\n")
