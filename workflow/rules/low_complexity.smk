lc_inter_dir = intermediate_dir / "LowComplexity"
lc_final_dir = final_dir / "LowComplexity"


rule download_ref:
    output:
        resources_dir / "ref.fna.gz",
    params:
        url=config["ref_url"],
    conda:
        envs_path("utils.yml")
    shell:
        "curl -sS -L -o {output} {params.url}"


use rule download_ref as download_rmsk with:
    output:
        resources_dir / "rmsk.txt.gz",
    params:
        url=config["rmsk_url"],


use rule download_ref as download_simreps with:
    output:
        resources_dir / "simreps.txt.gz",
    params:
        url=config["simreps_url"],


rule unzip_ref:
    input:
        rules.download_ref.output,
    output:
        resources_dir / "ref.fna",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        grep -Ev '^chr[0-9XYM]_|^chr[0-9][0-9]_|^chrUn_|^chrEBV' \
        > {output}
        """


rule index_ref:
    input:
        rules.download_ref.output,
    output:
        resources_dir / "ref.fna.fai",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        samtools faidx - -o - | \
        grep -Ev '^chr[0-9XYM]_|^chr[0-9][0-9]_|^chrUn_|^chrEBV' > \
        {output}
        """


use rule download_ref as download_repseq with:
    output:
        resources_dir / "tools" / "repseq.tar.gz",
    params:
        url=config["tools"]["repseq"],


hp_grid = {
    1: [4, 6, 12, 20],
    2: [10, 50, 200],
    3: [14, 50, 200],
    4: [19, 50, 200],
}


rule find_repeats:
    input:
        ref=rules.unzip_ref.output,
        bin=rules.build_repseq.output,
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{len}.bed",
    log:
        log_dir / "low_complexity" / "repeats_R{rep}_L{len}.log",
    shell:
        "{input.bin} {rep} {len} {input.ref} > {output} 2> {log}"


rule subtract_repeats:
    input:
        a=results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}.bed",
        b=results_dir / "low_complexity" / "repeats_R{rep}_L{lenB}.bed",
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}to{lenB}.bed",
    shell:
        "subtractBed -a {input.a} -b {input.b} > {output}"


rule add_repeat_slop:
    input:
        bed=rules.subtract_repeats.output,
        # genome=TODO
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}to{lenB}.bed",


# TODO get multi-inter of this
rule get_imperfect_homopolymers:
    input:
        results_dir / "low_complexity" / "repeats_R1_L4.bed",
    output:
        results_dir / "low_complexity" / "repeats_imp_hp_{len}_{base}.bed",
    shell:
        """
        grep 'unit=C' {input} | \
        mergeBed -i stdin -d 1 | \
        awk '$3-$2>10' > \
        {output}
        """


rule get_simreps_inter:
    input:
        rules.download_simreps.output,
    output:
        lc_inter_dir / "simpleRepeats.bed.gz",
    shell:
        """
        gunzip -c {input} | \
        cut -f2-4  | \
        sed 's/^chr//' | \
        grep "^[0-9XY]" | grep -v '_' | \
        sed 's/^X/23/;s/^Y/24/' |  \
        sort -k1,1n -k2,2n -k3,3n |  \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' |  \
        mergeBed -i stdin |  \
        bgzip -c > {output}
        """


simreps = {
    "lt51": {"lower": 0, "upper": 51},
    "51to200": {"lower": 50, "upper": 201},
    "51to200": {"lower": 200, "upper": 1000000},
}


# TODO this might be easier/cleaner with pandas
rule get_simreps:
    input:
        rules.get_simreps_inter,
    output:
        lc_final_dir / "GRCh38_trf_simpleRepeat_{simrep_key}.bed.gz",
    params:
        lower=lambda wildcards: simreps[wildcards.simrep_key]["lower"],
        upper=lambda wildcards: simreps[wildcards.simrep_key]["upper"],
    shell:
        """
        gunzip -c {input} | \
        awk '$3-$2>{params.lower} && $3-$2<{params.upper}' | \
        bgzip -c > {output}
        """


rule get_satellites_inter:
    input:
        rmsk=rules.download_rmsk.output,
        idx=rules.index_ref.output,
    output:
        lc_inter_dir / "satellites.bed.gz",
    shell:
        """
        gzcat {input.rmsk} | \
        grep "Satellite" | \
        awk '{ print $6 "\t" $7 "\t" $8 ; }' | \
        grep -Ev '^chr[0-9XYM]_|^chr[0-9][0-9]_|^chrUn_' | \
        sortBed -faidx {input.idx} -i stdin | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


rule get_satellites:
    input:
        bed=rules.get_satellites_inter.rules,
        genome=TODO,
    output:
        lc_final_dir / "GRCh38_satellites_slop5.bed.gz",
    shell:
        """
        slopBed -i {input.bed} -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        bgzip -c  > {output}
        """


rule get_all_hp:
    input:
        p6=TODO,
        imp10=TODO,
        genome=TODO,
    output:
        lc_final_dir / "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    shell:
        """
        grep -Ev '_|^chrEBV' {input.p6} | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        multiIntersectBed -i stdin {input.imp10} | \
        sed 's/^chr//' | \
        cut -f1-3 | grep "^[0-9XY]" | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n -k2,2n -k3,3n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """


rule get_all_tr_inter:
    input:
        hp=[],
        genome=TODO,
    output:
        lc_inter_dir / "AllTandemRepeats_intermediate.bed",
    shell:
        """
        multiIntersectBed -i {input.hp} \
        sed 's/^chr//' | 
        cut -f1-3 | grep "^[0-9XY]" | grep -v '_' | 
        sed 's/^/chr/' | 
        slopBed -i stdin -b 5 -g {input.genome} | 
        sed 's/^chr//' | 
        sed 's/^X/23/;s/^Y/24/' | 
        sort -k1,1n -k2,2n -k3,3n | 
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' | 
        mergeBed -i stdin > {output}
        """


# generalize this so I don't need to repeat awk alot
rule get_all_tr:
    input:
        tr=TODO,
        hp=TODO,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeats_gt10000bp_slop5.bed.gz",
    shell:
        """
        awk '$3-$2>10010' {input.tr} | \
        subtractBed -a stdin -b {input.hp} | \
        bgzip -c > {output}
        """


rule get_all_tr_and_hp:
    input:
        rules.get_all_tr.output,
        rules.get_all_hp_p6_imp10.output,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
    shell:
        """
        multiIntersectBed -i {input} | \
        sed 's/^chr//' | \
        cut -f1-3 | grep "^[0-9XY]" | grep -v '_' | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n -k2,2n -k3,3n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' | \
        mergeBed -i stdin | \
        bgzip -c > {output}
        """
