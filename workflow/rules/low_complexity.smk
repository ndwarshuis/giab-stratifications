lc_inter_dir = intermediate_dir / "LowComplexity"
lc_final_dir = final_dir / "LowComplexity"
lc_log_dir = log_dir / "LowComplexity"


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
        "gunzip -c {input} > {output}"


rule index_ref:
    input:
        rules.download_ref.output,
    output:
        ref_dir / "ref.fna.fai",
    conda:
        envs_path("utils.yml")
    shell:
        """
        gunzip -c {input} | \
        samtools faidx - -o - | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        sed 's/^chr//' | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n -k2,2n -k3,3n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' > \
        {output}
        """


rule get_genome:
    input:
        rules.index_ref.output,
    output:
        ref_dir / "genome.txt",
    shell:
        "cut -f 1,2 {input} > {output}"


rule filter_sort_rmsk:
    input:
        bed=rules.download_rmsk.output,
        genome=rules.get_genome.output,
    output:
        ref_dir / "rmsk.txt.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input.bed} | \
        cut -f 6,7,8,12 | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        bedtools sort -i stdin -g {input.genome} | \
        gzip -c > {output}
        """


rule filter_sort_simreps:
    input:
        bed=rules.download_simreps.output,
        genome=rules.get_genome.output,
    output:
        ref_dir / "simreps.txt.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input.bed} | \
        cut -f 2,3,4 | \
        grep -Pv '^\S+_|^\S+EBV\s|^\S+M\s' | \
        bedtools sort -i stdin -g {input.genome} | \
        gzip -c > {output}
        """


rule filter_sort_ref:
    input:
        fa=rules.unzip_ref.output,
        genome=rules.get_genome.output,
    output:
        ref_dir / "ref_filtered.fa",
    conda:
        envs_path("utils.yml")
    shell:
        """
        samtools faidx {input.fa} $(cut -f1 {input.genome} | tr '\n' ' ') > \
        {output}
        """


hp_grid = {
    1: [4, 6, 12, 20],
    2: [10, 50, 200],
    3: [14, 50, 200],
    4: [19, 50, 200],
}


rule get_perfect_homopolymers:
    input:
        ref=rules.filter_sort_ref.output,
        bin=rules.build_repseq.output,
        genome=rules.get_genome.output,
    output:
        lc_inter_dir / "repeats_R{rep}_L{len}.bed",
    log:
        lc_log_dir / "low_complexity" / "repeats_R{rep}_L{len}.log",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        {input.bin} {wildcards.rep} {wildcards.len} {input.ref} \
        > {output} 2> {log}
        """


rule subtract_repeats:
    input:
        a=results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}.bed",
        b=results_dir / "low_complexity" / "repeats_R{rep}_L{lenB}.bed",
    output:
        lc_inter_dir / "repeats_R{rep}_L{lenA}to{lenB}.bed",
    shell:
        "subtractBed -a {input.a} -b {input.b} > {output}"


rule add_repeat_slop:
    input:
        bed=rules.subtract_repeats.output,
        genome=rules.get_genome.output,
    output:
        results_dir / "low_complexity" / "repeats_R{rep}_L{lenA}to{lenB}.bed",


rule get_imperfect_homopolymers:
    input:
        expand(rules.get_perfect_homopolymers.output, rep=1, len=4),
    output:
        lc_inter_dir / "repeats_imp_hp_{len}_{base}.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        grep 'unit={wildcards.base}' {input} | \
        mergeBed -i stdin -d 1 | \
        awk '$3-$2>{wildcards.len}' > \
        {output}
        """


rule get_final_imperfect_hp:
    input:
        perfect_hp=expand(
            rules.get_imperfect_homopolymers.output,
            base=["A", "C", "G", "T"],
            len=10,
        ),
        genome=rules.get_genome.output,
    output:
        lc_final_dir="GRCh38_SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        multiIntersectBed -i {input.perfect_hp} | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


simreps = {
    "lt51": {"lower": 0, "upper": 51},
    "51to200": {"lower": 50, "upper": 201},
    "gt200": {"lower": 200, "upper": 1000000},
}


rule get_simple_repeats:
    input:
        rules.filter_sort_rmsk.output,
    output:
        lc_inter_dir / "GRCh38_rmsk_Simple_repeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        zgrep Simple_repeat {input} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


rule get_low_complexity:
    input:
        rules.filter_sort_rmsk.output,
    output:
        lc_inter_dir / "GRCh38_rmsk_Low_complexity.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        zgrep Low_complexity {input} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


# rule get_low_complexity:
#     input:
#         rules.get_pre_low_complexity.output,
#     output:
#         lc_inter_dir / "GRCh38_rmsk_Low_complexity_{key}.bed.gz",
#     params:
#         lower=lambda wildcards: simreps[wildcards.simrep_key]["lower"],
#         upper=lambda wildcards: simreps[wildcards.simrep_key]["upper"],
#     shell:
#         """
#         gunzip -c {input} | \
#         awk '$3-$2>{params.lower} && $3-$2<{params.upper}' | \
#         bgzip -c > {output}
#         """


# TODO make sure this is sorted
rule get_pre_simreps:
    input:
        rules.filter_sort_simreps.output,
    output:
        lc_inter_dir / "simpleRepeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input} | \
        mergeBed -i stdin |  \
        gzip -c > {output}
        """


# TODO this might be easier/cleaner with pandas
rule get_simreps:
    input:
        rules.get_pre_simreps.output,
    output:
        lc_final_dir / "GRCh38_trf_simpleRepeat_{simrep_key}.bed.gz",
    params:
        lower=lambda wildcards: simreps[wildcards.simrep_key]["lower"],
        upper=lambda wildcards: simreps[wildcards.simrep_key]["upper"],
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        gunzip -c {input} | \
        awk '$3-$2>{params.lower} && $3-$2<{params.upper}' | \
        gzip -c > {output}
        """


rule get_pre_satellites:
    input:
        rmsk=rules.filter_sort_rmsk.output,
        idx=rules.index_ref.output,
    output:
        lc_inter_dir / "satellites.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        zgrep Satellite {input.rmsk} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


rule get_satellites:
    input:
        bed=rules.get_pre_satellites.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_satellites_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        slopBed -i {input.bed} -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


# TODO ensure this is sorted
rule get_final_imperfect_homopolymers:
    input:
        p6=expand(rules.get_perfect_homopolymers.output, rep=1, len=5),
        imp10=rules.get_final_imperfect_hp.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        grep -Ev '_|^chrEBV' {input.p6} | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        mergeBed -i stdin | \
        multiIntersectBed -i stdin {input.imp10} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


# TODO add di11 repeats here
# TODO add tri15 repeats here
# TODO add quad20 repeats here
rule get_all_tr_pre:
    input:
        hp=[
            rules.get_low_complexity.output,
            rules.get_simple_repeats.output,
            rules.get_pre_simreps.output,
            rules.get_satellites.output,
        ],
        genome=rules.get_genome.output,
    output:
        lc_inter_dir / "AllTandemRepeats_intermediate.bed",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        multiIntersectBed -i {input.hp} | \
        sed 's/^chr//' | \
        cut -f1-3 | \
        grep "^[0-9XY]" | \
        grep -v '_' | \
        sed 's/^/chr/' | \
        slopBed -i stdin -b 5 -g {input.genome} | \
        sed 's/^chr//' | \
        sed 's/^X/23/;s/^Y/24/' | \
        sort -k1,1n -k2,2n -k3,3n | \
        sed 's/^23/X/;s/^24/Y/;s/^/chr/' | \
        mergeBed -i stdin > {output}
        """


# NOTE these look funny because they have slop added
tr_bounds = {
    "lt51": {"lower": 0, "upper": 61},
    "51to200": {"lower": 60, "upper": 211},
    "201to10000": {"lower": 210, "upper": 10011},
    "gt10000": {"lower": 10010, "upper": 1000000},
}


rule get_final_tandem_repeats:
    input:
        tr=rules.get_all_tr_pre.output,
        hp=rules.get_final_imperfect_homopolymers.output,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeats_{tr_bound}_slop5.bed.gz",
    params:
        lower=lambda wildcards: tr_bounds[wildcards.tr_bound]["lower"],
        upper=lambda wildcards: tr_bounds[wildcards.tr_bound]["upper"],
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        awk '$3-$2>{params.lower} && $3-$2<{params.upper}' {input.tr} | \
        subtractBed -a stdin -b {input.hp} | \
        gzip -c > {output}
        """


# TODO this needs to be expanded for 4 different ranges
rule get_final_all_tr:
    input:
        beds=expand(rules.get_final_tandem_repeats.output, tr_bound=tr_bounds),
    output:
        lc_final_dir / "GRCh38_allTandemRepeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        multiIntersectBed -i {input.beds} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


rule get_final_all_tr_and_hp:
    input:
        rules.get_final_all_tr.output,
        rules.get_final_imperfect_homopolymers.output,
    output:
        lc_final_dir / "GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        """
        multiIntersectBed -i {input} | \
        mergeBed -i stdin | \
        gzip -c > {output}
        """


# TODO use complement for this
rule get_final_all_tr_compl:
    input:
        tr=rules.get_final_all_tr.output,
        genome=rules.get_genome.output,
    output:
        lc_final_dir / "GRCh38_notinallTandemRepeats.bed.gz",
    conda:
        envs_path("bedtools.yml")
    shell:
        "bedtools complement -i {input.tr} -g {input.genome} | gzip -c > {output}"


rule all_low_complexity:
    input:
        rules.get_final_all_tr_and_hp.output,
        rules.get_final_all_tr.output,
        rules.get_final_all_tr_compl.output,
