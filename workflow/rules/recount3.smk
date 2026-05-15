# recount3 ingestion rules.
#
# Active when monorail.backend is "recount3". Instead of running recount-pump
# on raw reads, these rules download the Monorail-processed outputs that
# recount3 already hosts: per-sample coverage BigWigs and the per-study
# junction matrix. The downloaded files are then reshaped per sample group
# into the lean MM/RR layout that extract_fastder_inputs feeds to fastder.
#
# This file is included by the main Snakefile after the path constants and
# the R3_* constants and _r3_shard helper are defined. It relies on op,
# DATA_DIR, LOG_DIR, BENCH_DIR, WORKFLOW_DIR, config, the R3_* constants and
# the _r3_shard helper.


# recount3 study metadata. Maps recount3 rail_id to SRA run accession, which
# subset_recount3_junctions.py needs to pick the right MM columns.
rule recount3_fetch_metadata:
    output:
        tsv=op.join(R3_DIR, "{study}.recount_project.tsv"),
    log:
        op.join(LOG_DIR, "recount3", "fetch_metadata_{study}.log"),
    benchmark:
        op.join(BENCH_DIR, "recount3", "fetch_metadata_{study}.tsv")
    params:
        url=lambda wc: (
            f"{R3_METADATA_URL}/human/data_sources/sra/metadata/"
            f"{_r3_shard(wc.study)}/{wc.study}/"
            f"sra.recount_project.{wc.study}.MD.gz"
        ),
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        curl -fSL --retry 3 '{params.url}' 2> {log} \
            | gunzip -c > {output.tsv} 2>> {log}
        """


# Per-study junction matrix: coordinates (RR), counts (MM) and the rail_id of
# each matrix column (ID). Downloaded once and shared by all groups.
rule recount3_fetch_junctions:
    output:
        rr=op.join(R3_DIR, "junctions", "{study}.ALL.RR"),
        mm=op.join(R3_DIR, "junctions", "{study}.ALL.MM"),
        idf=op.join(R3_DIR, "junctions", "{study}.ALL.ID"),
    log:
        op.join(LOG_DIR, "recount3", "fetch_junctions_{study}.log"),
    benchmark:
        op.join(BENCH_DIR, "recount3", "fetch_junctions_{study}.tsv")
    params:
        base=lambda wc: (
            f"{R3_BASE_URL}/human/data_sources/sra/junctions/"
            f"{_r3_shard(wc.study)}/{wc.study}/sra.junctions.{wc.study}.ALL"
        ),
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p $(dirname {output.rr})
        curl -fSL --retry 3 '{params.base}.RR.gz' 2> {log} \
            | gunzip -c > {output.rr} 2>> {log}
        curl -fSL --retry 3 '{params.base}.MM.gz' 2>> {log} \
            | gunzip -c > {output.mm} 2>> {log}
        curl -fSL --retry 3 '{params.base}.ID.gz' 2>> {log} \
            | gunzip -c > {output.idf} 2>> {log}
        """


# Per-sample coverage BigWig. recount3 names these "base_sums": the per-base
# sum of read coverage, which is the coverage track fastder consumes.
rule recount3_fetch_bigwig:
    output:
        bw=op.join(R3_DIR, "bw", "{sample}.all.bw"),
    log:
        op.join(LOG_DIR, "recount3", "fetch_bigwig_{sample}.log"),
    benchmark:
        op.join(BENCH_DIR, "recount3", "fetch_bigwig_{sample}.tsv")
    params:
        url=lambda wc: (
            f"{R3_BASE_URL}/human/data_sources/sra/base_sums/"
            f"{_r3_shard(R3_STUDY)}/{R3_STUDY}/{_r3_shard(wc.sample)}/"
            f"sra.base_sums.{R3_STUDY}_{wc.sample}.ALL.bw"
        ),
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bw})
        curl -fSL --retry 3 -o {output.bw} '{params.url}' > {log} 2>&1
        """


# Reshape the study junction matrix into one lean MM/RR pair per sample
# group, keeping only that group's columns and the analysed chromosomes.
# The {scenario} wildcard is a group name (see SCENARIOS in the Snakefile).
rule recount3_group_junctions:
    input:
        rr=op.join(R3_DIR, "junctions", f"{R3_STUDY}.ALL.RR"),
        mm=op.join(R3_DIR, "junctions", f"{R3_STUDY}.ALL.MM"),
        idf=op.join(R3_DIR, "junctions", f"{R3_STUDY}.ALL.ID"),
        metadata=op.join(R3_DIR, f"{R3_STUDY}.recount_project.tsv"),
    output:
        rr=op.join(R3_DIR, "{scenario}", "junctions.ALL.RR"),
        mm=op.join(R3_DIR, "{scenario}", "junctions.ALL.MM"),
        samples_tsv=op.join(R3_DIR, "{scenario}", "junctions.ALL.samples.tsv"),
    log:
        op.join(LOG_DIR, "recount3", "group_junctions_{scenario}.log"),
    benchmark:
        op.join(BENCH_DIR, "recount3", "group_junctions_{scenario}.tsv")
    params:
        out_prefix=lambda wc: op.join(R3_DIR, wc.scenario, "junctions.ALL"),
        study=R3_STUDY,
        sample_args=lambda wc: " ".join(
            f"--sample {s}" for s in R3_GROUPS[wc.scenario]
        ),
        chrom_args=lambda wc: " ".join(
            config.get("fastder", {}).get("chromosomes")
            or [f"chr{i}" for i in range(1, 23)] + ["chrX"]
        ),
        script=op.join(WORKFLOW_DIR, "scripts", "subset_recount3_junctions.py"),
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p $(dirname {output.rr})
        python3 {params.script} \
            --rr {input.rr} --mm {input.mm} --id {input.idf} \
            --metadata {input.metadata} --study {params.study} \
            {params.sample_args} \
            --chromosomes {params.chrom_args} \
            --out-prefix {params.out_prefix} > {log} 2>&1
        """
