# ASimulatoR simulation rules.
#
# Included by the main Snakefile after the path constants are defined.


# 2. Run ASimulatoR
rule run_asimulator:
    input:
        gtf=REF_GTF,
        fastas=REF_FASTAS,
    output:
        # Explicit file outputs (rather than the parent directory) so snakemake
        # can chain make_scenario back to this rule via the FASTQ + GFF inputs.
        gff=op.join(ASIM_DIR, "{sample}", "splicing_variants.gff3"),
        fq1=simulated_reads(op.join(READS_DIR, "{sample}", "sample_01_1.fastq.gz")),
        fq2=simulated_reads(op.join(READS_DIR, "{sample}", "sample_01_2.fastq.gz")),
        meta=op.join(ASIM_DIR, "{sample}", "simulation_metadata.yaml"),
    benchmark:
        op.join(BENCH_DIR, "run_asimulator", "{sample}.tsv")
    log:
        op.join(LOG_DIR, "asimulator", "{sample}.log")
    params:
        # These read config["asimulator"] lazily, so a config without an
        # asimulator block still parses. This rule only runs when
        # pump_source is asimulator, where the block is present.
        outdir=lambda wc: op.join(READS_DIR, wc.sample),
        events=lambda wc: config["asimulator"]["samples"][wc.sample],
        seq_depth=lambda wc: config["asimulator"]["seq_depth"],
        multi_events_per_exon=lambda wc: config["asimulator"]["multi_events_per_exon"],
        strand_specific=lambda wc: config["asimulator"]["strand_specific"],
        probs_as_freq=lambda wc: config["asimulator"]["probs_as_freq"],
        seed=lambda wc: data_layout.sample_seed(config["seed"], ASIM_SAMPLES, wc.sample),
    threads: config["cores"]
    resources:
        mem_mb=32000,
        # Generous: a simulation killed near the end costs more.
        runtime=1440,
    container:
        ASIMULATOR_IMAGE
    script:
        "../scripts/runASimulatoR.R"


# 2b. Materialise the per-scenario asimulator outputs.
# template_and_variant keeps the original ASimulatoR output unchanged.
# variant_only drops the transcripts with template=TRUE in
# splicing_variants.gff3 from the truth, and their reads from the FASTQ, so
# the truth set used by gffcompare contains exactly the transcripts that
# produced the reads downstream rules will see.
# Truth and reads are separate rules: the reads are temp, and building them
# again must leave the truth file and its date alone. The input is ancient
# because run_asimulator writes it again, with the same seed, each time it
# rebuilds its reads.
rule make_scenario_truth:
    input:
        gff=ancient(op.join(ASIM_DIR, "{sample}", "splicing_variants.gff3")),
    output:
        gff=op.join(ASIM_DIR, "{sample}", "{scenario}", "splicing_variants.gff3"),
    log:
        op.join(LOG_DIR, "make_scenario_truth", "{sample}_{scenario}.log"),
    params:
        script=op.join(WORKFLOW_DIR, "scripts", "make_scenario.py"),
    resources:
        mem_mb=2000,
        runtime=10,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {params.script} --scenario {wildcards.scenario} \
            --gff-in {input.gff} --gff-out {output.gff} > {log} 2>&1
        """


rule make_scenario:
    input:
        gff=op.join(ASIM_DIR, "{sample}", "splicing_variants.gff3"),
        fq1=op.join(READS_DIR, "{sample}", "sample_01_1.fastq.gz"),
        fq2=op.join(READS_DIR, "{sample}", "sample_01_2.fastq.gz"),
    output:
        # Deleted once aligned. variant_only is a copy, hundreds of GB
        # over four depths. --notemp keeps them.
        fq1=temp(op.join(READS_DIR, "{sample}", "{scenario}", "sample_01_1.fastq.gz")),
        fq2=temp(op.join(READS_DIR, "{sample}", "{scenario}", "sample_01_2.fastq.gz")),
    benchmark:
        op.join(BENCH_DIR, "make_scenario", "{sample}_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "make_scenario", "{sample}_{scenario}.log"),
    params:
        script=op.join(WORKFLOW_DIR, "scripts", "make_scenario.py"),
    resources:
        mem_mb=4000,
        # About one minute per million reads.
        runtime=120,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {params.script} --scenario {wildcards.scenario} \
            --gff-in {input.gff} --fq1-in {input.fq1} --fq2-in {input.fq2} \
            --fq1-out {output.fq1} --fq2-out {output.fq2} \
            > {log} 2>&1
        """
