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
        gff=op.join(DATA_DIR, "asim", "{sample}", "splicing_variants.gff3"),
        fq1=op.join(DATA_DIR, "asim", "{sample}", "sample_01_1.fastq"),
        fq2=op.join(DATA_DIR, "asim", "{sample}", "sample_01_2.fastq"),
        meta=op.join(DATA_DIR, "asim", "{sample}", "simulation_metadata.yaml"),
    benchmark:
        op.join(BENCH_DIR, "run_asimulator", "{sample}.tsv")
    log:
        op.join(LOG_DIR, "asimulator", "{sample}.log")
    params:
        # These read config["asimulator"] lazily, so a config without an
        # asimulator block still parses. This rule only runs when
        # pump_source is asimulator, where the block is present.
        outdir=lambda wc: op.join(DATA_DIR, "asim", wc.sample),
        events=lambda wc: config["asimulator"]["samples"][wc.sample],
        seq_depth=lambda wc: config["asimulator"]["seq_depth"],
        multi_events_per_exon=lambda wc: config["asimulator"]["multi_events_per_exon"],
        strand_specific=lambda wc: config["asimulator"]["strand_specific"],
        probs_as_freq=lambda wc: config["asimulator"]["probs_as_freq"],
        seed=config["seed"],
        ncores=config["cores"],
    container:
        "docker://biomedbigdata/asimulator"
    script:
        "../scripts/runASimulatoR.R"


# 2b. Materialise the per-scenario asimulator outputs.
# template_and_variant symlinks the original ASimulatoR output unchanged.
# variant_only filters the FASTQ to drop reads whose source transcript carries
# template=TRUE in splicing_variants.gff3 and rewrites the GFF to keep only
# the alternative isoforms, so the truth set used by gffcompare contains
# exactly the transcripts that produced the reads downstream rules will see.
rule make_scenario:
    input:
        gff=op.join(DATA_DIR, "asim", "{sample}", "splicing_variants.gff3"),
        fq1=op.join(DATA_DIR, "asim", "{sample}", "sample_01_1.fastq"),
        fq2=op.join(DATA_DIR, "asim", "{sample}", "sample_01_2.fastq"),
    output:
        gff=op.join(DATA_DIR, "asim", "{sample}", "{scenario}", "splicing_variants.gff3"),
        fq1=op.join(DATA_DIR, "asim", "{sample}", "{scenario}", "sample_01_1.fastq"),
        fq2=op.join(DATA_DIR, "asim", "{sample}", "{scenario}", "sample_01_2.fastq"),
    benchmark:
        op.join(BENCH_DIR, "make_scenario", "{sample}_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "make_scenario", "{sample}_{scenario}.log"),
    params:
        script=op.join(WORKFLOW_DIR, "scripts", "make_scenario.py"),
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {params.script} --scenario {wildcards.scenario} \
            --gff-in {input.gff} --fq1-in {input.fq1} --fq2-in {input.fq2} \
            --gff-out {output.gff} --fq1-out {output.fq1} --fq2-out {output.fq2} \
            > {log} 2>&1
        """
