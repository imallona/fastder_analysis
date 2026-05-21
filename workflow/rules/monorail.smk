# Monorail (heavy) backend rules.
#
# Included by the main Snakefile after the path constants are defined.


def pump_input(wildcards):
    """Return file inputs for the pump rule.
    The asimulator source declares the per-sample ASimulatoR FASTQs as the
    dependency. The local source declares the user-provided FASTQ paths from
    monorail.local_samples. The SRA source has no file inputs because the
    reads are downloaded at run time.
    """
    if PUMP_SOURCE == "asimulator":
        return {
            "fq1": op.join(ASIM_DIR, wildcards.sample, "sample_01_1.fastq"),
            "fq2": op.join(ASIM_DIR, wildcards.sample, "sample_01_2.fastq"),
        }
    if PUMP_SOURCE == "local":
        sample_cfg = config["monorail"]["local_samples"][wildcards.sample]
        return {"fq1": sample_cfg["fq1"], "fq2": sample_cfg["fq2"]}
    return {}


# 3. Pull Singularity containers
rule pull_containers:
    output:
        pump=PUMP_SIF,
        unify=UNIFY_SIF,
    benchmark:
        op.join(BENCH_DIR, "pull_containers.tsv")
    log:
        op.join(LOG_DIR, "pull_containers.log")
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p {CONTAINERS_DIR}
        singularity pull --name {output.pump} docker://quay.io/broadsword/recount-pump:1.1.3 > {log} 2>&1
        singularity pull --name {output.unify} docker://quay.io/broadsword/recount-unify:1.1.1 >> {log} 2>&1
        """


# 4. Download Monorail reference indexes
rule download_monorail_refs:
    output:
        touch(op.join(MONORAIL_REF_DIR, "refs.DONE"))
    benchmark:
        op.join(BENCH_DIR, "download_monorail_refs.tsv")
    params:
        ref_version=config["monorail"]["ref_version"],
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p {MONORAIL_REF_DIR}
        cd {MONORAIL_REF_DIR}
        bash {MONORAIL_EXTERNAL}/get_human_ref_indexes.sh
        bash {MONORAIL_EXTERNAL}/get_unify_refs.sh {params.ref_version}
        """


# 5. Run recount-pump per sample
# Input source is controlled by monorail.pump_source in config.yaml:
#   asimulator: reads FASTQs from the run_asimulator output directories
#   sra:        downloads reads from SRA at runtime using the run/study accessions
#               defined in monorail.sra_samples
#   local:      reads FASTQs already on disk, at the paths given in
#               monorail.local_samples
rule pump:
    input:
        unpack(pump_input),
        container=PUMP_SIF,
        refs=op.join(MONORAIL_REF_DIR, "refs.DONE"),
    output:
        directory(op.join(DATA_DIR, "pump", "{sample}"))
    benchmark:
        op.join(BENCH_DIR, "pump", "{sample}.tsv")
    log:
        op.join(LOG_DIR, "pump", "{sample}.log")
    params:
        ref_version=config["monorail"]["ref_version"],
        ref_path=str(MONORAIL_REF_DIR),
        ncores=config["cores"],
        pump_script=op.join(MONORAIL_EXTERNAL, "singularity", "run_recount_pump.sh"),
        project_name=config["monorail"]["project_name"],
    run:
        os.makedirs(output[0], exist_ok=True)
        # asimulator and local both feed local FASTQ files to recount-pump;
        # they differ only in where the files are. sra downloads at run time.
        if PUMP_SOURCE in ("asimulator", "local"):
            if PUMP_SOURCE == "asimulator":
                fp1 = op.join(ASIM_DIR, wildcards.sample, "sample_01_1.fastq")
                fp2 = op.join(ASIM_DIR, wildcards.sample, "sample_01_2.fastq")
            else:
                sample_cfg = config["monorail"]["local_samples"][wildcards.sample]
                fp1 = sample_cfg["fq1"]
                fp2 = sample_cfg["fq2"]
            shell(
                "WORKING_DIR={output[0]}"
                " bash {params.pump_script}"
                " {input.container}"
                " {wildcards.sample}"
                " local"
                " {params.ref_version}"
                " {params.ncores}"
                " {params.ref_path}"
                f" {fp1}"
                f" {fp2}"
                " {params.project_name}"
                " > {log} 2>&1"
            )
        else:
            # TODO: test this part
            sample_cfg = config["monorail"]["sra_samples"][wildcards.sample]
            study_acc = sample_cfg["study_acc"]
            shell(
                "WORKING_DIR={output[0]}"
                " bash {params.pump_script}"
                " {input.container}"
                " {wildcards.sample}"
                f" {study_acc}"
                " {params.ref_version}"
                " {params.ncores}"
                " {params.ref_path}"
                " > {log} 2>&1"
            )


# 6. Recount-unify: aggregate all pump outputs
rule unify:
    input:
        pump_dirs=expand(op.join(DATA_DIR, "pump", "{sample}"), sample=PUMP_SAMPLES),
        container=UNIFY_SIF,
        refs=op.join(MONORAIL_REF_DIR, "refs.DONE"),
    output:
        directory(op.join(DATA_DIR, "unify")),
        op.join(DATA_DIR, "sample_metadata.tsv"),
    benchmark:
        op.join(BENCH_DIR, "unify.tsv")
    log:
        op.join(LOG_DIR, "unify.log")
    params:
        ref_version=config["monorail"]["ref_version"],
        ref_path=str(MONORAIL_REF_DIR),
        project="{short}:{pid}".format(
            short=config["monorail"]["project_name"],
            pid=config["monorail"]["project_id"],
        ),
        ncores=config["cores"],
        unify_script=op.join(MONORAIL_EXTERNAL, "singularity", "run_recount_unify.sh"),
        pump_parent=op.join(DATA_DIR, "pump"),
    run:
        os.makedirs(output[0], exist_ok=True)

        # Build the sample metadata (study_id<TAB>sample_id, with header).
        # study_id must match what pump encoded in its output filenames:
        #   asimulator and local: pump used project_name as the study
        #     -> filenames: sample!project_name!...
        #   sra: pump used the per-sample study_acc
        #     -> filenames: sample!study_acc!...
        if PUMP_SOURCE in ("asimulator", "local"):
            rows = [(config["monorail"]["project_name"], s) for s in PUMP_SAMPLES]
        else:
            rows = [(config["monorail"]["sra_samples"][s]["study_acc"], s) for s in PUMP_SAMPLES]
        with open(output[1], "w") as fh:
            fh.write("study_id\tsample_id\n")
            fh.writelines(f"{study}\t{sample}\n" for study, sample in rows)

        shell(
            "bash {params.unify_script}"
            " {input.container}"
            " {params.ref_version}"
            " {params.ref_path}"
            " {output[0]}"
            " {params.pump_parent}"
            " {output[1]}"
            " {params.ncores}"
            " {params.project}"
            " > {log} 2>&1"
        )


# 7a. Generate strand-specific BigWigs from pump BAM output.
# Uses bedtools genomecov to compute per-base coverage split by alignment strand,
# then converts to BigWig.  The -split flag ensures only aligned bases (not
# intron spans) are counted, matching megadepth's coverage semantics.
# -F 256 excludes secondary alignments (same as megadepth).
#
# TODO(deprecate): this rule only fires when fastder.stranded=true, a hack
# not used by the paper. Remove together with the STRANDED branches in
# fastder.smk if the path stays unused, subject to the original
# contributor's agreement.
rule generate_stranded_bigwigs:
    input:
        pump_dir=op.join(DATA_DIR, "pump", "{sample}"),
    output:
        plus_bw=op.join(DATA_DIR, "stranded_bigwigs", "{sample}.plus.bw"),
        minus_bw=op.join(DATA_DIR, "stranded_bigwigs", "{sample}.minus.bw"),
    benchmark:
        op.join(BENCH_DIR, "generate_stranded_bigwigs", "{sample}.tsv")
    log:
        op.join(LOG_DIR, "stranded_bigwigs", "{sample}.log")
    params:
        bam=lambda wc: op.join(DATA_DIR, "pump", wc.sample, "output", f"{wc.sample}_att0", "{}!{}!{}!local~sorted.bam".format(
                wc.sample,
                config["monorail"]["project_name"]
                    if PUMP_SOURCE in ("asimulator", "local")
                    else config["monorail"]["sra_samples"][wc.sample]["study_acc"],
                config["monorail"]["ref_version"],
            )),
        chrom_sizes=lambda wc: op.join(DATA_DIR, "stranded_bigwigs", f"{wc.sample}.chrom.sizes"),
    conda:
        "../envs/stranded_bigwig.yaml"
    shell:
        """
        bam="{params.bam}"

        # Index BAM if needed
        if [ ! -f "${{bam}}.bai" ]; then
            samtools index "$bam" >> {log} 2>&1
        fi

        # Get chromosome sizes from BAM header
        samtools idxstats "$bam" | awk '$2 > 0 {{print $1"\\t"$2}}' > {params.chrom_sizes} 2>> {log}

        # Plus strand: -F 256 excludes secondary alignments, -split respects CIGAR
        # -bga includes zero-coverage regions (required by fastder's Averager which
        # expects all samples to have identical per-chromosome vector lengths)
        samtools view -u -F 256 "$bam" \
            | bedtools genomecov -ibam stdin -bga -split -strand '+' \
            > {output.plus_bw}.unsorted.bg 2>> {log}
        sort -k1,1 -k2,2n {output.plus_bw}.unsorted.bg > {output.plus_bw}.bg 2>> {log}
        bedGraphToBigWig {output.plus_bw}.bg {params.chrom_sizes} {output.plus_bw} >> {log} 2>&1
        rm -f {output.plus_bw}.unsorted.bg {output.plus_bw}.bg

        # Minus strand
        samtools view -u -F 256 "$bam" \
            | bedtools genomecov -ibam stdin -bga -split -strand '-' \
            > {output.minus_bw}.unsorted.bg 2>> {log}
        sort -k1,1 -k2,2n {output.minus_bw}.unsorted.bg > {output.minus_bw}.bg 2>> {log}
        bedGraphToBigWig {output.minus_bw}.bg {params.chrom_sizes} {output.minus_bw} >> {log} 2>&1
        rm -f {output.minus_bw}.unsorted.bg {output.minus_bw}.bg
        """
