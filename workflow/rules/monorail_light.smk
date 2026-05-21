# =====================================================================
# monorail_light backend (BACKEND == "monorail_light")
# =====================================================================
# These rules replace pull_containers, download_monorail_refs, pump and unify
# with a chromosome-restricted STAR alignment plus a Python aggregator that
# emits the lean MM/RR files fastder consumes. They are gated on BACKEND so
# they never enter the DAG when the user asks for the heavy monorail backend.
#
# Included by the main Snakefile after the path constants are defined.


# Build a STAR genome index restricted to the chromosomes fastder will analyse
rule ml_star_index:
    input:
        gtf=REF_GTF,
        fastas=REF_FASTAS,
    output:
        idx=[op.join(LIGHT_STAR_IDX, f) for f in STAR_IDX_FILES],
    benchmark:
        op.join(BENCH_DIR, "ml_star_index.tsv")
    log:
        op.join(LOG_DIR, "ml_star_index.log"),
    params:
        idx_dir=LIGHT_STAR_IDX,
    threads: config["cores"]
    conda:
        "../envs/star.yaml"
    shell:
        """
        mkdir -p {params.idx_dir}
        # Normalize chromosome names to the UCSC-style 'chrN' the rest of the
        # pipeline expects. The Ensembl reference download stores FASTAs as
        # '21.fa' with header '>21' and a GTF whose first column is '21';
        # without this rewrite STAR's BAM, BigWig, SJ.out.tab and the lean
        # RR file would all use '21' while the config and gffcompare label
        # use 'chr21', leaving fastder with 0 junctions retained.
        tmpdir=$(mktemp -d)
        fastas=()
        for fa in {input.fastas}; do
            stripped=$(basename "$fa" .fa)
            stripped="${{stripped#chr}}"
            sed -E 's/^>(chr)?([^ \\t]*)/>chr\\2/' "$fa" \
                > "$tmpdir/chr${{stripped}}.fa"
            fastas+=("$tmpdir/chr${{stripped}}.fa")
        done
        # Prefix the GTF first column with 'chr' (skip header lines starting with '#').
        awk 'BEGIN{{FS=OFS="\\t"}}
             /^#/ {{print; next}}
             $1 !~ /^chr/ {{$1 = "chr" $1; print; next}}
             {{print}}' \
            {input.gtf} > "$tmpdir/annotation.gtf"
        # genomeSAindexNbases must be tuned down for small genomes (STAR manual).
        genome_size=$(cat "${{fastas[@]}}" | awk '!/^>/{{tot+=length($0)}} END{{print tot}}')
        sa=$(python3 -c "import math; print(min(14, int(math.log2($genome_size)/2 - 1)))")
        STAR --runMode genomeGenerate \
            --genomeDir {params.idx_dir} \
            --genomeFastaFiles "${{fastas[@]}}" \
            --sjdbGTFfile "$tmpdir/annotation.gtf" \
            --sjdbOverhang 100 \
            --genomeSAindexNbases "$sa" \
            --runThreadN {threads} > {log} 2>&1
        rm -rf "$tmpdir"
        """


# Align one sample's paired FASTQs against the chr-restricted STAR index.
# STAR emits the BAM unsorted; samtools sorts it. STAR's own BAM sort is
# capped by --limitBAMsortRAM, whose default (0) ties the sort buffer to the
# genome index size, which is too small for a chr-restricted index.
# Outputs a coordinate-sorted BAM and STAR's SJ.out.tab (used by ml_emit_mm_rr).
def ml_star_fastq_input(wc):
    """Paired FASTQs for ml_star_align. ASimulatoR input uses the per-scenario
    simulated reads; local input uses the paths from monorail.local_samples."""
    if PUMP_SOURCE == "local":
        sample_cfg = config["monorail"]["local_samples"][wc.sample]
        return {"fq1": sample_cfg["fq1"], "fq2": sample_cfg["fq2"]}
    return {
        "fq1": op.join(DATA_DIR, "asim", wc.sample, wc.scenario, "sample_01_1.fastq"),
        "fq2": op.join(DATA_DIR, "asim", wc.sample, wc.scenario, "sample_01_2.fastq"),
    }


rule ml_star_align:
    input:
        unpack(ml_star_fastq_input),
        idx=[op.join(LIGHT_STAR_IDX, f) for f in STAR_IDX_FILES],
    output:
        bam=op.join(LIGHT_DIR, "{scenario}", "{sample}", "Aligned.sortedByCoord.out.bam"),
        sj=op.join(LIGHT_DIR, "{scenario}", "{sample}", "SJ.out.tab"),
    benchmark:
        op.join(BENCH_DIR, "ml_star_align", "{sample}_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "ml_star_align", "{sample}_{scenario}.log"),
    params:
        outprefix=lambda wc: op.join(LIGHT_DIR, wc.scenario, wc.sample) + "/",
        idx_dir=LIGHT_STAR_IDX,
    threads: config["cores"]
    conda:
        "../envs/star.yaml"
    shell:
        """
        mkdir -p {params.outprefix}
        STAR --runMode alignReads \
            --genomeDir {params.idx_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --runThreadN {threads} \
            --outSAMtype BAM Unsorted \
            --outSAMstrandField intronMotif \
            --outFileNamePrefix {params.outprefix} > {log} 2>&1
        samtools sort -@ {threads} -o {output.bam} \
            {params.outprefix}Aligned.out.bam >> {log} 2>&1
        rm {params.outprefix}Aligned.out.bam
        samtools index {output.bam} >> {log} 2>&1
        """


# BigWig generation from STAR's BAM. Stranded/unstranded mirrors the heavy
# backend: stranded -> plus.bw + minus.bw, unstranded -> a single all.bw.
# -split / -F 256 match the megadepth coverage semantics used by recount-pump.
#
# TODO(deprecate): the stranded branch below (plus.bw + minus.bw) is a
# hack not used by the paper. Remove if the path stays unused, subject
# to the original contributor's agreement.
rule ml_bam_to_bigwig:
    input:
        bam=op.join(LIGHT_DIR, "{scenario}", "{sample}", "Aligned.sortedByCoord.out.bam"),
    output:
        bws=(
            [op.join(LIGHT_DIR, "{scenario}", "{sample}.plus.bw"),
             op.join(LIGHT_DIR, "{scenario}", "{sample}.minus.bw")]
            if STRANDED else
            [op.join(LIGHT_DIR, "{scenario}", "{sample}.all.bw")]
        ),
    benchmark:
        op.join(BENCH_DIR, "ml_bam_to_bigwig", "{sample}_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "ml_bam_to_bigwig", "{sample}_{scenario}.log"),
    params:
        chrom_sizes=lambda wc: op.join(LIGHT_DIR, wc.scenario, f"{wc.sample}.chrom.sizes"),
        stranded=STRANDED,
        outdir=lambda wc: op.join(LIGHT_DIR, wc.scenario),
    conda:
        "../envs/stranded_bigwig.yaml"
    shell:
        """
        bam="{input.bam}"
        samtools idxstats "$bam" | awk '$2 > 0 {{print $1"\\t"$2}}' > {params.chrom_sizes} 2>> {log}

        emit_bw() {{
            # $1 = strand_flag (+, -, or "" for unstranded); $2 = output bigwig path
            local strand_arg=""
            if [ -n "$1" ]; then strand_arg="-strand $1"; fi
            samtools view -u -F 256 "$bam" \
                | bedtools genomecov -ibam stdin -bga -split $strand_arg \
                > "$2.unsorted.bg" 2>> {log}
            sort -k1,1 -k2,2n "$2.unsorted.bg" > "$2.bg" 2>> {log}
            bedGraphToBigWig "$2.bg" {params.chrom_sizes} "$2" >> {log} 2>&1
            rm -f "$2.unsorted.bg" "$2.bg"
        }}

        if [ "{params.stranded}" = "True" ]; then
            emit_bw "+" "{params.outdir}/{wildcards.sample}.plus.bw"
            emit_bw "-" "{params.outdir}/{wildcards.sample}.minus.bw"
        else
            emit_bw "" "{params.outdir}/{wildcards.sample}.all.bw"
        fi
        """


# Aggregate STAR SJ.out.tab across all samples into a single MatrixMarket MM
# and a TSV RR. Junctions are filtered to fastder.chromosomes (so RR has only
# the rows fastder will analyse, not whole-genome) and the RR annotation
# columns are emitted as "." since fastder reads them but never accesses them
# downstream; see emit_lean_mm_rr.py for the rationale.
rule ml_emit_mm_rr:
    input:
        sj_files=expand(op.join(LIGHT_DIR, "{{scenario}}", "{sample}", "SJ.out.tab"),
                        sample=PUMP_SAMPLES),
    output:
        rr=op.join(LIGHT_DIR, "{scenario}", "junctions.ALL.RR"),
        mm=op.join(LIGHT_DIR, "{scenario}", "junctions.ALL.MM"),
        samples_tsv=op.join(LIGHT_DIR, "{scenario}", "samples.tsv"),
    benchmark:
        op.join(BENCH_DIR, "ml_emit_mm_rr_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "ml_emit_mm_rr_{scenario}.log"),
    params:
        out_prefix=lambda wc: op.join(LIGHT_DIR, wc.scenario, "junctions.ALL"),
        chroms=lambda wc: FASTDER_CFG.get("chromosomes") or [f"chr{i}" for i in range(1, 23)] + ["chrX"],
        samples=PUMP_SAMPLES,
        project=config["monorail"]["project_name"],
        emit_script=op.join(WORKFLOW_DIR, "scripts", "emit_lean_mm_rr.py"),
    run:
        # Build paired --sample / --sj args
        sample_args = []
        for s, sj in zip(params.samples, input.sj_files):
            sample_args.extend(["--sample", s, "--sj", sj])
        chr_args = list(params.chroms)
        shell(
            "python3 {params.emit_script} "
            + " ".join(f"'{a}'" for a in sample_args)
            + " --chromosomes " + " ".join(f"'{c}'" for c in chr_args)
            + " --out-prefix {params.out_prefix}"
            + " > {log} 2>&1"
        )
        # Emit a 3-column samples.tsv (rail_id, sample_id, study_id) so the
        # existing create_bigwig_list.sh can build the BigWig URL CSV from it.
        with open(output.samples_tsv, "w") as fh:
            fh.write("rail_id\tsample_id\tstudy_id\n")
            for rail_id, sample in enumerate(params.samples, start=1):
                fh.write(f"{rail_id}\t{sample}\t{params.project}\n")
