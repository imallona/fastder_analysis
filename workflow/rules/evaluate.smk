# Evaluation rules: gffcompare, fuzzy metrics, and result aggregation.
#
# Included by the main Snakefile after the path constants and the parameter
# grid block are defined.


# 12. Run gffcompare for each (tool, scenario, sample, param_id) combination
rule run_gffcompare:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
        gtf=op.join(DATA_DIR, "tools", "{tool}", "{scenario}", "{param_id}", "output.gtf"),
    output:
        stats=op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "gffcompare.stats"),
    benchmark:
        op.join(BENCH_DIR, "run_gffcompare", "{tool}_{scenario}_{sample}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "gffcompare", "{tool}_{scenario}_{sample}_{param_id}.log")
    params:
        label_gff=lambda wc: op.join(FASTDER_DIR, wc.scenario, f"{wc.sample}_label{LABEL_EXT}"),
        out_prefix=lambda wc: os.path.join(
            str(RESULTS_DIR), wc.tool, wc.scenario, wc.sample, wc.param_id, "gffcompare",
        ),
    conda:
        "../envs/gffcompare.yaml"
    shell:
        """
        mkdir -p $(dirname {params.out_prefix})
        gffcompare -r {params.label_gff} -o {params.out_prefix} {input.gtf} > {log} 2>&1
        """


# 13. Collect all gffcompare stats into a summary CSV. Parses every level
# (Base, Exon, Intron, Intron chain, Transcript, Locus), the matching
# counts, the missed and novel ratios, and the header summary numbers.
rule collect_results:
    input:
        lambda wc: expand_eval_paths(
            op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "gffcompare.stats"),
        ),
    output:
        op.join(RESULTS_DIR, "summary.csv"),
    benchmark:
        op.join(BENCH_DIR, "collect_results.tsv")
    params:
        parser=op.join(WORKFLOW_DIR, "scripts", "parse_gffcompare.py"),
    run:
        import csv
        import sys
        sys.path.insert(0, op.dirname(params.parser))
        import parse_gffcompare

        rows = []
        for stats_file in input:
            parts = Path(stats_file).parts
            pid = parts[-2]
            sample = parts[-3]
            scenario = parts[-4]
            tool = parts[-5]
            row = {"tool": tool, "scenario": scenario, "sample": sample, "param_id": pid}
            row.update(PARAM_DICT.get(pid, {}))
            row.update(parse_gffcompare.parse(stats_file))
            rows.append(row)

        # Stable column order: tool, scenario, sample, param_id, any grid
        # params, then all parsed metric keys in sorted order.
        meta_keys = {"tool", "scenario", "sample", "param_id"}
        param_cols = sorted({
            k for r in rows for k in r
            if k in PARAM_DICT.get(r["param_id"], {}) or k in (
                set(PARAM_DICT.get(r["param_id"], {}).keys()) - meta_keys
            )
        })
        metric_cols = sorted({k for r in rows for k in r if k not in meta_keys and k not in param_cols})
        fieldnames = ["tool", "scenario", "sample", "param_id"] + param_cols + metric_cols

        with open(output[0], "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in sorted(rows, key=lambda r: (r["tool"], r["scenario"], r["sample"], r["param_id"])):
                writer.writerow(row)


# 14. Per-transcript stats parsed from each fastder GTF.
rule collect_chain_stats:
    input:
        gtf_paths=expand(op.join(FASTDER_DIR, "{scenario}", "run_fastder_{param_id}.gtf_path"),
                         scenario=SCENARIOS, param_id=PARAM_IDS),
    output:
        op.join(RESULTS_DIR, "chain_stats.csv"),
    benchmark:
        op.join(BENCH_DIR, "collect_chain_stats.tsv")
    log:
        op.join(LOG_DIR, "collect_chain_stats.log"),
    params:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_chain_stats.py"),
        scenarios=SCENARIOS,
        param_ids=PARAM_IDS,
    run:
        # Each gtf_path file is at FASTDER_DIR/<scenario>/run_fastder_<pid>.gtf_path
        # so we recover (scenario, param_id) by matching the path components.
        args = []
        for gtf_path_file in input.gtf_paths:
            parts = Path(gtf_path_file).parts
            scenario = parts[-2]
            fname = parts[-1]
            pid = fname.removeprefix("run_fastder_").removesuffix(".gtf_path")
            args.extend(["--gtf-path-file", gtf_path_file,
                         "--param-id", pid, "--scenario", scenario])
        shell("python3 {params.script} " + " ".join(f"'{a}'" for a in args)
              + " --out {output} > {log} 2>&1")


# 15a. Fuzzy fastder-vs-reference comparison without gffcompare's strict
# exon-boundary requirement. Three per-(sample, param_id) CSVs:
# reciprocal-best Jaccard, fastder-exon-boundary distance to nearest ref
# boundary, and locus-recall curve at fixed thresholds.
rule eval_fuzzy_metrics:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
        gtf=op.join(DATA_DIR, "tools", "{tool}", "{scenario}", "{param_id}", "output.gtf"),
    output:
        jaccard=op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_jaccard.csv"),
        distances=op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_distances.csv"),
        recall=op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_locus_recall.csv"),
        strand=op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_strand.csv"),
    benchmark:
        op.join(BENCH_DIR, "eval_fuzzy_metrics", "{tool}_{scenario}_{sample}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "eval_fuzzy", "{tool}_{scenario}_{sample}_{param_id}.log")
    params:
        ref_gff=lambda wc: op.join(FASTDER_DIR, wc.scenario, f"{wc.sample}_label{LABEL_EXT}"),
        script=op.join(WORKFLOW_DIR, "scripts", "eval_fuzzy.py"),
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        python3 {params.script} \
            --ref-gff {params.ref_gff} \
            --fastder-gtf {input.gtf} \
            --sample {wildcards.sample} --param-id {wildcards.param_id} \
            --scenario {wildcards.scenario} \
            --out-jaccard {output.jaccard} \
            --out-distances {output.distances} \
            --out-recall {output.recall} \
            --out-strand {output.strand} > {log} 2>&1
        """


# 15b. Concatenate the per-(tool, scenario, sample, param_id) fuzzy CSVs.
rule collect_fuzzy_metrics:
    input:
        jaccard=lambda wc: expand_eval_paths(
            op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_jaccard.csv")),
        distances=lambda wc: expand_eval_paths(
            op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_distances.csv")),
        recall=lambda wc: expand_eval_paths(
            op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_locus_recall.csv")),
        strand=lambda wc: expand_eval_paths(
            op.join(RESULTS_DIR, "{tool}", "{scenario}", "{sample}", "{param_id}", "fuzzy_strand.csv")),
    output:
        jaccard=op.join(RESULTS_DIR, "fuzzy_jaccard.csv"),
        distances=op.join(RESULTS_DIR, "fuzzy_distances.csv"),
        recall=op.join(RESULTS_DIR, "fuzzy_locus_recall.csv"),
        strand=op.join(RESULTS_DIR, "fuzzy_strand.csv"),
    benchmark:
        op.join(BENCH_DIR, "collect_fuzzy_metrics.tsv")
    run:
        import csv as _csv
        # The per-(tool, scenario, sample, param_id) CSVs do not carry a
        # tool column; we recover it from the input file path:
        #   results/{tool}/{scenario}/{sample}/{param_id}/fuzzy_*.csv
        for src_files, outfile in [
            (input.jaccard, output.jaccard),
            (input.distances, output.distances),
            (input.recall, output.recall),
            (input.strand, output.strand),
        ]:
            header = None
            with open(outfile, "w", newline="") as out:
                writer = None
                for path in src_files:
                    parts = Path(path).parts
                    tool = parts[-5]
                    with open(path) as src:
                        reader = _csv.reader(src)
                        rows = list(reader)
                        if not rows:
                            continue
                        if header is None:
                            header = ["tool"] + rows[0]
                            writer = _csv.writer(out)
                            writer.writerow(header)
                        for row in rows[1:]:
                            writer.writerow([tool] + row)


# 16. Render the Rmarkdown summary report.
rule collect_truth_stats:
    input:
        # Only ASimulatoR input has simulated truth GFFs. recount3, SRA and
        # local input have nothing to aggregate; the rule then writes a
        # header-only CSV.
        gffs=lambda wc: (expand(
            op.join(ASIM_DIR, "{sample}", "{scenario}", "splicing_variants.gff3"),
            sample=PUMP_SAMPLES, scenario=SCENARIOS) if HAS_SIM_TRUTH else []),
    output:
        op.join(RESULTS_DIR, "truth_chain_stats.csv"),
    benchmark:
        op.join(BENCH_DIR, "collect_truth_stats.tsv")
    log:
        op.join(LOG_DIR, "collect_truth_stats.log"),
    params:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_truth_stats.py"),
    run:
        if not input.gffs:
            with open(output[0], "w") as fh:
                fh.write("scenario,sample,transcript_id,chrom,"
                         "n_exons,total_exon_length,strand\n")
        else:
            # input.gffs is in the same expand order: sample changes fastest,
            # scenario slowest.
            args = []
            for gff in input.gffs:
                parts = Path(gff).parts
                scenario = parts[-2]
                sample = parts[-3]
                args.extend(["--gff", gff, "--scenario", scenario, "--sample", sample])
            shell("python3 {params.script} " + " ".join(f"'{a}'" for a in args)
                  + " --out {output} > {log} 2>&1")
