# Report rendering and archiving rules.
#
# Included by the main Snakefile after the path constants and the parameter
# grid block are defined.


# Snapshot the rendered reports into archive/<timestamp>/ so previous runs
# of the same config are not lost when summary.html / benchmarks.html get
# re-rendered. Fires whenever the input HTMLs are newer than the marker.
rule archive_reports:
    input:
        summary=op.join(RESULTS_DIR, "summary.html"),
        benchmarks=op.join(RESULTS_DIR, "benchmarks.html"),
    output:
        marker=op.join(RESULTS_DIR, "archive.DONE"),
    log:
        op.join(LOG_DIR, "archive_reports.log"),
    params:
        archive_root=op.join(RESULTS_DIR, "archive"),
    conda:
        "../envs/base.yaml"
    shell:
        """
        ts=$(date +%Y%m%d_%H%M%S)
        dest={params.archive_root}/$ts
        mkdir -p $dest
        cp -f {input.summary} $dest/summary.html
        cp -f {input.benchmarks} $dest/benchmarks.html
        echo "archived to $dest" | tee {output.marker} > {log}
        """


# One scenario and one parameter combination per tool feed the browser-style
# track view in the report. The first scenario and each tool's first param_id
# are representative; those GTFs are already built by the time summary.csv is.
_REPORT_SCENARIO = SCENARIOS[0]
_REPORT_TRUTH = op.join(FASTDER_DIR, _REPORT_SCENARIO,
                        SAMPLES_BY_SCENARIO[_REPORT_SCENARIO][0] + "_label" + LABEL_EXT)


rule render_summary_report:
    input:
        summary=op.join(RESULTS_DIR, "summary.csv"),
        chain_stats=op.join(RESULTS_DIR, "chain_stats.csv"),
        truth_chain_stats=op.join(RESULTS_DIR, "truth_chain_stats.csv"),
        jaccard=op.join(RESULTS_DIR, "fuzzy_jaccard.csv"),
        distances=op.join(RESULTS_DIR, "fuzzy_distances.csv"),
        recall=op.join(RESULTS_DIR, "fuzzy_locus_recall.csv"),
        strand=op.join(RESULTS_DIR, "fuzzy_strand.csv"),
        fastder_gtf=op.join(DATA_DIR, "tools", "fastder", _REPORT_SCENARIO,
                            PARAM_IDS_BY_TOOL["fastder"][0], "output.gtf"),
        derfinder_gtf=op.join(DATA_DIR, "tools", "derfinder", _REPORT_SCENARIO,
                              PARAM_IDS_BY_TOOL["derfinder"][0], "output.gtf"),
        megadepth_gtf=op.join(DATA_DIR, "tools", "megadepth_baseline", _REPORT_SCENARIO,
                              PARAM_IDS_BY_TOOL["megadepth_baseline"][0], "output.gtf"),
        # The truth label GFF is a side effect of extract_fastder_inputs;
        # depend on match_chr_prefix.DONE so it is in place before rendering.
        chr_prefix=op.join(FASTDER_DIR, _REPORT_SCENARIO, "match_chr_prefix.DONE"),
        # summary.Rmd grades the tools against a simulated truth set. Runs
        # without one (recount3, sra, local) get the descriptive
        # summary_custom.Rmd, which reports what the tools called and how
        # the calls agree, with no precision or recall.
        rmd=op.join(WORKFLOW_DIR, "reports",
                    "summary.Rmd" if HAS_SIM_TRUTH else "summary_custom.Rmd"),
    output:
        op.join(RESULTS_DIR, "summary.html"),
    log:
        op.join(LOG_DIR, "render_summary_report.log"),
    params:
        truth_gff=_REPORT_TRUTH,
        track_scenario=_REPORT_SCENARIO,
    conda:
        "../envs/rmarkdown.yaml"
    shell:
        """
        Rscript -e "rmarkdown::render(
            input = '{input.rmd}',
            output_file = '$(realpath -m {output})',
            params = list(summary_csv = '$(realpath {input.summary})',
                          chain_stats_csv = '$(realpath {input.chain_stats})',
                          truth_chain_stats_csv = '$(realpath {input.truth_chain_stats})',
                          fuzzy_jaccard_csv = '$(realpath {input.jaccard})',
                          fuzzy_distances_csv = '$(realpath {input.distances})',
                          fuzzy_recall_csv = '$(realpath {input.recall})',
                          fuzzy_strand_csv = '$(realpath {input.strand})',
                          fastder_gtf = '$(realpath {input.fastder_gtf})',
                          derfinder_gtf = '$(realpath {input.derfinder_gtf})',
                          megadepth_gtf = '$(realpath {input.megadepth_gtf})',
                          truth_gff = '$(realpath -m {params.truth_gff})',
                          track_scenario = '{params.track_scenario}'),
            quiet = TRUE)" > {log} 2>&1
        """


# 15. Render the Rmarkdown benchmarks report from logs/benchmarks/.
rule render_benchmarks_report:
    input:
        summary=op.join(RESULTS_DIR, "summary.csv"),
        rmd=op.join(WORKFLOW_DIR, "reports", "benchmarks.Rmd"),
    output:
        op.join(RESULTS_DIR, "benchmarks.html"),
    log:
        op.join(LOG_DIR, "render_benchmarks_report.log"),
    params:
        bench_dir=BENCH_DIR,
    conda:
        "../envs/rmarkdown.yaml"
    shell:
        """
        Rscript -e "rmarkdown::render(
            input = '{input.rmd}',
            output_file = '$(realpath -m {output})',
            params = list(bench_dir = '$(realpath {params.bench_dir})'),
            quiet = TRUE)" > {log} 2>&1
        """


# Build the manifest the recount3 report reads: one row per sample, with the
# group it belongs to, its coverage BigWig, and its group's called-region GTF
# for each tool (fastder, derfinder, megadepth_baseline).
rule recount3_report_manifest:
    input:
        fastder_gtfs=expand(
            op.join(DATA_DIR, "tools", "fastder", "{scenario}",
                    PARAM_IDS_BY_TOOL["fastder"][0], "output.gtf"),
            scenario=SCENARIOS),
        derfinder_gtfs=expand(
            op.join(DATA_DIR, "tools", "derfinder", "{scenario}",
                    PARAM_IDS_BY_TOOL["derfinder"][0], "output.gtf"),
            scenario=SCENARIOS),
        megadepth_gtfs=expand(
            op.join(DATA_DIR, "tools", "megadepth_baseline", "{scenario}",
                    PARAM_IDS_BY_TOOL["megadepth_baseline"][0], "output.gtf"),
            scenario=SCENARIOS),
        bigwigs=expand(op.join(R3_DIR, "bw", "{sample}.all.bw"),
                       sample=R3_ALL_SAMPLES),
    output:
        manifest=op.join(RESULTS_DIR, "recount3_manifest.csv"),
    params:
        groups=R3_GROUPS,
    run:
        import csv as _csv

        def by_group(gtfs):
            return {Path(g).parts[-3]: os.path.abspath(g) for g in gtfs}

        fastder_by_group = by_group(input.fastder_gtfs)
        derfinder_by_group = by_group(input.derfinder_gtfs)
        megadepth_by_group = by_group(input.megadepth_gtfs)
        bw_by_sample = {}
        for bw in input.bigwigs:
            name = Path(bw).name
            if name.endswith(".all.bw"):
                bw_by_sample[name[: -len(".all.bw")]] = os.path.abspath(bw)
        with open(output.manifest, "w", newline="") as fh:
            writer = _csv.writer(fh)
            writer.writerow(["group", "sample", "bigwig",
                             "fastder_gtf", "derfinder_gtf", "megadepth_gtf"])
            for group, samples in params.groups.items():
                for sample in samples:
                    writer.writerow([
                        group, sample, bw_by_sample[sample],
                        fastder_by_group.get(group, ""),
                        derfinder_by_group.get(group, ""),
                        megadepth_by_group.get(group, ""),
                    ])


# Render the recount3 comparison report. Only requested by rule all when the
# backend is recount3. Shows the coverage view at the TDP-43 cryptic exons and
# the sample dissimilarity summary, knockdown group versus control group.
rule render_recount3_report:
    input:
        manifest=op.join(RESULTS_DIR, "recount3_manifest.csv"),
        summary=op.join(RESULTS_DIR, "summary.csv"),
        reference_gtf=(REF_GTF if BACKEND == "recount3" else []),
        rmd=op.join(WORKFLOW_DIR, "reports", "recount3.Rmd"),
    output:
        op.join(RESULTS_DIR, "recount3.html"),
    log:
        op.join(LOG_DIR, "render_recount3_report.log"),
    params:
        study=R3_STUDY,
        reference_gtf=REF_ANNOTATION,
    conda:
        "../envs/rmarkdown.yaml"
    shell:
        """
        Rscript -e "rmarkdown::render(
            input = '{input.rmd}',
            output_file = '$(realpath -m {output})',
            params = list(manifest_csv = '$(realpath {input.manifest})',
                          summary_csv = '$(realpath {input.summary})',
                          reference_gtf = '{params.reference_gtf}',
                          study = '{params.study}'),
            quiet = TRUE)" > {log} 2>&1
        """
