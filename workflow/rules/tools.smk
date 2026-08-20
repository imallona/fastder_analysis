# Competing-tool runner rules (derfinder, megadepth_baseline, grohmm).
#
# Included by the main Snakefile after the path constants and the parameter
# grid block are defined.


# Per-sample library sizes, shared by every tool.
#
# The library size is the whole-file sum of length times value, the same
# quantity fastder accumulates internally and recount3 calls the AUC. Its own
# rule keeps it out of the benchmarked ones. fastder takes the total from the
# BigWig summary header, which reads no intervals, whereas derfinder would have
# to import every chromosome to reach the same number. Charged inside its own
# rule, that would inflate derfinder's wall time by the ratio of the genome to
# the analysed subset.
rule compute_library_sizes:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
    output:
        tsv=op.join(FASTDER_DIR, "{scenario}", "library_sizes.tsv"),
    log:
        op.join(LOG_DIR, "compute_library_sizes", "{scenario}.log"),
    params:
        bigwig_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
        script=op.join(WORKFLOW_DIR, "scripts", "compute_library_sizes.py"),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
    conda:
        "../envs/megadepth_baseline.yaml"
    shell:
        """
        python3 {params.script} \
            --bigwig-dir {params.bigwig_dir} \
            --out {output.tsv} > {log} 2>&1
        """


# derfinder: Bioconductor coverage-based ER detection. CPM-normalised
# per-base mean coverage thresholded at --cutoff, post-filtered by
# --min-length, with optional gap-bridging via --maxregiongap. Sweeps the
# (min_coverage, position_tolerance) plane shared with fastder. Wildcard
# param_id has the form mc<v>_pt<v>; we parse it to recover the values.
rule run_derfinder:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
        library_sizes=op.join(FASTDER_DIR, "{scenario}", "library_sizes.tsv"),
    output:
        gtf=op.join(DATA_DIR, "tools", "derfinder", "{scenario}", "{param_id}", "output.gtf"),
    benchmark:
        op.join(BENCH_DIR, "run_derfinder", "{scenario}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "run_derfinder", "{scenario}_{param_id}.log"),
    params:
        bigwig_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
        chroms=lambda wc: FASTDER_CFG.get("chromosomes") or [f"chr{i}" for i in range(1, 23)] + ["chrX"],
        cutoff=lambda wc: _parse_mc_pt(wc.param_id)[0],
        maxregiongap=lambda wc: int(_parse_mc_pt(wc.param_id)[1]),
        min_length=lambda wc: (FASTDER_CFG.get("min_length") or [10])[0],
        script=op.join(WORKFLOW_DIR, "scripts", "run_derfinder.R"),
    threads: 1
    resources:
        # derfinder imports the coverage of a whole chromosome into R. Peak
        # resident set was 3.8 GB in the first submission's benchmarks.
        mem_mb=8000,
        runtime=120,
    conda:
        "../envs/derfinder.yaml"
    shell:
        """
        Rscript --vanilla {params.script} \
            --bigwig-dir {params.bigwig_dir} \
            --out-gtf {output.gtf} \
            --cutoff {params.cutoff} \
            --min-length {params.min_length} \
            --maxregiongap {params.maxregiongap} \
            --library-sizes {input.library_sizes} \
            --chromosomes {params.chroms} > {log} 2>&1
        """


# megadepth_baseline: thresholded mean-CPM segmenter. Same coverage
# pipeline as derfinder, but no gap-bridging or SJ stitching: one
# transcript per maximal run of bases at or above --cutoff. Sweeps the
# min_coverage axis shared with fastder. Wildcard param_id is mc<v>.
rule run_megadepth_baseline:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
        library_sizes=op.join(FASTDER_DIR, "{scenario}", "library_sizes.tsv"),
    output:
        gtf=op.join(DATA_DIR, "tools", "megadepth_baseline", "{scenario}", "{param_id}", "output.gtf"),
    benchmark:
        op.join(BENCH_DIR, "run_megadepth_baseline", "{scenario}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "run_megadepth_baseline", "{scenario}_{param_id}.log"),
    params:
        bigwig_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
        chroms=lambda wc: FASTDER_CFG.get("chromosomes") or [f"chr{i}" for i in range(1, 23)] + ["chrX"],
        cutoff=lambda wc: _parse_mc_pt(wc.param_id)[0],
        min_length=lambda wc: (FASTDER_CFG.get("min_length") or [10])[0],
        script=op.join(WORKFLOW_DIR, "scripts", "run_megadepth_baseline.py"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    conda:
        "../envs/megadepth_baseline.yaml"
    shell:
        """
        python3 {params.script} \
            --bigwig-dir {params.bigwig_dir} \
            --out-gtf {output.gtf} \
            --cutoff {params.cutoff} \
            --min-length {params.min_length} \
            --library-sizes {input.library_sizes} \
            --chromosomes {params.chroms} > {log} 2>&1
        """


# grohmm: HMM-based segmenter (Chae et al. 2015), originally for GRO-seq,
# adopted on RNA-seq coverage by TAR-scRNA-seq (Wang et al. 2021). Run here
# on the same unstranded coverage bigWigs as the other two tools. The runner
# averages per-50bp-window mean coverage across samples via kent's
# bigWigAverageOverBed, CPM-normalises with the shared library_size formula,
# integer-scales the result and feeds it to groHMM::detectTranscripts via the
# Fp (plus-strand) interface. Output is unstranded, matching derfinder and
# the megadepth baseline on recount3 input. Wildcard param_id is
# lp<v>_uts<v>; we parse it to recover (LtProbB, UTS).
rule run_grohmm:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
        library_sizes=op.join(FASTDER_DIR, "{scenario}", "library_sizes.tsv"),
    output:
        gtf=op.join(DATA_DIR, "tools", "grohmm", "{scenario}", "{param_id}", "output.gtf"),
    benchmark:
        op.join(BENCH_DIR, "run_grohmm", "{scenario}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "run_grohmm", "{scenario}_{param_id}.log"),
    params:
        bigwig_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
        chroms=lambda wc: FASTDER_CFG.get("chromosomes") or [f"chr{i}" for i in range(1, 23)] + ["chrX"],
        ltprobb=lambda wc: _parse_lp_uts(wc.param_id)[0],
        uts=lambda wc: _parse_lp_uts(wc.param_id)[1],
        min_length=lambda wc: (FASTDER_CFG.get("min_length") or [10])[0],
        window_size=lambda wc: (config.get("grohmm") or {}).get("window_size", 50),
        count_scale=lambda wc: (config.get("grohmm") or {}).get("count_scale", 100),
        script=op.join(WORKFLOW_DIR, "scripts", "run_grohmm.R"),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120,
    conda:
        "../envs/grohmm.yaml"
    shell:
        """
        Rscript --vanilla {params.script} \
            --bigwig-dir {params.bigwig_dir} \
            --out-gtf {output.gtf} \
            --ltprobb={params.ltprobb} \
            --uts {params.uts} \
            --window-size {params.window_size} \
            --min-length {params.min_length} \
            --count-scale {params.count_scale} \
            --library-sizes {input.library_sizes} \
            --chromosomes {params.chroms} > {log} 2>&1
        """
