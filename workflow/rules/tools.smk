# Competing-tool runner rules (derfinder, megadepth_baseline).
#
# Included by the main Snakefile after the path constants and the parameter
# grid block are defined.


# derfinder: Bioconductor coverage-based ER detection. CPM-normalised
# per-base mean coverage thresholded at --cutoff, post-filtered by
# --min-length, with optional gap-bridging via --maxregiongap. Sweeps the
# (min_coverage, position_tolerance) plane shared with fastder. Wildcard
# param_id has the form mc<v>_pt<v>; we parse it to recover the values.
rule run_derfinder:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
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
            --chromosomes {params.chroms} > {log} 2>&1
        """


# megadepth_baseline: thresholded mean-CPM segmenter. Same coverage
# pipeline as derfinder, but no gap-bridging or SJ stitching: one
# transcript per maximal run of bases at or above --cutoff. Sweeps the
# min_coverage axis shared with fastder. Wildcard param_id is mc<v>.
rule run_megadepth_baseline:
    input:
        chr_prefix_done=op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"),
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
    conda:
        "../envs/megadepth_baseline.yaml"
    shell:
        """
        python3 {params.script} \
            --bigwig-dir {params.bigwig_dir} \
            --out-gtf {output.gtf} \
            --cutoff {params.cutoff} \
            --min-length {params.min_length} \
            --chromosomes {params.chroms} > {log} 2>&1
        """
