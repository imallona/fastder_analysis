# Manuscript figure assembly.
#
# These figures aggregate results across configs (the depth sweep,
# config_gtex_comparison, config_gtex_concordance, config_klim_2019_tdp43_recount3),
# so they run after the per-config pipelines, like the depth-sweep meta-report.
# The panels are built by the scripts under scripts/figures/; the environment
# variables point them at the results tree and the figure directory.
#
# Two panel families need inputs from outside the per-config CSV results: the
# genome-browser track and the similarity and concordance heatmaps are the
# rendered report figures, expected as PNGs in FASTDER_FIG_DIR; the troponin
# marker loci read the per-sub-group GTFs of config_gtex_concordance.

FIG_SCRIPTS = op.join(WORKFLOW_DIR, "scripts", "figures")
FIG_DIR = config.get("figures_dir", op.join(WORKFLOW_DIR, "results", "figures"))
FIG_RESULTS = op.join(WORKFLOW_DIR, "results")
FIG_ENV = {
    "FASTDER_RESULTS_ROOT": FIG_RESULTS,
    "FASTDER_FIG_DIR": FIG_DIR,
    "FASTDER_BENCH_DIR": op.join(WORKFLOW_DIR, "logs", "benchmarks", "config_full_simulation"),
}
_fig_exports = " ".join(f"{k}={v}" for k, v in FIG_ENV.items())


# Schematics drawn from scratch (matplotlib): the pipeline, the simulation
# design, and the two worked-example sample designs.
rule figure_schematics:
    input:
        sim=op.join(FIG_SCRIPTS, "make_sim_schematic.py"),
        samples=op.join(FIG_SCRIPTS, "make_sample_schematics.py"),
    output:
        op.join(FIG_DIR, "fig_sim_schematic.pdf"),
        op.join(FIG_DIR, "fig_tdp43_scheme.pdf"),
        op.join(FIG_DIR, "fig_gtex_scheme.pdf"),
    log:
        op.join(LOG_DIR, "figure_schematics.log"),
    conda:
        "../envs/figures.yaml"
    shell:
        "mkdir -p {FIG_DIR} && cd {FIG_DIR} && "
        "python {input.sim} > {log} 2>&1 && python {input.samples} >> {log} 2>&1"


# Troponin ER exons in the marker windows, pulled from the per-sub-group GTFs.
rule figure_marker_loci:
    input:
        gtfs=op.join(FIG_RESULTS, "config_gtex_concordance", "archive.DONE"),
        script=op.join(FIG_SCRIPTS, "extract_marker_loci.sh"),
    output:
        op.join(FIG_DIR, "marker_loci.csv"),
    log:
        op.join(LOG_DIR, "figure_marker_loci.log"),
    shell:
        "bash {input.script} {FIG_RESULTS}/config_gtex_concordance/fastder {output} > {log} 2>&1"


# Novel ER exons per tissue, from the per-sub-group GTFs against REF_GTF, which
# is genome-wide only under a genome-wide config (config_gtex_concordance).
rule figure_novel_exons:
    input:
        gtfs=op.join(FIG_RESULTS, "config_gtex_concordance", "archive.DONE"),
        script=op.join(FIG_SCRIPTS, "extract_novel_exons.R"),
        reference=REF_GTF,
    output:
        op.join(FIG_DIR, "novel_exons.csv"),
    log:
        op.join(LOG_DIR, "figure_novel_exons.log"),
    conda:
        "../envs/figures.yaml"
    shell:
        "Rscript {input.script} {FIG_RESULTS}/config_gtex_concordance/fastder "
        "{input.reference} {output} > {log} 2>&1"


rule figure_main_1:
    input:
        helpers=op.join(FIG_SCRIPTS, "helpers.R"),
        script=op.join(FIG_SCRIPTS, "figure_main_1.R"),
        schematics=op.join(FIG_DIR, "fig_sim_schematic.pdf"),
    output:
        op.join(FIG_DIR, "figure_main_1.pdf"),
    log:
        op.join(LOG_DIR, "figure_main_1.log"),
    conda:
        "../envs/figures.yaml"
    shell:
        "{_fig_exports} Rscript {input.script} {output} > {log} 2>&1"


rule figure_main_2:
    input:
        helpers=op.join(FIG_SCRIPTS, "helpers.R"),
        script=op.join(FIG_SCRIPTS, "figure_main_2.R"),
        schematics=op.join(FIG_DIR, "fig_tdp43_scheme.pdf"),
        markers=op.join(FIG_DIR, "marker_loci.csv"),
        novel=op.join(FIG_DIR, "novel_exons.csv"),
    output:
        op.join(FIG_DIR, "figure_main_2.pdf"),
    log:
        op.join(LOG_DIR, "figure_main_2.log"),
    conda:
        "../envs/figures.yaml"
    shell:
        "{_fig_exports} Rscript {input.script} {output} > {log} 2>&1"


rule manuscript_figures:
    input:
        op.join(FIG_DIR, "figure_main_1.pdf"),
        op.join(FIG_DIR, "figure_main_2.pdf"),
