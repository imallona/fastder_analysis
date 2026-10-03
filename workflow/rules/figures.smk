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

import yaml

FIG_SCRIPTS = op.join(WORKFLOW_DIR, "scripts", "figures")
FIG_DIR = config.get("figures_dir", op.join(WORKFLOW_DIR, "results", "figures"))
FIG_RESULTS = op.join(WORKFLOW_DIR, "results")
FIG_ENV = {
    "FASTDER_RESULTS_ROOT": FIG_RESULTS,
    "FASTDER_FIG_DIR": FIG_DIR,
    "FASTDER_BENCH_DIR": op.join(LOGS_ROOT, "benchmarks", "config_full_simulation"),
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
    resources:
        mem_mb=8000,
        runtime=60,
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
        reference=REF_GTF,
    output:
        op.join(FIG_DIR, "marker_loci.csv"),
    log:
        op.join(LOG_DIR, "figure_marker_loci.log"),
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        "bash {input.script} {FIG_RESULTS}/config_gtex_concordance/fastder "
        "{output} {input.reference} > {log} 2>&1"


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
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/figures.yaml"
    shell:
        "Rscript {input.script} {FIG_RESULTS}/config_gtex_concordance/fastder "
        "{input.reference} {output} > {log} 2>&1"


# TDP-43 inputs of figure 2, from the showcase and panel runs.
TDP43_SHOWCASE = "config_klim_2019_tdp43_recount3"
TDP43_PANEL = "config_klim_2019_tdp43_recount3_panel"


def tdp43_manifest(config_name):
    return op.join(FIG_RESULTS, config_name, "recount3_manifest.csv")


def tdp43_threshold(config_name):
    """The single min_coverage a TDP-43 config calls regions at."""
    path = op.join(WORKFLOW_DIR, "..", "config", config_name + ".yaml")
    with open(path) as handle:
        return yaml.safe_load(handle)["fastder"]["min_coverage"][0]


rule figure_tdp43_novel_exons:
    input:
        manifest=tdp43_manifest(TDP43_PANEL),
        script=op.join(FIG_SCRIPTS, "extract_novel_exons_tdp43.R"),
        reference=REF_GTF,
    output:
        op.join(FIG_DIR, "tdp43_novel_exons.csv"),
    log:
        op.join(LOG_DIR, "figure_tdp43_novel_exons.log"),
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/figures.yaml"
    shell:
        "Rscript {input.script} {input.manifest} {input.reference} {output} > {log} 2>&1"


rule figure_tdp43_jaccard:
    input:
        showcase=tdp43_manifest(TDP43_SHOWCASE),
        panel=tdp43_manifest(TDP43_PANEL),
        script=op.join(FIG_SCRIPTS, "extract_tdp43_jaccard.R"),
    output:
        op.join(FIG_DIR, "tdp43_jaccard.csv"),
    params:
        showcase=tdp43_threshold(TDP43_SHOWCASE),
        panel=tdp43_threshold(TDP43_PANEL),
    log:
        op.join(LOG_DIR, "figure_tdp43_jaccard.log"),
    resources:
        mem_mb=4000,
        runtime=60,
    conda:
        "../envs/figures.yaml"
    shell:
        "Rscript {input.script} {output} {params.showcase}={input.showcase} "
        "{params.panel}={input.panel} > {log} 2>&1"


rule figure_main_1:
    input:
        helpers=op.join(FIG_SCRIPTS, "helpers.R"),
        script=op.join(FIG_SCRIPTS, "figure_main_1.R"),
        schematics=op.join(FIG_DIR, "fig_sim_schematic.pdf"),
    output:
        op.join(FIG_DIR, "figure_main_1.pdf"),
    log:
        op.join(LOG_DIR, "figure_main_1.log"),
    resources:
        mem_mb=8000,
        runtime=60,
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
        tdp43_novel=op.join(FIG_DIR, "tdp43_novel_exons.csv"),
        tdp43_jaccard=op.join(FIG_DIR, "tdp43_jaccard.csv"),
    output:
        composite=op.join(FIG_DIR, "figure_main_2.pdf"),
        # Demoted from the composite, kept visible.
        supp=op.join(FIG_DIR, "supp_gtex_transcript_precision.pdf"),
    log:
        op.join(LOG_DIR, "figure_main_2.log"),
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/figures.yaml"
    shell:
        "{_fig_exports} Rscript {input.script} {output.composite} > {log} 2>&1"


# Tidy tables for the three new panels: what is plotted, on disk.
rule collect_ablation_table:
    input:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_param_sweeps.py"),
    output:
        csv=op.join(FIG_DIR, "ablation.csv"),
    log:
        op.join(LOG_DIR, "collect_ablation_table.log"),
    params:
        results_root=FIG_RESULTS,
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {input.script} --axis no_stitch \
            --results-root {params.results_root} \
            --out {output.csv} > {log} 2>&1
        """


rule collect_min_junction_reads_table:
    input:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_param_sweeps.py"),
    output:
        csv=op.join(FIG_DIR, "min_junction_reads.csv"),
    log:
        op.join(LOG_DIR, "collect_min_junction_reads_table.log"),
    params:
        results_root=FIG_RESULTS,
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {input.script} --axis min_junction_reads \
            --results-root {params.results_root} \
            --config-prefix config_min_junction_reads_sweep \
            --out {output.csv} > {log} 2>&1
        """


rule collect_annotation_table:
    input:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_param_sweeps.py"),
    output:
        csv=op.join(FIG_DIR, "annotation.csv"),
    log:
        op.join(LOG_DIR, "collect_annotation_table.log"),
    params:
        results_root=FIG_RESULTS,
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {input.script} --axis annotated_index \
            --results-root {params.results_root} \
            --out {output.csv} > {log} 2>&1
        """


# The sweep runs under whichever config declares fastder.scaling_cores.
rule collect_scaling_table:
    input:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_scaling.py"),
    output:
        csv=op.join(FIG_DIR, "scaling.csv"),
    log:
        op.join(LOG_DIR, "collect_scaling_table.log"),
    params:
        bench_dir=op.join(LOGS_ROOT, "benchmarks",
                          config.get("scaling_bench_config", "config_full_simulation")),
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {input.script} --bench-dir {params.bench_dir} \
            --out {output.csv} > {log} 2>&1
        """


rule figure_supp_revision:
    input:
        helpers=op.join(FIG_SCRIPTS, "helpers.R"),
        script=op.join(FIG_SCRIPTS, "figure_supp_revision.R"),
        ablation=op.join(FIG_DIR, "ablation.csv"),
        junction_filter=op.join(FIG_DIR, "min_junction_reads.csv"),
        scaling=op.join(FIG_DIR, "scaling.csv"),
    output:
        ablation=op.join(FIG_DIR, "supp_ablation.pdf"),
        junction_filter=op.join(FIG_DIR, "supp_min_junction_reads.pdf"),
        scaling=op.join(FIG_DIR, "supp_scaling.pdf"),
    log:
        op.join(LOG_DIR, "figure_supp_revision.log"),
    params:
        out_dir=FIG_DIR,
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/figures.yaml"
    shell:
        "{_fig_exports} Rscript {input.script} {params.out_dir} > {log} 2>&1"


# Capability table, replacing the two zero-bar panels. Reads no results.
rule capability_table:
    input:
        script=op.join(FIG_SCRIPTS, "make_capability_table.py"),
    output:
        csv=op.join(FIG_DIR, "tool_capabilities.csv"),
        tex=op.join(FIG_DIR, "tool_capabilities.tex"),
    log:
        op.join(LOG_DIR, "capability_table.log"),
    params:
        out_dir=FIG_DIR,
    resources:
        mem_mb=1000,
        runtime=10,
    conda:
        "../envs/base.yaml"
    shell:
        "python3 {input.script} {params.out_dir} > {log} 2>&1"


# Every value the text quotes, from the configs behind the paper. A config
# listed here whose results are missing fails the rule.
REPORTED_SIMULATIONS = ["config_full_simulation", "config_full_simulation_rep2",
                        "config_full_simulation_rep3", "config_full_simulation_5M",
                        "config_full_simulation_30M", "config_full_simulation_40M"]
REPORTED_RUNTIMES = ["config_gtex_comparison", "config_gtex_concordance",
                     "config_klim_2019_tdp43_recount3",
                     "config_klim_2019_tdp43_recount3_panel"]


rule collect_reported_numbers:
    input:
        script=op.join(WORKFLOW_DIR, "scripts", "collect_reported_numbers.py"),
    output:
        csv=op.join(FIG_DIR, "reported_numbers.csv"),
        tex=op.join(FIG_DIR, "reported_numbers.tex"),
    log:
        op.join(LOG_DIR, "collect_reported_numbers.log"),
    params:
        results_root=FIG_RESULTS,
        bench_root=op.join(LOGS_ROOT, "benchmarks"),
        simulations=" ".join(REPORTED_SIMULATIONS),
        runtimes=" ".join(REPORTED_RUNTIMES),
    resources:
        mem_mb=4000,
        runtime=60,
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 {input.script} \
            --results-root {params.results_root} --bench-root {params.bench_root} \
            --simulation {params.simulations} \
            --comparison config_gtex_comparison \
            --runtime {params.runtimes} \
            --threshold-range config_klim_2019_tdp43_recount3_ladder \
            --junction-filter config_min_junction_reads_sweep \
            --unannotated config_unannotated_alignment \
            --out-csv {output.csv} --out-tex {output.tex} > {log} 2>&1
        """


rule manuscript_figures:
    input:
        op.join(FIG_DIR, "figure_main_1.pdf"),
        op.join(FIG_DIR, "figure_main_2.pdf"),
        op.join(FIG_DIR, "supp_gtex_transcript_precision.pdf"),
        op.join(FIG_DIR, "tool_capabilities.tex"),
        op.join(FIG_DIR, "supp_ablation.pdf"),
        op.join(FIG_DIR, "supp_min_junction_reads.pdf"),
        op.join(FIG_DIR, "supp_scaling.pdf"),
        op.join(FIG_DIR, "annotation.csv"),
        op.join(FIG_DIR, "reported_numbers.tex"),
