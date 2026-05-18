# Reference download rule.
#
# Included by the main Snakefile after the path constants are defined.


# 1. Download reference
# Uses fastder.chromosomes from config to restrict which chromosomes are
# downloaded.  When chromosomes is omitted (or empty), all (1-22, X) are
# fetched, matching fastder's own --chr default. The output folder is keyed
# by chromosome scope (REF_DIR) and the outputs are explicit files, so
# different chromosome subsets do not clobber each other.
rule download_reference:
    output:
        gtf=REF_GTF,
        fastas=REF_FASTAS,
    benchmark:
        op.join(BENCH_DIR, "download_reference.tsv")
    log:
        op.join(LOG_DIR, "download_reference.log")
    params:
        outdir=REF_DIR,
        chr_args=" ".join(str(c) for c in FASTDER_CFG["chromosomes"])
            if FASTDER_CFG.get("chromosomes")
            else "",
    conda:
        "../envs/reference_download.yaml"
    shell:
        "bash scripts/download_reference_ensembl.sh '{params.outdir}' {params.chr_args} > {log} 2>&1"
