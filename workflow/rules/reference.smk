# Reference download rule.
#
# Included by the main Snakefile after the path constants are defined.


# 1. Download reference
# Uses fastder.chromosomes from config to restrict which chromosomes are
# downloaded.  When chromosomes is omitted (or empty), all (1-22, X) are
# fetched, matching fastder's own --chr default.
rule download_reference:
    output:
        directory("data/reference/Homo_sapiens_GRCh38_Ensembl115")
    benchmark:
        op.join(BENCH_DIR, "download_reference.tsv")
    log:
        op.join(LOG_DIR, "download_reference.log")
    params:
        chr_args=" ".join(str(c) for c in FASTDER_CFG["chromosomes"])
            if FASTDER_CFG.get("chromosomes")
            else "",
    conda:
        "../envs/reference_download.yaml"
    shell:
        "bash scripts/download_reference_ensembl.sh '{output}' {params.chr_args} > {log} 2>&1"
