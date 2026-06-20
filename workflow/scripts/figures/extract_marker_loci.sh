#!/usr/bin/env bash
# Troponin ER exons in the marker gene windows, read from the per-sub-group
# gffcompare GTFs, plus the Ensembl reference gene models for the same genes.
# Usage: extract_marker_loci.sh <fastder_dir> <out.csv> [reference_gtf]
# fastder_dir holds one <tissue>_<n>/reference/mc1.0/gffcompare.annotated.gtf
# per sub-group. Windows are hg38 (TNNT2 chr1, TNNT3 chr11, TNNI3 chr19). When
# reference_gtf is given, the reference exons of the three genes are appended
# with tissue "Ensembl" so the panel can draw them as a reference row.
set -euo pipefail

fastder_dir=$1
out=$2
reference_gtf=${3:-}

echo "gene,tissue,start,end" > "$out"
for gtf in "$fastder_dir"/*/reference/mc1.0/gffcompare.annotated.gtf; do
    rel=${gtf#"$fastder_dir"/}
    sub=${rel%%/*}
    tissue=${sub%_*}
    awk -v t="$tissue" 'BEGIN { OFS = "," } $3 == "exon" {
        if ($1 == "chr1"  && $4 <= 201382000 && $5 >= 201355000) print "TNNT2", t, $4, $5;
        else if ($1 == "chr11" && $4 <= 1970000  && $5 >= 1925000)  print "TNNT3", t, $4, $5;
        else if ($1 == "chr19" && $4 <= 55162000 && $5 >= 55148000) print "TNNI3", t, $4, $5;
    }' "$gtf" >> "$out"
done

# Ensembl reference exons for the three genes. The Ensembl GTF names the genes
# in the gene_name attribute and drops the "chr" prefix, so match by gene_name.
if [ -n "$reference_gtf" ] && [ -f "$reference_gtf" ]; then
    awk 'BEGIN { OFS = "," } $3 == "exon" {
        if      ($0 ~ /gene_name "TNNT2"/) print "TNNT2", "Ensembl", $4, $5;
        else if ($0 ~ /gene_name "TNNT3"/) print "TNNT3", "Ensembl", $4, $5;
        else if ($0 ~ /gene_name "TNNI3"/) print "TNNI3", "Ensembl", $4, $5;
    }' "$reference_gtf" >> "$out"
fi
