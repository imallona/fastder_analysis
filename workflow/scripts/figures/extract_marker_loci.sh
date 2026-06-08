#!/usr/bin/env bash
# Troponin ER exons in the marker gene windows, read from the per-sub-group
# gffcompare GTFs. Usage: extract_marker_loci.sh <fastder_dir> <out.csv>
# fastder_dir holds one <tissue>_<n>/reference/mc1.0/gffcompare.annotated.gtf
# per sub-group. Windows are hg38 (TNNT2 chr1, TNNT3 chr11, TNNI3 chr19).
set -euo pipefail

fastder_dir=$1
out=$2

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
