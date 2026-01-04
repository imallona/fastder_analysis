#!/bin/bash
# Force all features to positive (+) strand in GFF3/GTF files
# Usage: bash force_strand.sh input.gff3 output.gff3

awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$7="+"; print}' "$1" > "$2"
