#!/bin/bash
# Add "chr" prefix to chromosome column in GFF3 files if not already present
# Usage: bash add_chr_prefix_to_gff3.sh input.gff3 output.gff3

input_file="$1"
output_file="$2"

awk 'BEGIN{FS=OFS="\t"}
     /^##/ {print; next}  # Print header lines as-is
     /^[^#]/ {
         if ($1 !~ /^chr/) {  # Only add chr prefix if not already present
             $1="chr"$1
         }
         print
     }
' "$input_file" > "$output_file"
