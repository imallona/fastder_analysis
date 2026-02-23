#!/usr/bin/env bash
# Add "chr" prefix to the chromosome column of a GFF3/GTF file if not already present.
#
# Usage: bash add_chr_prefix_to_gff3.sh <input> <output>

set -euo pipefail

INPUT="$1"
OUTPUT="$2"

awk 'BEGIN{FS=OFS="\t"}
    /^#/ { print; next }
    $1 !~ /^chr/ { $1 = "chr" $1 }
    { print }
' "$INPUT" > "$OUTPUT"
