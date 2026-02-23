#!/usr/bin/env bash
# Force the strand column (column 7) to "+" for all data lines in a GFF3/GTF file.
#
# Usage: bash force_strand.sh <input> <output>

set -euo pipefail

INPUT="$1"
OUTPUT="$2"

awk 'BEGIN{FS=OFS="\t"}
    /^#/ { print; next }
    { $7 = "+"; print }
' "$INPUT" > "$OUTPUT"
