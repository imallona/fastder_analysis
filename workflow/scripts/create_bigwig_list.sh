#!/usr/bin/env bash
# Generate the BigWig-list CSV that fastder expects as metadata input.
# Maps rail_id -> external_id (sample name) so fastder can link BigWig files
# to their splice-junction matrix columns.
#
# Usage: bash create_bigwig_list.sh <samples_tsv> <output_csv> <mode>
#
# <samples_tsv>  Path to the unify samples.tsv (columns: rail_id, sample_id, study_id, ...)
# <output_csv>   Path to write the output CSV
# <mode>         "stranded" or "unstranded" (default: unstranded)

set -euo pipefail

SAMPLES_TSV="$1"
OUTPUT_CSV="$2"
MODE="${3:-unstranded}"

[ ! -f "$SAMPLES_TSV" ] && echo "Error: $SAMPLES_TSV not found" && exit 1

if [ "$MODE" = "stranded" ]; then
    awk -F'\t' -v OFS="," '
        NR == 1 {
            print "rail_id","external_id","study","BigWigURL"
            next
        }
        {
            rail_id   = $1
            sample_id = $2
            study_id  = $3
            print rail_id, sample_id, study_id, "./" sample_id ".plus.bw"
        }
    ' "$SAMPLES_TSV" > "$OUTPUT_CSV"
else
    awk -F'\t' -v OFS="," '
        NR == 1 {
            print "rail_id","external_id","study","BigWigURL"
            next
        }
        {
            rail_id   = $1
            sample_id = $2
            study_id  = $3
            print rail_id, sample_id, study_id, "./" sample_id ".all.bw"
        }
    ' "$SAMPLES_TSV" > "$OUTPUT_CSV"
fi
