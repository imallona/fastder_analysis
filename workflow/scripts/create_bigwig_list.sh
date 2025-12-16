#!/usr/bin/env bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

INPUT="$WORKFLOW_DIR/data/unify_work/samples.tsv"
OUTPUT="$WORKFLOW_DIR/data/fastder_input/recount3_study-explorer_BigWig_list.csv"
FASTDER_DIR="$WORKFLOW_DIR/data/fastder_input"

[ ! -f "$INPUT" ] && echo "Error: $INPUT not found" && exit 1

awk -F'\t' -v OFS="," -v fastder_dir="$FASTDER_DIR" '
    NR == 1 {
        print "rail_id","external_id","study","BigWigURL"
        next
    }
    {
        rail_id   = $1
        sample_id = $2
        study_id  = $3

        # BigWig file is now in flat directory with pattern: sample_id.all.bw
        bigwig_file = fastder_dir "/" sample_id ".all.bw"

        # Check if file exists
        cmd = "test -f " bigwig_file " && echo yes || echo no"
        cmd | getline exists
        close(cmd)

        if (exists == "yes") {
            # Relative path from fastder_input directory
            rel_path = "./" sample_id ".all.bw"
            print rail_id, sample_id, study_id, rel_path
        }
    }
' "$INPUT" > "$OUTPUT"

echo "Created BigWig list CSV at: $OUTPUT"
