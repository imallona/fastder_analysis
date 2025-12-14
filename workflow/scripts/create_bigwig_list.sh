#!/usr/bin/env bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

INPUT="$WORKFLOW_DIR/data/unify_work/samples.tsv"
OUTPUT="$WORKFLOW_DIR/data/fastder_input/recount3_study-explorer_BigWig_list.csv"
SAMPLES_DIR="$WORKFLOW_DIR/data/fastder_input/samples"

[ ! -f "$INPUT" ] && echo "Error: $INPUT not found" && exit 1

awk -F'\t' -v OFS="," -v samples_dir="$SAMPLES_DIR" '
    NR == 1 {
        print "rail_id","external_id","study","BigWigURL"
        next
    }
    {
        rail_id   = $1
        sample_id = $2
        study_id  = $3
        
        # Find BigWig file
        cmd = "find " samples_dir "/" sample_id " -name \"*!*!*!local.all.bw\" -type f 2>/dev/null | head -n 1"
        cmd | getline bigwig_file
        close(cmd)
        
        if (bigwig_file) {
            # Extract just the filename
            split(bigwig_file, parts, "/")
            filename = parts[length(parts)]
            rel_path = "./samples/" sample_id "/" filename
            print rail_id, sample_id, study_id, rel_path
        }
    }
' "$INPUT" > "$OUTPUT"

echo "Created BigWig list CSV at: $OUTPUT"