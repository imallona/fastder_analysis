#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_gffcompare_minimal.sh <reference.gtf> <query.gtf> <out_prefix>
# Example:
#   ./run_gffcompare_minimal.sh sample_01_label_splicing_variants.gtf \
#     FASTDER_RESULT_POS_TOL_20_MIN_COV_0.050000_COV_TOL_0.800000_MIN_LENGTH_10.gtf \
#     cmp

REF_IN="${1:?reference.gtf required}"
QRY_IN="${2:?query.gtf required}"
OUT_PREFIX="${3:-cmp}"

REF_FINAL="${OUT_PREFIX}.ref.processed.gtf"
QRY_FINAL="${OUT_PREFIX}.qry.processed.gtf"

process_gtf () {
  local in="$1" out="$2"
  awk 'BEGIN{FS=OFS="\t"}
    /^#/ {print; next}
    NF<9 {next}
    {
      # Re-pack attributes into column 9 (fixes embedded TABs in attr field)
      attrs=$9
      for(i=10;i<=NF;i++) attrs=attrs " " $i
      gsub(/[ \t]+/, " ", attrs)
      sub(/[ ]+$/, "", attrs)
      if (attrs !~ /;[ ]*$/) attrs=attrs ";"

      # Strand-agnostic
      $7="+"

      # Optional contig harmonization: strip leading "chr"
      sub(/^chr/,"",$1)

      print $1,$2,$3,$4,$5,$6,$7,$8,attrs
    }' "$in" > "$out"
}

process_gtf "$REF_IN" "$REF_FINAL"
process_gtf "$QRY_IN" "$QRY_FINAL"

gffcompare -r "$REF_FINAL" -o "$OUT_PREFIX" "$QRY_FINAL"

echo
echo "==== Base level block from ${OUT_PREFIX}.stats ===="
grep -n -A6 "Base level" "${OUT_PREFIX}.stats" || true
