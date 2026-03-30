#!/usr/bin/env bash
set -euo pipefail

# basic settings
RELEASE=115
ASSEMBLY=GRCh38

# Usage: download_reference_ensembl.sh <outdir> [chr1 chr2 ...]
# If no chromosomes are specified, all (1-22, X) are downloaded.
# Chromosomes can be given with or without "chr" prefix (e.g. chr21 or 21).
OUTDIR="${1:-data/reference/Homo_sapiens_${ASSEMBLY}_Ensembl${RELEASE}}"
shift || true

# Build chromosome list: strip "chr" prefix, default to all if none given
if [ $# -gt 0 ]; then
  chroms=()
  for arg in "$@"; do
    chroms+=("${arg#chr}")
  done
else
  chroms=($(seq 1 22) X)
fi

BASE="https://ftp.ensembl.org/pub/release-${RELEASE}"
GTF_URL="${BASE}/gtf/homo_sapiens/Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf.gz"
DNA_BASE="${BASE}/fasta/homo_sapiens/dna"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

ABS_OUTDIR="$(realpath .)"
echo "Output directory: ${ABS_OUTDIR}"
echo "Downloading GTF..."
wget -nv -c -N "${GTF_URL}"

echo "Downloading chromosomes: ${chroms[*]}..."
for c in "${chroms[@]}"; do
  wget -nv -c -N "${DNA_BASE}/Homo_sapiens.${ASSEMBLY}.dna.chromosome.${c}.fa.gz"
done

echo "Extracting .gz files..."
shopt -s nullglob
for f in *.gz; do
  gunzip -f "${f}"
done
shopt -u nullglob

GTF_FILE="Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf"
# Build awk regex from chromosome list
chrom_pattern=$(IFS='|'; echo "${chroms[*]}")
echo "Filtering GTF to chromosomes: ${chroms[*]}..."
tmpfile=$(mktemp)
awk -v pat="^(${chrom_pattern})$" 'BEGIN{FS=OFS="\t"}
     /^#/ {print; next}
     $1 ~ pat {print}' \
    "${GTF_FILE}" > "${tmpfile}"
mv "${tmpfile}" "${GTF_FILE}"

echo "Renaming FASTA files to <chrom>.fa..."
shopt -s nullglob
for f in Homo_sapiens.${ASSEMBLY}.dna.chromosome.*.fa; do
  newname=$(echo "${f}" | sed -E "s/^Homo_sapiens\.${ASSEMBLY}\.dna\.chromosome\.([^.]+)\.fa$/\1.fa/")
  mv -f "${f}" "${newname}"
done
shopt -u nullglob

echo "Done. Files in ${ABS_OUTDIR}"
