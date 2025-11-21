#!/usr/bin/env bash
set -euo pipefail

# basic settings
RELEASE=115
ASSEMBLY=GRCh38

# root output directory fixed to ../../reference
OUTROOT=../../reference
OUTDIR_NAME="Homo_sapiens_${ASSEMBLY}_Ensembl${RELEASE}_chr1_22_X"
OUTDIR="${OUTROOT}/${OUTDIR_NAME}"

BASE="https://ftp.ensembl.org/pub/release-${RELEASE}"
GTF_URL="${BASE}/gtf/homo_sapiens/Homo_sapiens.${ASSEMBLY}.${RELEASE}.chr.gtf.gz"
DNA_BASE="${BASE}/fasta/homo_sapiens/dna"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

ABS_OUTDIR="$(realpath "${OUTDIR}")"
echo "Output directory: ${ABS_OUTDIR}"
echo "Downloading GTF..."
wget -nv -c -N "${GTF_URL}"

echo "Downloading chromosomes 1-22 and X..."
chroms=($(seq 1 22) X)
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
echo "Filtering GTF to chromosomes 1-22 and X only..."
tmpfile=$(mktemp)
awk 'BEGIN{FS=OFS="\t"}
     /^#/ {print; next}
     $1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X)$/ {print}' \
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
