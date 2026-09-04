#!/usr/bin/env bash
set -euo pipefail

#############################################
# Map PacBio Iso-Seq HQ isoforms to genome
# using minimap2, then sort + index with samtools
#############################################

# ---------- Location of executables ----------
MINIMAP2_PATH=$(which minimap2 2>/dev/null)
if [ -z "$MINIMAP2_PATH" ]; then
    echo "Error: minimap2 not found in the active conda environment. Please activate the appropriate environment."
    exit 1
fi

# ---------- User inputs ----------
THREADS=16

# Required: reference genome FASTA and HQ isoform FASTA
Mother_dir="$(cd "$(dirname "$0")/../../.." && pwd)"   # project root
working_dir=${Mother_dir}/Syn1_Transcriptomics/PacBio/PacBio_Processing

Ref_file=${Mother_dir}/Genomes_Input/syn1_genome.fasta
# Isoform_file=${Isoform_dir}2-hq_transcripts.fasta
FLNC_file=${working_dir}/intermediate_files/all.FLNC.fastq.gz

# ---------- Checks ----------
# for exe in minimap2 samtools; do
#   command -v "${exe}" >/dev/null 2>&1 || { echo "ERROR: ${exe} not found in PATH" >&2; exit 127; }
# done
[[ -s "${Ref_file}" ]] || { echo "ERROR: Ref_file not found or empty: ${Ref_file}" >&2; exit 2; }
[[ -s "${FLNC_file}" ]] || { echo "ERROR: FLNC_file not found or empty: ${FLNC_file}" >&2; exit 2; }

# ---------- Outputs ----------
OUT_PREFIX=syn1.PacBio.FLNC
BAM_SORTED=${OUT_PREFIX}.sorted.bam
BAM_SORTED_CSI=${BAM_SORTED}.csi
LOG=${OUT_PREFIX}.minimap2.log

# ---------- Optional: FASTA index (useful for downstream tooling) ----------
# if [[ ! -f "${Ref_file}.fai" ]]; then
#   echo "[1/4] Indexing reference FASTA with samtools faidx..."
#   samtools faidx "${Ref_file}"
# else
#   echo "[1/4] Reference FASTA index exists: ${Ref_file}.fai"
# fi

# ---------- Mapping + sorting ----------
echo "[2/4] Mapping with minimap2 and sorting with samtools..."
# Notes:
# -ax map-hifi : PacBio HiFi/CCS genomic reads
# -uf        : full-length cDNA (recommended for Iso-Seq)
# --secondary=no : suppress secondary alignments (helps interpretation)
$MINIMAP2_PATH -t "${THREADS}" -ax map-hifi --secondary=no "${Ref_file}" "${FLNC_file}" 2> "${LOG}" \
  | samtools sort -@ "${THREADS}" -o "${BAM_SORTED}"

# ---------- Indexing ----------
echo "[3/4] Indexing BAM..."
# Use CSI by default (works for large references too)
samtools index -c "${BAM_SORTED}"  # creates .csi

# ---------- Quick sanity checks ----------
echo "[4/4] Sanity checks:"
samtools flagstat "${BAM_SORTED}" | sed -n '1,12p' >> "${LOG}"
echo
echo "Minimap2 log: ${LOG}"
echo "BAM index   : ${BAM_SORTED_CSI}"
echo "Done."
 
