# Adapter and poly tail removed/trimmed by Guppy confirmed by quick check.
# We can proceed with alignment and quantification

#!/usr/bin/env bash
#############################################
# Map ONT Direct RNA reads (SQK-RNA002) to genome
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

# Required: reference genome FASTA and ONT direct RNA FASTQ
Mother_dir="$(cd "$(dirname "$0")/../../.." && pwd)"   # project root
working_dir=${Mother_dir}/Syn3A_Transcriptomics/ONT/ONT_Raw
Ref_file=${Mother_dir}/Genomes_Input/syn3A_genome.fasta
# Replicate index passed as first argument (e.g., 1, 2, 3); defaults to 1
FQ_file=${working_dir}/ONT.syn3A.rep1.fastq.gz

# ---------- Checks ----------
[[ -s "${Ref_file}" ]] || { echo "ERROR: Ref_file not found or empty: ${Ref_file}" >&2; exit 2; }
[[ -s "${FQ_file}"  ]] || { echo "ERROR: FQ_file not found or empty: ${FQ_file}"  >&2; exit 2; }

# ---------- Outputs ----------
OUT_PREFIX=syn3A.ONT.rep1
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
# Notes (ONT Direct RNA, SQK-RNA002 on JCVI-Syn3A — a prokaryote):
# -ax map-ont   : long-read, non-spliced ONT preset. JCVI-Syn3A is intron-less,
#                 so the splice-aware preset (-ax splice) is biologically wrong:
#                 it allowed minimap2 to emit alignments with very long `N`
#                 (skipped-reference) CIGAR operations — single reads ended up
#                 mapped across 10s-100s kb of the genome, producing "isoforms"
#                 wider than any transcript could be.
# -p 0.99       : secondary alignments must score >=99% of primary (strict)
# --MD          : output MD tag (mismatch encoding; useful for downstream tools)
# Notes on flags dropped vs the previous splice version:
#   -uf is splice-mode-only (controls GT-AG canonical-splice direction);
#       strandedness of direct RNA reads is preserved in BAM flag bits.
#   -k14 is the splice preset's recommended k-mer size; map-ont's default 15 is correct.
# (No --secondary=no here, to keep secondary alignments for downstream rRNA-locus
#  inspection as before.)

$MINIMAP2_PATH -t "${THREADS}" -ax map-ont -p 0.99 --MD \
  "${Ref_file}" "${FQ_file}" 2> "${LOG}" \
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