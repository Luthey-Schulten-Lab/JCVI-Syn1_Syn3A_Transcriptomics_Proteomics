#!/usr/bin/env bash
set -eu
#############################################
# Map Syn1 ONT direct-RNA reads (SQK-RNA002) to the syn1 genome with minimap2,
# sort + index with samtools, one BAM per rep, then a merged BAM.
#
# Two runs, both stored in 5'->3' sense orientation by 00_retrieve_fastq.sh
# (rep1 SRR36199726 kept as-is; rep2 SRR36199725 reversed; verified by
# 00b_orient_check.sh -> sense-dominant coverage on the syn1 genome).
#
# minimap2 flags identical to the Syn3A ONT pipeline:
#   -ax map-ont : long-read, NON-spliced preset (Mycoplasma is intron-less;
#                 the splice preset emits spurious long-N "isoforms").
#   -p 0.99     : secondary alignments must score >=99% of the primary.
#   --MD        : emit MD tag for downstream tools.
# (No --secondary=no, to keep secondary alignments for rRNA-locus inspection.)
#############################################

THREADS=16
Mother_dir="$(cd "$(dirname "$0")/../../.." && pwd)"      # project root
RAW_dir="${Mother_dir}/Syn1_Transcriptomics/ONT/ONT_Raw"
Ref_file="${Mother_dir}/Genomes_Input/syn1_genome.fasta"

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found (activate Omics env)" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found (activate Omics env)" >&2; exit 1; }
[[ -s "${Ref_file}" ]] || { echo "ERROR: Ref not found: ${Ref_file}" >&2; exit 2; }

cd "$(dirname "$0")"

REP_BAMS=()
for i in 1 2; do
  FQ="${RAW_dir}/ONT.syn1.rep${i}.fastq.gz"
  [[ -s "${FQ}" ]] || { echo "ERROR: FQ not found: ${FQ}" >&2; exit 2; }
  OUT="syn1.ONT.rep${i}"
  BAM="${OUT}.sorted.bam"
  LOG="${OUT}.minimap2.log"
  echo "[map] rep${i}: ${FQ} -> ${BAM}"
  minimap2 -t "${THREADS}" -ax map-ont -p 0.99 --MD "${Ref_file}" "${FQ}" 2> "${LOG}" \
    | samtools sort -@ "${THREADS}" -o "${BAM}"
  samtools index -c "${BAM}"
  samtools flagstat "${BAM}" | sed -n '1,12p' >> "${LOG}"
  REP_BAMS+=("${BAM}")
done

# ---------- Merged BAM (both runs, for the genome-browser depth tracks) ----------
MERGED="syn1.ONT.merged.sorted.bam"
echo "[merge] ${REP_BAMS[*]} -> ${MERGED}"
samtools merge -@ "${THREADS}" -f "${MERGED}" "${REP_BAMS[@]}"
samtools index -c "${MERGED}"

echo "Done."
