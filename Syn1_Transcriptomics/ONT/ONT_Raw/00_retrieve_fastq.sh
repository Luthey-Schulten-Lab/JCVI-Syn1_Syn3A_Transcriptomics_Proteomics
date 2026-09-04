#!/usr/bin/env bash
set -eu
#############################################
# Retrieve Syn1 ONT direct-RNA runs from NCBI/SRA
#
# Two ONT direct-RNA runs on JCVI-Syn1.0 (SQK-RNA002), per the same
# submission as the Syn3A run (SRR36199724):
#   SRR36199726 : Syn1, rep1
#   SRR36199725 : Syn1, rep2  (annotated 3'->5' direction on NCBI)
#
# Steps: fasterq-dump -> rename to ONT.syn1.rep{1,2}.fastq -> gzip.
#
# NOTE on read orientation (5'->3' vs native 3'->5'):
#   ONT direct-RNA reads can be deposited in native 3'->5' order, in which
#   case they must be flipped with `seqkit seq -r` (reverse base order only,
#   NOT reverse-complement, since RNA is single-stranded) to recover the
#   5'->3' mRNA sense before mapping. Which reps need this is DECIDED
#   EMPIRICALLY (see 00b_orient_check.sh) rather than trusted from the NCBI
#   annotation, then applied here via the REVERSE_REPS list below.
#############################################

# rep index -> SRA accession
declare -A SRR=( [1]=SRR36199726 [2]=SRR36199725 )

# Rep numbers whose fastq must be reversed to 5'->3' (decided by 00b_orient_check.sh).
#   rep1 SRR36199726 : correct AS-IS  (86.6% mapping, 97.7% sense) -> NOT reversed
#   rep2 SRR36199725 : native 3'->5'  (reverse -> 79.1% mapping, 99.2% sense) -> reverse
REVERSE_REPS=(2)

# ---------- Tool checks ----------
for tool in fasterq-dump seqkit pigz; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "Error: $tool not found in the active environment. Activate the Omics conda env." >&2
    exit 1
  }
done

cd "$(dirname "$0")"

# ---------- 1) Download + rename ----------
for i in 1 2; do
  ID=${SRR[$i]}
  REP="ONT.syn1.rep${i}"
  if [ -f "${REP}.fastq" ] || [ -f "${REP}.fastq.gz" ]; then
    echo "[rep${i}] ${REP} already present, skipping download."
    continue
  fi
  if [ -f "${ID}.fastq" ]; then
    echo "[rep${i}] ${ID}.fastq present, renaming to ${REP}.fastq."
    mv "${ID}.fastq" "${REP}.fastq"
  else
    echo "[rep${i}] fasterq-dump ${ID} ..."
    fasterq-dump --split-files --progress -e 16 "$ID"
    mv "${ID}.fastq" "${REP}.fastq"
  fi
done

# ---------- 2) Optional reversal (empirically decided) ----------
for i in "${REVERSE_REPS[@]:-}"; do
  [ -z "$i" ] && continue
  REP="ONT.syn1.rep${i}"
  MARK="${REP}.reversed"
  if [ -f "${MARK}" ]; then
    echo "[rep${i}] reversal marker present, skipping."
    continue
  fi
  if [ -f "${REP}.fastq" ]; then
    echo "[rep${i}] reversing ${REP}.fastq (seqkit seq -r, base order only) ..."
    seqkit seq -r "${REP}.fastq" -o "${REP}.reversed.fastq"
    mv "${REP}.reversed.fastq" "${REP}.fastq"
  elif [ -f "${REP}.fastq.gz" ]; then
    echo "[rep${i}] reversing ${REP}.fastq.gz ..."
    seqkit seq -r "${REP}.fastq.gz" -o "${REP}.reversed.fastq.gz"
    mv "${REP}.reversed.fastq.gz" "${REP}.fastq.gz"
  else
    echo "ERROR: no fastq for rep${i} to reverse." >&2; exit 3
  fi
  touch "${MARK}"
done

# ---------- 3) Compress ----------
for i in 1 2; do
  REP="ONT.syn1.rep${i}"
  if [ -f "${REP}.fastq" ] && [ ! -f "${REP}.fastq.gz" ]; then
    echo "[rep${i}] compressing ${REP}.fastq ..."
    pigz -p 16 "${REP}.fastq"
  fi
done

echo "Done."
