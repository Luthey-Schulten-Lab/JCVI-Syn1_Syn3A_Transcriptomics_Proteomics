#!/usr/bin/env bash
# 1) Define the list of accessions of ONT datasets
# SRR36199724: Syn3A, 3' to 5' direction on NCBI
# SRR36199726: Syn1, rep1, 
# SRR36199725: Syn1, rep2, 3' to 5' direction on NCBI

SRR_ACCESSIONS=("SRR36199724")

# ---------- Tool checks ----------
FASTERQ_PATH=$(which fasterq-dump 2>/dev/null)
if [ -z "$FASTERQ_PATH" ]; then
  echo "Error: fasterq-dump not found in the active conda environment. Please activate the appropriate environment."
  exit 1
fi
 
SEQKIT_PATH=$(which seqkit 2>/dev/null)
if [ -z "$SEQKIT_PATH" ]; then
  echo "Error: seqkit not found in the active conda environment. Please activate the appropriate environment."
  exit 1
fi
 
PIGZ_PATH=$(which pigz 2>/dev/null)
if [ -z "$PIGZ_PATH" ]; then
  echo "Error: pigz not found in the active conda environment. Please activate the appropriate environment."
  exit 1
fi
 
# ---------- Skip checks ----------
# Check if steps 2 and 3 can be skipped
SRR_PRESENT=true
for ID in "${SRR_ACCESSIONS[@]}"; do
    [ -f "${ID}.fastq" ] || { SRR_PRESENT=false; break; }
done
 
REP_PRESENT=true
for i in 1; do
    [ -f "ONT.syn3A.rep${i}.fastq" ] || [ -f "ONT.syn3A.rep${i}.fastq.gz" ] || { REP_PRESENT=false; break; }
done
 
if $REP_PRESENT; then
  echo "Rep files already present, skipping steps 2 and 3."
elif $SRR_PRESENT; then
  echo "SRR fastq files already present, skipping download (step 2). Running rename (step 3)."
  mv SRR36199724.fastq ONT.syn3A.rep1.fastq
else
  # 2) Loop through each accession
  for ID in "${SRR_ACCESSIONS[@]}"; do
    echo "Processing $ID..."
    # Run fasterq-dump (using -e 16 for multi-threading speed)
    "$FASTERQ_PATH" --split-files --progress -e 16 "$ID"
  done
  # 3) Rename the files
  mv SRR36199724.fastq ONT.syn3A.rep1.fastq
fi
 
# 4) Reverse reads from 3'->5' (native pore order) to 5'->3' (mRNA sense)
# ONT direct RNA reads from SRA are in native 3'->5' orientation. We use
# `seqkit seq -r` (reverse only, NOT reverse-complement) to flip base order,
# since RNA is single-stranded — only the read direction needs to change.
# The reversed file overwrites the rep file so downstream scripts can use it directly.
REV_PRESENT=true
for i in 1; do
    # Marker file confirms reversal step has already been run
    [ -f "ONT.syn3A.rep${i}.reversed" ] || { REV_PRESENT=false; break; }
done
 
if $REV_PRESENT; then
  echo "Reversal already done (marker file present), skipping step 4."
else
  for i in 1; do
    if [ -f "ONT.syn3A.rep${i}.fastq" ]; then
      echo "Reversing ONT.syn3A.rep${i}.fastq with seqkit seq -r..."
      "$SEQKIT_PATH" seq -r "ONT.syn3A.rep${i}.fastq" -o "ONT.syn3A.rep${i}.reversed.fastq"
      mv "ONT.syn3A.rep${i}.reversed.fastq" "ONT.syn3A.rep${i}.fastq"
    elif [ -f "ONT.syn3A.rep${i}.fastq.gz" ]; then
      echo "Reversing ONT.syn3A.rep${i}.fastq.gz with seqkit seq -r..."
      "$SEQKIT_PATH" seq -r "ONT.syn3A.rep${i}.fastq.gz" -o "ONT.syn3A.rep${i}.reversed.fastq.gz"
      mv "ONT.syn3A.rep${i}.reversed.fastq.gz" "ONT.syn3A.rep${i}.fastq.gz"
    else
      echo "ERROR: No fastq found for rep${i} to reverse." >&2
      exit 3
    fi
  done
fi
 
# 5) Compress them into fastq.gz
# -p 16: use 16 threads
GZ_PRESENT=true
for i in 1; do
    [ -f "ONT.syn3A.rep${i}.fastq.gz" ] || { GZ_PRESENT=false; break; }
done
 
if $GZ_PRESENT; then
  echo "Compressed rep files already present, skipping step 5."
else
  "$PIGZ_PATH" -p 16 ONT.syn3A.rep1.fastq
fi
 
echo "Done."
