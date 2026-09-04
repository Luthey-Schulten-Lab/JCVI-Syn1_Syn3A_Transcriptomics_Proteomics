#!/usr/bin/env bash
# Retrieve the syn3A Illumina RNA-seq data, wild type at log phase
# (Sandberg et al. 2023, iScience).
# Library protocol details (Ribo-Zero Gram-Negative + Kapa Stranded RNA-Seq,
# dUTP / fr-firststrand) are in:
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008596#sec009
#
# Note: The depositor uploaded the two mates of ONE paired-end run as two
# separate SRA accessions (SRR19432056 = R1, SRR19432057 = R2 — confirmed
# by identical Illumina cluster IDs across the two .fastq files). We
# download both then rename to the paired _1.fastq / _2.fastq layout that
# 01_quality_control.sh and 02_alignment_seqdepth.sh expect.

# Sample-name <-> SRR-pair mapping. Add more pairs as more replicates appear.
SAMPLE_NAME="syn3A_rep1"
SRR_R1="SRR19432056"   # Read 1 (mate 1) of the paired-end run
SRR_R2="SRR19432057"   # Read 2 (mate 2) of the same paired-end run

FASTERQ_PATH=$(which fasterq-dump 2>/dev/null)
if [ -z "$FASTERQ_PATH" ]; then
    echo "Error: fasterq-dump not found in the active conda environment. Please activate the appropriate environment."
    exit 1
fi

# 1) Download each SRA accession (skip if already on disk)
for ID in "${SRR_R1}" "${SRR_R2}"; do
    if [ -s "${ID}.fastq" ]; then
        echo "[SKIP] ${ID}.fastq already present."
        continue
    fi
    if [ -s "${SAMPLE_NAME}_1.fastq" ] && [ -s "${SAMPLE_NAME}_2.fastq" ]; then
        echo "[SKIP] ${SAMPLE_NAME}_{1,2}.fastq already present; nothing to download."
        break
    fi
    echo "[DOWNLOAD] ${ID} ..."
    "$FASTERQ_PATH" --split-files --progress -e 16 "${ID}"
done

# 2) Rename the two single-file downloads to a paired _1 / _2 layout.
if [ -s "${SAMPLE_NAME}_1.fastq" ] && [ -s "${SAMPLE_NAME}_2.fastq" ]; then
    echo "[RENAME] ${SAMPLE_NAME}_{1,2}.fastq already in place."
elif [ -s "${SRR_R1}.fastq" ] && [ -s "${SRR_R2}.fastq" ]; then
    echo "[RENAME] ${SRR_R1}.fastq  ->  ${SAMPLE_NAME}_1.fastq"
    mv -n "${SRR_R1}.fastq" "${SAMPLE_NAME}_1.fastq"
    echo "[RENAME] ${SRR_R2}.fastq  ->  ${SAMPLE_NAME}_2.fastq"
    mv -n "${SRR_R2}.fastq" "${SAMPLE_NAME}_2.fastq"
else
    echo "WARNING: expected ${SRR_R1}.fastq and ${SRR_R2}.fastq after download but did not find both." >&2
    exit 2
fi

echo "[DONE] Paired FASTQs ready: ${SAMPLE_NAME}_1.fastq, ${SAMPLE_NAME}_2.fastq"
