# Nice explanation of Illumina pair-ended sequencing https://www.youtube.com/watch?v=fCd6B5HRaZ8

#!/bin/bash

mkdir -p qc logs

THREADS=16

for R1 in ../Illumina_Raw/*_1.fastq; do
    R2=${R1/_1.fastq/_2.fastq}
    SAMPLE=$(basename $R1 _1.fastq)

    echo "[QC] $SAMPLE"

    fastqc -t $THREADS \
        -o qc \
        "$R1" "$R2" \
        > logs/${SAMPLE}_fastqc.log 2>&1
done

# Aggregate
multiqc qc -o qc