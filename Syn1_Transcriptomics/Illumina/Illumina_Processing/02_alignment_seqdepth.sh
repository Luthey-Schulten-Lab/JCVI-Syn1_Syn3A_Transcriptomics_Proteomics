#!/bin/bash
mkdir -p bam logs depth_bedgraph ref
THREADS=16
INDEX="ref/syn1"
REF_FASTA="../../../Genomes_Input/syn1_genome.fasta"

# ---------------------------------------------------------------
# Build the bowtie2 reference index from the syn1 genome FASTA if it
# doesn't already exist. Six files are produced under ref/:
#   ref/syn1.{1,2,3,4,rev.1,rev.2}.bt2
# The index is a derived artefact and is gitignored, so it must be
# rebuilt on a fresh clone.
# ---------------------------------------------------------------
if [ ! -s "${INDEX}.1.bt2" ]; then
    echo "[BUILD-REF] bowtie2-build ${REF_FASTA} -> ${INDEX}"
    bowtie2-build "${REF_FASTA}" "${INDEX}" > logs/bowtie2_build.log 2>&1
else
    echo "[BUILD-REF] ${INDEX}.1.bt2 exists, skipping bowtie2-build."
fi

# Fail fast: without an index every downstream samtools step fails on an
# empty stream, which buries the real error under a cascade.
if [ ! -s "${INDEX}.1.bt2" ]; then
    echo "ERROR: bowtie2 index ${INDEX} was not built. See logs/bowtie2_build.log" >&2
    exit 1
fi

for R1 in ../Illumina_Raw/*_1.fastq; do
    R2=${R1/_1.fastq/_2.fastq}
    SAMPLE=$(basename $R1 _1.fastq)
    echo "[ALIGN] $SAMPLE"

    # # Align with bowtie2 (default paired-end mode)
    bowtie2 \
        -x $INDEX \
        -1 $R1 -2 $R2 \
        --threads $THREADS \
        2> logs/${SAMPLE}_bowtie2.log \
    | samtools view -bS - \
    | samtools sort -@ $THREADS -o bam/${SAMPLE}.sorted.bam

    # samtools index bam/${SAMPLE}.sorted.bam

    # ---------------------------------------------------------------
    # Split into strand-specific BAMs
    # Assumes dUTP / fr-firststrand library (most common):
    #   Read2 maps to the transcript strand (sense)
    #   Read1 maps to the opposite strand (antisense)
    #
    # Plus-strand reads (reads whose alignment represents + strand transcription):
    #   - Read2 forward on + strand:  -f 128 -F 16   (read2, not reverse)
    #   - Read1 reverse on + strand:  -f 80           (read1, reverse) = -f 64 + -f 16
    #
    # Minus-strand reads (reads representing - strand transcription):
    #   - Read2 reverse on - strand:  -f 144          (read2, reverse) = -f 128 + -f 16
    #   - Read1 forward on - strand:  -f 64 -F 16     (read1, not reverse)
    # ---------------------------------------------------------------

    BAM="bam/${SAMPLE}.sorted.bam"

    echo "[STRAND SPLIT] $SAMPLE - plus strand"
    samtools view -b -f 128 -F 16 $BAM > bam/${SAMPLE}.fwd1.bam
    samtools view -b -f 80        $BAM > bam/${SAMPLE}.fwd2.bam
    samtools merge -f bam/${SAMPLE}.plus.bam bam/${SAMPLE}.fwd1.bam bam/${SAMPLE}.fwd2.bam
    samtools index bam/${SAMPLE}.plus.bam

    echo "[STRAND SPLIT] $SAMPLE - minus strand"
    samtools view -b -f 144       $BAM > bam/${SAMPLE}.rev1.bam
    samtools view -b -f 64 -F 16  $BAM > bam/${SAMPLE}.rev2.bam
    samtools merge -f bam/${SAMPLE}.minus.bam bam/${SAMPLE}.rev1.bam bam/${SAMPLE}.rev2.bam
    samtools index bam/${SAMPLE}.minus.bam

    # Generate bedGraph files for each strand
    echo "[DEPTH] $SAMPLE"
    samtools depth -a -d 0 bam/${SAMPLE}.plus.bam \
        | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' \
        > depth_bedgraph/${SAMPLE}.plus.bedGraph
 
    samtools depth -a -d 0 bam/${SAMPLE}.minus.bam \
        | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' \
        > depth_bedgraph/${SAMPLE}.minus.bedGraph

    # Clean up intermediate files
    rm -f bam/${SAMPLE}.fwd1.bam bam/${SAMPLE}.fwd2.bam
    rm -f bam/${SAMPLE}.rev1.bam bam/${SAMPLE}.rev2.bam

    echo "[DONE] $SAMPLE"
done