#!/bin/bash
mkdir -p depth_bedgraph

# Define threads
THREADS=16
BAM="./syn1.PacBio.FLNC.sorted.HQ.bam"   # BAM lives next to this script after the reorg

# 1) Total Genome-scale bedGraph
# Uses -a (all positions) and -d 0 (no depth cap)
samtools depth -@ $THREADS -a -d 0 "$BAM" | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ./depth_bedgraph/syn1.PacBio.FLNC.HQ.total.bedGraph

# 2) Forward-strand bedGraph (+ strand)
# -F 0x10 filters out reverse-mapped reads
samtools view -@ $THREADS -u -F 0x10 "$BAM" | \
samtools depth -a -d 0 - | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ./depth_bedgraph/syn1.PacBio.FLNC.HQ.plus.bedGraph

# 3) Reverse-strand bedGraph (- strand)
# -f 0x10 keeps only reverse-mapped reads
samtools view -@ $THREADS -u -f 0x10 "$BAM" | \
samtools depth -a -d 0 - | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ./depth_bedgraph/syn1.PacBio.FLNC.HQ.minus.bedGraph