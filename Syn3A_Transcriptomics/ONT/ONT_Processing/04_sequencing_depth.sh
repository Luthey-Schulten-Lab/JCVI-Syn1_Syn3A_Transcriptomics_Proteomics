#!/bin/bash
mkdir -p depth_bedgraph

# Define threads
THREADS=16
BAM="./syn3A.ONT.rep1.sorted.bam"

# Primary-only filter:
#   -F 0x904 excludes unmapped (0x4), secondary (0x100), and supplementary (0x800)
# Per qc_report.txt: 14579 secondary + 3812 supplementary alignments are dropped,
# leaving the 541272 primary mapped reads for depth calculation.
PRIMARY_FILTER="-F 0x904"

# 1) Total Genome-scale bedGraph (primary alignments only)
# Uses -a (all positions) and -d 0 (no depth cap)
samtools view -@ $THREADS -u $PRIMARY_FILTER "$BAM" | \
samtools depth -a -d 0 - | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ./depth_bedgraph/syn3A.ONT.rep1.total.bedGraph

# 2) Forward-strand bedGraph (+ strand, primary only)
# -F 0x914 = unmapped + secondary + supplementary + reverse-mapped (0x10)
samtools view -@ $THREADS -u -F 0x914 "$BAM" | \
samtools depth -a -d 0 - | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ./depth_bedgraph/syn3A.ONT.rep1.plus.bedGraph

# 3) Reverse-strand bedGraph (- strand, primary only)
# -f 0x10 keeps reverse-mapped reads; -F 0x904 excludes unmapped/secondary/supplementary
samtools view -@ $THREADS -u -f 0x10 -F 0x904 "$BAM" | \
samtools depth -a -d 0 - | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ./depth_bedgraph/syn3A.ONT.rep1.minus.bedGraph
