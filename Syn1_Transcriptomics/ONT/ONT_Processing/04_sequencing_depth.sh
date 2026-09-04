#!/usr/bin/env bash
set -eu
#############################################
# Strand-split per-base depth bedGraphs for Syn1 ONT direct-RNA BAMs.
# Reads are 5'->3' sense (see 00_retrieve_fastq.sh), so BAM strand == transcript
# strand: forward-mapped reads (FLAG & 0x10 == 0) are the + (plus) transcripts.
#
# Emits total/plus/minus for each rep and for the merged BAM.
#############################################
THREADS=16
mkdir -p depth_bedgraph
cd "$(dirname "$0")"
command -v samtools >/dev/null || { echo "ERROR: samtools not found (activate Omics env)" >&2; exit 1; }

# Primary-only filter: -F 0x904 drops unmapped(0x4)+secondary(0x100)+supplementary(0x800)
for OUT in syn1.ONT.rep1 syn1.ONT.rep2 syn1.ONT.merged; do
  BAM="${OUT}.sorted.bam"
  [[ -s "${BAM}" ]] || { echo "WARN: ${BAM} missing, skipping" >&2; continue; }
  echo "[depth] ${BAM}"

  # 1) total (primary alignments, both strands)
  samtools view -@ $THREADS -u -F 0x904 "$BAM" \
    | samtools depth -a -d 0 - \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' > "depth_bedgraph/${OUT}.total.bedGraph"

  # 2) plus strand: also drop reverse-mapped (0x10) -> -F 0x914
  samtools view -@ $THREADS -u -F 0x914 "$BAM" \
    | samtools depth -a -d 0 - \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' > "depth_bedgraph/${OUT}.plus.bedGraph"

  # 3) minus strand: keep reverse-mapped (0x10), drop unmapped/secondary/supplementary
  samtools view -@ $THREADS -u -f 0x10 -F 0x904 "$BAM" \
    | samtools depth -a -d 0 - \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' > "depth_bedgraph/${OUT}.minus.bedGraph"
done
echo "Done."
