#!/usr/bin/env bash
set -eu
#############################################
# QC report for ONT Direct RNA mapping output
#
# Reports:
#   1. Overall mapping statistics
#   2. Quality/length features of primary mapped reads
#   3. Features of unmapped reads
#   4. Contribution and locations of secondary alignments
#
# Usage: ./qc_ONT_mapping.sh <sorted.bam> [output_report.txt]
#############################################

# ---------- Tool checks ----------
SAMTOOLS_PATH=$(which samtools 2>/dev/null)
if [ -z "$SAMTOOLS_PATH" ]; then
  echo "Error: samtools not found in the active conda environment." >&2
  exit 1
fi

# ---------- Inputs ----------
# if [ $# -lt 1 ]; then
#   echo "Usage: $0 <sorted.bam> [output_report.txt]" >&2
#   exit 1
# fi

BAM="./syn3A.ONT.rep1.sorted.bam"
REPORT="${BAM}.qc_report.txt"

[[ -s "${BAM}" ]] || { echo "ERROR: BAM not found or empty: ${BAM}" >&2; exit 2; }

# ---------- Helper: write to both stdout and report file ----------
out() { echo "$@" | tee -a "${REPORT}"; }

# Wipe the report file
> "${REPORT}"

out "=========================================="
out "ONT Direct RNA Mapping QC Report"
out "BAM: ${BAM}"
out "Date: $(date)"
out "=========================================="

#############################################
# 1) Overall mapping statistics
#############################################
out ""
out "[1] OVERALL MAPPING STATISTICS"
out "----------------------------------------"
$SAMTOOLS_PATH flagstat "${BAM}" | tee -a "${REPORT}"

# Extract key numbers for downstream use and percentages
TOTAL=$($SAMTOOLS_PATH view -c "${BAM}")
PRIMARY=$($SAMTOOLS_PATH view -c -F 0x100 -F 0x800 "${BAM}")
PRIMARY_MAPPED=$($SAMTOOLS_PATH view -c -F 0x100 -F 0x800 -F 0x4 "${BAM}")
UNMAPPED=$($SAMTOOLS_PATH view -c -f 0x4 "${BAM}")
SECONDARY=$($SAMTOOLS_PATH view -c -f 0x100 "${BAM}")
SUPPLEMENTARY=$($SAMTOOLS_PATH view -c -f 0x800 "${BAM}")

out ""
out "Summary numbers:"
out "  Total alignment records:      ${TOTAL}"
out "  Primary reads (input reads):  ${PRIMARY}"
out "  Primary mapped:               ${PRIMARY_MAPPED}"
out "  Unmapped:                     ${UNMAPPED}"
out "  Secondary alignments:         ${SECONDARY}"
out "  Supplementary alignments:     ${SUPPLEMENTARY}"
if [ "${PRIMARY}" -gt 0 ]; then
  PCT_MAPPED=$(awk -v m=${PRIMARY_MAPPED} -v t=${PRIMARY} 'BEGIN {printf "%.2f", 100*m/t}')
  PCT_UNMAPPED=$(awk -v u=${UNMAPPED} -v t=${PRIMARY} 'BEGIN {printf "%.2f", 100*u/t}')
  out "  Mapping rate:                 ${PCT_MAPPED}%"
  out "  Unmapped rate:                ${PCT_UNMAPPED}%"
fi

#############################################
# 2) Quality and length features of primary mapped reads
#############################################
out ""
out "[2] PRIMARY MAPPED READS - Length and Quality Features"
out "----------------------------------------"

# -F 0x4   exclude unmapped
# -F 0x100 exclude secondary
# -F 0x800 exclude supplementary
$SAMTOOLS_PATH view -F 0x4 -F 0x100 -F 0x800 "${BAM}" | \
  awk '
  BEGIN { for(i=0; i<256; i++) _ord_[sprintf("%c", i)] = i }
  function ord(c) { return _ord_[c] }
  {
    seq = $10
    qual = $11
    L = length(seq)

    # length bins
    if(L<100) lbin="<100";
    else if(L<300) lbin="100-300";
    else if(L<500) lbin="300-500";
    else if(L<1000) lbin="500-1000";
    else if(L<2000) lbin="1000-2000";
    else lbin=">2000";
    lcount[lbin]++

    # mean quality (Phred, ASCII offset 33)
    qsum = 0
    n = length(qual)
    for(i=1; i<=n; i++) qsum += ord(substr(qual, i, 1)) - 33
    meanQ = (n > 0 ? qsum / n : 0)

    # quality bins
    if(meanQ<7) qbin="<7";
    else if(meanQ<10) qbin="7-10";
    else if(meanQ<15) qbin="10-15";
    else if(meanQ<20) qbin="15-20";
    else qbin=">=20";
    qcount[qbin]++

    # totals
    total++
    Lsum += L
    Qsum += meanQ
    if(L>maxL) maxL=L
    if(minL==0 || L<minL) minL=L
  }
  END {
    if(total == 0) { print "No primary mapped reads."; exit }
    printf "  Total primary mapped reads: %d\n", total
    printf "  Length:    mean=%.0f  min=%d  max=%d\n", Lsum/total, minL, maxL
    printf "  Mean qual: mean=%.2f\n", Qsum/total
    print ""
    print "  Length distribution:"
    n=split("<100 100-300 300-500 500-1000 1000-2000 >2000", order, " ")
    for(i=1; i<=n; i++) {
      b=order[i]
      printf "    %-12s %10d  %6.2f%%\n", b, lcount[b]+0, 100*(lcount[b]+0)/total
    }
    print ""
    print "  Mean per-read quality distribution (Phred):"
    n=split("<7 7-10 10-15 15-20 >=20", order2, " ")
    for(i=1; i<=n; i++) {
      b=order2[i]
      printf "    %-8s %10d  %6.2f%%\n", b, qcount[b]+0, 100*(qcount[b]+0)/total
    }
  }' | tee -a "${REPORT}"

#############################################
# 3) Features of UNMAPPED reads
#############################################
out ""
out "[3] UNMAPPED READS - Length and Quality Features"
out "----------------------------------------"

$SAMTOOLS_PATH view -f 0x4 "${BAM}" | \
  awk '
  BEGIN { for(i=0; i<256; i++) _ord_[sprintf("%c", i)] = i }
  function ord(c) { return _ord_[c] }
  {
    seq = $10
    qual = $11
    L = length(seq)

    if(L<100) lbin="<100";
    else if(L<300) lbin="100-300";
    else if(L<500) lbin="300-500";
    else if(L<1000) lbin="500-1000";
    else if(L<2000) lbin="1000-2000";
    else lbin=">2000";
    lcount[lbin]++

    qsum = 0
    n = length(qual)
    for(i=1; i<=n; i++) qsum += ord(substr(qual, i, 1)) - 33
    meanQ = (n > 0 ? qsum / n : 0)

    if(meanQ<7) qbin="<7";
    else if(meanQ<10) qbin="7-10";
    else if(meanQ<15) qbin="10-15";
    else if(meanQ<20) qbin="15-20";
    else qbin=">=20";
    qcount[qbin]++

    total++
    Lsum += L
    Qsum += meanQ
  }
  END {
    if(total == 0) { print "No unmapped reads."; exit }
    printf "  Total unmapped reads:  %d\n", total
    printf "  Length:    mean=%.0f\n", Lsum/total
    printf "  Mean qual: mean=%.2f\n", Qsum/total
    print ""
    print "  Length distribution:"
    n=split("<100 100-300 300-500 500-1000 1000-2000 >2000", order, " ")
    for(i=1; i<=n; i++) {
      b=order[i]
      printf "    %-12s %10d  %6.2f%%\n", b, lcount[b]+0, 100*(lcount[b]+0)/total
    }
    print ""
    print "  Mean per-read quality distribution (Phred):"
    n=split("<7 7-10 10-15 15-20 >=20", order2, " ")
    for(i=1; i<=n; i++) {
      b=order2[i]
      printf "    %-8s %10d  %6.2f%%\n", b, qcount[b]+0, 100*(qcount[b]+0)/total
    }
  }' | tee -a "${REPORT}"

#############################################
# 4) Secondary alignment contribution
#############################################
out ""
out "[4] SECONDARY ALIGNMENTS - Contribution and Locations"
out "----------------------------------------"

if [ "${SECONDARY}" -eq 0 ]; then
  out "  No secondary alignments."
else
  if [ "${PRIMARY_MAPPED}" -gt 0 ]; then
    PCT_SEC=$(awk -v s=${SECONDARY} -v p=${PRIMARY_MAPPED} 'BEGIN {printf "%.2f", 100*s/p}')
    out "  Secondary alignments: ${SECONDARY}"
    out "  As fraction of primary mapped reads: ${PCT_SEC}%"
  fi

  out ""
  out "  How many primary mapped reads have at least one secondary alignment?"
  # Read names that have a secondary alignment
  N_READS_WITH_SEC=$($SAMTOOLS_PATH view -f 0x100 "${BAM}" | awk '{print $1}' | sort -u | wc -l)
  out "    ${N_READS_WITH_SEC} unique reads have >=1 secondary alignment"
  if [ "${PRIMARY_MAPPED}" -gt 0 ]; then
    PCT_READS_SEC=$(awk -v r=${N_READS_WITH_SEC} -v p=${PRIMARY_MAPPED} 'BEGIN {printf "%.2f", 100*r/p}')
    out "    (${PCT_READS_SEC}% of primary mapped reads)"
  fi

  out ""
  out "  Top 20 secondary alignment positions (chrom, start) — likely repeated/duplicated loci:"
  $SAMTOOLS_PATH view -f 0x100 "${BAM}" | \
    awk '{print $3"\t"$4}' | sort | uniq -c | sort -rn | head -20 | \
    awk '{printf "    %8d  %s:%s\n", $1, $2, $3}' | tee -a "${REPORT}"
fi

#############################################
# Done
#############################################
out ""
out "=========================================="
out "Report written to: ${REPORT}"
out "=========================================="