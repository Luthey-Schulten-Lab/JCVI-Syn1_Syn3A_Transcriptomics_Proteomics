#!/usr/bin/env bash
set -eu
#############################################
# QC report for Syn1 ONT direct-RNA mapping (rep1, rep2, merged).
# Per BAM: overall mapping stats; length/quality of primary mapped and unmapped
# reads; secondary-alignment contribution + top repeated loci.
#############################################
command -v samtools >/dev/null || { echo "ERROR: samtools not found (activate Omics env)" >&2; exit 1; }
cd "$(dirname "$0")"

report_bam() {
  local BAM="$1"
  local REPORT="${BAM}.qc_report.txt"
  [[ -s "${BAM}" ]] || { echo "WARN: ${BAM} missing, skipping" >&2; return; }
  out() { echo "$@" | tee -a "${REPORT}"; }
  > "${REPORT}"

  out "=========================================="
  out "ONT Direct RNA Mapping QC Report"
  out "BAM: ${BAM}"
  out "Date: $(date)"
  out "=========================================="

  out ""
  out "[1] OVERALL MAPPING STATISTICS"
  out "----------------------------------------"
  samtools flagstat "${BAM}" | tee -a "${REPORT}"

  local TOTAL PRIMARY PRIMARY_MAPPED UNMAPPED SECONDARY SUPPLEMENTARY
  TOTAL=$(samtools view -c "${BAM}")
  PRIMARY=$(samtools view -c -F 0x100 -F 0x800 "${BAM}")
  PRIMARY_MAPPED=$(samtools view -c -F 0x100 -F 0x800 -F 0x4 "${BAM}")
  UNMAPPED=$(samtools view -c -f 0x4 "${BAM}")
  SECONDARY=$(samtools view -c -f 0x100 "${BAM}")
  SUPPLEMENTARY=$(samtools view -c -f 0x800 "${BAM}")
  out ""
  out "Summary numbers:"
  out "  Total alignment records:      ${TOTAL}"
  out "  Primary reads (input reads):  ${PRIMARY}"
  out "  Primary mapped:               ${PRIMARY_MAPPED}"
  out "  Unmapped:                     ${UNMAPPED}"
  out "  Secondary alignments:         ${SECONDARY}"
  out "  Supplementary alignments:     ${SUPPLEMENTARY}"
  if [ "${PRIMARY}" -gt 0 ]; then
    out "  Mapping rate:                 $(awk -v m=$PRIMARY_MAPPED -v t=$PRIMARY 'BEGIN{printf "%.2f", 100*m/t}')%"
    out "  Unmapped rate:                $(awk -v u=$UNMAPPED -v t=$PRIMARY 'BEGIN{printf "%.2f", 100*u/t}')%"
  fi

  out ""
  out "[2] PRIMARY MAPPED READS - Length and Quality Features"
  out "----------------------------------------"
  samtools view -F 0x4 -F 0x100 -F 0x800 "${BAM}" | _len_qual | tee -a "${REPORT}"

  out ""
  out "[3] UNMAPPED READS - Length and Quality Features"
  out "----------------------------------------"
  samtools view -f 0x4 "${BAM}" | _len_qual | tee -a "${REPORT}"

  out ""
  out "[4] SECONDARY ALIGNMENTS - Contribution and Locations"
  out "----------------------------------------"
  if [ "${SECONDARY}" -eq 0 ]; then
    out "  No secondary alignments."
  else
    if [ "${PRIMARY_MAPPED}" -gt 0 ]; then
      out "  Secondary alignments: ${SECONDARY}"
      out "  As fraction of primary mapped reads: $(awk -v s=$SECONDARY -v p=$PRIMARY_MAPPED 'BEGIN{printf "%.2f", 100*s/p}')%"
    fi
    local N_READS_WITH_SEC
    N_READS_WITH_SEC=$(samtools view -f 0x100 "${BAM}" | awk '{print $1}' | sort -u | wc -l)
    out ""
    out "  ${N_READS_WITH_SEC} unique reads have >=1 secondary alignment"
    out ""
    out "  Top 20 secondary alignment positions (chrom, start) - likely repeated/duplicated loci:"
    samtools view -f 0x100 "${BAM}" | awk '{print $3"\t"$4}' | sort | uniq -c | sort -rn | head -20 \
      | awk '{printf "    %8d  %s:%s\n", $1, $2, $3}' | tee -a "${REPORT}"
  fi
  out ""
  out "Report written to: ${REPORT}"
}

# shared length/quality summariser (reads SAM records on stdin)
_len_qual() {
  awk '
  BEGIN { for(i=0;i<256;i++) _ord_[sprintf("%c",i)]=i }
  function ord(c){ return _ord_[c] }
  {
    L=length($10); qual=$11
    if(L<100)lbin="<100"; else if(L<300)lbin="100-300"; else if(L<500)lbin="300-500";
    else if(L<1000)lbin="500-1000"; else if(L<2000)lbin="1000-2000"; else lbin=">2000";
    lcount[lbin]++
    qsum=0; n=length(qual); for(i=1;i<=n;i++) qsum+=ord(substr(qual,i,1))-33
    meanQ=(n>0?qsum/n:0)
    if(meanQ<7)qbin="<7"; else if(meanQ<10)qbin="7-10"; else if(meanQ<15)qbin="10-15";
    else if(meanQ<20)qbin="15-20"; else qbin=">=20";
    qcount[qbin]++
    total++; Lsum+=L; Qsum+=meanQ; if(L>maxL)maxL=L; if(minL==0||L<minL)minL=L
  }
  END {
    if(total==0){print "  No reads."; exit}
    printf "  Total reads: %d\n", total
    printf "  Length:    mean=%.0f  min=%d  max=%d\n", Lsum/total, minL, maxL
    printf "  Mean qual: mean=%.2f\n\n", Qsum/total
    print "  Length distribution:"
    n=split("<100 100-300 300-500 500-1000 1000-2000 >2000",order," ")
    for(i=1;i<=n;i++){b=order[i]; printf "    %-12s %10d  %6.2f%%\n", b, lcount[b]+0, 100*(lcount[b]+0)/total}
    print "\n  Mean per-read quality distribution (Phred):"
    n=split("<7 7-10 10-15 15-20 >=20",order2," ")
    for(i=1;i<=n;i++){b=order2[i]; printf "    %-8s %10d  %6.2f%%\n", b, qcount[b]+0, 100*(qcount[b]+0)/total}
  }'
}
export -f _len_qual

for BAM in syn1.ONT.rep1.sorted.bam syn1.ONT.rep2.sorted.bam syn1.ONT.merged.sorted.bam; do
  echo "########## ${BAM} ##########"
  report_bam "${BAM}"
  echo
done
