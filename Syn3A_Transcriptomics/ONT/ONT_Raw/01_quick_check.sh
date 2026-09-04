#!/bin/bash
# Check the adapter and polyA status of ONT Direct RNA reads before mapping

# Quick QC check for ONT Direct RNA (SQK-RNA002) reads from SRA
# Context: Guppy was run with --reverse_sequence yes, so reads are 5'->3' orientation
#   - Adapter (RMX) trimming: --trim_strategy rna should have removed 3' adapter
#   - PolyA tail: if present, would be at the 3' END of reads (post-reversal)
#   - Genome is ~70% AT-rich (Mycoplasma), so A-runs alone don't mean polyA
#
# Usage: ./ont_quick_check.sh <fastq_file>
#   Accepts .fastq or .fastq.gz

set -eu
# Note: not using pipefail because `head` closes pipes early (SIGPIPE),
# which would falsely abort the script on otherwise-successful commands.

# if [ $# -ne 1 ]; then
#     echo "Usage: $0 <fastq_file>"
#     exit 1
# fi

FQ="ONT.syn3A.rep1.fastq.gz" # hardcoded for now, can be made an argument later

# Pick the right reader
if [[ "$FQ" == *.gz ]]; then
    READER="zcat"
else
    READER="cat"
fi

echo "=========================================="
echo "ONT Direct RNA QC Check: $FQ"
echo "=========================================="

# ---- 1. Read length distribution ----
echo ""
echo "[1] Read length distribution:"
$READER "$FQ" | awk 'NR%4==2 {print length($0)}' | \
  awk '{
    if($1<100) bin="<100";
    else if($1<500) bin="100-500";
    else if($1<1000) bin="500-1000";
    else if($1<2000) bin="1000-2000";
    else bin=">2000";
    counts[bin]++; total++
  } END {
    printf "  %-12s %10s %8s\n", "Length", "Count", "Percent"
    for (b in counts) printf "  %-12s %10d %7.2f%%\n", b, counts[b], 100*counts[b]/total
    printf "  Total reads: %d\n", total
  }' | sort

# ---- 2. Check 3' end for polyA (post-reversal, polyA would be here) ----
echo ""
echo "[2] Last 60 nt of long reads (>300 nt) - checking for polyA tail at 3' end:"
echo "    (Look for clear A-runs of 20+ consecutive A's)"
$READER "$FQ" | awk 'NR%4==2 && length($0)>300' | head -10 | \
  awk '{print "  ..." substr($0, length($0)-59)}'

# ---- 3. Quantify polyA presence at 3' end ----
echo ""
echo "[3] PolyA quantification at 3' end (long reads only):"
echo "    Counting reads ending in >=15 consecutive A's within last 50 nt"
$READER "$FQ" | awk 'NR%4==2 && length($0)>300' | \
  awk '{
    tail = substr($0, length($0)-49)
    total++
    if (tail ~ /A{15,}/) polya++
  } END {
    if (total > 0) {
      printf "  Reads with polyA-like 3-prime end: %d / %d (%.2f%%)\n", polya, total, 100*polya/total
    } else {
      print "  No reads >300 nt found"
    }
  }'

# ---- 4. Check 5' end for adapter remnants ----
echo ""
echo "[4] First 60 nt of long reads (>300 nt) - checking for adapter remnants at 5' end:"
echo "    (Should look like normal mRNA sequence, not adapter motifs)"
$READER "$FQ" | awk 'NR%4==2 && length($0)>300' | head -10 | \
  awk '{print "  " substr($0, 1, 60) "..."}'

# ---- 5. Overall AT content (sanity check vs ~70% expected for Mycoplasma) ----
echo ""
echo "[5] Overall AT content (Mycoplasma syn1.0/syn3A expected ~70%):"
$READER "$FQ" | awk 'NR%4==2' | head -10000 | \
  awk '{
    seq = toupper($0)
    for (i=1; i<=length(seq); i++) {
      base = substr(seq, i, 1)
      counts[base]++
      total++
    }
  } END {
    at = (counts["A"] + counts["T"]) / total * 100
    gc = (counts["G"] + counts["C"]) / total * 100
    printf "  AT: %.2f%%  |  GC: %.2f%%  (sampled from first 10k reads)\n", at, gc
  }'

echo ""
echo "=========================================="
echo "Interpretation guide:"
echo "  - If [3] shows >50%% polyA-positive reads: polyA is INTACT, leave it for minimap2 soft-clipping"
echo "  - If [3] shows <10%% polyA-positive reads: polyA was TRIMMED before SRA upload"
echo "  - If [4] shows clean sequence (no repetitive motifs): adapter trimming worked"
echo "  - If [5] shows ~70%% AT: confirms Mycoplasma genome, A-richness is biological not polyA"
echo "=========================================="