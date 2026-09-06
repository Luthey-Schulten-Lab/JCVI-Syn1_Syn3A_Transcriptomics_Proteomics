#!/usr/bin/env bash
# ============================================================
# 0. Setup & sanity check
# ============================================================
set -euo pipefail

# --- Executables (resolved from PATH; activate your conda env first) ---
# e.g. `conda activate Omics` (env recipe in environment.yml)
NUCMER=$(which nucmer 2>/dev/null)
DELTA_FILTER=$(which delta-filter 2>/dev/null)
SHOW_COORDS=$(which show-coords 2>/dev/null)
SHOW_SNPS=$(which show-snps 2>/dev/null)
DNADIFF=$(which dnadiff 2>/dev/null)
SAMTOOLS=$(which samtools 2>/dev/null)
BEDTOOLS=$(which bedtools 2>/dev/null)

# Fail loudly if anything is missing rather than silently producing junk
for tool in NUCMER DELTA_FILTER SHOW_COORDS SHOW_SNPS DNADIFF SAMTOOLS BEDTOOLS; do
    if [ -z "${!tool}" ]; then
        echo "ERROR: $tool not found on PATH. Activate the conda env (e.g. 'conda activate Omics')." >&2
        exit 1
    fi
done

# --- Inputs / params -----------------------------------------
# Resolve to absolute paths so they survive the `cd aln` below.
# Run this script from the Genome_Reduction/ directory.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
SYN1_FA="$PROJECT_ROOT/Genomes_Input/syn1_genome.fasta"
SYN3A_FA="$PROJECT_ROOT/Genomes_Input/syn3A_genome.fasta"
SYN1_GFF="$PROJECT_ROOT/Genomes_Input/syn1.genes.gff3"   # for deletion gene annotation
THREADS=8
MIN_DEL_BP=50                       # drop deletions smaller than this
MERGE_GAP=50                        # close alignment gaps ≤ this before complement

for f in "$SYN1_FA" "$SYN3A_FA" "$SYN1_GFF"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done

mkdir -p "$SCRIPT_DIR/aln/raw" && cd "$SCRIPT_DIR/aln"
exec > >(tee -a run.log) 2>&1
echo "[$(date)] starting alignment pipeline"

# Confirm sizes (~1,078,809 bp syn1; ~543,379 bp syn3A)
for f in "$SYN1_FA" "$SYN3A_FA"; do
    echo -n "$f: "; grep -v '^>' "$f" | tr -d '\n' | wc -c
done

# ============================================================
# 1. nucmer (MUMmer 3.23; NUCmer 3.1) — primary aligner
# ============================================================
# --maxmatch reports ALL matches (not just unique) — important because
# syn3A has rearrangements and possibly duplicated regions
# -c 100 = minimum cluster length 100 bp (default for closely related)
# -l 20 = minimum exact match seed (default; lower = more sensitive)

"$NUCMER" --maxmatch -c 100 -l 20 \
       -p raw/syn1_vs_syn3A \
       "$SYN1_FA" "$SYN3A_FA"

# Filter: keep only 1-to-1 best alignments above 95% identity, min 100 bp
"$DELTA_FILTER" -1 -i 95 -l 100 \
    raw/syn1_vs_syn3A.delta > raw/syn1_vs_syn3A.1delta

# Generate the coordinate table (this is your alignment-block table).
# Header columns match the show-coords -rclTH output.
{
    printf "S1\tE1\tS2\tE2\tLEN1\tLEN2\t%%IDY\tLENR\tLENQ\tCOVR\tCOVQ\tTAGR\tTAGQ\n"
    "$SHOW_COORDS" -rclTH raw/syn1_vs_syn3A.1delta
} > raw/syn1_vs_syn3A.coords.tsv
echo "nucmer 1-to-1 blocks: $(($(wc -l < raw/syn1_vs_syn3A.coords.tsv) - 1))"

# Full diff report (gives you SNPs, indels, breakpoints, summary stats)
"$DNADIFF" -p raw/dnadiff_out -d raw/syn1_vs_syn3A.delta
# Look at: dnadiff_out.report (summary), dnadiff_out.snps (single-bp diffs),
#          dnadiff_out.1coords (1-to-1 blocks), dnadiff_out.mcoords (many-to-many)

# Per-base SNP / indel list — gives you the single-bp resolution you want at edges.
# Header columns match the show-snps -ClrTH output.
{
    printf "P1\tSUB_R\tSUB_Q\tP2\tBUFF\tDIST\tLEN_R\tLEN_Q\tFRM_R\tFRM_Q\tTAG_R\tTAG_Q\n"
    "$SHOW_SNPS" -ClrTH raw/syn1_vs_syn3A.1delta
} > raw/syn1_vs_syn3A.snps.tsv

# ============================================================
# 2. Deletion-breakpoint extraction
# ============================================================
# Deletions = stretches of syn1 NOT covered by any nucmer 1-to-1 block.
# Convert coords to BED, merge near-adjacent blocks, take the complement.

# nucmer coords -> BED (syn1 coordinates).
# `#`-prefixed header is ignored by bedtools but visible to humans / pandas.
{
    printf "#chrom\tstart0\tend\tblock_id\tpct_identity\tstrand\n"
    awk 'BEGIN{OFS="\t"} NR>1 {print $12, $1-1, $2, "nuc_"(NR-1), $7, "+"}' \
        raw/syn1_vs_syn3A.coords.tsv
} > raw/nucmer_syn1.bed

# Genome file (chrom, length) for bedtools complement
"$SAMTOOLS" faidx "$SYN1_FA"
cut -f1,2 "${SYN1_FA}.fai" > raw/syn1.genome

# Merge blocks closing gaps ≤ MERGE_GAP, complement against the syn1 genome,
# then drop sub-MIN_DEL_BP intervals (boundary noise). Add a `#` header line.
{
    printf "#chrom\tstart0\tend\tlength_bp\n"
    sort -k1,1 -k2,2n raw/nucmer_syn1.bed \
        | "$BEDTOOLS" merge -d "$MERGE_GAP" -i - \
        | "$BEDTOOLS" complement -i - -g raw/syn1.genome \
        | awk -v m="$MIN_DEL_BP" 'BEGIN{OFS="\t"} ($3-$2) >= m {print $1, $2, $3, $3-$2}'
} > raw/syn1_deleted_regions.bed

echo "Number of deletion intervals (≥${MIN_DEL_BP} bp): $(($(wc -l < raw/syn1_deleted_regions.bed) - 1))"
echo -n "Total syn1 bp deleted in syn3A: "
awk 'BEGIN{FS=OFS="\t"} !/^#/ {sum += $3-$2} END {print sum+0}' raw/syn1_deleted_regions.bed

# ============================================================
# 3. Annotate deletions with syn1 gene content
# ============================================================
# For each deleted interval, list the syn1 genes it overlaps (full or partial).
# This is the Phase-1.1 deliverable: which genes were lost during minimization.
if [ -f "$SYN1_GFF" ]; then
    {
        printf "del_chrom\tdel_start0\tdel_end\tdel_length_bp\tgff_chrom\tgff_source\tgff_feature\tgff_start1\tgff_end\tgff_score\tgff_strand\tgff_frame\tgff_attributes\n"
        awk '$3=="gene"' "$SYN1_GFF" \
          | "$BEDTOOLS" intersect -a raw/syn1_deleted_regions.bed \
                                  -b - -wa -wb
    } > raw/syn1_deleted_genes.tsv
    echo "Deleted intervals overlapping a syn1 gene: $(awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$2,$3}' raw/syn1_deleted_genes.tsv | sort -u | wc -l)"
else
    echo "WARN: $SYN1_GFF not found — skipping gene annotation step"
fi

echo "[$(date)] pipeline complete"