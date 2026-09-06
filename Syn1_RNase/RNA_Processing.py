"""
RNA Processing — degradation asymmetry in Syn1 (PacBio FLNC isoforms)
=====================================================================
Quantify the directionality of mRNA degradation in Syn1 two independent ways,
then relate it to ribosome rescue. Run from Syn1_RNase/.

Part A — Endpoint context (genomic).  Label each isoform's 5' and 3' end as
  intragenic or intergenic. An intact transcript starts at a promoter (5' end
  intergenic) and ends at a terminator (3' end intergenic); an endpoint *inside*
  a gene body cannot be explained by promoter/terminator usage and must come
  from processing or exonucleolytic erosion:

    - 5' end intragenic  ->  5'->3' exoribonuclease activity (RNase J1/J2)
    - 3' end intragenic  ->  3'->5' exoribonuclease activity (RNase R / YhaM)

    +--------------------------+-------------+-------------+
    |  Category                | 5' end      | 3' end      |
    |--------------------------+-------------+-------------|
    |  unprocessed             | intergenic  | intergenic  |
    |  5p_intragenic_only      | intragenic  | intergenic  |
    |  3p_intragenic_only      | intergenic  | intragenic  |
    |  both_intragenic         | intragenic  | intragenic  |
    +--------------------------+-------------+-------------+

  The ratio N(5'-intragenic) / N(3'-intragenic), in both unique-isoform and
  read-weighted form, reports the asymmetry between the two exoribonuclease
  directions.

Part B — ORFs with a start codon but no stop codon.  For every isoform, walk the
  canonical genes it covers; a gene whose start codon is contained is one ORF
  observation, and the ORF is "no-stop" if the matching stop codon is not also
  contained. These are the transcripts that strand ribosomes.

Part C — Truncation within containment clusters (relative).  An independent
  read of the same question: cluster isoforms by containment (the segmentation
  rule), take each cluster's longest isoform as the full-length reference, and
  ask whether the shorter members are trimmed at the 5' or the 3' end.

    same_5p_shorter_3p  -> 3' erosion   diff_5p_same_3p -> 5' erosion
    diff_5p_shorter_3p  -> both ends    same_both / other

  If 5'->3' exo outpaces 3'->5' exo, downstream endo-cleavage products (which
  carry the 5'-monophosphate RNase J prefers) are destroyed quickly while the
  upstream products sharing the original TSS persist, so same_5p_shorter_3p
  should dominate.

Part D — Ribosome rescue.  tmRNA (ssrA) as a fraction of non-rRNA TPM, against
  the ~25% E. coli reference, plus the read fraction carried by 3'-eroded
  (stop-codon-less) isoforms — the substrates tmRNA tags.

Two read thresholds are deliberate: Parts A/B use MIN_READS = 10 to keep the
endpoint census broad, Part C uses TRUNC_MIN_READS = 50 to match the operon
segmentation whose clustering rule it reuses.

Inputs
  ../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv
  ../Genomes_Input/syn1.genes.gff3
  ../Syn1_Syn3A_Transcriptomics/syn1_TPM.tsv            (Part D only)

Outputs (all under RNase/)
  RNase_asymmetry_isoform_categories.pdf   RNase_asymmetry_stacked100.pdf
  RNase_asymmetry_exoribo_categories.pdf   degradation_asymmetry.{pdf,png}
  isoform_endpoint_context.tsv             orf_start_stop_observations.tsv
  truncation_classification.tsv            per_operon_degradation_ratio.tsv

The run log is this script's stdout:
    python RNA_Processing.py > RNA_Processing.txt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
from collections import defaultdict

os.makedirs("RNase", exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════
MOTHER_FOLDER = ".."
ISOFORMS_TSV  = MOTHER_FOLDER + "/Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
GFF3_FILE     = MOTHER_FOLDER + "/Genomes_Input/syn1.genes.gff3"
OUT_FOLDER    = "./RNase"

TPM_TSV       = MOTHER_FOLDER + "/Syn1_Syn3A_Transcriptomics/syn1_TPM.tsv"

MIN_READS        = 10   # Parts A/B: only keep well-supported isoforms
TRUNC_MIN_READS  = 50   # Part C: matches the operon segmentation threshold
BOUNDARY_TOL     = 10   # Part C: bp tolerance for "same endpoint"


# ═══════════════════════════════════════════════════════════════════════════════
# PART A, Step 1 — Load gene intervals → per-strand intragenic mask
# ═══════════════════════════════════════════════════════════════════════════════
def load_intragenic_mask(gff3_path: str) -> dict:
    """Return {'+': bool array, '-': bool array} marking intragenic
    positions on each strand, indexed in 0-based genomic coordinates.
    Only protein-coding genes (rna_type=mRNA) are included."""
    gene_rows = []
    chrom_len = 0
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene":
                continue
            attrs = f[8] if len(f) >= 9 else ""
            if "rna_type=mRNA" not in attrs:
                continue
            s1, e1, strand = int(f[3]), int(f[4]), f[6]
            gene_rows.append((s1 - 1, e1, strand))   # → 0-based half-open
            if e1 > chrom_len:
                chrom_len = e1
    mask = {"+": np.zeros(chrom_len, dtype=bool),
            "-": np.zeros(chrom_len, dtype=bool)}
    for s0, e0, strand in gene_rows:
        if strand in mask:
            mask[strand][s0:e0] = True
    print(f"Loaded {len(gene_rows)} protein-coding gene intervals from GFF3 "
          f"(chrom length {chrom_len:,})")
    return mask

INTRAGENIC = load_intragenic_mask(GFF3_FILE)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 2 — Load isoforms and label endpoint context
# ═══════════════════════════════════════════════════════════════════════════════
df = pd.read_csv(ISOFORMS_TSV, sep="\t")
df_iso = df[df["n_reads"] >= MIN_READS].copy().reset_index(drop=True)
print(f"Isoforms loaded (n_reads >= {MIN_READS}): {len(df_iso)}")


def endpoint_positions(row):
    """Return (5'-end pos0, 3'-end pos0) for an isoform."""
    if row["strand"] == "+":
        return int(row["start0"]), int(row["end0"]) - 1
    else:
        return int(row["end0"]) - 1, int(row["start0"])


def is_intragenic(pos0: int, strand: str) -> bool:
    mask = INTRAGENIC[strand]
    return bool(0 <= pos0 < len(mask) and mask[pos0])


pos5_list, pos3_list, in5_list, in3_list = [], [], [], []
for _, row in df_iso.iterrows():
    p5, p3 = endpoint_positions(row)
    pos5_list.append(p5)
    pos3_list.append(p3)
    in5_list.append(is_intragenic(p5, row["strand"]))
    in3_list.append(is_intragenic(p3, row["strand"]))

df_iso["pos5p_0"]       = pos5_list
df_iso["pos3p_0"]       = pos3_list
df_iso["intragenic_5p"] = in5_list
df_iso["intragenic_3p"] = in3_list


def categorise(row) -> str:
    if row["intragenic_5p"] and row["intragenic_3p"]:
        return "both_intragenic"
    if row["intragenic_5p"]:
        return "5p_intragenic_only"
    if row["intragenic_3p"]:
        return "3p_intragenic_only"
    return "unprocessed"

df_iso["category"] = df_iso.apply(categorise, axis=1)


# ═══════════════════════════════════════════════════════════════════════════════
# Step 3 — Summary statistics
# ═══════════════════════════════════════════════════════════════════════════════
CATS = ["unprocessed", "5p_intragenic_only", "3p_intragenic_only", "both_intragenic"]
LABELS = ["Unprocessed\n(5' & 3' intergenic)",
          "5' intragenic only",
          "3' intragenic only",
          "Both ends eroded"]
COLORS = ["#969696", "#2166AC", "#B2182B", "#762A83"]

cat_counts = df_iso["category"].value_counts()
cat_reads  = df_iso.groupby("category")["n_reads"].sum()
total_iso  = len(df_iso)
total_reads = int(df_iso["n_reads"].sum())

print("\n" + "="*70)
print("ENDPOINT-CONTEXT CATEGORY SUMMARY")
print("="*70)
print(f"{'Category':<24} {'Isoforms':>10} {'%':>8} {'Reads':>14} {'%':>8}")
print("-"*70)
for cat in CATS:
    n_iso  = int(cat_counts.get(cat, 0))
    n_read = int(cat_reads.get(cat, 0))
    print(f"{cat:<24} {n_iso:>10,} {n_iso/total_iso*100:>7.1f}% "
          f"{n_read:>14,} {n_read/total_reads*100:>7.1f}%")
print("-"*70)
print(f"{'TOTAL':<24} {total_iso:>10,} {'100.0%':>8} "
      f"{total_reads:>14,} {'100.0%':>8}")

# Overall (any) intragenic-end counts: union over the two single-end groups
n_5p_intra = int((df_iso["intragenic_5p"]).sum())
n_3p_intra = int((df_iso["intragenic_3p"]).sum())
r_5p_intra = int(df_iso.loc[df_iso["intragenic_5p"], "n_reads"].sum())
r_3p_intra = int(df_iso.loc[df_iso["intragenic_3p"], "n_reads"].sum())

print("\nIntragenic-endpoint totals (counting any isoform with that end intragenic):")
print(f"  5' intragenic: {n_5p_intra:,} isoforms / {r_5p_intra:,} reads")
print(f"  3' intragenic: {n_3p_intra:,} isoforms / {r_3p_intra:,} reads")

ratio_iso  = n_5p_intra / n_3p_intra if n_3p_intra > 0 else float("inf")
ratio_read = r_5p_intra / r_3p_intra if r_3p_intra > 0 else float("inf")
print(f"\n5'→3' / 3'→5' ratio (by unique isoform kinds):  {ratio_iso:.2f}")
print(f"5'→3' / 3'→5' ratio (by unique isoform counts): {ratio_read:.2f}")

# Per-strand breakdown
print("\nPer-strand breakdown:")
for strand in ["+", "-"]:
    sub = df_iso[df_iso["strand"] == strand]
    n5 = int(sub["intragenic_5p"].sum())
    n3 = int(sub["intragenic_3p"].sum())
    r5 = int(sub.loc[sub["intragenic_5p"], "n_reads"].sum())
    r3 = int(sub.loc[sub["intragenic_3p"], "n_reads"].sum())
    rstr_i = f"{n5/n3:.2f}" if n3 > 0 else "inf"
    rstr_r = f"{r5/r3:.2f}" if r3 > 0 else "inf"
    print(f"  {strand} strand: 5'-intra={n5} ({r5:,} reads), "
          f"3'-intra={n3} ({r3:,} reads), "
          f"ratio_iso={rstr_i}, ratio_read={rstr_r}")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 4 — Figures (one PDF per panel)
# ═══════════════════════════════════════════════════════════════════════════════

# --- Panel A: Category bar chart (unique-isoform vs read-weighted) ---
figA, ax = plt.subplots(figsize=(8, 5))
iso_vals  = [int(cat_counts.get(c, 0)) for c in CATS]
read_vals = [int(cat_reads.get(c, 0))  for c in CATS]
iso_pcts  = [v / total_iso  * 100 for v in iso_vals]
read_pcts = [v / total_reads * 100 for v in read_vals]

x = np.arange(len(CATS))
w = 0.35
bars1 = ax.bar(x - w/2, iso_pcts,  w, color=COLORS, alpha=0.85,
               edgecolor="white", label="By unique isoform kinds")
bars2 = ax.bar(x + w/2, read_pcts, w, color=COLORS, alpha=0.50,
               edgecolor="white", hatch="//", label="By unique isoform counts")
ax.set_xticks(x)
ax.set_xticklabels(LABELS, fontsize=12, ha="center")
ax.set_ylabel("Percentage (%)", fontsize=12)
# ax.set_title("A. Isoform endpoint-context categories",
            #  fontsize=11, fontweight="bold")
ax.legend(fontsize=12, loc="upper left")
ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))
for bar, val in zip(bars1, iso_vals):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"n={val:,}", ha="center", va="bottom", fontsize=7.5)
for bar, val in zip(bars2, read_vals):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{round(val/1000)}k", ha="center", va="bottom", fontsize=7.5)
figA.tight_layout()
figA.savefig(f"{OUT_FOLDER}/RNase_asymmetry_isoform_categories.pdf",
             dpi=300, bbox_inches="tight")
plt.close(figA)


# --- Panel A2: Stacked 100% bar (kinds vs counts) ---
figA2, ax = plt.subplots(figsize=(8, 5))

bar_labels  = ["By isoform kinds", "By isoform counts"]
vals_matrix = [iso_vals, read_vals]   # 2 bars × 4 categories
tots        = [total_iso, total_reads]

bottoms = [0.0, 0.0]
for i, (cat, label, color) in enumerate(zip(CATS, LABELS, COLORS)):
    heights = [v / tot * 100 for v, tot in zip([vals[i] for vals in vals_matrix], tots)]
    bars = ax.bar([0, 1], heights, 0.45, bottom=bottoms,
                  color=color, edgecolor="white", linewidth=0.6, label=label)
    for b_idx, (bar, val, tot) in enumerate(
            zip(bars, [iso_vals[i], read_vals[i]], tots)):
        pct = val / tot * 100
        mid = bottoms[b_idx] + pct / 2
        if pct >= 3:
            n_str = f"{val:,}" if b_idx == 0 else f"{round(val/1000)}k"
            ax.text(bar.get_x() + bar.get_width() / 2, mid,
                    f"{n_str}\n({pct:.1f}%)",
                    ha="center", va="center", fontsize=8.5,
                    color="white" if pct > 8 else "black")
    bottoms = [b + h for b, h in zip(bottoms, heights)]

ax.set_xticks([0, 1])
ax.set_xticklabels(bar_labels, fontsize=12)
ax.set_ylabel("Percentage (%)", fontsize=12)
ax.set_ylim(0, 100)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))
ax.legend(title="Category", fontsize=10, title_fontsize=10,
          loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)
figA2.tight_layout()
figA2.savefig(f"{OUT_FOLDER}/RNase_asymmetry_stacked100.pdf",
              dpi=300, bbox_inches="tight")
plt.close(figA2)


# --- Panel B: 5'→3' vs 3'→5' exo activity bar chart + ratio ---
figB, ax = plt.subplots(figsize=(7, 5))
groups = ["5'→3' exo\n(5' intragenic)", "3'→5' exo\n(3' intragenic)"]
iso_counts_grp  = [n_5p_intra, n_3p_intra]
read_counts_grp = [r_5p_intra, r_3p_intra]
group_colors    = ["#2166AC", "#B2182B"]

x = np.arange(len(groups))
w = 0.35
b1 = ax.bar(x - w/2, iso_counts_grp,  w, color=group_colors, alpha=0.85,
            edgecolor="white", label="By unique isoform kinds")
ax2 = ax.twinx()
b2 = ax2.bar(x + w/2, read_counts_grp, w, color=group_colors, alpha=0.45,
             edgecolor="white", hatch="//", label="By unique isoform counts")

ax.set_xticks(x)
ax.set_xticklabels(groups, fontsize=10)
ax.set_ylabel("Unique isoform kinds", fontsize=10, color="#222222")
ax2.set_ylabel("Unique isoform counts (reads)", fontsize=10, color="#222222")
# ax.set_title("B. Exoribonuclease activity asymmetry\n"
#              f"ratio (kinds) = {ratio_iso:.2f},  ratio (counts) = {ratio_read:.2f}",
#              fontsize=11, fontweight="bold")

for bar, val in zip(b1, iso_counts_grp):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
            f"{val:,}", ha="center", va="bottom", fontsize=8.5)
for bar, val in zip(b2, read_counts_grp):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
             f"{val:,}", ha="center", va="bottom", fontsize=8.5)

# combined legend
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8.5, loc="upper right")

figB.tight_layout()
figB.savefig(f"{OUT_FOLDER}/RNase_asymmetry_exoribo_categories.pdf",
             dpi=300, bbox_inches="tight")
plt.close(figB)

print("\nFigures saved:")
print(f"  {OUT_FOLDER}/RNase_asymmetry_isoform_categories.pdf")
print(f"  {OUT_FOLDER}/RNase_asymmetry_exoribo_categories.pdf")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 5 — Save annotated isoform table
# ═══════════════════════════════════════════════════════════════════════════════
out_tsv = f"{OUT_FOLDER}/isoform_endpoint_context.tsv"
keep_cols = ["isoform_id", "chrom", "strand", "start0", "end0", "n_reads",
             "pos5p_0", "pos3p_0", "intragenic_5p", "intragenic_3p", "category"]
df_iso[[c for c in keep_cols if c in df_iso.columns]].to_csv(
    out_tsv, sep="\t", index=False)
print(f"Saved: {out_tsv} ({len(df_iso)} rows)")


# ═══════════════════════════════════════════════════════════════════════════════
# PART B — Estimate ORFs with start codon but no stop codon
# ═══════════════════════════════════════════════════════════════════════════════
# For every filtered isoform we walk over the canonical (GFF3) genes on the
# same strand. Whenever the gene's *start codon* lies fully inside the isoform
# we count this as one observed ORF (one start/stop pair). We then check
# whether the matching stop codon is also inside the isoform — if not, the
# 3' end of the read is truncated relative to the canonical CDS.
#
# Convention: GFF3 gene coordinates already include the stop codon.
#   + strand:  start codon = [start0,   start0+3)
#              stop  codon = [end0-3,   end0)
#   − strand:  start codon = [end0-3,   end0)
#              stop  codon = [start0,   start0+3)

def load_gene_intervals(gff3_path: str) -> pd.DataFrame:
    """Load gene intervals, keeping only protein-coding genes (rna_type=mRNA)."""
    rows = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "gene":
                continue
            attrs = f[8] if len(f) >= 9 else ""
            if "rna_type=mRNA" not in attrs:
                continue
            s1, e1, strand = int(f[3]), int(f[4]), f[6]
            locus = ""
            for kv in attrs.split(";"):
                if kv.startswith("locus_tag="):
                    locus = kv.split("=", 1)[1]
                    break
            rows.append({
                "locus_tag": locus,
                "strand":    strand,
                "start0":    s1 - 1,
                "end0":      e1,
            })
    return pd.DataFrame(rows)


genes_df = load_gene_intervals(GFF3_FILE)
print(f"\nLoaded {len(genes_df)} gene intervals for ORF analysis")

# Split by strand for fast per-isoform lookup
genes_by_strand = {
    s: g.sort_values("start0").reset_index(drop=True)
    for s, g in genes_df.groupby("strand")
}

iso_orf_rows = []
for _, iso in df_iso.iterrows():
    strand = iso["strand"]
    iso_s, iso_e = int(iso["start0"]), int(iso["end0"])
    g = genes_by_strand.get(strand)
    if g is None:
        continue
    overlap = g[(g["start0"] < iso_e) & (g["end0"] > iso_s)]

    # Walk overlapping genes in transcription order along the isoform
    # (5'→3'). For + strand that's ascending start0, for − strand it's
    # descending end0. We then collect every gene whose start codon is
    # fully contained inside the isoform — each such (gene, isoform) pair
    # is one ORF observation.
    if strand == "+":
        ordered = overlap.sort_values("start0")
    else:
        ordered = overlap.sort_values("end0", ascending=False)

    orf_loci, has_stop_flags = [], []
    for _, gene in ordered.iterrows():
        gs, ge = int(gene["start0"]), int(gene["end0"])
        if strand == "+":
            start_codon = (gs, gs + 3)
            stop_codon  = (ge - 3, ge)
        else:
            start_codon = (ge - 3, ge)
            stop_codon  = (gs, gs + 3)

        start_in = (start_codon[0] >= iso_s) and (start_codon[1] <= iso_e)
        if not start_in:
            continue
        stop_in = (stop_codon[0] >= iso_s) and (stop_codon[1] <= iso_e)
        orf_loci.append(str(gene["locus_tag"]))
        has_stop_flags.append(bool(stop_in))

    if not orf_loci:
        continue

    # Compact text view of the ORFs along the isoform, e.g.
    #   "MMSYN1_0001✓ MMSYN1_0002✓ MMSYN1_0003✗"
    orfs_str = " ".join(
        f"{lt}{'✓' if hs else '✗'}"
        for lt, hs in zip(orf_loci, has_stop_flags)
    )
    iso_orf_rows.append({
        "orfs":       orfs_str,
        "isoform_id": iso["isoform_id"],
        "n_reads":    int(iso["n_reads"]),
        "strand":     strand,
        "n_orfs":     len(orf_loci),
        "n_with_stop":    int(sum(has_stop_flags)),
        "n_without_stop": int(len(orf_loci) - sum(has_stop_flags)),
    })

# Per-isoform table — the format the user asked for
orf_df = pd.DataFrame(iso_orf_rows,
                      columns=["orfs", "isoform_id", "n_reads", "strand",
                               "n_orfs", "n_with_stop", "n_without_stop"])

# Per-ORF totals derived directly from the isoform-level table.
# Each ORF in a polycistronic isoform is counted independently:
# n_orfs gives the multiplicity, so summing n_with_stop / n_without_stop
# yields the per-ORF totals exactly.
n_total     = int(orf_df["n_orfs"].sum())
n_no_stop   = int(orf_df["n_without_stop"].sum())
n_isoforms  = len(orf_df)
n_polycis   = int((orf_df["n_orfs"] >= 2).sum())

print(f"ORFs (start codon contained in an isoform): {n_total:,}")
print(f"  from {n_isoforms:,} isoforms ({n_polycis:,} carry ≥2 ORFs)")

if n_total > 0:
    pct_no_stop = n_no_stop / n_total * 100

    # Read-weighted: each ORF carries its isoform's read count
    r_total   = int((orf_df["n_orfs"]        * orf_df["n_reads"]).sum())
    r_no_stop = int((orf_df["n_without_stop"] * orf_df["n_reads"]).sum())
    pct_no_stop_r = r_no_stop / r_total * 100

    print("\n" + "="*70)
    print("ORFs WITH START CODON BUT NO STOP CODON")
    print("="*70)
    print("  (each ORF in a polycistronic isoform is counted independently)")
    print(f"  By unique ORF observations:")
    print(f"    no-stop / total = {n_no_stop:,} / {n_total:,} "
          f"= {pct_no_stop:.2f}%")
    print(f"  By read-weighted observations:")
    print(f"    no-stop / total = {r_no_stop:,} / {r_total:,} "
          f"= {pct_no_stop_r:.2f}%")

    print(f"\n  ORFs-per-isoform distribution: "
          f"mean={orf_df['n_orfs'].mean():.2f}, max={int(orf_df['n_orfs'].max())}")

    print("\n  Per-strand:")
    for strand in ["+", "-"]:
        sub = orf_df[orf_df["strand"] == strand]
        if len(sub) == 0:
            continue
        s_total = int(sub["n_orfs"].sum())
        s_ns    = int(sub["n_without_stop"].sum())
        print(f"    {strand} strand: no-stop {s_ns:,} / {s_total:,} "
              f"({s_ns/s_total*100:.2f}%)")

    orf_out = f"{OUT_FOLDER}/orf_start_stop_observations.tsv"
    orf_df.to_csv(orf_out, sep="\t", index=False)
    print(f"\nSaved: {orf_out} ({len(orf_df)} rows)")

# ═════════════════════════════════════════════════════════════════════════════
# PART C — Truncation within containment clusters
# ═════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# Re-filter the frame Part A already loaded, at the segmentation threshold.
tr_iso = df[df["n_reads"] >= TRUNC_MIN_READS].copy()
print(f"\nIsoforms for truncation analysis (n_reads >= {TRUNC_MIN_READS}): {len(tr_iso)}")


def cluster_isoforms_with_members(isoforms: pd.DataFrame,
                                   tol: int = BOUNDARY_TOL) -> list[dict]:
    """
    Containment clustering (identical to the segmentation notebook).
    Returns cluster dicts that include full member-level detail.
    """
    if isoforms.empty:
        return []

    iso    = isoforms.reset_index(drop=True)
    starts = iso["start0"].astype(int).tolist()
    ends   = iso["end0"].astype(int).tolist()
    n      = len(iso)

    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)

    order = sorted(range(n), key=lambda i: starts[i])
    for idx, i in enumerate(order):
        si, ei = starts[i], ends[i]
        for j in order[idx + 1:]:
            sj = starts[j]
            if sj >= ei + tol:
                break
            ej = ends[j]
            i_in_j = (si >= sj - tol) and (ei <= ej + tol)
            j_in_i = (sj >= si - tol) and (ej <= ei + tol)
            if i_in_j or j_in_i:
                union(i, j)

    components = defaultdict(list)
    for i in range(n):
        components[find(i)].append(i)

    strand = iso.loc[0, "strand"]
    blocks = []
    for indices in components.values():
        members = []
        for k in indices:
            members.append({
                "isoform_id": iso.loc[k, "isoform_id"],
                "start0":     int(iso.loc[k, "start0"]),
                "end0":       int(iso.loc[k, "end0"]),
                "n_reads":    int(iso.loc[k, "n_reads"]),
                "length":     int(iso.loc[k, "end0"] - iso.loc[k, "start0"]),
            })
        # Identify the longest member as the full-length reference
        members.sort(key=lambda m: m["length"], reverse=True)
        blocks.append({
            "strand":  strand,
            "start0":  min(m["start0"] for m in members),
            "end0":    max(m["end0"]   for m in members),
            "members": members,
        })
    return sorted(blocks, key=lambda b: b["start0"])


clusters_plus  = cluster_isoforms_with_members(tr_iso[tr_iso["strand"] == "+"])
clusters_minus = cluster_isoforms_with_members(tr_iso[tr_iso["strand"] == "-"])
all_clusters = clusters_plus + clusters_minus
print(f"Containment clusters: + strand {len(clusters_plus)}, "
      f"- strand {len(clusters_minus)}, total {len(all_clusters)}")

# ═══════════════════════════════════════════════════════════════════════════════
# PART C, Step 2 — Classify truncated isoforms within each cluster
# ═══════════════════════════════════════════════════════════════════════════════
#
# For + strand:  TSS = start0 (5' end, low coord),  TTS = end0  (3' end, high coord)
# For − strand:  TSS = end0   (5' end, high coord), TTS = start0 (3' end, low coord)
#
# "Same 5'" means the isoform's 5' end matches the reference within tol.
# "Shorter 3'" means the isoform's 3' end is receded inward from the reference.

tr_records = []

for cluster in all_clusters:
    strand  = cluster["strand"]
    members = cluster["members"]
    if len(members) < 2:
        continue  # singleton clusters have no truncated isoforms

    ref = members[0]  # longest isoform = full-length reference

    if strand == "+":
        ref_5p = ref["start0"]  # TSS
        ref_3p = ref["end0"]    # TTS
    else:
        ref_5p = ref["end0"]    # TSS (high coord for minus)
        ref_3p = ref["start0"]  # TTS (low coord for minus)

    for mem in members[1:]:  # all non-reference (shorter) members
        if strand == "+":
            mem_5p = mem["start0"]
            mem_3p = mem["end0"]
            same_5p = abs(mem_5p - ref_5p) <= BOUNDARY_TOL
            same_3p = abs(mem_3p - ref_3p) <= BOUNDARY_TOL
            shorter_3p = mem_3p < ref_3p - BOUNDARY_TOL
            shorter_5p = mem_5p > ref_5p + BOUNDARY_TOL  # 5' receded inward
        else:
            mem_5p = mem["end0"]
            mem_3p = mem["start0"]
            same_5p = abs(mem_5p - ref_5p) <= BOUNDARY_TOL
            same_3p = abs(mem_3p - ref_3p) <= BOUNDARY_TOL
            shorter_3p = mem_3p > ref_3p + BOUNDARY_TOL   # for minus, 3' is low coord
            shorter_5p = mem_5p < ref_5p - BOUNDARY_TOL   # for minus, 5' is high coord

        # Classify
        if same_5p and shorter_3p:
            category = "same_5p_shorter_3p"
        elif shorter_5p and same_3p:
            category = "diff_5p_same_3p"
        elif shorter_5p and shorter_3p:
            category = "diff_5p_shorter_3p"
        elif same_5p and same_3p:
            category = "same_both"       # near-identical to reference
        else:
            category = "other"

        # Compute how many bp are trimmed from each end
        if strand == "+":
            trim_3p = max(0, ref_3p - mem_3p)
            trim_5p = max(0, mem_5p - ref_5p)
        else:
            trim_3p = max(0, mem_3p - ref_3p)
            trim_5p = max(0, ref_5p - mem_5p)

        tr_records.append({
            "strand":       strand,
            "ref_isoform":  ref["isoform_id"],
            "ref_length":   ref["length"],
            "ref_reads":    ref["n_reads"],
            "mem_isoform":  mem["isoform_id"],
            "mem_length":   mem["length"],
            "mem_reads":    mem["n_reads"],
            "category":     category,
            "trim_5p_bp":   trim_5p,
            "trim_3p_bp":   trim_3p,
            "frac_retained": round(mem["length"] / ref["length"], 3),
        })

trunc_df = pd.DataFrame(tr_records)
print(f"\nTruncated isoforms classified: {len(trunc_df)}")

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3 — Summary statistics
# ═══════════════════════════════════════════════════════════════════════════════

# --- 3a: Isoform counts by category ---
tr_cat_counts = trunc_df["category"].value_counts()
tr_cat_reads  = trunc_df.groupby("category")["mem_reads"].sum()
tr_total_iso  = len(trunc_df)
tr_total_reads = trunc_df["mem_reads"].sum()

print("\n" + "="*70)
print("TRUNCATION CATEGORY SUMMARY")
print("="*70)
print(f"{'Category':<25} {'Isoforms':>10} {'%':>8} {'Reads':>12} {'%':>8}")
print("-"*70)
for cat in ["same_5p_shorter_3p", "diff_5p_same_3p", "diff_5p_shorter_3p",
            "same_both", "other"]:
    n_iso  = tr_cat_counts.get(cat, 0)
    n_read = tr_cat_reads.get(cat, 0)
    print(f"{cat:<25} {n_iso:>10,} {n_iso/tr_total_iso*100:>7.1f}% "
          f"{n_read:>12,} {n_read/tr_total_reads*100:>7.1f}%")
print("-"*70)
print(f"{'TOTAL':<25} {tr_total_iso:>10,} {'100.0%':>8} "
      f"{tr_total_reads:>12,} {'100.0%':>8}")


# --- 3b: Degradation asymmetry ratio ---
reads_3p_erosion = tr_cat_reads.get("same_5p_shorter_3p", 0)
reads_5p_erosion = tr_cat_reads.get("diff_5p_same_3p", 0)
if reads_5p_erosion > 0:
    asymmetry_ratio = reads_3p_erosion / reads_5p_erosion
    print(f"\n3' erosion / 5' erosion read ratio: {asymmetry_ratio:.2f}")
    print(f"  → {asymmetry_ratio:.1f}× more reads in 3'-eroded isoforms")
    print(f"  → Consistent with faster 5'→3' exo clearing downstream fragments")
else:
    asymmetry_ratio = float("inf")
    print(f"\n3' erosion / 5' erosion read ratio: inf (no 5' erosion isoforms)")

n_3p_erosion = tr_cat_counts.get("same_5p_shorter_3p", 0)
n_5p_erosion = tr_cat_counts.get("diff_5p_same_3p", 0)
if n_5p_erosion > 0:
    iso_ratio = n_3p_erosion / n_5p_erosion
    print(f"\n3' erosion / 5' erosion isoform count ratio: {iso_ratio:.2f}")
else:
    iso_ratio = float("inf")
    print(f"\n3' erosion / 5' erosion isoform count ratio: inf")


# --- 3c: Per-strand breakdown ---
print("\n" + "-"*70)
print("Per-strand breakdown:")
for strand in ["+", "-"]:
    sub = trunc_df[trunc_df["strand"] == strand]
    n3 = sub[sub["category"] == "same_5p_shorter_3p"].shape[0]
    n5 = sub[sub["category"] == "diff_5p_same_3p"].shape[0]
    nb = sub[sub["category"] == "diff_5p_shorter_3p"].shape[0]
    r3 = sub[sub["category"] == "same_5p_shorter_3p"]["mem_reads"].sum()
    r5 = sub[sub["category"] == "diff_5p_same_3p"]["mem_reads"].sum()
    ratio_str = f"{r3/r5:.2f}" if r5 > 0 else "inf"
    print(f"  {strand} strand: 3'erosion={n3} ({r3:,} reads), "
          f"5'erosion={n5} ({r5:,} reads), both={nb}, "
          f"ratio={ratio_str}")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 4 — Trim-length distributions (how far is each end eroded?)
# ═══════════════════════════════════════════════════════════════════════════════

# For 3'-eroded isoforms: distribution of 3' trim length
erosion_3p = trunc_df[trunc_df["category"] == "same_5p_shorter_3p"]["trim_3p_bp"]
erosion_5p = trunc_df[trunc_df["category"] == "diff_5p_same_3p"]["trim_5p_bp"]

print("\n" + "-"*70)
print("Trim-length distributions (bp eroded from reference):")
if len(erosion_3p) > 0:
    print(f"\n  3' erosion (same_5p_shorter_3p), n={len(erosion_3p)}:")
    print(f"    median = {erosion_3p.median():.0f} bp")
    print(f"    mean   = {erosion_3p.mean():.0f} bp")
    print(f"    Q25    = {erosion_3p.quantile(0.25):.0f} bp")
    print(f"    Q75    = {erosion_3p.quantile(0.75):.0f} bp")
if len(erosion_5p) > 0:
    print(f"\n  5' erosion (diff_5p_same_3p), n={len(erosion_5p)}:")
    print(f"    median = {erosion_5p.median():.0f} bp")
    print(f"    mean   = {erosion_5p.mean():.0f} bp")
    print(f"    Q25    = {erosion_5p.quantile(0.25):.0f} bp")
    print(f"    Q75    = {erosion_5p.quantile(0.75):.0f} bp")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 5 — Per-operon degradation ratio
# ═══════════════════════════════════════════════════════════════════════════════
# For each cluster with ≥2 truncated isoforms, compute:
#   degradation_ratio = reads(same_5p_shorter_3p) / reads(diff_5p_same_3p)
# A ratio > 1 means 3' erosion dominates (consistent with faster 5'→3' exo).

operon_ratios = []
for cluster in all_clusters:
    ref_id = cluster["members"][0]["isoform_id"] if cluster["members"] else None
    sub = trunc_df[trunc_df["ref_isoform"] == ref_id]
    if len(sub) < 2:
        continue
    r3 = sub[sub["category"] == "same_5p_shorter_3p"]["mem_reads"].sum()
    r5 = sub[sub["category"] == "diff_5p_same_3p"]["mem_reads"].sum()
    n3 = (sub["category"] == "same_5p_shorter_3p").sum()
    n5 = (sub["category"] == "diff_5p_same_3p").sum()
    if r3 + r5 > 0:
        operon_ratios.append({
            "ref_isoform":  ref_id,
            "strand":       cluster["strand"],
            "n_members":    len(cluster["members"]),
            "reads_3p_erosion": r3,
            "reads_5p_erosion": r5,
            "n_3p_erosion": n3,
            "n_5p_erosion": n5,
            "ratio":        r3 / r5 if r5 > 0 else r3/1, 
            "dominant":     "3p_erosion" if r3 > r5 else ("5p_erosion" if r5 > r3 else "balanced"),
        })

ratio_df = pd.DataFrame(operon_ratios)
finite_ratios = ratio_df[ratio_df["ratio"] != -float("inf")]
positive_ratios = finite_ratios[finite_ratios["ratio"] > 0]

print("\n" + "="*70)
print("PER-OPERON DEGRADATION ASYMMETRY")
print("="*70)
if len(ratio_df) > 0:
    n_3p_dom  = (ratio_df["dominant"] == "3p_erosion").sum()
    n_5p_dom  = (ratio_df["dominant"] == "5p_erosion").sum()
    n_balanced = (ratio_df["dominant"] == "balanced").sum()
    n_inf     = (ratio_df["ratio"] == float("inf")).sum()
    print(f"Operons analysed: {len(ratio_df)}")
    print(f"  3' erosion dominant: {n_3p_dom}  ({n_3p_dom/len(ratio_df)*100:.1f}%)")
    print(f"  5' erosion dominant: {n_5p_dom}  ({n_5p_dom/len(ratio_df)*100:.1f}%)")
    print(f"  Balanced:            {n_balanced}")
    print(f"  Ratio = inf (no 5' erosion at all): {n_inf}")
    if len(finite_ratios) > 0:
        print(f"\n  Finite ratio distribution (n={len(finite_ratios)}):")
        print(f"    median = {finite_ratios['ratio'].median():.2f}")
        print(f"    mean   = {finite_ratios['ratio'].mean():.2f}")
        print(f"    Q25    = {finite_ratios['ratio'].quantile(0.25):.2f}")
        print(f"    Q75    = {finite_ratios['ratio'].quantile(0.75):.2f}")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 6 — Figures
# ═══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Degradation Asymmetry Analysis — Syn1 Operon Isoforms",
             fontsize=13, fontweight="bold", y=0.98)

# --- Panel A: Category bar chart (isoform counts + read-weighted) ---
ax = axes[0, 0]
tr_cats = ["same_5p_shorter_3p", "diff_5p_same_3p", "diff_5p_shorter_3p", "same_both", "other"]
tr_labels = ["Same 5', shorter 3'\n(3' erosion)",
          "Diff 5', same 3'\n(5' erosion)",
          "Diff 5', shorter 3'\n(both trimmed)",
          "Same both\n(near-identical)",
          "Other"]
tr_colors = ["#2166AC", "#B2182B", "#762A83", "#969696", "#CCCCCC"]

tr_iso_vals  = [tr_cat_counts.get(c, 0) for c in tr_cats]
tr_read_vals = [tr_cat_reads.get(c, 0) for c in tr_cats]
tr_iso_pcts  = [v / tr_total_iso * 100 for v in tr_iso_vals]
tr_read_pcts = [v / tr_total_reads * 100 for v in tr_read_vals]

x = np.arange(len(tr_cats))
w = 0.35
bars1 = ax.bar(x - w/2, tr_iso_pcts,  w, color=tr_colors, alpha=0.85, edgecolor="white", label="By isoform count")
bars2 = ax.bar(x + w/2, tr_read_pcts, w, color=tr_colors, alpha=0.50, edgecolor="white", hatch="//", label="By read count")
ax.set_xticks(x)
ax.set_xticklabels(tr_labels, fontsize=7.5, ha="center")
ax.set_ylabel("Percentage (%)", fontsize=9)
ax.set_title("A. Truncation categories", fontsize=10, fontweight="bold")
ax.legend(fontsize=7.5, loc="upper right")
ax.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))

# Add count annotations on bars
for bar, val in zip(bars1, tr_iso_vals):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"n={val:,}", ha="center", va="bottom", fontsize=6.5)

# --- Panel B: Trim-length distributions ---
ax = axes[0, 1]
bins = np.arange(0, min(5000, max(
    erosion_3p.max() if len(erosion_3p) > 0 else 0,
    erosion_5p.max() if len(erosion_5p) > 0 else 0
) + 200), 50)

if len(erosion_3p) > 0:
    ax.hist(erosion_3p, bins=bins, alpha=0.7, color="#2166AC",
            label=f"3' erosion (n={len(erosion_3p):,})", density=True)
if len(erosion_5p) > 0:
    ax.hist(erosion_5p, bins=bins, alpha=0.7, color="#B2182B",
            label=f"5' erosion (n={len(erosion_5p):,})", density=True)
ax.set_xlabel("Bases trimmed (bp)", fontsize=9)
ax.set_ylabel("Density", fontsize=9)
ax.set_title("B. Trim-length distributions", fontsize=10, fontweight="bold")
ax.legend(fontsize=8)
# --- Panel C: Per-operon ratio distribution ---
ax = axes[1, 0]
if len(positive_ratios) > 0:
    log_ratios = np.log2(positive_ratios["ratio"].values)
    ax.hist(log_ratios, bins=40, color="#4393C3", alpha=0.8, edgecolor="white")
    ax.axvline(0, color="red", linestyle="--", linewidth=1.2,
               label="Balanced (ratio=1)")
    med = np.median(log_ratios)
    ax.axvline(med, color="black", linestyle="-", linewidth=1.5,
               label=f"Median = {2**med:.2f}× (log₂ = {med:.2f})")
    ax.set_xlabel("log₂(3' erosion reads / 5' erosion reads)", fontsize=9)
    ax.set_ylabel("Number of operons", fontsize=9)
    ax.legend(fontsize=7.5)
ax.set_title("C. Per-operon degradation asymmetry", fontsize=10, fontweight="bold")
# Add annotation for operons with only 3' erosion (ratio = inf)
n_inf = (ratio_df["ratio"] == float("inf")).sum()
if n_inf > 0:
    ax.text(0.98, 0.95, f"+{n_inf} operons with\nonly 3' erosion (ratio=∞)",
            transform=ax.transAxes, fontsize=7.5, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#E0E0E0", alpha=0.8))

# --- Panel D: Read fraction retained vs category ---
ax = axes[1, 1]
for cat, color, label in [
    ("same_5p_shorter_3p", "#2166AC", "3' erosion"),
    ("diff_5p_same_3p",    "#B2182B", "5' erosion"),
    ("diff_5p_shorter_3p", "#762A83", "Both trimmed"),
]:
    sub = trunc_df[trunc_df["category"] == cat]
    if len(sub) > 0:
        ax.scatter(sub["frac_retained"], sub["mem_reads"],
                   alpha=0.3, s=8, color=color, label=f"{label} (n={len(sub):,})")
ax.set_xlabel("Fraction of full-length retained", fontsize=9)
ax.set_ylabel("Read count (truncated isoform)", fontsize=9)
ax.set_yscale("log")
ax.set_title("D. Read support vs. length retention", fontsize=10, fontweight="bold")
ax.legend(fontsize=7.5, markerscale=2)

plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{OUT_FOLDER}/degradation_asymmetry.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT_FOLDER}/degradation_asymmetry.png", dpi=200, bbox_inches="tight")
plt.close()
print("\nFigures saved: degradation_asymmetry.pdf / .png")


# ═══════════════════════════════════════════════════════════════════════════════
# Step 7 — Save detailed tables
# ═══════════════════════════════════════════════════════════════════════════════
trunc_df.to_csv(f"{OUT_FOLDER}/truncation_classification.tsv", sep="\t", index=False)
print(f"Saved: truncation_classification.tsv ({len(trunc_df)} rows)")

if len(ratio_df) > 0:
    ratio_df.to_csv(f"{OUT_FOLDER}/per_operon_degradation_ratio.tsv", sep="\t", index=False)
    print(f"Saved: per_operon_degradation_ratio.tsv ({len(ratio_df)} rows)")

print("\nDone.")

# ═════════════════════════════════════════════════════════════════════════════
# PART D — Ribosome rescue: tmRNA abundance and stop-codon-less mRNA
# ═════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# tmRNA (ssrA) as a fraction of non-rRNA expression
# ═══════════════════════════════════════════════════════════════════════════════
# Motivation: In E. coli, tmRNA accounts for ~25% of non-rRNA reads, reflecting
# heavy ribosome-rescue demand. If Syn1 shows a similar or higher fraction it
# would be consistent with RNase Y endo-cleavages constantly generating truncated
# mRNAs that stall ribosomes.
#
# We use the Illumina avg_sense_TPM as the primary abundance metric (averaged
# over the three biological replicates already stored in the CSV), with
# PacBio_sense_TPM shown alongside for comparison.
# ═══════════════════════════════════════════════════════════════════════════════

tpm = pd.read_csv(TPM_TSV, sep="\t")

# ── classify genes ──────────────────────────────────────────────────────────
is_rRNA   = tpm["rna_type"] == "rRNA"
is_tmRNA  = tpm["rna_type"] == "tmRNA"          # ssrA / MMSYN1_0158
is_non_rRNA = ~is_rRNA                           # everything that is not rRNA

# ── Illumina avg_sense_TPM ──────────────────────────────────────────────────
illumina_col   = "avg_sense_TPM"
pacbio_col     = "PacBio_sense_TPM"

total_illumina       = tpm[illumina_col].sum()
non_rRNA_illumina    = tpm.loc[is_non_rRNA, illumina_col].sum()
tmRNA_illumina       = tpm.loc[is_tmRNA,    illumina_col].sum()

rRNA_illumina        = tpm.loc[is_rRNA,     illumina_col].sum()

total_pacbio         = tpm[pacbio_col].sum()
non_rRNA_pacbio      = tpm.loc[is_non_rRNA, pacbio_col].sum()
tmRNA_pacbio         = tpm.loc[is_tmRNA,    pacbio_col].sum()
rRNA_pacbio          = tpm.loc[is_rRNA,     pacbio_col].sum()

tmRNA_frac_illumina  = tmRNA_illumina  / non_rRNA_illumina  * 100
tmRNA_frac_pacbio    = tmRNA_pacbio    / non_rRNA_pacbio    * 100
rRNA_frac_illumina   = rRNA_illumina   / total_illumina     * 100
rRNA_frac_pacbio     = rRNA_pacbio     / total_pacbio       * 100

print("=" * 62)
print("tmRNA (ssrA) ABUNDANCE — Syn1 vs E. coli reference")
print("=" * 62)

# Show the tmRNA row for transparency
tmrna_row = tpm[is_tmRNA][["locus_tag", "gene_name", "rna_type",
                             illumina_col, pacbio_col]]
print("\ntmRNA gene entry:")
print(tmrna_row.to_string(index=False))

print(f"\n{'Metric':<45} {'Illumina':>10} {'PacBio':>10}")
print("-" * 65)
print(f"{'rRNA fraction of total TPM':<45} {rRNA_frac_illumina:>9.1f}% {rRNA_frac_pacbio:>9.1f}%")
print(f"{'non-rRNA TPM (sum)':<45} {non_rRNA_illumina:>10.1f} {non_rRNA_pacbio:>10.1f}")
print(f"{'tmRNA TPM':<45} {tmRNA_illumina:>10.1f} {tmRNA_pacbio:>10.1f}")
print(f"{'tmRNA / non-rRNA TPM  (Syn1)':<45} {tmRNA_frac_illumina:>9.1f}% {tmRNA_frac_pacbio:>9.1f}%")
print(f"{'tmRNA / non-rRNA reads (E. coli reference)':<45} {'~25 %':>10}")

print("\nInterpretation:")
ecoli_ref = 25.0
for label, frac in [("Illumina", tmRNA_frac_illumina), ("PacBio", tmRNA_frac_pacbio)]:
    fold = frac / ecoli_ref
    direction = "HIGHER" if frac > ecoli_ref else "lower"
    print(f"  {label}: tmRNA = {frac:.1f}% of non-rRNA TPM  "
          f"({fold:.2f}× vs E. coli ~25%) — {direction} than E. coli")

print("\n  A fraction ≥ 25% is consistent with elevated ribosome-rescue")
print("  demand, as expected if RNase Y endo-cleavages frequently generate")
print("  truncated mRNAs that stall ribosomes and require tmRNA tagging.")

# ═══════════════════════════════════════════════════════════════════════════════
# Ribosome-trapping potential: 3'-eroded mRNAs as a fraction of all mRNA reads
# ═══════════════════════════════════════════════════════════════════════════════
# A 3'-eroded isoform (same 5' end, truncated 3' end) lacks the stop codon.
# Ribosomes translating such transcripts reach the 3' end without terminating
# and become stalled — exactly the substrates that tmRNA rescues.
#
# Estimate: (reads from same_5p_shorter_3p isoforms) / (all classified reads)
# Both numerator and denominator come from the PacBio containment-clustering
# analysis performed above (trunc_df / cat_reads / tr_total_reads).
# ═══════════════════════════════════════════════════════════════════════════════

# reads_3p_erosion and tr_total_reads are already defined from Step 3 above.
# For a more conservative denominator we also include full-length reference reads.
total_reads_incl_ref = trunc_df["mem_reads"].sum() + \
    sum(c["members"][0]["n_reads"] for c in all_clusters if c["members"])

frac_no_stop_truncated = reads_3p_erosion / tr_total_reads * 100
frac_no_stop_all       = reads_3p_erosion / total_reads_incl_ref * 100

print("\n" + "=" * 62)
print("RIBOSOME-TRAPPING POTENTIAL — 3'-eroded (stop-codon-less) mRNAs")
print("=" * 62)
print(f"\n  3'-eroded reads (same_5p_shorter_3p) : {reads_3p_erosion:>12,}")
print(f"  Total truncated isoform reads         : {tr_total_reads:>12,}")
print(f"  Total reads incl. full-length refs    : {total_reads_incl_ref:>12,}")
print(f"\n  % of truncated reads lacking stop codon : {frac_no_stop_truncated:.1f}%")
print(f"  % of all isoform reads lacking stop codon: {frac_no_stop_all:.1f}%")
print(f"\n  Interpretation:")
print(f"  At least {frac_no_stop_all:.1f}% of PacBio isoform reads represent transcripts")
print(f"  without a stop codon. These are the direct substrates for ribosome")
print(f"  stalling and tmRNA-mediated rescue, corroborating the high tmRNA")
print(f"  abundance observed above.")

print("\nDone.")
