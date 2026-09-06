"""
Combined syn1 vs syn3A comparison at RNA (TPM) and protein (iPM) levels.

Background
----------
syn3A is a reduced descendant of syn1: ~430 genes were deleted, so syn3A is
annotated with ~490 genes vs syn1's ~918. TPM and iPM are *relative* abundance
units (each sums to ~1e6 within an organism). When the gene complement shrinks,
the same fixed pool is divided among fewer genes, so every *retained* gene's
relative value is mechanically inflated in syn3A. A naive fold change
(syn3A / syn1) therefore has a baseline well above 1 even for genes whose true
abundance did not change.

Normalization (FC baseline correction)
--------------------------------------
We express each gene relative to the *per-gene mean* of its own organism and
subset (coding vs non-coding):

    rel_i = value_i / mean(value over genes DETECTED in that subset)

This makes the average gene equal to 1 in both organisms, so an unchanged gene
gives FC ~= 1 regardless of the gene-count difference.

Platform policy for TPM
-----------------------
ONT direct-RNA is poor for quantification, so the *standard* TPM comparison and
the combined CSV use Illumina for both organisms:
    TPM_fold_change = relTPM(Illumina syn3A) / relTPM(Illumina syn1)
ONT is kept only for QC correlation scatters:
    - ONT syn3A vs Illumina syn1   (cross-organism, cross-platform)
    - ONT syn3A vs Illumina syn3A  (within-organism platform agreement)
Exception: rRNA/tRNA (non-coding) are absent from the Ribo-Zero Illumina library
(syn3A Illumina detects 0 of them), so the NON-coding CSV falls back to ONT for
the syn3A side — it is the only quantification available there.

Two complementary change metrics
--------------------------------
  - fold change       : rel_syn3A / rel_syn1            (relative; baseline ~1)
  - absolute change   : rel_syn3A - rel_syn1            (robust to tiny denoms)

Deleted-gene occupancy
----------------------
Quantifies how much of syn1's transcriptome (total + mRNA-only) and proteome was
carried by genes that syn3A deleted (raw syn1 TPM / iPM shares), and classifies
the deleted loci by RNA type. Deleted = syn1 loci absent from syn3a_genome.gff3.

Function-category TPM change
----------------------------
Generalizes the ribosomal-protein observation (r-proteins occupy more of the
mRNA pool in syn3A) to every curated function group. Coding genes are joined to
the Secondary/Tertiary function from syn3A_proteome.tsv (by locus_num)
and per category we report mRNA-pool share (syn1 vs syn3A) + its change, median
TPM fold/absolute change, and Mann-Whitney U of log10 TPM FC vs the rest.

Outputs (all under Compare_RNA_Protein/)
----------------------------------------
  - distribution + correlation PDFs (standard TPM = Illumina; iPM; 2 ONT QC scatters)
  - operon summary txt for extreme-FC groups
  - deleted_gene_occupancy.txt
  - Compare_RNA_Protein.txt   : occupancy report + FC / abs-change outliers
  - syn1_vs_syn3a_RNA_protein.tsv      : coding (mRNA + pseudo) + protein layer
  - syn1_vs_syn3a_noncoding_RNA.tsv    : rRNA / tRNA / ncRNA / tmRNA (RNA only)
  - TPM_change_by_{secondary,tertiary}.tsv   : per-category pool share + FC stats
                                               (retained-pool / deletion-corrected)
  - TPM_FC_by_{secondary,tertiary}_function.pdf       : log10 TPM-FC boxplots
  - mRNA_pool_share_change_by_{secondary,tertiary}.pdf: share-change diverging bars
  - mRNA_pool_composition_by_secondary.pdf   : tall (x:2x) two full-pool stacked bars
        (syn1 w/ hatched deleted block | syn3A); blocks = secondary (Unclear & Cellular
        collapsed to Primary; 'other'/unannotated dropped & renormalized), Roman-indexed,
        colored by Primary; no spines/title/y-axis; legend below.
  - tertiary_share_change_dumbbell.pdf       : tall (x:2x) retained-pool deletion-corrected
        tertiary dumbbell (|delta|>0.5pp); y-labels = Secondary over Tertiary; y-grid, no
        left/right spines; per-dot share values. (Designed to sit beside the bars in Illustrator.)
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from adjustText import adjust_text
from scipy import stats

# Nature print-ready rcParams (per OUTPUT.md). Used as a local rc_context wrapper
# around figures we've explicitly retrofitted to Nature spec; the other figures
# in this script keep their original styling (forward-only retrofit).
NATURE_RC = {
    "font.size": 7,
    "font.family": "sans-serif",
    # Arial first (Mac/Win); Liberation Sans (Arial-metric-clone) is the Linux
    # fallback. With pdf.fonttype=42 the PDF embeds whichever font was actually
    # used — install ttf-mscorefonts-installer to get true Arial on Linux.
    "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


# ── Inputs ───────────────────────────────────────────────────────────────────
SYN1_TPM_CSV          = "../Syn1_Syn3A_Transcriptomics/syn1_TPM.tsv"
SYN3A_TPM_TSV         = "../Syn1_Syn3A_Transcriptomics/syn3A_TPM.tsv"
SYN3A_ILLUMINA_TPM_CSV = "../Syn1_Syn3A_Transcriptomics/Processed_TPM_Palsson/GSM6204176_3A.csv"
SYN1_PTN_CSV   = "../Syn1_Syn3A_Proteomics/syn1_proteome.tsv"
SYN3A_PTN_CSV  = "../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv"
OPERON_COV     = "../Syn1_Operon/segmentation/gene_operon_coverage.tsv"
SYN3A_GFF      = "../Genomes_Input/syn3a_genome.gff3"
IMPACT_CTX     = "delete_gene/retained_gene_context.tsv"  # gene_impact_class (08)
SYN3A_FUNC_TSV = SYN3A_PTN_CSV   # one table now: abundance + curated function

# The syn3A proteome table is snake_case; this script uses the curated Title-Case
# names throughout, so rename at the door and leave the rest of the script alone.
SYN3A_FUNC_RENAME = {"locus_tag": "Locus Tag", "gene_name": "Gene Name",
                     "gene_product": "Gene Product", "essentiality": "Essentiality",
                     "localization": "Localization",
                     "primary_function": "Primary Function",
                     "secondary_function": "Secondary Function",
                     "tertiary_function": "Tertiary Function",
                     "exp_ptn_cnt_2019": "Exp. Ptn Cnt 2019",
                     "exp_ptn_cnt_2026": "Exp. Ptn Cnt 2026",
                     "sim_initial_ptn_cnt": "Sim. Initial Ptn Cnt",
                     "mRNA_copies_per_cell": "mRNA Copies/Cell",
                     "protein_length": "Protein Length",
                     "protein_sequence": "Protein Sequence"}


def load_syn3a_function():
    """Curated syn3A proteome, with the Title-Case column names this script expects."""
    return pd.read_csv(SYN3A_FUNC_TSV, sep="\t").rename(columns=SYN3A_FUNC_RENAME)


# gene_impact_class ordered by predicted transcriptional effect: down -> neutral -> up.
CLASS_ORDER = ["promoter_lost", "promoter_disconnected", "promoter_proximity_changed",
               "context_only", "unaffected", "readthrough_exposed", "new_promoter_fusion"]

SYN1_TPM_COL   = "avg_sense_TPM"   # Illumina merged across bio samples
SYN3A_TPM_COL  = "ONT_sense_TPM"
SYN3A_ILLUMINA_TPM_COL = "Illumina_TPM"

# ── Outputs ──────────────────────────────────────────────────────────────────
OUTDIR         = "Compare_RNA_Protein"
OUT_CODING     = f"{OUTDIR}/syn1_vs_syn3a_RNA_protein.tsv"
OUT_NONCODING  = f"{OUTDIR}/syn1_vs_syn3a_noncoding_RNA.tsv"
OUT_REPORT     = f"{OUTDIR}/Compare_RNA_Protein.txt"
os.makedirs(OUTDIR, exist_ok=True)

CODING_RNA_TYPES    = {"mRNA", "pseudo"}
NONCODING_RNA_TYPES = {"rRNA", "tRNA", "ncRNA", "tmRNA"}

# Scatter-highlight cutoff in mean-normalized units: rel > ABUNDANT_THRESH means
# the gene is at least as abundant as the average syn1 gene.
ABUNDANT_THRESH = 1.0

# ── Operon coverage (used by both layers) ────────────────────────────────────
operon_cov = pd.read_csv(OPERON_COV, sep="\t")
operon_cov["locus_num"] = operon_cov["locus_tag"].str.extract(r"(\d+)$").astype(int)


def _extract_locus_num(s: pd.Series) -> pd.Series:
    return s.str.extract(r"(\d+)$").astype(int)


def _load_gff_locus_nums(path: str) -> set:
    """Return the set of locus_num integers for gene + pseudogene records."""
    nums = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = re.search(r"locus_tag=([^;]+)", f[8])
            if not m:
                continue
            mm = re.search(r"(\d+)$", m.group(1))
            if mm:
                nums.add(int(mm.group(1)))
    return nums


def _fc_plots(result: pd.DataFrame, value_col: str, fc_col: str, log10fc_col: str,
              tag: str, ylabel_prefix: str, abundant_thresh: float = ABUNDANT_THRESH) -> None:
    """FC + log10FC distribution histograms and the syn1-vs-syn3a scatter for a
    syn3A/syn1 comparison. Values are mean-normalized, so FC~1 is the baseline."""
    fc = result[fc_col]
    fc_stats = fc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"\n[{tag}] fold change statistics (n={int(fc_stats['count'])}):")
    print(f"  mean={fc_stats['mean']:.3f}  median={fc_stats['50%']:.3f}  std={fc_stats['std']:.3f}")
    print(f"  5th={fc_stats['5%']:.3f}  25th={fc_stats['25%']:.3f}  "
          f"75th={fc_stats['75%']:.3f}  95th={fc_stats['95%']:.3f}")

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(1,   color="black",  linestyle="-",  linewidth=1.0, label="FC = 1 (baseline)")
    ax.axvline(10,  color="red",    linestyle="--", linewidth=1.2, label="FC = 10")
    ax.axvline(0.1, color="orange", linestyle="--", linewidth=1.2, label="FC = 0.1")
    ax.set_xlabel(f"{tag} fold change  (syn3A / syn1)", fontsize=11)
    ax.set_ylabel(f"Number of {ylabel_prefix}", fontsize=11)
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{tag}_fold_change_distribution.pdf")
    plt.close(fig)

    log10fc = result[log10fc_col]
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(log10fc, bins=50, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(0,             color="black",  linestyle="-",  linewidth=1.0, label="log10FC = 0")
    ax.axvline(np.log10(10),  color="red",    linestyle="--", linewidth=1.2, label="log10(10) = 1")
    ax.axvline(np.log10(0.1), color="orange", linestyle="--", linewidth=1.2, label="log10(0.1) = -1")
    ax.set_xlabel(f"log10({tag} fold change)  (syn3A / syn1)", fontsize=12)
    ax.set_ylabel(f"Number of {ylabel_prefix}", fontsize=11)
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{tag}_log10FC_distribution.pdf")
    plt.close(fig)

    # Restrict to rows where both syn1 and syn3a have positive, non-NaN values.
    plot_df = result[
        result[f"{value_col}_syn1"].notna() & (result[f"{value_col}_syn1"] > 0) &
        result[f"{value_col}_syn3a"].notna() & (result[f"{value_col}_syn3a"] > 0)
    ]
    x = plot_df[f"{value_col}_syn1"]
    y = plot_df[f"{value_col}_syn3a"]
    pearson_r, pearson_p   = stats.pearsonr(np.log10(x), np.log10(y))
    spearman_r, spearman_p = stats.spearmanr(x, y)

    abundant  = plot_df[f"{value_col}_syn1"] > abundant_thresh
    mask_high = (plot_df[fc_col] > 10)  & abundant
    mask_low  = (plot_df[fc_col] < 0.1) & abundant
    mask_mid  = ~mask_high & ~mask_low

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x[mask_mid],  y[mask_mid],  s=12, alpha=0.55, color="steelblue", edgecolors="none", label="other")
    ax.scatter(x[mask_high], y[mask_high], s=25, alpha=0.90, color="red",       edgecolors="none",
               label=f"FC > 10  & syn1 rel-{tag} > {abundant_thresh:g}")
    ax.scatter(x[mask_low],  y[mask_low],  s=25, alpha=0.90, color="orange",    edgecolors="none",
               label=f"FC < 0.1 & syn1 rel-{tag} > {abundant_thresh:g}")


    lims = [min(x.min(), y.min()) * 0.5, max(x.max(), y.max()) * 2]
    ax.plot(lims, lims,                    color="black",  linewidth=0.8, linestyle="--")
    ax.plot(lims, [v * 10  for v in lims], color="red",    linewidth=0.6, linestyle=":")
    ax.plot(lims, [v * 0.1 for v in lims], color="orange", linewidth=0.6, linestyle=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"syn1  rel-{tag} (log)",  fontsize=12)
    ax.set_ylabel(f"syn3A rel-{tag} (log)",  fontsize=12)
    ax.text(0.04, 0.96,
            f"Pearson r = {pearson_r:.3f}  (log10)\nSpearman r = {spearman_r:.3f}\nn = {len(plot_df)}",
            transform=ax.transAxes, fontsize=12, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
    ax.legend(fontsize=10, loc="lower right")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{tag}_correlation_syn1_vs_syn3a.pdf")
    plt.close(fig)
    print(f"[{tag}] Pearson r (log10) = {pearson_r:.3f} (p={pearson_p:.2e}); "
          f"Spearman r = {spearman_r:.3f} (p={spearman_p:.2e})")


def _corr_scatter(merged: pd.DataFrame, xcol: str, ycol: str,
                  xlabel: str, ylabel: str, outname: str, title: str) -> None:
    """Pure log-log correlation scatter for an arbitrary x vs y (no FC framing).
    Used for the ONT QC comparisons."""
    sub = merged[merged[xcol].notna() & (merged[xcol] > 0) &
                 merged[ycol].notna() & (merged[ycol] > 0)]
    x, y = sub[xcol], sub[ycol]
    pr, pp = stats.pearsonr(np.log10(x), np.log10(y))
    sr, sp = stats.spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(x, y, s=14, alpha=0.55, color="steelblue", edgecolors="none")
    lims = [min(x.min(), y.min()) * 0.5, max(x.max(), y.max()) * 2]
    ax.plot(lims, lims, color="black", linewidth=0.8, linestyle="--", label="y = x")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.text(0.04, 0.96,
            f"Pearson r = {pr:.3f}  (log10)\nSpearman r = {sr:.3f}\nn = {len(sub)}",
            transform=ax.transAxes, fontsize=11, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
    ax.legend(fontsize=10, loc="lower right")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{outname}.pdf")
    plt.close(fig)
    print(f"[{outname}] Pearson r (log10) = {pr:.3f} (p={pp:.2e}); "
          f"Spearman r = {sr:.3f} (p={sp:.2e}); n={len(sub)}")


def _fc_vs_abs_plot(df: pd.DataFrame, fc_col: str, abs_col: str, base_col: str,
                    val3_col: str, tag: str, outname: str,
                    highlight_col: str = None, highlight_label: str = "highlighted") -> None:
    """Fold change (x, log) vs absolute change (y, symlog). syn1 baseline is
    encoded by marker size + opacity (single crimson hue), so abundant genes
    -- the ones that dominate the absolute change -- stand out. Because
    abs_change = baseline * (FC - 1), points fan out by baseline abundance.
    If highlight_col (a bool column) is given, those genes are recolored green."""
    sub = df[df[fc_col].notna() & (df[fc_col] > 0) &
             df[abs_col].notna() &
             df[base_col].notna() & (df[base_col] > 0)].copy()
    x, y = sub[fc_col], sub[abs_col]

    nz = y[y != 0].abs()
    linthresh = max(1e-3, float(np.nanmedian(nz))) if len(nz) else 1e-3

    # Encode baseline by alpha only (single hue, fixed dot size).
    DOT_SIZE = 28
    logb = np.log10(sub[base_col])
    lo, hi = float(logb.min()), float(logb.max())
    norm = (logb - lo) / (hi - lo) if hi > lo else pd.Series(0.5, index=logb.index)
    hue = (165 / 255, 15 / 255, 21 / 255)    # crimson #a50f15
    green = (16 / 255, 130 / 255, 60 / 255)  # highlight
    rgba = np.zeros((len(sub), 4))
    rgba[:, 0], rgba[:, 1], rgba[:, 2] = hue
    rgba[:, 3] = 0.12 + 0.78 * norm.to_numpy()

    hlmask = None
    if highlight_col is not None and highlight_col in sub:
        hlmask = sub[highlight_col].fillna(False).to_numpy().astype(bool)
        rgba[hlmask, 0], rgba[hlmask, 1], rgba[hlmask, 2] = green

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(x, y, s=DOT_SIZE, facecolors=rgba, edgecolors="none")
    if hlmask is not None:
        ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=green,
                                  markersize=7, label=f"{highlight_label} (n={int(hlmask.sum())})")],
                  loc="lower right", fontsize=9)

    ax.axvline(1, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.set_xlabel(rf"{tag} fold change  ($\mathrm{{rel}}_{{\mathrm{{syn3A}}}}/\mathrm{{rel}}_{{\mathrm{{syn1}}}}$)", fontsize=12)
    ax.set_ylabel(rf"{tag} absolute change  ($\mathrm{{rel}}_{{\mathrm{{syn3A}}}}-\mathrm{{rel}}_{{\mathrm{{syn1}}}}$)", fontsize=12)
    # ax.set_title(f"{tag}: fold change vs absolute change", fontsize=12)

    # Symmetric limits around FC=1 (x) and 0 (y) so loss-vs-gain asymmetry shows.
    mlog = float(np.log10(x).abs().max()) * 1.05
    ax.set_xlim(10 ** (-mlog), 10 ** mlog)
    ymax = float(y.abs().max()) * 1.15
    ax.set_ylim(-ymax, ymax)

    # Stats box (upper-left): linearity of syn1 vs syn3A abundances (log10).
    both = sub[(sub[base_col] > 0) & (sub[val3_col] > 0)]
    pr, _ = stats.pearsonr(np.log10(both[base_col]), np.log10(both[val3_col]))
    ax.text(0.03, 0.97, f"n = {len(sub)}\nPearson r \n({tag} syn1 vs syn3A) = {pr:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85))

    # Single-hue alpha-ramp colorbar (transparent -> solid steelblue), mirroring
    # the dot opacity encoding of the syn1 baseline.
    cmap = LinearSegmentedColormap.from_list("hue_alpha", [(*hue, 0.12), (*hue, 0.90)])
    sm = plt.cm.ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(f"log10 syn1 baseline (rel-{tag})", fontsize=10)

    # Annotate the top-5 absolute movers (both signs) AND the top-5 fold-change
    # outliers (both high and low), repelled. Union, de-duplicated by locus.
    top = pd.concat([sub[sub[abs_col] > 0].nlargest(5, abs_col),
                     sub[sub[abs_col] < 0].nsmallest(5, abs_col),
                     sub.nlargest(5, fc_col),
                     sub.nsmallest(5, fc_col)]).drop_duplicates(subset="locus_num")
    texts = [ax.text(r[fc_col], r[abs_col], f"{int(r['locus_num']):04d}", fontsize=7)
             for _, r in top.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{outname}.pdf")
    plt.close(fig)
    print(f"[{outname}] plotted n={len(sub)}; Pearson(syn1 vs syn3A)={pr:.3f}; "
          f"symlog linthresh={linthresh:.3g}")


def _operon_group_report(df_extreme: pd.DataFrame, folder: str, label: str,
                         tag: str, value_col: str, fc_col: str) -> None:
    """Write the per-operon summary txt for an extreme-FC group (no per-group CSVs)."""
    os.makedirs(folder, exist_ok=True)

    all_operons = set()
    for ops_str in df_extreme["sense_covering_ops"].dropna():
        for op in str(ops_str).split(","):
            all_operons.add(op.strip())

    lines = [f"# {tag} {label.upper()} genes  (n={len(df_extreme)})\n"]
    for opid in sorted(all_operons):
        op_genes = operon_cov[operon_cov["sense_covering_ops"].str.contains(opid, na=False)]
        extreme_in_op = df_extreme[df_extreme["sense_covering_ops"].str.contains(opid, na=False)]
        lines.append(f"\n## Operon {opid}  —  {len(op_genes)} total genes, "
                     f"{len(extreme_in_op)} with extreme FC\n")
        lines.append(f"{'locus_tag':<18} {'gene_name':<14} {f'{tag}_syn1':>12} "
                     f"{f'{tag}_syn3a':>12} {'FC':>8}\n")
        lines.append("-" * 70 + "\n")
        for _, row in op_genes.iterrows():
            r = df_extreme[df_extreme["locus_num"] == row["locus_num"]]
            if r.empty:
                r = result_lookup.get(tag, pd.DataFrame())
                r = r[r["locus_num"] == row["locus_num"]] if not r.empty else r
            if r.empty:
                v1 = v3 = fc_str = "N/A"
            else:
                rr = r.iloc[0]
                v1 = f"{rr[f'{value_col}_syn1']:.3f}"
                v3 = f"{rr[f'{value_col}_syn3a']:.3f}"
                fc_str = f"{rr[fc_col]:.3f}"
            marker = " <<<" if row["locus_num"] in extreme_in_op["locus_num"].values else ""
            gname  = str(row["gene_name"]) if pd.notna(row.get("gene_name")) else ""
            lines.append(f"{row['locus_tag']:<18} {gname:<14} {v1:>12} {v3:>12} {fc_str:>8}{marker}\n")

    with open(f"{folder}/operon_summary_{label}.txt", "w") as fh:
        fh.writelines(lines)
    print(f"[{tag} {label}] operon_summary written to {folder}/")


# ─────────────────────────────────────────────────────────────────────────────
# Load RNA tables. rna_type comes from the syn1 TPM CSV; we use syn1's
# classification to split syn3A genes (paired by locus_num) into subsets.
# ─────────────────────────────────────────────────────────────────────────────
syn1_tpm_full  = pd.read_csv(SYN1_TPM_CSV, sep="\t")
syn3a_tpm_full = pd.read_csv(SYN3A_TPM_TSV, sep="\t")
syn3a_ill_full = pd.read_csv(SYN3A_ILLUMINA_TPM_CSV).rename(
    columns={"Geneid": "locus_tag", SYN3A_ILLUMINA_TPM_COL: "TPM_mean"})

# Mask zero TPM as NaN so missing measurements stay NaN (not inf) through division.
syn1_tpm_full[SYN1_TPM_COL]   = syn1_tpm_full[SYN1_TPM_COL].mask(syn1_tpm_full[SYN1_TPM_COL]   <= 0)
syn3a_tpm_full[SYN3A_TPM_COL] = syn3a_tpm_full[SYN3A_TPM_COL].mask(syn3a_tpm_full[SYN3A_TPM_COL] <= 0)
syn3a_ill_full["TPM_mean"]    = syn3a_ill_full["TPM_mean"].mask(syn3a_ill_full["TPM_mean"]    <= 0)

for _df in (syn1_tpm_full, syn3a_tpm_full, syn3a_ill_full):
    _df["locus_num"] = _extract_locus_num(_df["locus_tag"])

syn1_tpm_full  = syn1_tpm_full.rename(columns={SYN1_TPM_COL:  "TPM_mean"})
syn3a_tpm_full = syn3a_tpm_full.rename(columns={SYN3A_TPM_COL: "TPM_mean"})

# Reference metadata: locus_num -> rna_type / gene_name / gene_product (from syn1).
syn1_meta = syn1_tpm_full[["locus_num", "rna_type", "gene_name", "gene_product"]].copy()
locus_to_rnatype = dict(zip(syn1_meta["locus_num"], syn1_meta["rna_type"]))


# ─────────────────────────────────────────────────────────────────────────────
# Load protein tables (raw iPM kept for the occupancy block before normalization).
# ─────────────────────────────────────────────────────────────────────────────
syn1_ptn  = pd.read_csv(SYN1_PTN_CSV, sep="\t")
syn3a_ptn = pd.read_csv(SYN3A_PTN_CSV, sep="\t")

syn1_ptn  = syn1_ptn[syn1_ptn["iPM_mean"].notna()   & (syn1_ptn["iPM_mean"]   > 0)].copy()
syn3a_ptn = syn3a_ptn[syn3a_ptn["iPM_mean"].notna() & (syn3a_ptn["iPM_mean"] > 0)].copy()
# The three syn3A pseudogenes carry MS intensity but no curated protein (no sequence,
# length, MW or localization), so they are kept out of the protein pool and its
# mean normaliser. They remain in the transcriptome tables, where all 496 loci count.
syn3a_ptn = syn3a_ptn[syn3a_ptn["primary_function"] != "Pseudogene"].copy()

n_syn1_ptn, n_syn3a_ptn = len(syn1_ptn), len(syn3a_ptn)
syn1_ptn["locus_num"]  = _extract_locus_num(syn1_ptn["locus_tag"])
syn3a_ptn["locus_num"] = _extract_locus_num(syn3a_ptn["locus_tag"])
syn1_ptn["rank_syn1_iPM"]   = syn1_ptn["iPM_mean"].rank(ascending=False, method="min")
syn3a_ptn["rank_syn3a_iPM"] = syn3a_ptn["iPM_mean"].rank(ascending=False, method="min")


# ─────────────────────────────────────────────────────────────────────────────
# #3  Deleted-gene occupancy of syn1's transcriptome and proteome.
#     Deleted gene = syn1 locus_num absent from the syn3A annotation.
#     Shares use RAW (un-normalized) syn1 TPM / iPM, i.e. fraction of the pool.
#     Returns the report text so it can be reused at the head of the outlier report.
# ─────────────────────────────────────────────────────────────────────────────
def _deleted_gene_occupancy() -> str:
    syn3a_anno = _load_gff_locus_nums(SYN3A_GFF)
    syn1_loci  = set(syn1_tpm_full["locus_num"])
    deleted    = syn1_loci - syn3a_anno

    coding_loci = {ln for ln, rt in locus_to_rnatype.items() if rt in CODING_RNA_TYPES}

    tpm = syn1_tpm_full[["locus_num", "TPM_mean", "rna_type"]].copy()
    tot_all   = tpm["TPM_mean"].sum(skipna=True)
    del_all   = tpm.loc[tpm["locus_num"].isin(deleted), "TPM_mean"].sum(skipna=True)
    tpm_cod   = tpm[tpm["locus_num"].isin(coding_loci)]
    tot_cod   = tpm_cod["TPM_mean"].sum(skipna=True)
    del_cod   = tpm_cod.loc[tpm_cod["locus_num"].isin(deleted), "TPM_mean"].sum(skipna=True)

    tot_ptn = syn1_ptn["iPM_mean"].sum(skipna=True)
    del_ptn = syn1_ptn.loc[syn1_ptn["locus_num"].isin(deleted), "iPM_mean"].sum(skipna=True)

    n_del_rna = int(tpm.loc[tpm["locus_num"].isin(deleted), "TPM_mean"].notna().sum())
    n_del_ptn = int(syn1_ptn["locus_num"].isin(deleted).sum())

    # RNA-type classification of the deleted loci
    del_types = (pd.Series([locus_to_rnatype.get(ln, "unknown") for ln in deleted])
                 .value_counts())

    lines = [
        "# Deleted-gene occupancy of syn1's transcriptome and proteome\n",
        f"# Deleted = syn1 loci absent from syn3A annotation ({SYN3A_GFF})\n\n",
        f"syn1 annotated loci (TPM table) : {len(syn1_loci)}\n",
        f"syn3A annotated loci (gff)      : {len(syn3a_anno)}\n",
        f"deleted loci                    : {len(deleted)}\n",
        f"  of which RNA-detected in syn1 : {n_del_rna}\n",
        f"  of which protein-detected     : {n_del_ptn}\n\n",
        "Deleted loci by RNA type\n",
    ]
    for rt, c in del_types.items():
        lines.append(f"  {str(rt):<8} : {c}\n")
    lines += [
        "\nTranscriptome (raw syn1 TPM share)\n",
        f"  total transcriptome (all RNA) : {del_all/tot_all*100:6.2f}%  "
        f"({del_all:,.0f} / {tot_all:,.0f} TPM)\n",
        f"  mRNA pool (coding only)       : {del_cod/tot_cod*100:6.2f}%  "
        f"({del_cod:,.0f} / {tot_cod:,.0f} TPM)\n\n",
        "Proteome (raw syn1 iPM share)\n",
        f"  total proteome                : {del_ptn/tot_ptn*100:6.2f}%  "
        f"({del_ptn:,.0f} / {tot_ptn:,.0f} iPM)\n",
    ]

    top_rna = (tpm[tpm["locus_num"].isin(deleted)]
               .merge(syn1_meta[["locus_num", "gene_name"]], on="locus_num", how="left")
               .dropna(subset=["TPM_mean"]).sort_values("TPM_mean", ascending=False).head(15))
    lines.append("\nTop 15 deleted genes by syn1 TPM\n")
    for _, r in top_rna.iterrows():
        gname = str(r["gene_name"]) if pd.notna(r["gene_name"]) else ""
        lines.append(f"  MMSYN1_{int(r['locus_num']):04d}  {gname:<12} {r['TPM_mean']:10.1f} TPM\n")

    top_ptn = (syn1_ptn.loc[syn1_ptn["locus_num"].isin(deleted), ["locus_num", "iPM_mean"]]
               .merge(syn1_meta[["locus_num", "gene_name"]], on="locus_num", how="left")
               .sort_values("iPM_mean", ascending=False).head(15))
    lines.append("\nTop 15 deleted genes by syn1 iPM\n")
    for _, r in top_ptn.iterrows():
        gname = str(r["gene_name"]) if pd.notna(r["gene_name"]) else ""
        lines.append(f"  MMSYN1_{int(r['locus_num']):04d}  {gname:<12} {r['iPM_mean']:10.1f} iPM\n")

    text = "".join(lines)
    with open(f"{OUTDIR}/deleted_gene_occupancy.txt", "w") as fh:
        fh.write(text)
    print(text)
    return text


OCCUPANCY_TEXT = _deleted_gene_occupancy()


# ─────────────────────────────────────────────────────────────────────────────
# Mean-normalize the protein layer (per-gene mean over detected proteins), then
# outer-merge syn1 <-> syn3a. After the occupancy block so raw shares are intact.
# ─────────────────────────────────────────────────────────────────────────────
syn1_ptn["iPM_mean"]  = syn1_ptn["iPM_mean"]  / syn1_ptn["iPM_mean"].mean()
syn3a_ptn["iPM_mean"] = syn3a_ptn["iPM_mean"] / syn3a_ptn["iPM_mean"].mean()

ptn = pd.merge(
    syn1_ptn[["locus_num", "locus_tag", "iPM_mean", "rank_syn1_iPM"]],
    syn3a_ptn[["locus_num", "locus_tag", "iPM_mean", "rank_syn3a_iPM"]],
    on="locus_num", how="outer", suffixes=("_syn1", "_syn3a"),
)
ptn["iPM_fold_change"] = ptn["iPM_mean_syn3a"] / ptn["iPM_mean_syn1"]
ptn["iPM_log10FC"]     = np.log10(ptn["iPM_fold_change"])
ptn["iPM_abs_change"]  = ptn["iPM_mean_syn3a"] - ptn["iPM_mean_syn1"]
ptn = ptn.merge(
    operon_cov[["locus_num", "sense_covering_ops", "antisense_covering_ops", "coverage_type"]],
    on="locus_num", how="left",
)
print(f"[iPM] syn1 detected={n_syn1_ptn}, syn3a detected={n_syn3a_ptn}, merged={len(ptn)}")


# ─────────────────────────────────────────────────────────────────────────────
# RNA layer: relativize each platform's subset, then build syn3A/syn1 comparison.
# ─────────────────────────────────────────────────────────────────────────────
def _relativize_subsets(rna_types: set):
    """Return mean-normalized (syn1 Illumina, syn3a ONT, syn3a Illumina) subsets."""
    keep = {ln for ln, rt in locus_to_rnatype.items() if rt in rna_types}

    def rel(df: pd.DataFrame) -> pd.DataFrame:
        sub = df[df["locus_num"].isin(keep)].copy()
        m = sub["TPM_mean"].mean(skipna=True)
        if pd.notna(m) and m > 0:
            sub["TPM_mean"] = sub["TPM_mean"] / m
        return sub

    return rel(syn1_tpm_full), rel(syn3a_tpm_full), rel(syn3a_ill_full), keep


def _compare(s_syn1: pd.DataFrame, s_syn3a: pd.DataFrame) -> pd.DataFrame:
    """Outer-merge a syn1 and a syn3a relativized subset into a wide comparison
    table with fold change + absolute change. Annotation comes from syn1_meta."""
    a = s_syn1[["locus_num", "locus_tag", "TPM_mean"]].copy()
    a["rank_syn1_TPM"] = a["TPM_mean"].rank(ascending=False, method="min")
    b = s_syn3a[["locus_num", "locus_tag", "TPM_mean"]].copy()
    b["rank_syn3a_TPM"] = b["TPM_mean"].rank(ascending=False, method="min")

    df = pd.merge(a, b, on="locus_num", how="outer", suffixes=("_syn1", "_syn3a"))
    df = df.merge(syn1_meta, on="locus_num", how="left")
    df["TPM_fold_change"] = df["TPM_mean_syn3a"] / df["TPM_mean_syn1"]
    df["TPM_log10FC"]     = np.log10(df["TPM_fold_change"].replace(0, np.nan))
    df["TPM_abs_change"]  = df["TPM_mean_syn3a"] - df["TPM_mean_syn1"]
    df = df.merge(
        operon_cov[["locus_num", "sense_covering_ops", "antisense_covering_ops", "coverage_type"]],
        on="locus_num", how="left",
    )
    return df


s1_c,  s3o_c,  s3i_c,  _ = _relativize_subsets(CODING_RNA_TYPES)
s1_nc, s3o_nc, s3i_nc, _ = _relativize_subsets(NONCODING_RNA_TYPES)

# Standard TPM comparison: Illumina for coding; ONT for non-coding (no Illumina rRNA/tRNA).
coding_tpm    = _compare(s1_c,  s3i_c)
noncoding_tpm = _compare(s1_nc, s3o_nc)
print(f"[coding]     Illumina standard: syn1 det={s1_c['TPM_mean'].notna().sum()}, "
      f"syn3a_Illumina det={s3i_c['TPM_mean'].notna().sum()}, syn3a_ONT det={s3o_c['TPM_mean'].notna().sum()}")
print(f"[non-coding] ONT standard:      syn1 det={s1_nc['TPM_mean'].notna().sum()}, "
      f"syn3a_ONT det={s3o_nc['TPM_mean'].notna().sum()}, syn3a_Illumina det={s3i_nc['TPM_mean'].notna().sum()}")

# Standard distribution + scatter PDFs (coding TPM = Illumina). iPM moved to 10.
_fc_plots(coding_tpm, value_col="TPM_mean", fc_col="TPM_fold_change",
          log10fc_col="TPM_log10FC", tag="TPM", ylabel_prefix="genes")

# ONT QC scatters (coding) — ONT is unreliable for quantification, shown only here.
ont_vs_syn1 = pd.merge(
    s1_c[["locus_num", "TPM_mean"]].rename(columns={"TPM_mean": "x"}),
    s3o_c[["locus_num", "TPM_mean"]].rename(columns={"TPM_mean": "y"}),
    on="locus_num", how="inner")
_corr_scatter(ont_vs_syn1, "x", "y",
              "syn1 Illumina rel-TPM (log)", "syn3A ONT rel-TPM (log)",
              "TPM_ONT_syn3A_vs_Illumina_syn1", "ONT syn3A vs Illumina syn1")

ont_vs_syn3i = pd.merge(
    s3i_c[["locus_num", "TPM_mean"]].rename(columns={"TPM_mean": "x"}),
    s3o_c[["locus_num", "TPM_mean"]].rename(columns={"TPM_mean": "y"}),
    on="locus_num", how="inner")
_corr_scatter(ont_vs_syn3i, "x", "y",
              "syn3A Illumina rel-TPM (log)", "syn3A ONT rel-TPM (log)",
              "TPM_ONT_syn3A_vs_Illumina_syn3A", "ONT syn3A vs Illumina syn3A (platform QC)")


# ─────────────────────────────────────────────────────────────────────────────
# Operon summaries for extreme-FC groups (standard TPM = Illumina; iPM).
# ─────────────────────────────────────────────────────────────────────────────
result_lookup = {"TPM": coding_tpm}

_operon_group_report(coding_tpm[coding_tpm["TPM_fold_change"] > 10].sort_values("TPM_fold_change", ascending=False),
                     f"{OUTDIR}/RNA_upgrade",   "FC_gt10",  "TPM", "TPM_mean", "TPM_fold_change")
_operon_group_report(coding_tpm[coding_tpm["TPM_fold_change"] < 0.1].sort_values("TPM_fold_change"),
                     f"{OUTDIR}/RNA_downgrade", "FC_lt0.1", "TPM", "TPM_mean", "TPM_fold_change")


# ─────────────────────────────────────────────────────────────────────────────
# Combined CSV outputs. Value columns are mean-normalized (rel* prefix). TPM is
# the standard platform per subset: coding = Illumina, non-coding = ONT.
# ─────────────────────────────────────────────────────────────────────────────
def _build_combined(tpm_df: pd.DataFrame, ptn_df=None) -> pd.DataFrame:
    out = tpm_df.rename(columns={
        "locus_tag_syn1":  "locus_syn1",
        "locus_tag_syn3a": "locus_syn3a",
        "TPM_mean_syn1":   "relTPM_syn1",
        "TPM_mean_syn3a":  "relTPM_syn3a",
    })[[
        "locus_num", "locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
        "relTPM_syn1", "rank_syn1_TPM",
        "relTPM_syn3a", "rank_syn3a_TPM",
        "TPM_fold_change", "TPM_abs_change",
        "sense_covering_ops",
    ]].copy()
    if ptn_df is not None:
        ptn_keep = ptn_df.rename(columns={
            "iPM_mean_syn1":  "relIPM_syn1",
            "iPM_mean_syn3a": "relIPM_syn3a",
        })[[
            "locus_num",
            "relIPM_syn1", "rank_syn1_iPM",
            "relIPM_syn3a", "rank_syn3a_iPM",
            "iPM_fold_change", "iPM_abs_change",
        ]]
        out = out.merge(ptn_keep, on="locus_num", how="outer")
    return out.sort_values("locus_num")


coding = _build_combined(coding_tpm, ptn)
# Prefer the syn3A proteome's gene_name for retained genes (syn1 as fallback);
# also brings in essentiality. syn3A genes carry the up-to-date annotation.
_syn3a_anno = pd.read_csv(SYN3A_PTN_CSV, sep="\t",
                          usecols=["locus_tag", "gene_name", "essentiality"])
_syn3a_anno["locus_num"] = _extract_locus_num(_syn3a_anno["locus_tag"])
coding = coding.merge(
    _syn3a_anno.rename(columns={"gene_name": "_gene_name_syn3a"})[
        ["locus_num", "_gene_name_syn3a", "essentiality"]],
    on="locus_num", how="left")
coding["gene_name"] = coding["_gene_name_syn3a"].fillna(coding["gene_name"])
coding = coding.drop(columns=["_gene_name_syn3a"])
# Protein-to-transcript ratio (PTR ~ translation-efficiency proxy; steady-state
# protein/mRNA, NOT Ribo-seq TE). PTR_fold_change == iPM_fold_change / TPM_fold_change.
coding["PTR_syn1"] = coding["relIPM_syn1"] / coding["relTPM_syn1"]
coding["PTR_syn3a"] = coding["relIPM_syn3a"] / coding["relTPM_syn3a"]
coding["PTR_fold_change"] = coding["PTR_syn3a"] / coding["PTR_syn1"]
coding_cols = [
    "locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
    "relTPM_syn1", "rank_syn1_TPM",
    "relTPM_syn3a", "rank_syn3a_TPM",
    "TPM_fold_change", "TPM_abs_change",
    "relIPM_syn1", "rank_syn1_iPM",
    "relIPM_syn3a", "rank_syn3a_iPM",
    "iPM_fold_change", "iPM_abs_change",
    "PTR_syn1", "PTR_syn3a", "PTR_fold_change",
    "sense_covering_ops",
]
coding[coding_cols].to_csv(OUT_CODING, sep="\t", index=False)
print(f"\nSaved: {OUT_CODING}  ({len(coding)} coding genes; TPM = Illumina)")

noncoding = _build_combined(noncoding_tpm, ptn_df=None)
noncoding_cols = [
    "locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
    "relTPM_syn1", "rank_syn1_TPM",
    "relTPM_syn3a", "rank_syn3a_TPM",
    "TPM_fold_change", "TPM_abs_change",
    "sense_covering_ops",
]
noncoding[noncoding_cols].to_csv(OUT_NONCODING, sep="\t", index=False)
print(f"Saved: {OUT_NONCODING}  ({len(noncoding)} non-coding genes; TPM = ONT)")


# Ribosomal protein flag (rps/rpl/rpm gene names or 'ribosomal protein' product).
_gn = coding["gene_name"].fillna("").astype(str)
_gp = coding["gene_product"].fillna("").astype(str).str.lower()
coding["is_rprotein"] = _gn.str.match(r"rp[slm]") | _gp.str.contains("ribosomal protein")

# Truncated product for the report tables (essentiality merged earlier).
coding["gene_product_disp"] = coding["gene_product"].fillna("").astype(str).str.slice(0, 45)

# Fold change vs absolute change diagnostic (TPM only; iPM moved to 10).
_fc_vs_abs_plot(coding, "TPM_fold_change", "TPM_abs_change", "relTPM_syn1",
                "relTPM_syn3a", "TPM", "TPM_FC_vs_absChange")
# Same TPM plot, ribosomal proteins highlighted green.
_fc_vs_abs_plot(coding, "TPM_fold_change", "TPM_abs_change", "relTPM_syn1",
                "relTPM_syn3a", "TPM", "TPM_FC_vs_absChange_rprotein",
                highlight_col="is_rprotein", highlight_label="ribosomal protein")


def _rprotein_stacked_bar() -> None:
    """Two stacked bars (syn1, syn3A) of the coding mRNA pool composition.
    syn1 carries a hatched red top block = mRNA share of the genes deleted in
    syn3A; both bars show the ribosomal-protein share (green). Visualizes that
    the r-protein share rises both because deleted genes vacate the pool and
    because r-proteins are genuinely upregulated."""
    syn3a_anno = _load_gff_locus_nums(SYN3A_GFF)
    c = coding.copy()
    c["is_deleted"] = ~c["locus_num"].isin(syn3a_anno)
    s1 = c[c["relTPM_syn1"].notna() & (c["relTPM_syn1"] > 0)]
    s3 = c[c["relTPM_syn3a"].notna() & (c["relTPM_syn3a"] > 0)]
    t1, t3 = s1["relTPM_syn1"].sum(), s3["relTPM_syn3a"].sum()
    rp1 = s1.loc[s1["is_rprotein"], "relTPM_syn1"].sum() / t1 * 100
    del1 = s1.loc[s1["is_deleted"], "relTPM_syn1"].sum() / t1 * 100
    other1 = 100 - rp1 - del1
    rp3 = s3.loc[s3["is_rprotein"], "relTPM_syn3a"].sum() / t3 * 100
    other3 = 100 - rp3

    green = (16 / 255, 130 / 255, 60 / 255)
    red = (0.84, 0.19, 0.15)
    gray = "lightgray"
    W = 0.6

    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    # syn1 (x=0): r-proteins, other retained, deleted (hatched red, top)
    ax.bar(0, rp1, W, color=green)
    ax.bar(0, other1, W, bottom=rp1, color=gray)
    ax.bar(0, del1, W, bottom=rp1 + other1, color=red, hatch="///", edgecolor="white", linewidth=0)
    # syn3A (x=1): r-proteins, other
    ax.bar(1, rp3, W, color=green)
    ax.bar(1, other3, W, bottom=rp3, color=gray)

    ax.text(0, rp1 / 2, f"r-proteins\n{rp1:.0f}%", ha="center", va="center",
            color="white", fontweight="bold", fontsize=9)
    ax.text(1, rp3 / 2, f"r-proteins\n{rp3:.0f}%", ha="center", va="center",
            color="white", fontweight="bold", fontsize=9)
    ax.text(0, rp1 + other1 + del1 / 2, f"deleted\n{del1:.0f}%", ha="center", va="center",
            color="white", fontweight="bold", fontsize=8)
    ax.text(0, rp1 + other1 / 2, "other\nmRNA", ha="center", va="center", color="dimgray", fontsize=8)
    ax.text(1, rp3 + other3 / 2, "other\nmRNA", ha="center", va="center", color="dimgray", fontsize=8)

    ax.set_xticks([0, 1]); ax.set_xticklabels(["syn1", "syn3A"], fontsize=13)
    ax.set_ylim(0, 100); ax.set_xlim(-0.6, 1.6)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/rprotein_share_stackedbar.pdf")
    plt.close(fig)
    print(f"[rprotein_share] syn1: rprot={rp1:.1f}% deleted={del1:.1f}% other={other1:.1f}% | "
          f"syn3A: rprot={rp3:.1f}% other={other3:.1f}%")


_rprotein_stacked_bar()


# ─────────────────────────────────────────────────────────────────────────────
# TPM fold change by gene_impact_class (08). Mechanism = promoter-source change,
# which acts on transcription -> tested on TPM only (not protein). Each class is
# compared to the `unaffected` baseline (Mann-Whitney U).
# ─────────────────────────────────────────────────────────────────────────────
def _impact_class_boxplot() -> None:
    ctx = pd.read_csv(IMPACT_CTX, sep="\t")
    ctx["locus_num"] = _extract_locus_num(ctx["locus_tag"])
    m = coding.merge(ctx[["locus_num", "gene_impact_class"]], on="locus_num", how="inner")
    m = m[m["TPM_fold_change"].notna() & (m["TPM_fold_change"] > 0)].copy()
    m["log10FC"] = np.log10(m["TPM_fold_change"])

    classes = [c for c in CLASS_ORDER if c in set(m["gene_impact_class"])]
    data = [m.loc[m["gene_impact_class"] == c, "log10FC"].values for c in classes]
    base = m.loc[m["gene_impact_class"] == "unaffected", "log10FC"].values

    fig, ax = plt.subplots(figsize=(10, 5.5))
    bp = ax.boxplot(data, positions=range(len(classes)), widths=0.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch in bp["boxes"]:
        patch.set_facecolor("#cfe2f3"); patch.set_alpha(0.85)

    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        ax.scatter(i + rng.uniform(-0.18, 0.18, size=len(d)), d,
                   s=10, color="#08519c", alpha=0.5, edgecolors="none", zorder=3)

    ax.axhline(0, color="red", linestyle="--", linewidth=1, label="FC = 1 (no change)")

    labels = []
    print("\n[gene_impact_class -> TPM fold change]")
    for c, d in zip(classes, data):
        if c == "unaffected" or len(d) < 2 or len(base) < 2:
            p, star = np.nan, ""
        else:
            p = stats.mannwhitneyu(d, base, alternative="two-sided").pvalue
            star = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
        labels.append(f"{c}\nn={len(d)}" + (f"\n{star}" if star else ""))
        med = float(np.median(d)) if len(d) else float("nan")
        print(f"  {c:<28} n={len(d):3d}  median log10FC={med:+.3f}  "
              f"median FC={10**med:.3f}  p(vs unaffected)={p:.2e}")

    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("log10 TPM fold change  (syn3A / syn1)", fontsize=12)
    ax.set_title("TPM fold change by gene_impact_class (transcriptional mechanism)", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/TPM_FC_by_impact_class.pdf")
    plt.close(fig)
    print(f"Saved: {OUTDIR}/TPM_FC_by_impact_class.pdf")


_impact_class_boxplot()


# ─────────────────────────────────────────────────────────────────────────────
# TPM change by curated FUNCTION CATEGORY (Secondary + Tertiary from
# syn3A_proteome.tsv). Generalizes the r-protein finding (a category
# can occupy more/less of the mRNA pool) to every functional group. For each
# category we report:
#   - mRNA-pool share in syn1 vs syn3A (% of the coding relTPM pool) + its change
#     (== which functional groups expand/shrink their slice of the mRNA pool),
#   - median TPM fold change and absolute change,
#   - Mann-Whitney U of log10 TPM FC vs all other categories.
# Outputs per level: a TSV table, a log10-FC boxplot, and a pool-share-change bar.
# Function categories annotate protein-coding genes; mapped to the mRNA layer by
# locus_num (MMSYN1_NNNN <-> JCVISYN3A_NNNN).
# ─────────────────────────────────────────────────────────────────────────────
def _function_category_tpm_analysis(min_n: int = 3) -> None:
    fa = load_syn3a_function()
    fa["locus_num"] = _extract_locus_num(fa["Locus Tag"])
    funccols = ["Primary Function", "Secondary Function", "Tertiary Function"]
    cm = coding.merge(fa[["locus_num"] + funccols], on="locus_num", how="left")

    # Deletion correction: normalize the syn1 side to the RETAINED-gene pool (loci
    # kept in syn3A), not the full syn1 pool, so genes deleted in syn3A no longer
    # inflate the retained genes' relative values. syn3A is entirely retained so
    # its own mean-normalization already references the retained pool.
    #   relTPM_syn1_ret = relTPM_syn1 / mean(relTPM_syn1 over retained, detected)
    # -> retained syn1 genes have mean 1, matched to syn3A; FC/abs use this basis.
    retained = _load_gff_locus_nums(SYN3A_GFF)
    ret = cm["locus_num"].isin(retained)
    r = cm.loc[ret & (cm["relTPM_syn1"] > 0), "relTPM_syn1"].mean()
    cm["relTPM_syn1_ret"] = cm["relTPM_syn1"] / r
    cm["TPM_FC_corr"]  = cm["relTPM_syn3a"] / cm["relTPM_syn1_ret"]
    cm["TPM_abs_corr"] = cm["relTPM_syn3a"] - cm["relTPM_syn1_ret"]
    print(f"\n[function categories] {int(cm['Secondary Function'].notna().sum())}/{len(cm)} "
          f"coding genes annotated; retained-pool renorm factor r={r:.3f}")

    # mRNA-pool totals over the RETAINED pool (deleted genes excluded from syn1)
    tot1 = cm.loc[ret & (cm["relTPM_syn1"] > 0), "relTPM_syn1"].sum()
    tot3 = cm.loc[cm["relTPM_syn3a"] > 0, "relTPM_syn3a"].sum()
    parent = {"Secondary Function": ["Primary Function"],
              "Tertiary Function":  ["Primary Function", "Secondary Function"]}

    def per_level(level_col: str, level_name: str) -> None:
        rows = []
        for cat, g in cm.dropna(subset=[level_col]).groupby(level_col):
            fc = g.loc[g["TPM_FC_corr"].notna() & (g["TPM_FC_corr"] > 0), "TPM_FC_corr"]
            ab = g["TPM_abs_corr"].dropna()
            s1 = g.loc[g["relTPM_syn1"]  > 0, "relTPM_syn1"].sum()
            s3 = g.loc[g["relTPM_syn3a"] > 0, "relTPM_syn3a"].sum()
            rest = cm[cm[level_col].notna() & (cm[level_col] != cat)]
            rfc = rest.loc[rest["TPM_FC_corr"].notna() & (rest["TPM_FC_corr"] > 0), "TPM_FC_corr"]
            p = (stats.mannwhitneyu(np.log10(fc), np.log10(rfc), alternative="two-sided").pvalue
                 if len(fc) >= min_n and len(rfc) >= min_n else np.nan)
            rec = {pc: (g[pc].dropna().iloc[0] if g[pc].notna().any() else "") for pc in parent[level_col]}
            rec.update(category=cat, n_genes=len(g), n_detected_FC=len(fc),
                       syn1_pool_share_pct=round(s1 / tot1 * 100, 3),
                       syn3a_pool_share_pct=round(s3 / tot3 * 100, 3),
                       pool_share_change=round((s3 / tot3 - s1 / tot1) * 100, 3),
                       pool_share_FC=round((s3 / tot3) / (s1 / tot1), 3) if s1 > 0 else np.nan,
                       median_TPM_FC_corr=round(float(fc.median()), 3) if len(fc) else np.nan,
                       median_TPM_abs_corr=round(float(ab.median()), 4) if len(ab) else np.nan,
                       mwu_p_vs_rest=p)
            rows.append(rec)
        cols = (parent[level_col] + ["category", "n_genes", "n_detected_FC",
                "syn1_pool_share_pct", "syn3a_pool_share_pct", "pool_share_change",
                "pool_share_FC", "median_TPM_FC_corr", "median_TPM_abs_corr", "mwu_p_vs_rest"])
        tbl = pd.DataFrame(rows)[cols].sort_values("syn1_pool_share_pct", ascending=False)
        out_tsv = f"{OUTDIR}/TPM_change_by_{level_name}.tsv"
        tbl.to_csv(out_tsv, sep="\t", index=False)
        print(f"Saved: {out_tsv}  ({len(tbl)} {level_name} categories; retained-pool normalized)")

        # log10 corrected-TPM-FC boxplot (categories with >= min_n; by median)
        def _fc(c):
            return np.log10(cm.loc[(cm[level_col] == c) & (cm["TPM_FC_corr"] > 0),
                                   "TPM_FC_corr"].dropna().values)
        cats = sorted([c for c in tbl["category"] if (tbl.loc[tbl["category"] == c, "n_detected_FC"] >= min_n).any()],
                      key=lambda c: np.median(_fc(c)))
        data = [_fc(c) for c in cats]
        fig, ax = plt.subplots(figsize=(max(7, len(cats) * 0.55), 6))
        bp = ax.boxplot(data, positions=range(len(cats)), widths=0.6, showfliers=False,
                        patch_artist=True, medianprops=dict(color="black"))
        for patch in bp["boxes"]:
            patch.set_facecolor("#cfe2f3"); patch.set_alpha(0.85)
        rng = np.random.default_rng(0)
        for i, d in enumerate(data):
            ax.scatter(i + rng.uniform(-0.18, 0.18, size=len(d)), d, s=10,
                       color="#08519c", alpha=0.5, edgecolors="none", zorder=3)
        ax.axhline(0, color="red", linestyle="--", linewidth=1, label="FC = 1 (no change)")
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels([f"{c}\nn={len(d)}" for c, d in zip(cats, data)],
                           rotation=40, ha="right", fontsize=8)
        ax.set_ylabel("log10 TPM fold change  (syn3A / syn1, retained-pool)", fontsize=12)
        ax.set_title(f"TPM fold change by {level_name} function (deletion-corrected)", fontsize=12)
        ax.legend(fontsize=9, loc="upper left")
        plt.tight_layout()
        fig.savefig(f"{OUTDIR}/TPM_FC_by_{level_name}_function.pdf")
        plt.close(fig)
        print(f"Saved: {OUTDIR}/TPM_FC_by_{level_name}_function.pdf ({len(cats)} categories, n>={min_n})")

        # mRNA-pool share change diverging bar (retained-pool shares; the
        # r-protein observation generalized: who expands/shrinks their slice)
        t = tbl.sort_values("pool_share_change")
        colors = ["#2c7fb8" if v >= 0 else "#d95f0e" for v in t["pool_share_change"]]
        fig, ax = plt.subplots(figsize=(7.5, max(4, len(t) * 0.3)))
        ax.barh(range(len(t)), t["pool_share_change"], color=colors,
                edgecolor="white", linewidth=0.4)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels(t["category"], fontsize=8)
        ax.set_xlabel("retained-pool share change  (syn3A% − syn1%)", fontsize=11)
        ax.set_title(f"Retained-pool share change by {level_name} function", fontsize=12)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        plt.tight_layout()
        fig.savefig(f"{OUTDIR}/mRNA_pool_share_change_by_{level_name}.pdf")
        plt.close(fig)
        print(f"Saved: {OUTDIR}/mRNA_pool_share_change_by_{level_name}.pdf")

    per_level("Secondary Function", "secondary")
    per_level("Tertiary Function", "tertiary")


_function_category_tpm_analysis()


# ─────────────────────────────────────────────────────────────────────────────
# Two story plots of how the mRNA pool is reapportioned across functions.
#   Plot 1 (secondary): two stacked bars of the FULL mRNA pool (syn1 incl. a
#     hatched 'deleted' block; syn3A). Segments = secondary function, colored by
#     Primary family (shades of one hue). Shows occupancy + the deleted slice.
#   Plot 2 (tertiary): dumbbell of the deletion-corrected RETAINED-pool share,
#     syn1 vs syn3A, for tertiary functions whose share changes by > 1 pp.
# ─────────────────────────────────────────────────────────────────────────────
PRIM_COLORS = {
    "Genetic Information Processing":       "#3b6db3",  # blue
    "Metabolism":                          "#3f9e5a",  # green
    "Unclear":                             "#9aa0a6",  # grey
    "Cellular Processes":                  "#8e6bb1",  # purple
    "Environmental Information Processing": "#2aa6a0",  # teal
    "Exogenous":                           "#c0654e",  # terracotta
}


def _shade(base: str, frac: float):
    """Lighten a base color toward white by frac in [0, 1]."""
    return tuple(c + (1 - c) * frac for c in mcolors.to_rgb(base))


def _text_color(rgb):
    lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return "white" if lum < 0.55 else "black"


def _roman(n: int) -> str:
    out = ""
    for v, s in [(10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]:
        while n >= v:
            out += s
            n -= v
    return out


def _load_function_for_plots():
    fa = load_syn3a_function()
    fa["locus_num"] = _extract_locus_num(fa["Locus Tag"])
    cm = coding.merge(fa[["locus_num", "Primary Function", "Secondary Function",
                          "Tertiary Function"]], on="locus_num", how="left")
    cm["is_deleted"] = ~cm["locus_num"].isin(_load_gff_locus_nums(SYN3A_GFF))
    return cm


def _pool_composition_secondary_bar() -> None:
    """LEFT story plot (tall, x:2x). Two full-pool stacked bars (syn1 w/ hatched
    'deleted' block | syn3A) of the mRNA pool. Blocks = secondary function, with
    Unclear & Cellular Processes collapsed to Primary; retained-but-unannotated
    ('other') is dropped and the shown blocks renormalized to 100%. Roman-indexed;
    'index %' inside / side-labeled if thin; no spines / title / y-axis; legend below."""
    from collections import defaultdict
    cm = _load_function_for_plots()
    cm["block"] = cm["Secondary Function"]
    coll = cm["Primary Function"].isin(["Unclear", "Cellular Processes"])
    cm.loc[coll, "block"] = cm.loc[coll, "Primary Function"]
    det1, det3 = cm["relTPM_syn1"] > 0, cm["relTPM_syn3a"] > 0
    annm = cm["block"].notna()

    # shares exclude retained-unannotated; syn1 = secondary blocks + deleted, syn3A
    # = secondary blocks; renormalized to 100% so 'other' is gone, not just hidden.
    s_sec1 = cm.loc[det1 & annm].groupby("block")["relTPM_syn1"].sum()
    s_del1 = cm.loc[det1 & cm["is_deleted"], "relTPM_syn1"].sum()
    tot1 = s_sec1.sum() + s_del1
    g1, deleted1 = s_sec1 / tot1 * 100, s_del1 / tot1 * 100
    s_sec3 = cm.loc[det3 & annm].groupby("block")["relTPM_syn3a"].sum()
    g3 = s_sec3 / s_sec3.sum() * 100

    block_prim = cm.dropna(subset=["block"]).groupby("block")["Primary Function"].first().to_dict()
    prim_order = list(PRIM_COLORS)
    blocks = sorted(set(g1.index) | set(g3.index),
                    key=lambda b: (prim_order.index(block_prim.get(b, "Unclear"))
                                   if block_prim.get(b) in prim_order else 99, -float(g1.get(b, 0))))
    prim_blocks = defaultdict(list)
    for b in blocks:
        prim_blocks[block_prim.get(b, "Unclear")].append(b)
    block_color = {}
    for p, bl in prim_blocks.items():
        base, n = PRIM_COLORS.get(p, "#9aa0a6"), len(bl)
        for i, b in enumerate(bl):
            block_color[b] = _shade(base, 0 if n == 1 else (i / (n - 1)) * 0.5)
    roman = {b: _roman(i + 1) for i, b in enumerate(blocks)}

    # Nature print-ready (OUTPUT.md): 3.5x7 in, Arial, 5-7pt, pdf.fonttype=42,
    # dpi=300, manual subplots_adjust (no bbox_inches='tight').
    W, MIN_INSIDE = 0.84, 2.0
    with plt.rc_context(NATURE_RC):
        fig, ax = plt.subplots(figsize=(7 / 3, 7 / 2))
        side = {0: [], 1: []}
        for xi, (g, deleted) in enumerate([(g1, deleted1), (g3, 0.0)]):
            bottom = 0.0
            for b in blocks:
                h = float(g.get(b, 0))
                if h <= 0:
                    continue
                ax.bar(xi, h, W, bottom=bottom, color=block_color[b], edgecolor="white", linewidth=0.4)
                if h >= MIN_INSIDE:
                    ax.text(xi, bottom + h / 2, f"{roman[b]} {h:.0f}%", ha="center", va="center",
                            fontsize=5, color=_text_color(block_color[b]))
                else:
                    side[xi].append((bottom + h / 2, f"{roman[b]} {h:.1f}%"))
                bottom += h
            if deleted > 0:
                ax.bar(xi, deleted, W, bottom=bottom, color=(0.84, 0.19, 0.15),
                       hatch="///", edgecolor="white", linewidth=0)
                ax.text(xi, bottom + deleted / 2, f"deleted {deleted:.0f}%", ha="center",
                        va="center", color="white", fontweight="bold", fontsize=5)

        def _place(items, x_text, ha, x_edge):
            items = sorted(items)
            ys = [y for y, _ in items]
            for i in range(1, len(ys)):                           # push up to maintain gap
                ys[i] = max(ys[i], ys[i - 1] + 2.4)
            if ys and ys[-1] > 99:                                # cap inside axes; propagate down
                ys[-1] = 99
                for i in range(len(ys) - 2, -1, -1):
                    ys[i] = min(ys[i], ys[i + 1] - 2.4)
            for (y0, txt), y in zip(items, ys):
                ax.annotate(txt, xy=(x_edge, y0), xytext=(x_text, y), ha=ha, va="center",
                            fontsize=5, color="#333",
                            arrowprops=dict(arrowstyle="-", color="gray", lw=0.3))
        _place(side[0], -0.60, "right", -W / 2)        # tighter for the narrow canvas
        _place(side[1],  1.60, "left",   1 + W / 2)

        ax.set_xticks([0, 1]); ax.set_xticklabels(["syn1", "syn3A"])   # xtick.labelsize=6
        ax.set_ylim(0, 100); ax.set_xlim(-1.35, 2.1)
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)

        handles = [Patch(facecolor=block_color[b], edgecolor="white", label=f"{roman[b]} — {b}")
                   for b in blocks]
        handles.append(Patch(facecolor=(0.84, 0.19, 0.15), hatch="///", edgecolor="white",
                             label="deleted (syn1 only)"))
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
                  ncol=2, frameon=False, fontsize=5,         # shrunk to fit the narrow canvas
                  handlelength=0.9, handletextpad=0.3, columnspacing=0.5)
        plt.subplots_adjust(left=0.10, right=0.95, top=0.97, bottom=0.24)
        fig.savefig(f"{OUTDIR}/mRNA_pool_composition_by_secondary.pdf", dpi=300)
        plt.close(fig)
    print(f"Saved: {OUTDIR}/mRNA_pool_composition_by_secondary.pdf "
          f"(deleted={deleted1:.1f}%, Translation {g1.get('Translation',0):.0f}->{g3.get('Translation',0):.0f}%)")


def _tertiary_share_dumbbell(min_change: float = 0.5) -> None:
    """RIGHT story plot (tall, x:2x). Retained-pool (deletion-corrected) tertiary
    dumbbell, syn1 vs syn3A, for tertiaries changing by > min_change pp. y-labels
    (Secondary above Tertiary) and dots are colored by the LEFT plot's block color;
    broken x-axis hides the empty mid band; share values sit above/below each dot;
    horizontal y-grid, no left/right/top spines."""
    import math
    from collections import defaultdict
    cm = _load_function_for_plots()
    det1, det3 = cm["relTPM_syn1"] > 0, cm["relTPM_syn3a"] > 0
    ret = ~cm["is_deleted"]

    # replicate the LEFT plot's block palette (secondary; Unclear/Cellular->Primary)
    cm["block"] = cm["Secondary Function"]
    coll = cm["Primary Function"].isin(["Unclear", "Cellular Processes"])
    cm.loc[coll, "block"] = cm.loc[coll, "Primary Function"]
    annm = cm["block"].notna()
    g1b = cm.loc[det1 & annm].groupby("block")["relTPM_syn1"].sum()
    g3k = set(cm.loc[det3 & annm, "block"])
    block_prim = cm.dropna(subset=["block"]).groupby("block")["Primary Function"].first().to_dict()
    prim_order = list(PRIM_COLORS)
    blocks = sorted(set(g1b.index) | g3k,
                    key=lambda b: (prim_order.index(block_prim.get(b, "Unclear"))
                                   if block_prim.get(b) in prim_order else 99, -float(g1b.get(b, 0))))
    prim_blocks = defaultdict(list)
    for b in blocks:
        prim_blocks[block_prim.get(b, "Unclear")].append(b)
    block_color = {}
    for p, bl in prim_blocks.items():
        base, nb = PRIM_COLORS.get(p, "#9aa0a6"), len(bl)
        for i, b in enumerate(bl):
            block_color[b] = _shade(base, 0 if nb == 1 else (i / (nb - 1)) * 0.5)

    # dumbbell data (retained pool, corrected)
    tot1 = cm.loc[ret & det1, "relTPM_syn1"].sum()
    tot3 = cm.loc[det3, "relTPM_syn3a"].sum()
    ann = cm.dropna(subset=["Tertiary Function"])
    s1 = ann.loc[det1.loc[ann.index]].groupby("Tertiary Function")["relTPM_syn1"].sum() / tot1 * 100
    s3 = ann.loc[det3.loc[ann.index]].groupby("Tertiary Function")["relTPM_syn3a"].sum() / tot3 * 100
    ter_prim = ann.groupby("Tertiary Function")["Primary Function"].first().to_dict()
    ter_sec = ann.groupby("Tertiary Function")["Secondary Function"].first().to_dict()
    d = pd.DataFrame({"syn1": s1, "syn3a": s3}).fillna(0.0)
    d["change"] = d["syn3a"] - d["syn1"]
    d = d[d["change"].abs() > min_change].sort_values("change")
    n = len(d)

    def _color(name):
        p, s = ter_prim.get(name), ter_sec.get(name)
        return block_color.get(p if p in ("Unclear", "Cellular Processes") else s, "#555555")

    # broken x-axis: hide the empty band between low cluster and high (ribosome)
    hi_row = d[["syn1", "syn3a"]].min(axis=1) >= 20
    lo_max = float(d.loc[~hi_row, ["syn1", "syn3a"]].values.max())
    hi_min = float(d.loc[hi_row, ["syn1", "syn3a"]].values.min())
    hi_max = float(d.loc[hi_row, ["syn1", "syn3a"]].values.max())
    lo_xlim = (0, math.ceil(lo_max) + 1)
    hi_xlim = (math.floor(hi_min) - 2, math.ceil(hi_max) + 2)
    dlo_w, dhi_w = lo_xlim[1] - lo_xlim[0], hi_xlim[1] - hi_xlim[0]

    # Nature print-ready (OUTPUT.md): 3.5x7 in, Arial, 5-7pt, pdf.fonttype=42,
    # dpi=300, manual subplots_adjust (no bbox_inches='tight').
    with plt.rc_context(NATURE_RC):
        import textwrap
        fig = plt.figure(figsize=(7 / 3, 7 / 2))
        gs = fig.add_gridspec(1, 2, width_ratios=[dlo_w, dhi_w], wspace=0.06)
        # full-width background axis: draws the y-grid continuously across the x break
        axbg = fig.add_subplot(gs[0, :])
        axbg.set_xlim(0, 1); axbg.set_ylim(-0.6, n - 0.4)
        for i in range(n):
            axbg.axhline(i, color="0.9", lw=0.5, zorder=0)
        axbg.set_xticks([]); axbg.set_yticks([])
        for sp in axbg.spines.values():
            sp.set_visible(False)
        axbg.patch.set_visible(False)
        axl = fig.add_subplot(gs[0]); axr = fig.add_subplot(gs[1], sharey=axl)
        axl.patch.set_visible(False); axr.patch.set_visible(False)
        axl.set_yticks(range(n))
        for i, (name, row) in enumerate(d.iterrows()):
            ax = axr if max(row["syn1"], row["syn3a"]) > lo_xlim[1] else axl
            col = _color(name)
            ax.plot([row["syn1"], row["syn3a"]], [i, i], color="lightgray", lw=0.8, zorder=2)
            ax.scatter(row["syn1"], i, s=22, facecolor="white", edgecolor="gray", linewidth=0.5, zorder=3)
            ax.scatter(row["syn3a"], i, s=30, color=col, edgecolor="white", linewidth=0.3, zorder=4)
            ax.text(row["syn3a"], i + 0.30, f"{row['syn3a']:.1f}", ha="center", va="center", fontsize=5, color="#333")
            ax.text(row["syn1"], i - 0.30, f"{row['syn1']:.1f}", ha="center", va="center", fontsize=5, color="#333")

        axl.set_xlim(*lo_xlim); axr.set_xlim(*hi_xlim)
        axl.set_ylim(-0.6, n - 0.4)

        # y labels (Secondary small italic over Tertiary), colored like the left blocks
        axl.set_yticklabels([])
        tr = axl.get_yaxis_transform()
        WRAP_W = 20            # wrap long Secondary/Tertiary names to fit the narrow panel
        for i, name in enumerate(d.index):
            c = _color(name)
            sec_txt = textwrap.fill(ter_sec.get(name, ""), WRAP_W)
            ter_txt = textwrap.fill(name, WRAP_W)
            axl.text(-0.04, i + 0.24, sec_txt, transform=tr, ha="right", va="center",
                     fontsize=5, color=c, style="italic", linespacing=0.9)
            axl.text(-0.04, i - 0.20, ter_txt, transform=tr, ha="right", va="center",
                     fontsize=5, color=c, linespacing=0.9)

        for sp in ("top", "left", "right"):
            axl.spines[sp].set_visible(False); axr.spines[sp].set_visible(False)
        axl.tick_params(left=False)        # xtick.labelsize=6 from rcParams
        axr.tick_params(left=False)
        plt.setp(axr.get_yticklabels(), visible=False)

        dd = 0.012                         # break marks at the x-axis only (bottom)
        kw = dict(transform=axl.transAxes, color="k", clip_on=False, lw=0.6)
        axl.plot((1 - dd, 1 + dd), (-dd, dd), **kw)
        kw.update(transform=axr.transAxes)
        axr.plot((-dd, dd), (-dd, dd), **kw)

        plt.subplots_adjust(left=0.32, right=0.97, top=0.97, bottom=0.14)
        dumb_mid = (axl.get_position().x0 + axr.get_position().x1) / 2
        fig.text(dumb_mid, 0.02, "mRNA Pool Share Change (%)", ha="center", fontsize=7)
        mk = [Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor="white",
                     markeredgecolor="gray", label="syn1"),
              Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor="#555",
                     markeredgecolor="white", label="syn3A")]
        axr.legend(handles=mk, loc="lower right", frameon=False)        # legend.fontsize=6
        fig.savefig(f"{OUTDIR}/tertiary_share_change_dumbbell.pdf", dpi=300)
        plt.close(fig)
    print(f"Saved: {OUTDIR}/tertiary_share_change_dumbbell.pdf "
          f"({n} tertiary, |Δ|>{min_change}pp; x-break {lo_xlim[1]:.0f}-{hi_xlim[0]:.0f}% hidden)")


if os.path.exists(f"{OUTDIR}/mRNA_pool_reallocation.pdf"):
    os.remove(f"{OUTDIR}/mRNA_pool_reallocation.pdf")
_pool_composition_secondary_bar()
_tertiary_share_dumbbell()


# ─────────────────────────────────────────────────────────────────────────────
# Outlier report -> Compare_RNA_Protein.txt (occupancy report merged at the top).
# Outliers via the Tukey 1.5*IQR rule: on log10 for fold change (multiplicative),
# on the raw delta for absolute change. Computed on the coding table.
# ─────────────────────────────────────────────────────────────────────────────
def _tukey_outliers(df: pd.DataFrame, col: str, log: bool):
    s = df[col].dropna()
    if log:
        s = s[s > 0]
        vals = np.log10(s)
    else:
        vals = s
    q1, q3 = vals.quantile(0.25), vals.quantile(0.75)
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    if log:
        lo, hi = 10 ** lo, 10 ** hi
    up = df[df[col] > hi].sort_values(col, ascending=False)
    dn = df[df[col] < lo].sort_values(col, ascending=True)
    return up, dn, lo, hi


def _outlier_section(title: str, col: str, log: bool, show_cols: list) -> str:
    up, dn, lo, hi = _tukey_outliers(coding, col, log)
    out = [f"\n{'=' * 78}\n## {title}\n",
           f"   Tukey fences ({'log10 ' if log else ''}1.5*IQR): low < {lo:.3g}, high > {hi:.3g}\n",
           f"\n-- HIGH outliers (n={len(up)}) --\n",
           (up[show_cols].to_string(index=False) + "\n") if len(up) else "   (none)\n",
           f"\n-- LOW outliers (n={len(dn)}) --\n",
           (dn[show_cols].to_string(index=False) + "\n") if len(dn) else "   (none)\n"]
    return "".join(out)


TPM_SHOW = ["locus_syn1", "gene_name", "essentiality", "gene_product_disp",
            "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change", "TPM_abs_change"]

report = [OCCUPANCY_TEXT,
          "\n\n", "#" * 78, "\n# TPM FOLD-CHANGE AND ABSOLUTE-CHANGE OUTLIERS (coding genes)\n",
          "#" * 78, "\n",
          "# TPM = Illumina (syn1 & syn3A). rel* = mean-normalized. iPM is in 10_Compare_Ptn.\n"]
# --- cross-organism conservation of the retained-gene expression landscape ---
# How well the relative expression hierarchy of genes kept in BOTH cells is
# preserved through minimization (Pearson on log10 of mean-normalized rel units).
def _conservation_r(df, a, b):
    d = df[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
    d = d[(d[a] > 0) & (d[b] > 0)]
    pr = float(np.corrcoef(np.log10(d[a]), np.log10(d[b]))[0, 1])
    sr = float(d[a].corr(d[b], method="spearman"))
    return pr, sr, len(d)

_tpr, _tsr, _tn = _conservation_r(coding[coding["rna_type"] == "mRNA"], "relTPM_syn1", "relTPM_syn3a")
_ipr, _isr, _in = _conservation_r(coding, "relIPM_syn1", "relIPM_syn3a")
print(f"[conservation] transcriptome (mRNA) r={_tpr:.3f} Spearman={_tsr:.3f} n={_tn}; "
      f"proteome r={_ipr:.3f} Spearman={_isr:.3f} n={_in}")
report += ["\n\n", "#" * 78,
           "\n# CROSS-ORGANISM CONSERVATION OF THE RETAINED-GENE EXPRESSION LANDSCAPE\n",
           "#" * 78, "\n",
           "# Pearson on log10 of mean-normalized relative units, genes retained in both cells.\n",
           f"transcriptome  relTPM syn1 vs syn3A (mRNA): Pearson(log10) r={_tpr:.3f}  Spearman={_tsr:.3f}  n={_tn}\n",
           f"proteome       relIPM syn1 vs syn3A:        Pearson(log10) r={_ipr:.3f}  Spearman={_isr:.3f}  n={_in}\n"]

report.append(_outlier_section("TPM fold change (Illumina)", "TPM_fold_change", True,  TPM_SHOW))
report.append(_outlier_section("TPM absolute change (relTPM delta)", "TPM_abs_change", False, TPM_SHOW))

# All ribosomal-protein genes (TPM + iPM), ascending locus number.
RPROT_COLS = ["locus_syn1", "gene_name", "essentiality", "gene_product_disp",
              "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change", "TPM_abs_change",
              "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change"]
rprot_tbl = coding[coding["is_rprotein"]].sort_values("locus_num")
report += ["\n\n", "#" * 78,
           f"\n# RIBOSOMAL PROTEIN GENES (n={len(rprot_tbl)}) -- TPM + iPM, ascending locus\n",
           "#" * 78, "\n",
           rprot_tbl[RPROT_COLS].to_string(index=False) + "\n"]

with open(OUT_REPORT, "w") as fh:
    fh.write("".join(report))
print(f"\nSaved: {OUT_REPORT}")

# Brief console summary of outlier counts
for title, col, log in [("TPM fold change", "TPM_fold_change", True),
                        ("TPM abs change",  "TPM_abs_change",  False)]:
    up, dn, lo, hi = _tukey_outliers(coding, col, log)
    print(f"  [{title}] high={len(up)}, low={len(dn)}  (fences {lo:.3g} .. {hi:.3g})")
