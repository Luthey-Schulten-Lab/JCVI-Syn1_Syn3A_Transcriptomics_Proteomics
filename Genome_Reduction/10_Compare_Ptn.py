"""
Protein (iPM) syn1 vs syn3A comparison — the proteome counterpart of 09.

Split out of 09_Compare_RNA_Protein.py: 09 owns the RNA/TPM story (and still
writes the joined TPM+iPM CSV); this script owns the protein-specific figures and
outlier report. It READS 09's combined table
(Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv), which already carries the
mean-normalized protein columns relIPM_syn1, relIPM_syn3a, iPM_fold_change,
iPM_abs_change. Outputs go to the same Compare_RNA_Protein/ folder.

Normalization / metrics are identical to 09 (see its docstring): rel = value /
per-gene mean over detected proteins; fold change = ratio, absolute change =
rel_syn3A - rel_syn1.

Proteome reallocation story plots (iPM) mirror 09's mRNA pair, using the curated
Secondary/Tertiary function annotation and the deletion-corrected (retained-pool)
normalization. NOTE: ribosomal-protein iPM is digestion-biased — read Ribosome /
Translation at the protein level with caution.

Outputs (Compare_RNA_Protein/)
  - iPM_fold_change_distribution.pdf, iPM_log10FC_distribution.pdf, iPM_correlation_syn1_vs_syn3a.pdf
  - iPM_FC_vs_absChange.pdf, iPM_FC_vs_absChange_rprotein.pdf
  - PTR_TPMfc_vs_iPMfc.pdf, PTR_by_category_boxplot.pdf
  - iPM_pool_composition_by_secondary.pdf : tall (x:2x) two full-pool stacked bars
        (syn1 w/ hatched deleted block | syn3A), blocks = secondary (Unclear & Cellular
        collapsed to Primary), Roman-indexed, colored by Primary.
  - iPM_tertiary_share_change_dumbbell.pdf : tall (x:2x) retained-pool deletion-corrected
        tertiary dumbbell (|delta|>0.5pp), broken x-axis; labels/dots colored like the bars.
  - protein_upgrade/, protein_downgrade/ (per-operon summaries)
  - Compare_Ptn.txt (iPM FC / abs-change outliers + ribosomal-protein iPM table)
"""

import os
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


# ── Inputs / outputs ─────────────────────────────────────────────────────────
IN_CSV       = "Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv"   # written by 09
SYN3A_PTN_CSV = "../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv"
OPERON_COV   = "../Syn1_Operon/segmentation/gene_operon_coverage.tsv"
OUTDIR       = "Compare_RNA_Protein"
OUT_REPORT   = f"{OUTDIR}/Compare_Ptn.txt"
os.makedirs(OUTDIR, exist_ok=True)

ABUNDANT_THRESH = 1.0

operon_cov = pd.read_csv(OPERON_COV, sep="\t")
operon_cov["locus_num"] = operon_cov["locus_tag"].str.extract(r"(\d+)$").astype(int)


def _extract_locus_num(s: pd.Series) -> pd.Series:
    return s.str.extract(r"(\d+)$").astype(int)


def _fc_plots(result, value_col, fc_col, log10fc_col, tag, ylabel_prefix,
              abundant_thresh=ABUNDANT_THRESH):
    """FC + log10FC histograms and the syn1-vs-syn3a scatter (mean-normalized)."""
    fc = result[fc_col]
    fc_stats = fc.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    print(f"\n[{tag}] fold change statistics (n={int(fc_stats['count'])}):")
    print(f"  mean={fc_stats['mean']:.3f}  median={fc_stats['50%']:.3f}  std={fc_stats['std']:.3f}")

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
    ax.plot(lims, lims, color="black", linewidth=0.8, linestyle="--")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"syn1  rel-{tag} (log)", fontsize=12)
    ax.set_ylabel(f"syn3A rel-{tag} (log)", fontsize=12)
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


def _fc_vs_abs_plot(df, fc_col, abs_col, base_col, val3_col, tag, outname,
                    highlight_col=None, highlight_label="highlighted"):
    """FC (x, log) vs absolute change (y, symlog); baseline by alpha; optional
    green highlight. Mirrors 09's _fc_vs_abs_plot."""
    sub = df[df[fc_col].notna() & (df[fc_col] > 0) &
             df[abs_col].notna() &
             df[base_col].notna() & (df[base_col] > 0)].copy()
    x, y = sub[fc_col], sub[abs_col]

    nz = y[y != 0].abs()
    linthresh = max(1e-3, float(np.nanmedian(nz))) if len(nz) else 1e-3

    DOT_SIZE = 28
    logb = np.log10(sub[base_col])
    lo, hi = float(logb.min()), float(logb.max())
    norm = (logb - lo) / (hi - lo) if hi > lo else pd.Series(0.5, index=logb.index)
    hue = (165 / 255, 15 / 255, 21 / 255)
    green = (16 / 255, 130 / 255, 60 / 255)
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
                  loc="lower left", fontsize=9)

    ax.axvline(1, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.set_xlabel(rf"{tag} fold change  ($\mathrm{{rel}}_{{\mathrm{{syn3A}}}}/\mathrm{{rel}}_{{\mathrm{{syn1}}}}$)", fontsize=12)
    ax.set_ylabel(rf"{tag} absolute change  ($\mathrm{{rel}}_{{\mathrm{{syn3A}}}}-\mathrm{{rel}}_{{\mathrm{{syn1}}}}$)", fontsize=12)

    mlog = float(np.log10(x).abs().max()) * 1.05
    ax.set_xlim(10 ** (-mlog), 10 ** mlog)
    ymax = float(y.abs().max()) * 1.15
    ax.set_ylim(-ymax, ymax)

    both = sub[(sub[base_col] > 0) & (sub[val3_col] > 0)]
    pr, _ = stats.pearsonr(np.log10(both[base_col]), np.log10(both[val3_col]))
    ax.text(0.03, 0.97, f"n = {len(sub)}\nPearson r ({tag} syn1 vs syn3A) = {pr:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85))

    cmap = LinearSegmentedColormap.from_list("hue_alpha", [(*hue, 0.12), (*hue, 0.90)])
    sm = plt.cm.ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(f"log10 syn1 baseline (rel-{tag})", fontsize=10)

    top = pd.concat([sub[sub[abs_col] > 0].nlargest(5, abs_col),
                     sub[sub[abs_col] < 0].nsmallest(5, abs_col)])
    texts = [ax.text(r[fc_col], r[abs_col], f"{int(r['locus_num']):04d}", fontsize=7)
             for _, r in top.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/{outname}.pdf")
    plt.close(fig)
    print(f"[{outname}] plotted n={len(sub)}; Pearson(syn1 vs syn3A)={pr:.3f}")


def _operon_group_report(df_extreme, folder, label, tag, value_col, fc_col):
    """Per-operon summary txt for an extreme-FC group."""
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
                v1 = f"{rr[f'{value_col}_syn1']:.3f}" if pd.notna(rr[f'{value_col}_syn1']) else "N/A"
                v3 = f"{rr[f'{value_col}_syn3a']:.3f}" if pd.notna(rr[f'{value_col}_syn3a']) else "N/A"
                fc_str = f"{rr[fc_col]:.3f}" if pd.notna(rr[fc_col]) else "N/A"
            marker = " <<<" if row["locus_num"] in extreme_in_op["locus_num"].values else ""
            gname = str(row["gene_name"]) if pd.notna(row.get("gene_name")) else ""
            lines.append(f"{row['locus_tag']:<18} {gname:<14} {v1:>12} {v3:>12} {fc_str:>8}{marker}\n")

    with open(f"{folder}/operon_summary_{label}.txt", "w") as fh:
        fh.writelines(lines)
    print(f"[{tag} {label}] operon_summary written to {folder}/")


def _tukey_outliers(df, col, log):
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


def _outlier_section(df, title, col, log, show_cols):
    up, dn, lo, hi = _tukey_outliers(df, col, log)
    out = [f"\n{'=' * 78}\n## {title}\n",
           f"   Tukey fences ({'log10 ' if log else ''}1.5*IQR): low < {lo:.3g}, high > {hi:.3g}\n",
           f"\n-- HIGH outliers (n={len(up)}) --\n",
           (up[show_cols].to_string(index=False) + "\n") if len(up) else "   (none)\n",
           f"\n-- LOW outliers (n={len(dn)}) --\n",
           (dn[show_cols].to_string(index=False) + "\n") if len(dn) else "   (none)\n"]
    return "".join(out)


# ─────────────────────────────────────────────────────────────────────────────
# Load 09's combined table and build the protein frame.
# ─────────────────────────────────────────────────────────────────────────────
coding = pd.read_csv(IN_CSV, sep="\t")
_ln = coding["locus_syn1"].fillna(coding["locus_syn3a"]).astype(str).str.extract(r"(\d+)$")[0]
coding["locus_num"] = pd.to_numeric(_ln, errors="coerce").astype("Int64")
_gn = coding["gene_name"].fillna("").astype(str)
_gp = coding["gene_product"].fillna("").astype(str).str.lower()
coding["is_rprotein"] = _gn.str.match(r"rp[slm]") | _gp.str.contains("ribosomal protein")

# Essentiality (syn3A design) + truncated product for the report tables.
_ess = pd.read_csv(SYN3A_PTN_CSV, sep="\t", usecols=["locus_tag", "essentiality"])
_ess["locus_num"] = _extract_locus_num(_ess["locus_tag"]).astype("Int64")
coding = coding.merge(_ess[["locus_num", "essentiality"]], on="locus_num", how="left")
coding["gene_product_disp"] = coding["gene_product"].fillna("").astype(str).str.slice(0, 45)

# Curated function annotation (Primary / Secondary / Tertiary) from the manually
# reviewed syn3A_proteome.tsv — replaces the earlier crude keyword-
# based func_category. PTR analysis groups dots by Primary (scatter) and Secondary
# (boxplot); both exclude ribosomal proteins (digestion-biased iPM).
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

PRIM_COLORS = {
    "Genetic Information Processing":       "#3b6db3",
    "Metabolism":                          "#3f9e5a",
    "Unclear":                             "#9aa0a6",
    "Cellular Processes":                  "#8e6bb1",
    "Environmental Information Processing": "#2aa6a0",
    "Exogenous":                           "#c0654e",
}
_fa = load_syn3a_function()
_fa["locus_num"] = _extract_locus_num(_fa["Locus Tag"]).astype("Int64")
coding = coding.merge(_fa[["locus_num", "Primary Function", "Secondary Function",
                           "Tertiary Function"]], on="locus_num", how="left")

ptn = coding.rename(columns={"relIPM_syn1": "iPM_mean_syn1",
                             "relIPM_syn3a": "iPM_mean_syn3a"}).copy()
ptn["iPM_log10FC"] = np.log10(ptn["iPM_fold_change"].replace(0, np.nan))
n_det = int((ptn["iPM_mean_syn1"].notna() & ptn["iPM_mean_syn3a"].notna()).sum())
print(f"[iPM] proteins with syn1&syn3A iPM: {n_det}; rows: {len(ptn)}")

# ── Distribution / correlation + FC-vs-abs figures ───────────────────────────
_fc_plots(ptn, value_col="iPM_mean", fc_col="iPM_fold_change",
          log10fc_col="iPM_log10FC", tag="iPM", ylabel_prefix="proteins")
_fc_vs_abs_plot(ptn, "iPM_fold_change", "iPM_abs_change", "iPM_mean_syn1",
                "iPM_mean_syn3a", "iPM", "iPM_FC_vs_absChange")
_fc_vs_abs_plot(ptn, "iPM_fold_change", "iPM_abs_change", "iPM_mean_syn1",
                "iPM_mean_syn3a", "iPM", "iPM_FC_vs_absChange_rprotein",
                highlight_col="is_rprotein", highlight_label="ribosomal protein")

# ── Operon summaries for extreme-FC protein groups ───────────────────────────
result_lookup = {"iPM": ptn}
_operon_group_report(ptn[ptn["iPM_fold_change"] > 10].sort_values("iPM_fold_change", ascending=False),
                     f"{OUTDIR}/protein_upgrade",   "FC_gt10",  "iPM", "iPM_mean", "iPM_fold_change")
_operon_group_report(ptn[ptn["iPM_fold_change"] < 0.1].sort_values("iPM_fold_change"),
                     f"{OUTDIR}/protein_downgrade", "FC_lt0.1", "iPM", "iPM_mean", "iPM_fold_change")

# ── PTR (protein-to-transcript ratio) = iPM_FC / TPM_FC ──────────────────────
# PTR > baseline => protein rose more than its mRNA (post-transcriptional gain).
# Categories driven by the curated function annotation (replaces the old keyword
# scheme). Ribosomal proteins are excluded throughout — their iPM is digestion-
# biased and would distort PTR estimates for the Translation/Ribosome group.
from matplotlib.patches import Patch as _Patch


def _ptr_scatter():
    d = coding[coding["TPM_fold_change"].notna() & (coding["TPM_fold_change"] > 0) &
               coding["iPM_fold_change"].notna() & (coding["iPM_fold_change"] > 0) &
               ~coding["is_rprotein"].fillna(False)].copy()
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    prim_order = list(PRIM_COLORS)
    for prim in prim_order:
        sub = d[d["Primary Function"] == prim]
        if sub.empty:
            continue
        ax.scatter(sub["TPM_fold_change"], sub["iPM_fold_change"],
                   s=22, alpha=0.7, color=PRIM_COLORS[prim], edgecolors="none",
                   label=f"{prim} (n={len(sub)})")
    unann = d[d["Primary Function"].isna()]
    if not unann.empty:
        ax.scatter(unann["TPM_fold_change"], unann["iPM_fold_change"],
                   s=14, alpha=0.5, color="lightgray", edgecolors="none",
                   label=f"unannotated (n={len(unann)})")
    # annotate top |log10(PTR)| outliers (any category)
    d["_dev"] = (np.log10(d["iPM_fold_change"]) - np.log10(d["TPM_fold_change"])).abs()
    top = d.nlargest(15, "_dev")
    texts = []
    for _, r in top.iterrows():
        lab = str(r.get("gene_name")) if pd.notna(r.get("gene_name")) and r.get("gene_name") else f"{int(r['locus_num']):04d}"
        texts.append(ax.text(r["TPM_fold_change"], r["iPM_fold_change"], lab, fontsize=6))
    lims = [min(d["TPM_fold_change"].min(), d["iPM_fold_change"].min()) * 0.5,
            max(d["TPM_fold_change"].max(), d["iPM_fold_change"].max()) * 2]
    ax.plot(lims, lims, "k--", lw=0.9, label="y = x (PTR unchanged)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(r"TPM fold change  ($\mathrm{rel}_{\mathrm{syn3A}}/\mathrm{rel}_{\mathrm{syn1}}$)", fontsize=12)
    ax.set_ylabel(r"iPM fold change  ($\mathrm{rel}_{\mathrm{syn3A}}/\mathrm{rel}_{\mathrm{syn1}}$)", fontsize=12)
    ax.set_title("Transcript vs protein fold change (PTR; r-proteins excluded; color = Primary)", fontsize=11)
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))
    ax.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/PTR_TPMfc_vs_iPMfc.pdf")
    plt.close(fig)
    print(f"[PTR scatter] n={len(d)} (r-proteins excluded)")


def _ptr_boxplot(min_n: int = 3):
    """Boxplot of log10 PTR-FC per Secondary function (curated annotation).
    r-proteins excluded; only Secondaries with >= min_n proteins shown;
    boxes colored by Primary family. Mann-Whitney p tests each category vs the
    overall (non-r-protein) baseline."""
    d = coding[coding["PTR_fold_change"].notna() & (coding["PTR_fold_change"] > 0)
               & ~coding["is_rprotein"].fillna(False)
               & coding["Secondary Function"].notna()].copy()
    d["log10PTR"] = np.log10(d["PTR_fold_change"])
    base = d["log10PTR"].values
    base_med = float(np.median(base)) if len(base) else 0.0

    sec_prim = d.groupby("Secondary Function")["Primary Function"].first().to_dict()
    counts = d.groupby("Secondary Function").size()
    prim_order = list(PRIM_COLORS)
    secs = sorted([s for s, n in counts.items() if n >= min_n],
                  key=lambda s: (prim_order.index(sec_prim.get(s, "Unclear"))
                                 if sec_prim.get(s) in prim_order else 99, s))
    data = [d.loc[d["Secondary Function"] == s, "log10PTR"].values for s in secs]

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    bp = ax.boxplot(data, positions=range(len(secs)), widths=0.65, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch, s in zip(bp["boxes"], secs):
        col = PRIM_COLORS.get(sec_prim.get(s, "Unclear"), "#777")
        patch.set_facecolor(col); patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, dd in enumerate(data):
        ax.scatter(i + rng.uniform(-0.20, 0.20, size=len(dd)), dd,
                   s=10, color="black", alpha=0.45, edgecolors="none", zorder=3)
    ax.axhline(base_med, color="red", linestyle="--", linewidth=1,
               label=f"all-genes median ({base_med:+.2f})")
    ax.axhline(0, color="black", linestyle=":", linewidth=0.7, alpha=0.6)

    labels = []
    print(f"\n[PTR by Secondary function (r-proteins excluded, n>={min_n})]")
    for s, dd in zip(secs, data):
        if len(dd) < 2 or len(base) < 2:
            p, sig = np.nan, "n/a"
        else:
            p = stats.mannwhitneyu(dd, base, alternative="two-sided").pvalue
            sig = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."
        labels.append(f"{s}\nn={len(dd)}  {sig}")
        med = float(np.median(dd))
        print(f"  {s:<35} n={len(dd):3d}  median log10PTR={med:+.3f}  p(vs all)={p:.2e}")
    ax.set_xticks(range(len(secs)))
    ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("log10 PTR fold change  (iPM_FC / TPM_FC)", fontsize=12)
    ax.set_title("Protein-to-transcript ratio change by Secondary function "
                 "(r-proteins excluded; color = Primary)", fontsize=11)
    used = []
    for p in prim_order:
        if any(sec_prim.get(s) == p for s in secs):
            used.append(p)
    handles = [_Patch(facecolor=PRIM_COLORS[p], alpha=0.55, label=p) for p in used]
    handles.append(plt.Line2D([0], [0], color="red", linestyle="--", label="all-genes median"))
    ax.legend(handles=handles, fontsize=8, loc="upper right", frameon=False)
    plt.tight_layout()
    fig.savefig(f"{OUTDIR}/PTR_by_category_boxplot.pdf")
    plt.close(fig)


_ptr_scatter()
_ptr_boxplot()


# ─────────────────────────────────────────────────────────────────────────────
# Proteome reallocation story plots (iPM) — the protein counterpart of 09's two
# mRNA plots. Same design: deletion-corrected (retained-pool) normalization, the
# curated Secondary/Tertiary function annotation, Primary-family colors.
#   1) iPM_pool_composition_by_secondary.pdf — two FULL-pool stacked bars (syn1
#      with hatched 'deleted' block | syn3A); blocks = secondary (Unclear &
#      Cellular Processes collapsed to Primary; retained-unannotated dropped),
#      Roman-indexed, colored by Primary; no spines/title/y-axis; legend below.
#   2) iPM_tertiary_share_change_dumbbell.pdf — retained-pool tertiary dumbbell
#      (|Δ|>0.5pp), broken x-axis, labels/dots colored like (1), share at dots.
# NOTE: ribosomal-protein iPM is digestion-biased; read Ribosome/Translation at
# the protein level with caution.
# ─────────────────────────────────────────────────────────────────────────────
import re as _re
SYN3A_GFF = "../Genomes_Input/syn3a_genome.gff3"
PRIM_COLORS = {
    "Genetic Information Processing":       "#3b6db3",
    "Metabolism":                          "#3f9e5a",
    "Unclear":                             "#9aa0a6",
    "Cellular Processes":                  "#8e6bb1",
    "Environmental Information Processing": "#2aa6a0",
    "Exogenous":                           "#c0654e",
}


def _shade(base, frac):
    return tuple(c + (1 - c) * frac for c in mcolors.to_rgb(base))


def _text_color(rgb):
    lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return "white" if lum < 0.55 else "black"


def _roman(num):
    out = ""
    for v, s in [(10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]:
        while num >= v:
            out += s; num -= v
    return out


def _load_gff_locus_nums(path):
    nums = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = _re.search(r"locus_tag=([^;]+)", f[8])
            mm = _re.search(r"(\d+)$", m.group(1)) if m else None
            if mm:
                nums.add(int(mm.group(1)))
    return nums


def _load_function_for_plots():
    cm = coding.copy()
    if "Primary Function" not in cm.columns:                # already merged at top now
        fa = load_syn3a_function()
        fa["locus_num"] = _extract_locus_num(fa["Locus Tag"]).astype("Int64")
        cm = cm.merge(fa[["locus_num", "Primary Function", "Secondary Function",
                          "Tertiary Function"]], on="locus_num", how="left")
    cm["is_deleted"] = ~cm["locus_num"].isin(_load_gff_locus_nums(SYN3A_GFF))
    return cm


def _block_palette(cm):
    """Add cm['block'] (secondary; Unclear/Cellular -> Primary) and return the
    Primary-family-shaded color per block. Ordering/shading use the syn1 mRNA
    (TPM) share so the Roman index + colors match 09's composition bar exactly
    (left-bar consistency across 09/10); bar heights here still use iPM."""
    from collections import defaultdict
    cm["block"] = cm["Secondary Function"]
    coll = cm["Primary Function"].isin(["Unclear", "Cellular Processes"])
    cm.loc[coll, "block"] = cm.loc[coll, "Primary Function"]
    annm = cm["block"].notna()
    t1, t3 = cm["relTPM_syn1"] > 0, cm["relTPM_syn3a"] > 0
    g1b = cm.loc[t1 & annm].groupby("block")["relTPM_syn1"].sum()
    g3k = set(cm.loc[t3 & annm, "block"])
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
    return blocks, block_color


def _ipm_composition_secondary_bar():
    """Plot 1: two full-pool proteome (iPM) stacked bars by secondary function."""
    cm = _load_function_for_plots()
    blocks, block_color = _block_palette(cm)
    det1, det3 = cm["relIPM_syn1"] > 0, cm["relIPM_syn3a"] > 0
    annm = cm["block"].notna()
    s_sec1 = cm.loc[det1 & annm].groupby("block")["relIPM_syn1"].sum()
    s_del1 = cm.loc[det1 & cm["is_deleted"], "relIPM_syn1"].sum()
    tot1 = s_sec1.sum() + s_del1
    g1, deleted1 = s_sec1 / tot1 * 100, s_del1 / tot1 * 100
    s_sec3 = cm.loc[det3 & annm].groupby("block")["relIPM_syn3a"].sum()
    g3 = s_sec3 / s_sec3.sum() * 100
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
        _place(side[0], -0.60, "right", -W / 2)       # tighter for the narrow canvas
        _place(side[1],  1.60, "left",  1 + W / 2)

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
        fig.savefig(f"{OUTDIR}/iPM_pool_composition_by_secondary.pdf", dpi=300)
        plt.close(fig)
    print(f"Saved: {OUTDIR}/iPM_pool_composition_by_secondary.pdf "
          f"(deleted={deleted1:.1f}%, Translation {g1.get('Translation',0):.0f}->{g3.get('Translation',0):.0f}%)")


def _ipm_tertiary_share_dumbbell(min_change=0.5):
    """Plot 2: retained-pool (deletion-corrected) tertiary iPM dumbbell."""
    import math
    cm = _load_function_for_plots()
    _, block_color = _block_palette(cm)
    det1, det3 = cm["relIPM_syn1"] > 0, cm["relIPM_syn3a"] > 0
    ret = ~cm["is_deleted"]
    tot1 = cm.loc[ret & det1, "relIPM_syn1"].sum()
    tot3 = cm.loc[det3, "relIPM_syn3a"].sum()
    ann = cm.dropna(subset=["Tertiary Function"])
    s1 = ann.loc[det1.loc[ann.index]].groupby("Tertiary Function")["relIPM_syn1"].sum() / tot1 * 100
    s3 = ann.loc[det3.loc[ann.index]].groupby("Tertiary Function")["relIPM_syn3a"].sum() / tot3 * 100
    ter_prim = ann.groupby("Tertiary Function")["Primary Function"].first().to_dict()
    ter_sec = ann.groupby("Tertiary Function")["Secondary Function"].first().to_dict()
    d = pd.DataFrame({"syn1": s1, "syn3a": s3}).fillna(0.0)
    d["change"] = d["syn3a"] - d["syn1"]
    d = d[d["change"].abs() > min_change].sort_values("change")
    n = len(d)

    def _color(name):
        p, s = ter_prim.get(name), ter_sec.get(name)
        return block_color.get(p if p in ("Unclear", "Cellular Processes") else s, "#555555")

    hi_row = d[["syn1", "syn3a"]].min(axis=1) >= 20
    assert hi_row.any() and (~hi_row).any(), "no clear high/low cluster for the broken axis"
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

        dd = 0.012
        kw = dict(transform=axl.transAxes, color="k", clip_on=False, lw=0.6)
        axl.plot((1 - dd, 1 + dd), (-dd, dd), **kw)
        kw.update(transform=axr.transAxes)
        axr.plot((-dd, dd), (-dd, dd), **kw)

        plt.subplots_adjust(left=0.32, right=0.97, top=0.97, bottom=0.14)
        dumb_mid = (axl.get_position().x0 + axr.get_position().x1) / 2
        fig.text(dumb_mid, 0.02, "Proteome Pool Share Change (%)", ha="center", fontsize=7)
        mk = [Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor="white",
                     markeredgecolor="gray", label="syn1"),
              Line2D([0], [0], marker="o", linestyle="", markersize=5, markerfacecolor="#555",
                     markeredgecolor="white", label="syn3A")]
        axr.legend(handles=mk, loc="lower right", frameon=False)        # legend.fontsize=6
        fig.savefig(f"{OUTDIR}/iPM_tertiary_share_change_dumbbell.pdf", dpi=300)
        plt.close(fig)
    print(f"Saved: {OUTDIR}/iPM_tertiary_share_change_dumbbell.pdf "
          f"({n} tertiary, |Δ|>{min_change}pp; x-break {lo_xlim[1]:.0f}-{hi_xlim[0]:.0f}% hidden)")


_ipm_composition_secondary_bar()
_ipm_tertiary_share_dumbbell()


# ─────────────────────────────────────────────────────────────────────────────
# Macromolecule complex abundance — limiting-subunit (stoichiometry-adjusted MIN)
# estimate for the transcription engine and the RNA-turnover machinery.
#
# Motivation: relate the syn1->syn3A change to syn3A's longer cell cycle
# (105 vs 60 min). Ribosomes can't be estimated here (no rRNA quantification),
# so we probe RNA polymerase (transcription engine) and the degradosome (RNA
# turnover) — both have only protein-coding subunits measured. Excess subunits
# get degraded, so the assembled complex copy number is set by the lowest
# stoichiometry-adjusted subunit (MIN). We compute it for relTPM and relIPM in
# syn1 and syn3A. rel values are per-organism mean-normalized; the syn3A/syn1
# ratio is the interpretable cross-organism quantity.
#
# Subunit identities (curated):
#   RNAP core α2ββ': rpoA=0645 (x2 per complex), rpoC=0803, rpoB=0804.
#   Degradosome: rny=0359, rnjA=0600, plus a 3'->5' exoribonuclease pair summed
#                (yhaM=0437 + rnr=0775).
# σ factors are not included (this is core RNAP, not holoenzyme).
# ─────────────────────────────────────────────────────────────────────────────
COMPLEXES = [
    # (name, printable formula, [ [ (locus, divisor, gene), ... ], ... ])
    ("RNA polymerase", "MIN(rpoA/2, rpoC, rpoB)",
     [[(645, 2, "rpoA")], [(803, 1, "rpoC")], [(804, 1, "rpoB")]]),
    ("Degradosome", "MIN(rny, rnjA, yhaM+rnr)",
     [[(359, 1, "rny")], [(600, 1, "rnjA")], [(437, 1, "yhaM"), (775, 1, "rnr")]]),
]
OUT_COMPLEX = f"{OUTDIR}/macromolecule_complex_abundance.tsv"
_LAYERS = ("relTPM_syn1", "relTPM_syn3a", "relIPM_syn1", "relIPM_syn3a")


def _gene_value(locus, col):
    row = coding[coding["locus_num"] == locus]
    if row.empty:
        return None
    v = row[col].iloc[0]
    return float(v) if pd.notna(v) else None


def _term_label(term):
    return "+".join(f"{gn}/{div}" if div > 1 else gn for _, div, gn in term)


def _complex_min(terms, col):
    vals = []
    for term in terms:
        s = 0.0
        for locus, div, _ in term:
            v = _gene_value(locus, col)
            if v is None:
                s = None; break
            s += v / div
        vals.append(s)
    if any(v is None for v in vals):
        return None, None
    i = min(range(len(vals)), key=lambda k: vals[k])
    return vals[i], _term_label(terms[i])


_cx_rows, _cx_detail = [], []
for cx_name, formula, terms in COMPLEXES:
    rec = {"complex": cx_name, "formula": formula}
    for col in _LAYERS:
        v, lab = _complex_min(terms, col)
        rec[col] = v
        rec[f"limiter_{col}"] = lab
    rec["TPM_fold_change"] = (rec["relTPM_syn3a"] / rec["relTPM_syn1"]
                              if rec["relTPM_syn1"] and rec["relTPM_syn3a"] else None)
    rec["iPM_fold_change"] = (rec["relIPM_syn3a"] / rec["relIPM_syn1"]
                              if rec["relIPM_syn1"] and rec["relIPM_syn3a"] else None)
    _cx_rows.append(rec)
    _cx_detail.append(f"\n  {cx_name}   [{formula}]")
    for term in terms:
        for locus, div, gn in term:
            vals = {col: _gene_value(locus, col) for col in _LAYERS}
            f = lambda v: f"{v:.3f}" if v is not None else "  N/A"
            _cx_detail.append(
                f"    {gn:<6} JCVISYN3A_{locus:04d} (/{div}) | "
                f"relTPM syn1={f(vals['relTPM_syn1'])} syn3A={f(vals['relTPM_syn3a'])} | "
                f"relIPM syn1={f(vals['relIPM_syn1'])} syn3A={f(vals['relIPM_syn3a'])}"
            )

complex_tbl = pd.DataFrame(_cx_rows)
_CX_COLS = ["complex", "formula",
            "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change",
            "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change",
            "limiter_relTPM_syn1", "limiter_relTPM_syn3a",
            "limiter_relIPM_syn1", "limiter_relIPM_syn3a"]
complex_tbl[_CX_COLS].to_csv(OUT_COMPLEX, sep="\t", index=False, float_format="%.4f")
print(f"\nSaved: {OUT_COMPLEX}")

print("\n[macromolecule complex abundance — limiting-subunit MIN]")
for rec in _cx_rows:
    f = lambda v: f"{v:.3f}" if v is not None else "N/A"
    fc = lambda v: f"{v:.2f}" if v is not None else "N/A"
    print(f"  {rec['complex']:<16} {rec['formula']}")
    print(f"    relTPM: syn1={f(rec['relTPM_syn1'])} (limiter:{rec['limiter_relTPM_syn1']})  "
          f"syn3A={f(rec['relTPM_syn3a'])} (limiter:{rec['limiter_relTPM_syn3a']})  "
          f"FC={fc(rec['TPM_fold_change'])}")
    print(f"    relIPM: syn1={f(rec['relIPM_syn1'])} (limiter:{rec['limiter_relIPM_syn1']})  "
          f"syn3A={f(rec['relIPM_syn3a'])} (limiter:{rec['limiter_relIPM_syn3a']})  "
          f"FC={fc(rec['iPM_fold_change'])}")

# Text section appended to Compare_Ptn.txt (joined into the report list below).
COMPLEX_REPORT_TEXT = ("\n\n" + "#" * 78
    + "\n# MACROMOLECULE COMPLEX ABUNDANCE (limiting-subunit estimate)\n"
    + "#" * 78
    + "\n# Each complex's assembled copy number is set by its lowest stoichiometry-"
    + "\n# adjusted subunit (excess subunits are degraded). rel* are per-organism"
    + "\n# mean-normalized; the syn3A/syn1 ratio is the cross-organism comparison."
    + "\n# RNAP is core (no sigma factors). Linked to syn3A's longer cell cycle"
    + "\n# (105 vs 60 min): RNAP = transcription engine, degradosome = RNA turnover.\n\n"
    + complex_tbl[_CX_COLS].to_string(index=False, float_format=lambda v: f"{v:.4f}")
    + "\n\n  Per-subunit values (transparency):"
    + "\n".join(_cx_detail) + "\n")


# ── Outlier report (Compare_Ptn.txt) ─────────────────────────────────────────
IPM_SHOW = ["locus_syn1", "gene_name", "essentiality", "gene_product_disp",
            "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change"]
report = ["#" * 78, "\n# PROTEIN (iPM) FOLD-CHANGE AND ABSOLUTE-CHANGE OUTLIERS (coding genes)\n",
          "#" * 78, "\n",
          "# iPM = proteomics; rel* = mean-normalized. Source: " + IN_CSV + "\n"]
report.append(_outlier_section(coding, "iPM fold change", "iPM_fold_change", True,  IPM_SHOW))
report.append(_outlier_section(coding, "iPM absolute change (relIPM delta)", "iPM_abs_change", False, IPM_SHOW))

# PTR (protein-to-transcript ratio) outliers — r-proteins excluded (unreliable iPM).
PTR_SHOW = ["locus_syn1", "gene_name", "Primary Function", "Secondary Function",
            "Tertiary Function", "essentiality", "gene_product_disp",
            "TPM_fold_change", "iPM_fold_change", "PTR_fold_change"]
report.append(_outlier_section(coding[~coding["is_rprotein"].fillna(False)],
                               "PTR fold change = iPM_FC / TPM_FC  (r-proteins excluded)",
                               "PTR_fold_change", True, PTR_SHOW))

RPROT_COLS = ["locus_syn1", "gene_name", "essentiality", "gene_product_disp",
              "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change", "TPM_abs_change",
              "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change"]
rprot_tbl = coding[coding["is_rprotein"]].sort_values("locus_num")
report += ["\n\n", "#" * 78,
           f"\n# RIBOSOMAL PROTEIN GENES (n={len(rprot_tbl)}) -- TPM + iPM, ascending locus\n",
           "#" * 78, "\n",
           rprot_tbl[RPROT_COLS].to_string(index=False) + "\n"]

report.append(COMPLEX_REPORT_TEXT)
with open(OUT_REPORT, "w") as fh:
    fh.write("".join(report))
print(f"\nSaved: {OUT_REPORT}")
for title, col, log in [("iPM fold change", "iPM_fold_change", True),
                        ("iPM abs change",  "iPM_abs_change",  False)]:
    up, dn, lo, hi = _tukey_outliers(coding, col, log)
    print(f"  [{title}] high={len(up)}, low={len(dn)}  (fences {lo:.3g} .. {hi:.3g})")
