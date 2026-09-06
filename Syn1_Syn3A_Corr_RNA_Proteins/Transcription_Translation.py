#!/usr/bin/env python
# coding: utf-8
"""
Transcription–Translation Correlation for Syn1
===============================================
Correlates Syn1 proteomics (iPM, intensity-based protein molecules) with
transcriptomics (Illumina / PacBio sense TPM) at the gene level.

Key findings:
  - Pearson r ~ 0.6 between Illumina TPM and iPM.
  - Pearson r ~ 0.5 between PacBio TPM and iPM.
  - Cytosolic-only subset: r ~ 0.67 (slightly higher).

iPM (intensity-based Protein Molecules) is derived from iBAQ and scales
more linearly with true protein copy number than raw iBAQ.
"""

# =============================================================================
# 1. Imports
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text

# =============================================================================
# 2. Paths
# =============================================================================

HOME_DIR = ".."   # project root; this script lives at Syn1_Syn3A_Corr_RNA_Proteins/
PROTEOME_FOLDER   = HOME_DIR + "/Syn1_Syn3A_Proteomics"
TRANSCRIPTOME_CSV = HOME_DIR + "/Syn1_Syn3A_Transcriptomics/syn1_TPM.tsv"
PROTEOME_CSV      = PROTEOME_FOLDER + "/syn1_proteome.tsv"
PROTEOME_SYN3A_CSV = PROTEOME_FOLDER + "/syn3A_proteome.tsv"

# Which transcriptome column to correlate against iPM:
#   "avg_sense_TPM"   → merged Illumina TPM
#   "PacBio_sense_TPM" → PacBio Iso-Seq TPM
TRANSCRIPTOME_COL = "avg_sense_TPM"

# =============================================================================
# 3. Load data
# =============================================================================

# --- 3a. Transcriptomics (Illumina + PacBio TPM profiles) ---
genes = pd.read_csv(TRANSCRIPTOME_CSV, sep="\t")
print(f"Transcriptomics loaded: {len(genes)} genes")

# --- 3b. Proteomics (iPM + localization) ---
prot = pd.read_csv(PROTEOME_CSV, sep="\t")
print(f"Proteomics loaded: {len(prot)} entries")
print(f"  Columns: {prot.columns.tolist()}")
print(f"  Localization categories: {sorted(prot['ptn_localization'].dropna().unique())}")

# --- 3c. Syn3A proteomics (for Syn3A gene_name mapping) ---
prot_syn3a = pd.read_csv(PROTEOME_SYN3A_CSV, sep="\t")
print(f"Syn3A proteomics loaded: {len(prot_syn3a)} entries")

# Map Syn3A gene_name onto Syn1 genes by numeric locus suffix
# (MMSYN1_0001 ↔ JCVISYN3A_0001)
syn3a_name_map = (
    prot_syn3a.assign(locus_num=prot_syn3a['locus_tag'].str.split('_').str[-1])
              .dropna(subset=['gene_name'])
              .drop_duplicates(subset=['locus_num'])
              .set_index('locus_num')['gene_name']
)
genes['gene_name_syn3A'] = (
    genes['locus_tag'].str.split('_').str[-1].map(syn3a_name_map)
                      .fillna('deleted')
)

# =============================================================================
# 4. Merge proteomics onto genes dataframe
# =============================================================================

# Map iPM replicates, mean, CV, and localization onto genes by locus_tag
iPM_cols = ['iPM_rep1', 'iPM_rep2', 'iPM_rep3', 'iPM_mean', 'iPM_CV', 'ptn_localization','ptn_copy_number']
iPM_map  = prot.set_index('locus_tag')[iPM_cols]

for col in iPM_cols:
    genes[col] = genes['locus_tag'].map(iPM_map[col])

# Coverage summary
n_mrna    = (genes['rna_type'] == 'mRNA').sum()
n_detected = genes.loc[genes['rna_type'] == 'mRNA', 'iPM_mean'].notna().sum()
print(f"\nProtein detection coverage (mRNA genes only):")
print(f"  Total mRNA genes : {n_mrna}")
print(f"  With iPM data    : {n_detected} ({n_detected / n_mrna * 100:.1f}%)")
print(f"  Without iPM data : {n_mrna - n_detected}")

# Save merged table
genes.to_csv("syn1_genes_transcriptomics_proteomics.csv", index=False)
print("\nSaved: syn1_genes_transcriptomics_proteomics.csv")

# =============================================================================
# 5. Correlation: iPM vs TPM
# =============================================================================

# Filter to genes with both iPM and TPM > 0 (all localization categories)
loc_categories = set(genes['ptn_localization'].dropna().unique())
mask = (
    (genes['iPM_mean'] > 0) &
    (genes[TRANSCRIPTOME_COL] > 0) &
    genes['ptn_localization'].isin(loc_categories)
)

corr_df = genes.loc[mask].copy()
corr_df['log10_iPM'] = np.log10(corr_df['iPM_mean'])
corr_df['log10_TPM'] = np.log10(corr_df[TRANSCRIPTOME_COL])

rho, p_spearman = spearmanr(corr_df['log10_iPM'], corr_df['log10_TPM'])
r,   p_pearson  = pearsonr(corr_df['log10_iPM'],  corr_df['log10_TPM'])

print(f"\nCorrelation: iPM vs {TRANSCRIPTOME_COL}")
print(f"  Localization filter : {', '.join(sorted(loc_categories))}")
print(f"  N genes             : {len(corr_df)}")
print(f"  Spearman ρ          : {rho:.4f}  (p = {p_spearman:.2e})")
print(f"  Pearson  r          : {r:.4f}  (p = {p_pearson:.2e})")
print(f"  Pearson r²         : {r**2:.4f}")


# =============================================================================
# 6. Scatter plot — log10(iPM) vs log10(TPM)
# =============================================================================
# Two panels side by side (same 6×6 per panel):
#   Left  — all detected proteins
#   Right — cytosolic proteins only (ptn_localization == 'cytoplasmic')

def _scatter_panel(ax, df, rho, r, p_rho, p_r):
    """Draw one correlation scatter panel with OLS fit, colored outliers, and non-overlapping labels."""
    df = df.copy()
    z      = np.polyfit(df['log10_TPM'], df['log10_iPM'], 1)
    x_line = np.linspace(df['log10_TPM'].min(), df['log10_TPM'].max(), 100)

    # Compute residuals to identify outliers
    df['residual'] = df['log10_iPM'] - np.polyval(z, df['log10_TPM'])
    top_high = df.nlargest(5,  'residual').index
    top_low  = df.nsmallest(5, 'residual').index
    outlier_idx = top_high.union(top_low)

    # Assign per-point colors: red = high residual, green = low residual, blue = regular
    colors = pd.Series('#4C8BB5', index=df.index)
    colors[top_high] = '#D94F4F'   # red
    colors[top_low]  = '#3A9E5F'   # green

    ax.scatter(df['log10_TPM'], df['log10_iPM'],
               alpha=0.55, s=25, c=colors,
               edgecolors='white', linewidths=0.3, zorder=2)

    # Re-draw outlier points on top with full opacity and larger size
    ax.scatter(df.loc[top_high, 'log10_TPM'], df.loc[top_high, 'log10_iPM'],
               s=40, color='#D94F4F', edgecolors='white', linewidths=0.4, zorder=4)
    ax.scatter(df.loc[top_low,  'log10_TPM'], df.loc[top_low,  'log10_iPM'],
               s=40, color='#3A9E5F', edgecolors='white', linewidths=0.4, zorder=4)

    ax.plot(x_line, np.polyval(z, x_line),
            color='#555555', linewidth=1.5, alpha=0.7, zorder=3)

    ax.set_xlabel('log$_{10}$(TPM)  —  Transcriptomics', fontsize=12)
    ax.set_ylabel('log$_{10}$(iPM)  —  Proteomics', fontsize=12)

    # Stats annotation box
    ax.text(0.04, 0.97,
            # f'Spearman ρ = {rho:.3f}  (p = {p_rho:.1e})\n'
            f'Pearson  r = {r:.3f}, r² = {r**2:.3f} (p = {p_r:.1e})\n'
            f'n = {len(df)} genes',
            transform=ax.transAxes, fontsize=9.5, verticalalignment='top')
            # bbox=dict(boxstyle='round,pad=0.45', facecolor='#FFF8E7',
                    #   edgecolor='#CCBBAA', alpha=0.85))

    # Annotate outliers — use adjustText to avoid label overlaps
    # Label with the numeric suffix of locus_tag (e.g. MMSYN1_0042 → 0042)
    texts = []
    for idx, row in df.loc[outlier_idx].iterrows():
        locus_num  = row['locus_tag'].split('_')[-1]
        txt_color  = '#D94F4F' if idx in top_high else '#3A9E5F'
        texts.append(ax.text(row['log10_TPM'], row['log10_iPM'], locus_num,
                             fontsize=7.5, color=txt_color, fontweight='bold', zorder=5))

    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.6),
                expand=(1.4, 1.6), force_text=(0.5, 0.8))

    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, alpha=0.25, linestyle='--')
    return df  # return with residual column attached


# --- All proteins ---
rho_all, p_rho_all = spearmanr(corr_df['log10_iPM'], corr_df['log10_TPM'])
r_all,   p_r_all   = pearsonr( corr_df['log10_iPM'], corr_df['log10_TPM'])

# --- Cytosolic proteins only ---
cyto_df = corr_df[corr_df['ptn_localization'] == 'cytoplasmic'].copy()
rho_cyto, p_rho_cyto = spearmanr(cyto_df['log10_iPM'], cyto_df['log10_TPM'])
r_cyto,   p_r_cyto   = pearsonr( cyto_df['log10_iPM'], cyto_df['log10_TPM'])

print(f"\nCorrelation — all proteins    : Spearman ρ={rho_all:.4f}, Pearson r={r_all:.4f}, Pearson r²={r_all**2:.4f}, n={len(corr_df)}")
print(f"Correlation — cytosolic only  : Spearman ρ={rho_cyto:.4f}, Pearson r={r_cyto:.4f}, Pearson r²={r_cyto**2:.4f}, n={len(cyto_df)}")

# All proteins
fig, ax = plt.subplots(figsize=(6, 6))
corr_df = _scatter_panel(ax, corr_df, rho_all, r_all, p_rho_all, p_r_all)
plt.tight_layout()
plt.savefig('./protein_vs_mRNA_correlation_all.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# Cytosolic proteins only
fig, ax = plt.subplots(figsize=(6, 6))
cyto_df = _scatter_panel(ax, cyto_df, rho_cyto, r_cyto, p_rho_cyto, p_r_cyto)
plt.tight_layout()
plt.savefig('./protein_vs_mRNA_correlation_cytosolic.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# =============================================================================
# 7. Residual analysis — genes deviating most from the mRNA–protein fit
# =============================================================================
# Positive residual: more protein than expected from mRNA
#   → high translational efficiency, or stable/long-lived protein
# Negative residual: less protein than expected from mRNA
#   → poor translation, rapid protein degradation, or low proteomics coverage

cols_show = ['locus_tag', 'gene_name', 'gene_name_syn3A', 'gene_product', 'ptn_localization', 'ptn_copy_number', 'log10_TPM', 'log10_iPM', 'residual']

print("\n=== Top 20: MORE protein than expected from mRNA ===")
print("(high translational efficiency or stable protein)\n")
print(corr_df.nlargest(20, 'residual')[cols_show].reset_index(drop=True).to_string())

print("\n=== Top 20: LESS protein than expected from mRNA ===")
print("(poor translation, rapid degradation, or low proteomics coverage)\n")
print(corr_df.nsmallest(20, 'residual')[cols_show].reset_index(drop=True).to_string())

print("\n=== Top 20: LESS protein in CYTOSOLIC than expected from mRNA ===")
print("(poor translation, rapid degradation, or low proteomics coverage)\n")
print(corr_df[corr_df['ptn_localization'] == 'cytoplasmic'].nsmallest(20, 'residual')[cols_show].reset_index(drop=True).to_string())

# =============================================================================
# 8. Protein copy number distribution — linear scale
# =============================================================================
# ptn_copy_number: absolute protein molecules per cell estimated from iPM.
# Excludes Non_proteins (tRNAs, rRNAs, etc.) and missing/zero values.

LOC_CATEGORIES = ['cytoplasmic', 'membrane', 'lipoprotein', 'extracellular']

plot_df = prot[
    prot['ptn_localization'].isin(LOC_CATEGORIES) &
    prot['ptn_copy_number'].notna() &
    (prot['ptn_copy_number'] > 0)
].copy()

median_copy = plot_df['ptn_copy_number'].median()

fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(plot_df['ptn_copy_number'], bins=50, color='steelblue',
        edgecolor='white', linewidth=0.3)
ax.axvline(median_copy, color='red', linewidth=1.5, linestyle='--',
           label=f'Median = {median_copy:.1f}')
ax.set_xlabel('Protein copy number per cell', fontsize=12)
ax.set_ylabel('Number of unique proteins', fontsize=12)
# ax.set_title('Syn1 Protein Copy Number Distribution (linear)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('./protein_copy_number_linear.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# =============================================================================
# 9. Protein copy number distribution — log10 scale
# =============================================================================

median_log10 = np.log10(median_copy)

fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(np.log10(plot_df['ptn_copy_number']), bins=50, color='steelblue',
        edgecolor='white', linewidth=0.3)
ax.axvline(median_log10, color='red', linewidth=1.5, linestyle='--',
           label=f'Median = {median_copy:.1f} ({median_log10:.2f} in log$_{{10}}$)')
ax.set_xlabel('log$_{10}$(protein copy number per cell)', fontsize=12)
ax.set_ylabel('Number of uniqueproteins', fontsize=12)
# ax.set_title('Syn1 Protein Copy Number Distribution (log$_{10}$)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('./protein_copy_number_log10.pdf', dpi=300, bbox_inches='tight')
# plt.show()

# =============================================================================
# 10. Copy number statistics by subcellular localization
# =============================================================================

LOC_COLORS = {
    'cytoplasmic' : '#4C8BB5',
    'membrane'    : '#E05A5A',
    'lipoprotein' : '#E8A838',
    'extracellular': '#7B5EA7',
}

stats_rows = []
for loc in LOC_CATEGORIES:
    subset = plot_df.loc[plot_df['ptn_localization'] == loc, 'ptn_copy_number']
    if subset.empty:
        continue
    stats_rows.append({
        'localization' : loc,
        'n'            : len(subset),
        'min'          : subset.min(),
        'Q1'           : subset.quantile(0.25),
        'median'       : subset.median(),
        'mean'         : subset.mean(),
        'Q3'           : subset.quantile(0.75),
        'max'          : subset.max(),
        'CV'           : subset.std() / subset.mean(),
    })

stats_df = pd.DataFrame(stats_rows)
print("\nProtein copy number statistics by localization (molecules per cell):")
print(stats_df.to_string(index=False, float_format=lambda x: f'{x:.2f}'))

# stats_df.to_csv('./protein_copy_number_stats_by_localization.csv', index=False, float_format='%.4f')
# print("\nSaved: protein_copy_number_stats_by_localization.csv")

# --- Histogram by localization (log10 scale) ---
fig, ax = plt.subplots(figsize=(6, 6))

for loc in LOC_CATEGORIES:
    subset = plot_df.loc[plot_df['ptn_localization'] == loc, 'ptn_copy_number']
    if subset.empty:
        continue
    median_val = subset.median()
    color      = LOC_COLORS[loc]

    ax.hist(np.log10(subset), bins=30, alpha=0.55, color=color,
            edgecolor='white', linewidth=0.3,
            label=f'{loc} (n={len(subset)})')

    # Vertical median line; label at top of axes in data-x / axes-y coords
    ax.axvline(np.log10(median_val), color=color, linewidth=1.5,
               linestyle='--', alpha=0.9)
    ax.text(np.log10(median_val), 1.01, f' {median_val:.0f}',
            transform=ax.get_xaxis_transform(),
            color=color, fontsize=12, va='bottom', ha='left', fontweight='bold')

ax.set_xlabel('log$_{10}$(protein copy number per cell)', fontsize=12)
ax.set_ylabel('Number of unique proteins', fontsize=12)
ax.legend(fontsize=10, framealpha=0.85)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(True, alpha=0.25, linestyle='--')
plt.tight_layout()
plt.savefig('./protein_copy_number_by_localization.pdf', dpi=300, bbox_inches='tight')
# plt.show()
