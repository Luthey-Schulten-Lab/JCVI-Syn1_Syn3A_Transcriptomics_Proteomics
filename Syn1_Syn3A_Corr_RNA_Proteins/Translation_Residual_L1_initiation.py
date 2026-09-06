"""
Residual analysis of differential translation — Level 1: Translation Initiation Rate

Part of the effort to explain the Pearson r ≈ 0.6 between transcriptome
(avg_sense_TPM) and proteome (iPM_mean) in JCVI-Syn1.0. Level 2 (elongation, via
codon adaptation) lives in Translation_Residual_L2_elongation.py.

Goal
----
Quantify per-gene translation initiation efficiency via thermodynamic RBS modelling
(OSTIR), then test how much of the log10(protein) ~ log10(mRNA) residual is
explained by initiation rate.

Inputs
------
- Genome FASTA : ../Genomes_Input/syn1_genome.fasta  (loaded as one circular string)
- Isoforms     : ../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv
                 Filtered to n_reads >= MIN_READS (= 10).
- Omics table  : ./syn1_genes_transcriptomics_proteomics.csv
                 Keep rna_type == "mRNA" with avg_sense_TPM > 0 and iPM_mean > 0.

Pipeline
--------
1. For every mRNA gene, select PacBio isoforms whose transcribed span (pos5p0..pos3p0)
   fully covers the ORF (start codon to stop codon), using strand-aware half-open
   coordinate logic.
2. For each covering isoform, carve the initiation window from the genome:
      [start_codon − min(UTR_WINDOW, actual 5' UTR) , start_codon + CDS_WINDOW]
   with UTR_WINDOW = CDS_WINDOW = 30 nt. Circular wraparound is handled via
   `_fetch_circular`; dnaA (MMSYN1_0001) at position 0 gets its effective 5'
   position shifted by −UTR_WINDOW so the upstream bases wrap from the genome end.
3. Feed each window to OSTIR with aSD = "ACCUCCUUU" (Syn1 16S rRNA 3' tail).
   Select the hit at the annotated start codon (1-based index = utr_len + 1),
   falling back to the nearest start_position if absent.
4. Per isoform: record TIR (`expression`), dG_total, dG_rRNA:mRNA, dG_standby,
   dG_spacing, dG_mRNA, dG_start_codon, RBS_distance_bp, utr_len.
   Per gene: read-count–weighted mean TIR across valid isoforms and the median
   of each dG sub-feature (restricted to isoforms with a valid TIR so the median
   set matches the weighted-average TIR set).

Statistics
----------
- 5' UTR length distribution across isoforms covering a start codon.
- Baseline OLS   : log10(iPM) ~ log10(TPM)                → r_base, R²_base
- Augmented OLS  : log10(iPM) ~ log10(TPM) + log10(TIR)   via sklearn
                   → R²_full, ΔR² = R²_full − R²_base
- Pearson r between log10(TIR) and the baseline residual.
- Pearson r between each dG sub-feature (median) and the residual.

Outputs (under ./residual_analysis/)
------------------------------------
- isoform_TIR.csv                   per (gene, isoform) OSTIR features
- gene_TIR.csv                      per-gene weighted TIR + median sub-features
- gene_TIR_omics_merged.csv         omics + TIR + residuals for mRNAs with valid TIR
- utr_length_distribution.pdf       histogram of UTR lengths > 0
- TIR_residual_analysis.pdf         (A) mRNA vs protein coloured by TIR with OLS line
                                    (B) residuals vs log10(TIR) with regression line
- TIR_subfeatures_residuals.pdf     4-panel: each dG sub-feature vs residual

Conventions
-----------
- Genome: JCVI-syn1.0 (CP002027.1, 1,078,809 bp, circular).
- Genetic code: Mycoplasma (UGA = Trp); relevant to OSTIR's start-codon handling.
- Coordinates: gene start0/end0 are 0-based half-open; isoform pos5p0/pos3p0 are
  inclusive 0-based transcribed ends.
"""

# =============================================================================
# 1. Imports
# =============================================================================

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import pearsonr
from scipy import stats
from sklearn.linear_model import LinearRegression

from ostir import run_ostir

# =============================================================================
# 2. Paths
# =============================================================================

HOME_DIR         = ".."
GENOME_FASTA     = HOME_DIR + "/Genomes_Input/syn1_genome.fasta"
ISOFORMS_TSV     = HOME_DIR + "/Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
OMICS_CSV        = "./syn1_genes_transcriptomics_proteomics.csv"
OUT_DIR          = "./residual_analysis"

os.makedirs(OUT_DIR, exist_ok=True)

# Cache: skip OSTIR if per-isoform and per-gene TIR tables already exist.
ISO_TIR_CSV_PATH  = f"{OUT_DIR}/isoform_TIR.csv"
GENE_TIR_CSV_PATH = f"{OUT_DIR}/gene_TIR.csv"
SKIP_OSTIR = os.path.exists(ISO_TIR_CSV_PATH) and os.path.exists(GENE_TIR_CSV_PATH)
if SKIP_OSTIR:
    print(f"Cached TIR files found ({ISO_TIR_CSV_PATH}, {GENE_TIR_CSV_PATH}) — "
          "skipping sections 5–7 (genome load, isoform match, OSTIR).")

# OSTIR parameters (Syn1-specific)
ANTI_SD          = "ACCUCCUUU"   # Syn1 16S rRNA 3' tail
UTR_WINDOW       = 30            # nt upstream of start codon to include
CDS_WINDOW       = 30            # nt downstream of start codon to include

# Isoform filtering
MIN_READS        = 10            # minimum PacBio read count per isoform

# Transcriptome column to use for mRNA abundance
TRANSCRIPTOME_COL = "avg_sense_TPM"

# =============================================================================
# 3. Helper: load genome as a single string
# =============================================================================

def load_genome(fasta_path):
    """Return genome sequence as a single uppercase string (no newlines)."""
    seq_parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


_RC_TABLE = str.maketrans("ACGTN", "TGCAN")

def revcomp(seq):
    """Return reverse complement of a DNA sequence."""
    return seq.translate(_RC_TABLE)[::-1]


# =============================================================================
# 4. Helper: extract translation initiation sequence for one gene / isoform
# =============================================================================

def _fetch_circular(genome, g_start, g_end):
    """
    Fetch genome[g_start:g_end] with circular wraparound.
    Handles g_start < 0 (wraps from genome end) and
    g_end > genome_len (wraps to genome start).
    """
    genome_len = len(genome)
    g_start = g_start % genome_len
    g_end   = g_end   % genome_len if g_end % genome_len != 0 else genome_len

    if g_start < g_end:
        return genome[g_start:g_end]
    else:
        # wraparound: straddles the origin
        return genome[g_start:] + genome[:g_end]


def extract_initiation_seq(genome, gene_start0, gene_end0, strand,
                            iso_pos5p0, utr_window=UTR_WINDOW, cds_window=CDS_WINDOW):
    """
    Extract the translation initiation window from the genome.

    For + strand: start codon begins at gene_start0.
      mRNA 5' end is at iso_pos5p0 (genomic, 0-based).
      UTR length = min(utr_window, gene_start0 - iso_pos5p0)
      Sequence extracted: genome[start0 - utr_len : start0 + cds_window]

    For - strand: start codon begins at gene_end0 - 3 (genomic).
      mRNA 5' end is at iso_pos5p0 (genomic, higher coordinate for minus strand).
      UTR length = min(utr_window, iso_pos5p0 - gene_end0)
      Extracted region (genomic): genome[gene_end0 - cds_window : gene_end0 + utr_len]
      Then take reverse complement.

    Handles circular genome wraparound at both ends.

    Returns
    -------
    seq : str  — DNA sequence (sense strand, 5'→3'), start codon at position utr_len
    utr_len : int — actual 5' UTR length in this sequence
    """
    if strand == "+":
        utr_len = min(utr_window, gene_start0 - iso_pos5p0)
        utr_len = max(0, utr_len)
        g_start = gene_start0 - utr_len
        g_end   = gene_start0 + cds_window
        seq = _fetch_circular(genome, g_start, g_end)
    else:
        # minus strand: iso_pos5p0 is the inclusive genomic 5' end (higher coordinate).
        # UTR bases span [gene_end0, iso_pos5p0] inclusive → length = iso_pos5p0 - gene_end0 + 1
        utr_len = min(utr_window, iso_pos5p0 - gene_end0 + 1)
        utr_len = max(0, utr_len)
        g_start = gene_end0 - cds_window
        g_end   = gene_end0 + utr_len
        seq = revcomp(_fetch_circular(genome, g_start, g_end))

    return seq, utr_len


# =============================================================================
# 5. Load data
# =============================================================================

print("Loading transcriptomics + proteomics table ...")
omics = pd.read_csv(OMICS_CSV)

mrna_genes = omics[
    (omics["rna_type"] == "mRNA") &
    (omics[TRANSCRIPTOME_COL] > 0) &
    (omics["iPM_mean"].notna()) &
    (omics["iPM_mean"] > 0)
].copy()
print(f"  mRNA genes with both TPM and iPM > 0: {len(mrna_genes)}")

if not SKIP_OSTIR:
    print("Loading genome ...")
    genome = load_genome(GENOME_FASTA)
    print(f"  Genome length: {len(genome):,} bp")

    print("Loading isoforms ...")
    isoforms = pd.read_csv(ISOFORMS_TSV, sep="\t")
    # Keep only sense isoforms with sufficient read support
    isoforms = isoforms[
        (isoforms["n_reads"] >= MIN_READS)
    ].copy()
    print(f"  Isoforms with >= {MIN_READS} reads: {len(isoforms):,}")

# =============================================================================
# 6. For each gene: find isoforms that cover the start codon
# =============================================================================

def isoforms_covering_orf(gene_row, iso_df):
    """
    Return subset of iso_df that span the entire ORF (start codon to stop codon).

    For + strand: iso.pos5p0 <= gene.start0  AND  iso.pos3p0 >= gene.end0
    For - strand: iso.pos5p0 >= gene.end0    AND  iso.pos3p0 <= gene.start0

    This ensures the isoform provides both the 5' UTR context for TIR calculation
    and confirms the full coding sequence is transcribed.
    """
    strand = gene_row["strand"]
    chrom  = gene_row["chrom"]

    sub = iso_df[(iso_df["chrom"] == chrom) & (iso_df["strand"] == strand)]

    # pos5p0/pos3p0 are inclusive 0-based coordinates of the transcribed ends.
    # gene start0/end0 are 0-based half-open, so the last transcribed base of
    # the ORF on + strand is end0-1, and on - strand the first transcribed base
    # (start codon A) is end0-1.
    if strand == "+":
        mask = (sub["pos5p0"] <= gene_row["start0"]) & (sub["pos3p0"] >= gene_row["end0"] - 1)
    else:
        mask = (sub["pos5p0"] >= gene_row["end0"] - 1) & (sub["pos3p0"] <= gene_row["start0"])

    return sub[mask]


# =============================================================================
# 7. Run OSTIR for each gene
# =============================================================================

if SKIP_OSTIR:
    iso_tir_df  = pd.read_csv(ISO_TIR_CSV_PATH)
    gene_tir_df = pd.read_csv(GENE_TIR_CSV_PATH)
    print(f"Loaded cached TIR tables: "
          f"{len(iso_tir_df)} isoform rows, {len(gene_tir_df)} gene rows")
else:
    print("\nRunning OSTIR for each gene ...")

    isoform_tir_records = []   # one row per (gene, isoform)
    gene_tir_records    = []   # one row per gene (aggregated)

    for _, gene in mrna_genes.iterrows():
        covering = isoforms_covering_orf(gene, isoforms)

        if covering.empty:
            gene_tir_records.append({
                "locus_tag":       gene["locus_tag"],
                "gene_name":       gene["gene_name"],
                "n_covering_iso":  0,
                "TIR_weighted":    np.nan,
                "dG_mRNA_rRNA_med": np.nan,
                "dG_standby_med":  np.nan,
                "dG_spacing_med":  np.nan,
                "dG_mRNA_med":     np.nan,
            })
            print(f"  {gene['locus_tag']} ({gene['gene_name']}): no isoforms covering ORF")
            continue

        print(f"  {gene['locus_tag']} ({gene['gene_name']}): {len(covering)} covering isoforms, running OSTIR ...")

        iso_tirs = []

        for _, iso in covering.iterrows():
            # Special case: dnaA (MMSYN1_0001) starts at position 0 on a circular genome.
            # Its isoforms also begin at pos5p0 = 0, yielding zero UTR by arithmetic.
            # Force the upstream window to wrap around from the genome end by shifting
            # the effective 5' position UTR_WINDOW bases upstream (negative → circular wrap).
            if gene["locus_tag"] == "MMSYN1_0001":
                effective_pos5p0 = iso["pos5p0"] - UTR_WINDOW
            else:
                effective_pos5p0 = iso["pos5p0"]

            seq, utr_len = extract_initiation_seq(
                genome,
                gene["start0"], gene["end0"],
                gene["strand"],
                effective_pos5p0,
            )

            # start codon A is at 0-based index utr_len → 1-based position for OSTIR
            start_codon_1based = utr_len + 1

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ostir_out = run_ostir(
                        seq,
                        aSD=ANTI_SD,
                        threads=16,
                    )
            except Exception as e:
                print(f"  OSTIR error ({gene['locus_tag']}): {e}")
                continue

            if not ostir_out or len(ostir_out) == 0:
                continue

            # Pick the hit at the annotated start position; fall back to nearest
            best = None
            for hit in ostir_out:
                if hit.get("start_position") == start_codon_1based:
                    best = hit
                    break
            if best is None:
                best = min(ostir_out,
                           key=lambda x: abs(x.get("start_position", 0) - start_codon_1based))

            # OSTIR column names: expression = TIR, dG_rRNA:mRNA = SD strength
            tir_val = float(best.get("expression", np.nan))

            isoform_tir_records.append({
                "locus_tag":      gene["locus_tag"],
                "isoform_id":     iso["isoform_id"],
                "n_reads":        iso["n_reads"],
                "utr_len":        utr_len,
                "start_position": best.get("start_position", np.nan),
                "TIR":            tir_val,
                "dG_total":       float(best.get("dG_total",       np.nan)),
                "dG_mRNA_rRNA":   float(best.get("dG_rRNA:mRNA",   np.nan)),
                "dG_standby":     float(best.get("dG_standby",     np.nan)),
                "dG_spacing":     float(best.get("dG_spacing",     np.nan)),
                "dG_mRNA":        float(best.get("dG_mRNA",        np.nan)),
                "dG_start_codon": float(best.get("dG_start_codon", np.nan)),
                "RBS_distance_bp": best.get("RBS_distance_bp",    np.nan),
            })
            iso_tirs.append((iso["n_reads"], tir_val, best))

        if not iso_tirs:
            gene_tir_records.append({
                "locus_tag":        gene["locus_tag"],
                "gene_name":        gene["gene_name"],
                "n_covering_iso":   len(covering),
                "TIR_weighted":     np.nan,
                "dG_mRNA_rRNA_med": np.nan,
                "dG_standby_med":   np.nan,
                "dG_spacing_med":   np.nan,
                "dG_mRNA_med":      np.nan,
            })
            print(f"  {gene['locus_tag']} ({gene['gene_name']}): no valid OSTIR hits for covering isoforms")
            continue

        # Weighted average TIR
        reads  = np.array([x[0] for x in iso_tirs], dtype=float)
        tirs   = np.array([x[1] for x in iso_tirs], dtype=float)
        valid  = ~np.isnan(tirs)
        if valid.sum() == 0:
            tir_w = np.nan
        else:
            tir_w = np.average(tirs[valid], weights=reads[valid])

        # Median sub-features — restrict to isoforms that produced a valid TIR
        # so the median set matches the weighted-average TIR set
        rows_df = pd.DataFrame([
            {
                "dG_mRNA_rRNA": float(x[2].get("dG_rRNA:mRNA", np.nan)),
                "dG_standby":   float(x[2].get("dG_standby",   np.nan)),
                "dG_spacing":   float(x[2].get("dG_spacing",   np.nan)),
                "dG_mRNA":      float(x[2].get("dG_mRNA",      np.nan)),
            }
            for x in iso_tirs if not np.isnan(x[1])   # x[1] is tir_val; x[2] is the OSTIR hit dict
        ])

        gene_tir_records.append({
            "locus_tag":        gene["locus_tag"],
            "gene_name":        gene["gene_name"],
            "n_covering_iso":   len(covering),
            "TIR_weighted":     tir_w,
            "dG_mRNA_rRNA_med": rows_df["dG_mRNA_rRNA"].median(),
            "dG_standby_med":   rows_df["dG_standby"].median(),
            "dG_spacing_med":   rows_df["dG_spacing"].median(),
            "dG_mRNA_med":      rows_df["dG_mRNA"].median(),
        })

    # Save per-isoform and per-gene TIR tables
    iso_tir_df  = pd.DataFrame(isoform_tir_records)
    gene_tir_df = pd.DataFrame(gene_tir_records)
    iso_tir_df.to_csv(ISO_TIR_CSV_PATH,  index=False)
    gene_tir_df.to_csv(GENE_TIR_CSV_PATH, index=False)
    print(f"  Genes with TIR computed: {gene_tir_df['TIR_weighted'].notna().sum()} / {len(gene_tir_df)}")

# =============================================================================
# 8. 5' UTR statistics
# =============================================================================

if not iso_tir_df.empty and "utr_len" in iso_tir_df.columns:
    utr_lengths = iso_tir_df["utr_len"]
    n_valid_utr = (utr_lengths > 0).sum()
    print(f"\n5' UTR statistics (per isoform):")
    print(f"  Isoforms with UTR > 0 bp : {n_valid_utr} / {len(utr_lengths)}")
    print(f"  Median UTR length        : {utr_lengths.median():.1f} bp")
    print(f"  Mean UTR length          : {utr_lengths.mean():.1f} bp")

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(utr_lengths[utr_lengths > 0], bins=30, color="#4C8BB5", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("5' UTR length (nt)", fontsize=11)
    ax.set_ylabel("Number of isoforms",  fontsize=11)
    ax.set_title("Distribution of 5' UTR lengths\n(isoforms covering start codon)", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/utr_length_distribution.pdf", dpi=300)
    plt.close(fig)
    print("Saved: utr_length_distribution.pdf")

# =============================================================================
# 9. Merge TIR onto omics table & compute regression residuals
# =============================================================================

# Merge per-gene TIR into the omics table
analysis = mrna_genes.merge(
    gene_tir_df[["locus_tag", "TIR_weighted", "n_covering_iso",
                 "dG_mRNA_rRNA_med", "dG_standby_med", "dG_spacing_med", "dG_mRNA_med"]],
    on="locus_tag", how="left"
)

# Work on genes with both omics and a valid TIR
valid_mask = (
    analysis["TIR_weighted"].notna() &
    (analysis["TIR_weighted"] > 0)
)
df = analysis[valid_mask].copy()
df["log10_TPM"] = np.log10(df[TRANSCRIPTOME_COL])
df["log10_iPM"] = np.log10(df["iPM_mean"])
df["log10_TIR"] = np.log10(df["TIR_weighted"])

print(f"\nGenes for residual analysis: {len(df)}")
if len(df) == 0:
    print("No genes with valid TIR — check OSTIR output above. Exiting residual analysis.")
    import sys; sys.exit(1)

# --- Baseline regression: log10(iPM) ~ log10(TPM) ---
slope, intercept, r_base, p_base, _ = stats.linregress(df["log10_TPM"], df["log10_iPM"])
df["predicted_iPM"] = slope * df["log10_TPM"] + intercept
df["residual"]      = df["log10_iPM"] - df["predicted_iPM"]

r_sq_base = r_base ** 2
print(f"\nBaseline regression (iPM ~ TPM):")
print(f"  Pearson r  = {r_base:.4f}   R² = {r_sq_base:.4f}")

# --- Multiple regression: log10(iPM) ~ log10(TPM) + log10(TIR) ---
X = df[["log10_TPM", "log10_TIR"]].values
y = df["log10_iPM"].values
lr = LinearRegression().fit(X, y)
y_pred_full = lr.predict(X)
ss_res = np.sum((y - y_pred_full) ** 2)
ss_tot = np.sum((y - y.mean()) ** 2)
r_sq_full = 1 - ss_res / ss_tot
delta_r_sq = r_sq_full - r_sq_base

print(f"\nMultiple regression (iPM ~ TPM + TIR):")
print(f"  R²_baseline   = {r_sq_base:.4f}")
print(f"  R²_full       = {r_sq_full:.4f}")
print(f"  ΔR² from TIR  = {delta_r_sq:.4f}  ({delta_r_sq * 100:.1f}% of total variance)")

# --- Residuals vs TIR ---
r_res, p_res = pearsonr(df["log10_TIR"], df["residual"])
print(f"\nResidual vs log10(TIR):")
print(f"  Pearson r = {r_res:.4f}  (p = {p_res:.2e})")

df.to_csv(f"{OUT_DIR}/gene_TIR_omics_merged.csv", index=False)
print("Saved: gene_TIR_omics_merged.csv")

# =============================================================================
# 10. Visualization
# =============================================================================

# --- Figure 1: log10(mRNA) vs log10(protein), colored by TIR ---
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Colormap: TIR on log scale
tir_vals = df["log10_TIR"]
norm = mcolors.Normalize(vmin=tir_vals.quantile(0.05), vmax=tir_vals.quantile(0.95))
cmap = plt.cm.RdYlBu_r   # red = high TIR, blue = low TIR

# Panel A: mRNA vs protein colored by TIR
ax = axes[0]
sc = ax.scatter(
    df["log10_TPM"], df["log10_iPM"],
    c=tir_vals, cmap=cmap, norm=norm,
    s=30, alpha=0.7, edgecolors="white", linewidths=0.3, zorder=2
)
x_line = np.linspace(df["log10_TPM"].min(), df["log10_TPM"].max(), 100)
ax.plot(x_line, slope * x_line + intercept,
        color="#444444", linewidth=1.5, alpha=0.7, zorder=3, label="OLS fit")
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label("log₁₀(TIR)", fontsize=10)
ax.set_xlabel("log₁₀(mRNA, avg_sense_TPM)", fontsize=11)
ax.set_ylabel("log₁₀(protein, iPM)", fontsize=11)
ax.set_title(f"mRNA vs Protein (colored by TIR)\nr = {r_base:.3f},  R² = {r_sq_base:.3f}", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)

# Panel B: Residuals vs log10(TIR)
ax = axes[1]
res_slope, res_intercept, _, _, _ = stats.linregress(df["log10_TIR"], df["residual"])
x_res = np.linspace(tir_vals.min(), tir_vals.max(), 100)

ax.axhline(0, color="#999999", linewidth=1, linestyle="--", zorder=1)
ax.scatter(
    df["log10_TIR"], df["residual"],
    c=tir_vals, cmap=cmap, norm=norm,
    s=30, alpha=0.7, edgecolors="white", linewidths=0.3, zorder=2
)
ax.plot(x_res, res_slope * x_res + res_intercept,
        color="#444444", linewidth=1.5, alpha=0.7, zorder=3)
ax.set_xlabel("log₁₀(TIR, read-count weighted)", fontsize=11)
ax.set_ylabel("Residual  [log₁₀(iPM) − predicted]", fontsize=11)
ax.set_title(
    f"Residuals vs Translation Initiation Rate\n"
    f"r = {r_res:.3f}  (p = {p_res:.1e})   ΔR² = {delta_r_sq:.3f}",
    fontsize=11
)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(f"{OUT_DIR}/TIR_residual_analysis.pdf", dpi=300)
plt.close(fig)
print("Saved: TIR_residual_analysis.pdf")

# Standalone TIR vs residual panel (SI) — same format as CAI_vs_residual
fig, ax = plt.subplots(figsize=(3, 3))
ax.scatter(df["log10_TIR"], df["residual"], s=12, alpha=0.55,
           color="#4C8BB5", edgecolors="none")
m, b = np.polyfit(df["log10_TIR"], df["residual"], 1)
xs = np.linspace(df["log10_TIR"].min(), df["log10_TIR"].max(), 50)
ax.plot(xs, m * xs + b, color="#c0392b", linewidth=1.2)
ax.set_xlabel("log₁₀(TIR)", fontsize=10)
ax.set_ylabel("log₁₀(iPM) − fit(log₁₀TPM)", fontsize=10)
ax.set_title(f"\nr = {r_res:.3f}  p = {p_res:.1e}  n = {len(df)}", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/TIR_vs_residual.pdf", dpi=300)
plt.close(fig)
print("Saved: TIR_vs_residual.pdf")

# --- Figure 2: Feature breakdown --- dG sub-components vs residuals ---
feature_cols = [
    ("dG_mRNA_rRNA_med", "Shine–Dalgarno ΔG (kcal/mol)"),
    ("dG_standby_med",   "Standby site ΔG (kcal/mol)"),
    ("dG_spacing_med",   "RBS spacing ΔG (kcal/mol)"),
    ("dG_mRNA_med",      "5' mRNA folding ΔG (kcal/mol)"),
]
valid_feat = df.dropna(subset=[c for c, _ in feature_cols])

if len(valid_feat) > 5:
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    for ax, (col, label) in zip(axes, feature_cols):
        r_f, p_f = pearsonr(valid_feat[col], valid_feat["residual"])
        ax.scatter(valid_feat[col], valid_feat["residual"],
                   s=20, alpha=0.6, color="#4C8BB5",
                   edgecolors="white", linewidths=0.3)
        m, b, *_ = stats.linregress(valid_feat[col], valid_feat["residual"])
        xs = np.linspace(valid_feat[col].min(), valid_feat[col].max(), 100)
        ax.plot(xs, m * xs + b, color="#444444", linewidth=1.5, alpha=0.7)
        ax.axhline(0, color="#999999", linewidth=1, linestyle="--")
        ax.set_xlabel(label, fontsize=9)
        ax.set_ylabel("Residual" if ax == axes[0] else "", fontsize=9)
        ax.set_title(f"r = {r_f:.3f}  (p = {p_f:.1e})", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("OSTIR sub-features vs translation residuals", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/TIR_subfeatures_residuals.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved: TIR_subfeatures_residuals.pdf")

print("\nLevel 1 analysis complete.")
print(f"Outputs in: {OUT_DIR}/")
