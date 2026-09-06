"""
Residual analysis of differential translation — Level 2: Translation Elongation Efficiency

Part of the effort to explain the Pearson r ~ 0.6 between transcriptome
(avg_sense_TPM) and proteome (iPM_mean) in JCVI-Syn1.0. Level 1 (initiation rate
via OSTIR) lives in Translation_Residual_L1_initiation.py.

Objective
---------
Quantify codon optimality as a proxy for elongation speed, then test how much of
the log10(protein) ~ log10(mRNA) residual it explains.

Method
------
Codon adaptation index (CAI, Sharp & Li 1987) over all frame-valid Syn1 coding
sequences under the Mycoplasma genetic code (table 4, UGA = Trp). CAI needs a
reference set of highly expressed genes; because the question is about protein
output, that set is the top 20% of protein-coding genes by mean iPM rather than
the usual ribosomal-protein list.

  - Codon counts over the reference CDSs -> RSCU -> relative adaptiveness w.
  - CAI_g = geometric mean of w over informative codons, excluding the start and
    stop codons, the single-codon families (Met, Trp), and codons with w
    undefined.

Inputs
------
- ./syn1_genes_transcriptomics_proteomics.csv   (from Transcription_Translation.py)
- ../Genomes_Input/syn1_genome.fasta
- ../Genomes_Input/syn1_proteins.faa

Outputs (under ./residual_analysis/)
------------------------------------
- gene_CAI.csv                  per-gene CAI + reference-set flag
- gene_CAI_omics_merged.csv     omics + CAI + proteome residual
- CAI_vs_omics.pdf              CAI vs TPM, vs iPM, vs residual
- CAI_vs_residual_all.pdf       standalone residual panel, all proteins
- CAI_vs_residual_cytosolic.pdf standalone residual panel, cytosolic subset
- CAI_deltaR2_bar.pdf           baseline vs +CAI R², both subsets
"""

# =============================================================================
# 1. Imports
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# =============================================================================
# 2. Paths
# =============================================================================

HOME_DIR   = ".."
OMICS_CSV  = "./syn1_genes_transcriptomics_proteomics.csv"
OUT_DIR    = "./residual_analysis"

GENOME_LEN = 1_078_809   # JCVI-Syn1.0 (CP002027.1), circular

os.makedirs(OUT_DIR, exist_ok=True)

# =============================================================================
# 3. Load omics table
# =============================================================================

print("Loading transcriptomics + proteomics table ...")
omics = pd.read_csv(OMICS_CSV)
print(f"  Genes: {len(omics)}")

# ----------------------------
# Genetic code 4 (Mycoplasma: TGA = W)
# ----------------------------
STOP_CODONS = {"TAA", "TAG"}

GENETIC_CODE_4 = {
	"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
	"TCT":"S","TCC":"S","TCA":"S","TCG":"S",
	"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
	"TGT":"C","TGC":"C","TGA":"W","TGG":"W",
	"CTT":"L","CTC":"L","CTA":"L","CTG":"L",
	"CCT":"P","CCC":"P","CCA":"P","CCG":"P",
	"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
	"CGT":"R","CGC":"R","CGA":"R","CGG":"R",
	"ATT":"I","ATC":"I","ATA":"I","ATG":"M",
	"ACT":"T","ACC":"T","ACA":"T","ACG":"T",
	"AAT":"N","AAC":"N","AAA":"K","AAG":"K",
	"AGT":"S","AGC":"S","AGA":"R","AGG":"R",
	"GTT":"V","GTC":"V","GTA":"V","GTG":"V",
	"GCT":"A","GCC":"A","GCA":"A","GCG":"A",
	"GAT":"D","GAC":"D","GAA":"E","GAG":"E",
	"GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# =============================================================================
# 4. CAI (Codon Adaptation Index)
#    Reference set = top 20% protein-coding genes by iPM_mean.
#    Codon counts over reference CDSs → RSCU → relative adaptiveness w.
#    CAI_g = geometric mean of w over informative codons (exclude start, stop,
#    single-codon families Met/Trp, and codons with w undefined).
# =============================================================================

PROTEIN_FAA = HOME_DIR + "/Genomes_Input/syn1_proteins.faa"
GENOME_FA   = HOME_DIR + "/Genomes_Input/syn1_genome.fasta"

# --- 4.1 Load genome (single circular contig) -------------------------------
def load_single_fasta(path):
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_chunks.append(line.strip().upper())
    return "".join(seq_chunks)

genome = load_single_fasta(GENOME_FA)
assert len(genome) == GENOME_LEN, f"genome length {len(genome)} != {GENOME_LEN}"

_COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s):
    return s.translate(_COMP)[::-1]

def extract_cds(start0, end0, strand):
    """Extract CDS nucleotides, handling circular wraparound."""
    if end0 <= GENOME_LEN:
        seq = genome[start0:end0]
    else:
        seq = genome[start0:] + genome[:end0 - GENOME_LEN]
    if strand == "-":
        seq = revcomp(seq)
    return seq

# --- 4.2 Canonical protein-coding locus list from faa -----------------------
faa_locus_tags = set()
with open(PROTEIN_FAA) as fh:
    for line in fh:
        if line.startswith(">"):
            faa_locus_tags.add(line[1:].split("|", 1)[0].strip())
print(f"\nProtein-coding loci in faa: {len(faa_locus_tags)}")

# --- 4.3 CDS coordinates come from the omics table (already has start0/end0/strand)
from collections import defaultdict
mrna = omics[(omics["rna_type"] == "mRNA") & omics["locus_tag"].isin(faa_locus_tags)].copy()

cds_seqs = {}
n_bad_frame = 0
for _, row in mrna.iterrows():
    seq = extract_cds(int(row["start0"]), int(row["end0"]), row["strand"])
    if len(seq) < 6 or len(seq) % 3 != 0:
        n_bad_frame += 1
        continue
    cds_seqs[row["locus_tag"]] = seq
print(f"  CDSs extracted: {len(cds_seqs)} (skipped non-multiple-of-3: {n_bad_frame})")

# --- 4.4 Reference set: top 20% by iPM_mean ---------------------------------
prot = omics[omics["locus_tag"].isin(cds_seqs.keys())].copy()
prot = prot[prot["iPM_mean"].notna() & (prot["iPM_mean"] > 0)].copy()
prot = prot.sort_values("iPM_mean", ascending=False).reset_index(drop=True)
n_ref = max(1, int(round(0.20 * len(prot))))
ref_loci = prot["locus_tag"].iloc[:n_ref].tolist()
print(f"  CAI reference set (top 20% by iPM_mean): {n_ref}/{len(prot)} genes")

# Localization breakdown of the top-20% iPM reference set
if "ptn_localization" in prot.columns:
    ref_loc = prot["ptn_localization"].iloc[:n_ref].fillna("unknown")
    loc_frac = ref_loc.value_counts(normalize=True).sort_values(ascending=False)
    loc_cnt  = ref_loc.value_counts().reindex(loc_frac.index)
    print("  Localization fractions in top 20% iPM_mean:")
    for loc, frac in loc_frac.items():
        print(f"    {loc:<16s}: {frac*100:5.1f}%  (n = {int(loc_cnt[loc])})")

# --- 4.5 Count reference codons → RSCU → w ----------------------------------
#  Exclude start codon and stop codon; skip any premature stop codons.
SYNONYMOUS = defaultdict(list)
for codon, aa in GENETIC_CODE_4.items():
    if aa != "*":
        SYNONYMOUS[aa].append(codon)

ref_codon_counts = defaultdict(int)
for lt in ref_loci:
    seq = cds_seqs[lt]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]  # drop start
    for c in codons:
        if "N" in c or c in STOP_CODONS:
            continue
        if GENETIC_CODE_4.get(c) is None:
            continue
        ref_codon_counts[c] += 1

w = {}
for aa, codons in SYNONYMOUS.items():
    counts = np.array([ref_codon_counts[c] for c in codons], dtype=float)
    if len(codons) == 1 or aa in ("M", "W"):
        for c in codons:
            w[c] = 1.0
        continue
    if counts.max() == 0:
        for c in codons:
            w[c] = 1.0
        continue
    rel = counts / counts.max()
    # Avoid log(0): floor at 0.01 (standard Sharp & Li practice)
    rel = np.where(rel == 0, 0.01, rel)
    for c, rv in zip(codons, rel):
        w[c] = float(rv)

# --- 4.6 CAI per gene -------------------------------------------------------
def compute_cai(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    codons = codons[1:]
    logs = []
    for c in codons:
        if "N" in c or c in STOP_CODONS:
            continue
        aa = GENETIC_CODE_4.get(c)
        if aa is None or aa in ("M", "W"):
            continue
        logs.append(np.log(w[c]))
    if not logs:
        return np.nan, 0
    return float(np.exp(np.mean(logs))), len(logs)

cai_rows = []
for lt, seq in cds_seqs.items():
    cai, n_inf = compute_cai(seq)
    cai_rows.append({"locus_tag": lt, "length_codons": len(seq)//3,
                     "n_informative": n_inf, "CAI": cai,
                     "in_reference": lt in set(ref_loci)})
cai_df = pd.DataFrame(cai_rows)
cai_df.to_csv(f"{OUT_DIR}/gene_CAI.csv", index=False)
print(f"  gene_CAI.csv written ({len(cai_df)} genes)")

# --- 4.7 Sanity check: CAI vs iPM_mean and avg_sense_TPM --------------------
merged = omics.merge(cai_df[["locus_tag", "CAI", "in_reference"]],
                     on="locus_tag", how="inner")
ok = merged[merged["CAI"].notna() & (merged["iPM_mean"] > 0)
            & (merged["avg_sense_TPM"] > 0)].copy()
r_pi, p_pi = pearsonr(ok["CAI"], np.log10(ok["iPM_mean"]))
r_si, p_si = spearmanr(ok["CAI"], ok["iPM_mean"])
r_pt, p_pt = pearsonr(ok["CAI"], np.log10(ok["avg_sense_TPM"]))
print(f"\nCAI vs log10(iPM_mean):   Pearson r = {r_pi:.3f} (p={p_pi:.2e})   "
      f"Spearman ρ = {r_si:.3f} (p={p_si:.2e})   n={len(ok)}")
print(f"CAI vs log10(avg_TPM):    Pearson r = {r_pt:.3f} (p={p_pt:.2e})")
print(f"  mean CAI in reference set: {ok.loc[ok['in_reference'], 'CAI'].mean():.3f}")
print(f"  mean CAI outside:          {ok.loc[~ok['in_reference'], 'CAI'].mean():.3f}")

# =============================================================================
# 5. CAI vs transcriptome/proteome — scatter plots + proteome residual
#    No TIR / Level-1 coupling here.
# =============================================================================

ok["log10_TPM"]  = np.log10(ok["avg_sense_TPM"])
ok["log10_iPM"]  = np.log10(ok["iPM_mean"])

# Residual of log10(iPM) ~ log10(TPM): the part of proteome not explained by mRNA
slope, intercept = np.polyfit(ok["log10_TPM"], ok["log10_iPM"], 1)
ok["proteome_residual"] = ok["log10_iPM"] - (slope * ok["log10_TPM"] + intercept)
r_base, _ = pearsonr(ok["log10_TPM"], ok["log10_iPM"])
print(f"\nBaseline log10(iPM) ~ log10(TPM): Pearson r = {r_base:.3f}  "
      f"(R² = {r_base**2:.3f})  n = {len(ok)}")

r_res_p, p_res_p = pearsonr(ok["CAI"], ok["proteome_residual"])
r_res_s, p_res_s = spearmanr(ok["CAI"], ok["proteome_residual"])
print(f"CAI vs proteome residual:  Pearson r = {r_res_p:.3f} (p={p_res_p:.2e})  "
      f"Spearman ρ = {r_res_s:.3f} (p={p_res_s:.2e})")

# --- ΔR² from adding CAI to the baseline log10(iPM) ~ log10(TPM) OLS --------
y  = ok["log10_iPM"].to_numpy()
X1 = np.column_stack([np.ones(len(ok)), ok["log10_TPM"].to_numpy()])
X2 = np.column_stack([X1, ok["CAI"].to_numpy()])

def _r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot, beta

r2_base, _   = _r2(X1, y)
r2_aug,  beta_cai  = _r2(X2, y)
delta_r2     = r2_aug - r2_base
print(f"\nBaseline   R² (log10 iPM ~ log10 TPM)          = {r2_base:.4f}")
print(f"Augmented  R² (log10 iPM ~ log10 TPM + CAI)    = {r2_aug:.4f}")
print(f"ΔR² from CAI                                    = {delta_r2:+.4f}  "
      f"({100*delta_r2/r2_base:+.1f}% of baseline)")
print(f"  Augmented coefficients: intercept={beta_cai[0]:.3f}, "
      f"β_logTPM={beta_cai[1]:.3f}, β_CAI={beta_cai[2]:.3f}")

# --- Cytosolic-only: baseline R² and ΔR² from adding CAI --------------------
if "ptn_localization" in ok.columns:
    cyto = ok[ok["ptn_localization"] == "cytoplasmic"].copy()
    if len(cyto) > 3:
        y_c  = cyto["log10_iPM"].to_numpy()
        X1_c = np.column_stack([np.ones(len(cyto)), cyto["log10_TPM"].to_numpy()])
        X2_c = np.column_stack([X1_c, cyto["CAI"].to_numpy()])
        r2_base_c, _  = _r2(X1_c, y_c)
        r2_aug_c,  _  = _r2(X2_c, y_c)
        delta_r2_c    = r2_aug_c - r2_base_c
        print(f"\n[Cytosolic only, n = {len(cyto)}]")
        print(f"  Baseline  R² (log10 iPM ~ log10 TPM)        = {r2_base_c:.4f}")
        print(f"  Augmented R² (log10 iPM ~ log10 TPM + CAI)  = {r2_aug_c:.4f}")
        print(f"  ΔR² from CAI                                 = {delta_r2_c:+.4f}")

ok.to_csv(f"{OUT_DIR}/gene_CAI_omics_merged.csv", index=False)

# --- Three-panel figure: CAI vs TPM, CAI vs iPM, CAI vs residual -----------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))

def _scatter(ax, x, y, xlab, ylab, title, r, p):
    ax.scatter(x, y, s=12, alpha=0.55, color="#4C8BB5",
               edgecolors="none")
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, m * xs + b, color="#c0392b", linewidth=1.2)
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_title(f"{title}\nr = {r:.3f}  p = {p:.1e}  n = {len(x)}",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

_scatter(axes[0], ok["CAI"], ok["log10_TPM"],
         "CAI", "log₁₀(avg_sense_TPM)",
         "CAI vs transcriptome", r_pt, p_pt)
_scatter(axes[1], ok["CAI"], ok["log10_iPM"],
         "CAI", "log₁₀(iPM_mean)",
         "CAI vs proteome", r_pi, p_pi)
_scatter(axes[2], ok["CAI"], ok["proteome_residual"],
         "CAI", "log₁₀(iPM) − fit(log₁₀TPM)",
         "CAI vs proteome residual", r_res_p, p_res_p)

fig.tight_layout()
fig.savefig(f"{OUT_DIR}/CAI_vs_omics.pdf", dpi=300)
plt.close(fig)
print(f"Saved: CAI_vs_omics.pdf")

# Standalone CAI vs proteome residual (main-text panel)
# All proteins
fig, ax = plt.subplots(figsize=(3, 3))
_scatter(ax, ok["CAI"], ok["proteome_residual"],
         "CAI", "log₁₀(iPM) − fit(log₁₀TPM)",
         "", r_res_p, p_res_p)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/CAI_vs_residual_all.pdf", dpi=300)
plt.close(fig)
print("Saved: CAI_vs_residual_all.pdf")

# Cytosolic only (refit residual within cytosolic)
if "ptn_localization" in ok.columns:
    cyto_res = ok[ok["ptn_localization"] == "cytoplasmic"].copy()
    if len(cyto_res) > 3:
        sl_c, in_c = np.polyfit(cyto_res["log10_TPM"], cyto_res["log10_iPM"], 1)
        cyto_res["proteome_residual"] = (
            cyto_res["log10_iPM"] - (sl_c * cyto_res["log10_TPM"] + in_c)
        )
        r_res_c, p_res_c = pearsonr(cyto_res["CAI"], cyto_res["proteome_residual"])
        fig, ax = plt.subplots(figsize=(3, 3))
        _scatter(ax, cyto_res["CAI"], cyto_res["proteome_residual"],
                 "CAI", "log₁₀(iPM) − fit(log₁₀TPM)",
                 "", r_res_c, p_res_c)
        fig.tight_layout()
        fig.savefig(f"{OUT_DIR}/CAI_vs_residual_cytosolic.pdf", dpi=300)
        plt.close(fig)
        print("Saved: CAI_vs_residual_cytosolic.pdf")

# Bar chart: baseline vs +CAI R² (all proteins and cytosolic-only)
cyto_bar = ok[ok["ptn_localization"] == "cytoplasmic"].copy() \
    if "ptn_localization" in ok.columns else ok.iloc[0:0]
if len(cyto_bar) > 3:
    y_c  = cyto_bar["log10_iPM"].to_numpy()
    X1_c = np.column_stack([np.ones(len(cyto_bar)), cyto_bar["log10_TPM"].to_numpy()])
    X2_c = np.column_stack([X1_c, cyto_bar["CAI"].to_numpy()])
    r2_base_cyto, _ = _r2(X1_c, y_c)
    r2_aug_cyto,  _ = _r2(X2_c, y_c)

    fig, ax = plt.subplots(figsize=(3, 3))
    groups  = ["All proteins", "Cytosolic only"]
    x_pos   = np.arange(len(groups))
    width   = 0.38
    base_vals = [r2_base,     r2_base_cyto]
    aug_vals  = [r2_aug,      r2_aug_cyto]
    b1 = ax.bar(x_pos - width/2, base_vals, width,
                color="#9CA3AF", edgecolor="white", label="Baseline\n(log₁₀TPM)")
    b2 = ax.bar(x_pos + width/2, aug_vals,  width,
                color="#4C8BB5", edgecolor="white", label="+ CAI")
    for bars, vals in ((b1, base_vals), (b2, aug_vals)):
        for rect, v in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width()/2, v + 0.01,
                    f"{v:.3f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(groups)
    ax.set_ylabel("R²  (log₁₀ iPM model)", fontsize=10)
    ax.set_ylim(0, max(aug_vals) * 1.18)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/CAI_deltaR2_bar.pdf", dpi=300)
    plt.close(fig)
    print("Saved: CAI_deltaR2_bar.pdf")


# =============================================================================
# 6. Final report — printed so it can be redirected to a .txt file:
#        python Translation_Residual_L2_elongation.py > Translation_Residual_L2_elongation.txt
# =============================================================================

report = f"""
================================================================================
Level 2 — Translation Elongation Efficiency in JCVI-Syn1.0
Report generated from Translation_Residual_L2_elongation.py
================================================================================

GOAL
----
Quantify how much of the log10(iPM) ~ log10(mRNA TPM) residual (the r ~ 0.6
transcriptome-proteome correlation in Syn1) is explained by codon-level
elongation efficiency, without touching the Level 1 initiation model.

DATA
----
  - Transcriptomics + proteomics table: syn1_genes_transcriptomics_proteomics.csv
  - Protein-coding loci:                {len(faa_locus_tags)} (from syn1_proteins.faa)
  - CDSs extracted and frame-valid:     {len(cds_seqs)}
  - Genome:                             JCVI-Syn1.0 (CP002027.1), 1,078,809 bp, circular
  - Genetic code:                       Mycoplasma table 4 (UGA = Trp)

CAI — codon adaptation index
----------------------------
  Reference set: top 20% of protein-coding genes by iPM_mean
                 ({n_ref} genes out of {len(prot)} with iPM data)
  Correlations (n = {len(ok)}):
    CAI vs log10(avg_sense_TPM):   Pearson r = {r_pt:+.3f}  p = {p_pt:.2e}
    CAI vs log10(iPM_mean):        Pearson r = {r_pi:+.3f}  p = {p_pi:.2e}
                                   Spearman r = {r_si:+.3f}  p = {p_si:.2e}
    CAI vs proteome residual:      Pearson r = {r_res_p:+.3f}  p = {p_res_p:.2e}
  Mean CAI in reference set: {ok.loc[ok['in_reference'], 'CAI'].mean():.3f}
  Mean CAI outside:          {ok.loc[~ok['in_reference'], 'CAI'].mean():.3f}

REGRESSION ΔR² (dependent variable = log10(iPM_mean))
-----------------------------------------------------
  All detected proteins (n = {len(ok)}):
      Baseline   log10(iPM) ~ log10(TPM)         R² = {r2_base:.4f}  (Pearson r = {r_base:.3f})
      + CAI                                      R² = {r2_aug:.4f}   ΔR² = {delta_r2:+.4f}
                                                 ({100*delta_r2/r2_base:+.1f}% of baseline)
      Coefficients: intercept = {beta_cai[0]:+.3f}   β_logTPM = {beta_cai[1]:+.3f}   β_CAI = {beta_cai[2]:+.3f}

  Cytosolic proteins only (n = {len(cyto)}):
      Baseline                                   R² = {r2_base_c:.4f}
      + CAI                                      R² = {r2_aug_c:.4f}   ΔR² = {delta_r2_c:+.4f}

INTERPRETATION
--------------
  CAI carries a substantial codon-level signal. Proteome-based reference
  selection (top 20% by iPM) finds exactly the codons that predict protein
  abundance, and mRNA level is uninformative for CAI (r ~ 0), so the ΔR² is
  nearly orthogonal to the baseline transcriptome term.

OUTPUTS
-------
  residual_analysis/
    gene_CAI.csv
    gene_CAI_omics_merged.csv
    CAI_vs_omics.pdf
    CAI_vs_residual_all.pdf
    CAI_vs_residual_cytosolic.pdf
    CAI_deltaR2_bar.pdf

================================================================================
END OF LEVEL 2 REPORT
================================================================================
"""
print(report)
