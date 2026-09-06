#!/usr/bin/env python
# coding: utf-8
"""
R3 syn3A panels — transcriptome-proteome correlation repeated on JCVI-syn3A.

Companion to R3_figure_panels.py in this folder. Builds the syn3A side of
the main correlation figure + its SI, born-at-size, OUTPUT.md default fonts (7 pt
labels / 6 pt ticks). Every panel carries a bold red "Syn3A" tag (the Fig 5/6
organism convention; the syn1 twin is blue). Same analysis as syn1, with two
organism-specific differences:

  * TIR  — regenerated FROM SCRATCH at the gene level. syn3A has no PacBio isoform
           set, so there is no read-weighting over covering isoforms; OSTIR runs
           on a fixed 30-nt window around each annotated start codon.
  * CAI  — recomputed on syn3A (genetic code 4, TGA=Trp) with syn3A's OWN top-20%-
           by-iPM reference set (the reference gene set differs from syn1).

MAIN figure (correlation.pdf), syn3A row:
  g  7/3 x 7/3  Illumina sense TPM vs iPM (log10), OPEN circles by localization
  f  7/3 x 7/6  copy-number distribution by localization (4 merged classes)
SI figure (si-correlation.pdf), syn3A row:
  h  7/4 x 7/4  proteome residual vs log10(TIR)   [gene-level OSTIR]
  i  7/4 x 7/4  proteome residual vs CAI           [syn3A reference set]
  l  7/4 x 7/4  model Pearson R (baseline vs +CAI) for all + cytoplasmic

Run from Syn1_Syn3A_Corr_RNA_Proteins/ in the Omics conda env (needs ostir):
  conda run -n Omics python Corr_RNA_Protein_Syn3A.py
Output: R3_panels_syn3A/panel_{f,g,h,i,l}.pdf + syn3A_genes_transcriptomics_proteomics.csv
        + R3_panels_syn3A/R3_syn3A.txt (+ gene_TIR_syn3A.csv cache)
"""

import os
import warnings
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, gaussian_kde
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size': 7,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Nimbus Sans', 'Helvetica', 'DejaVu Sans'],
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

HALF, QUART = 7 / 2, 7 / 4          # 3.5 in, 1.75 in (SI residual panels)
W    = 7 / 3                        # 2.333 in -- main-panel width
CORR = (W, W)                       # 7/3 x 7/3  (TPM vs iPM)
DIST = (W, 7 / 6)                   # 7/3 x 7/6  (copy-number distribution)
OUT = 'R3_panels_syn3A'
os.makedirs(OUT, exist_ok=True)

PROT_CSV = '../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv'
TPM_TSV  = '../Syn1_Syn3A_Transcriptomics/syn3A_TPM.tsv'
GENOME_FA = '../Genomes_Input/syn3A_genome.fasta'
TIR_CACHE = f'{OUT}/gene_TIR_syn3A.csv'

GENOME_LEN = 543_379   # JCVI-syn3A (CP016816.2), circular

# --- localization: collapse syn3A's finer scheme to syn1's four classes --------
LOC_MAP = {
    'cytoplasm': 'cytoplasmic',
    'peripheral membrane': 'cytoplasmic',
    'unidentified': 'cytoplasmic',
    'trans-membrane': 'membrane',
    'lipoprotein': 'lipoprotein',
    'extracellular': 'extracellular',
}
LOC_ORDER = ['cytoplasmic', 'lipoprotein', 'membrane', 'extracellular']
LOC_LABEL = {'cytoplasmic': 'Cyto', 'lipoprotein': 'Lipo',
             'membrane': 'Mem', 'extracellular': 'Extra'}
LOC_COLORS = {'cytoplasmic': '#0072B2', 'lipoprotein': '#009E73',
              'membrane': '#D55E00', 'extracellular': '#CC79A7'}

# Organism identity colors (match Fig 5/6): Syn1 blue, Syn3A red. Used only as a
# bold left-title tag on every panel -- never for the localization-coloured points.
SYN1_COL, SYN3A_COL = '#3182bd', '#c0392b'
def org_tag(ax, which='syn3a'):
    name, col = ('Syn1.0', SYN1_COL) if which == 'syn1' else ('Syn3A', SYN3A_COL)
    ax.set_title(name, loc='left', color=col, fontweight='bold', fontsize=7)

log = []
def say(s):
    print(s); log.append(s)

def fit_line(ax, x, y, color='crimson'):
    b, a = np.polyfit(x, y, 1)
    xs = np.array([x.min(), x.max()])
    ax.plot(xs, b * xs + a, color=color, lw=1.0, zorder=5)

def r2_lstsq(X, y):
    A = np.column_stack([np.ones(len(y))] + [X[:, i] for i in range(X.shape[1])])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    return 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)

# =============================================================== load & merge
prot = pd.read_csv(PROT_CSV, sep='\t')
# The three syn3A pseudogenes carry MS intensity but no curated protein (no sequence,
# length, MW or localization), so they are dropped from every protein-level statistic
# here. They stay in the transcriptome tables, where all 496 loci are quantified.
n_pseudo = int(prot['primary_function'].eq('Pseudogene').sum())
prot = prot[~prot['primary_function'].eq('Pseudogene')].reset_index(drop=True)
prot['loc4'] = prot['localization'].map(LOC_MAP)
tpm = pd.read_csv(TPM_TSV, sep='\t')

df = prot.merge(
    tpm[['locus_tag', 'start0', 'end0', 'strand', 'gene_len', 'Illumina_sense_TPM']],
    on='locus_tag', how='left')
say(f"syn3A proteins: {len(prot)} ({n_pseudo} pseudogenes dropped); "
    f"merged with TPM coords: {df['start0'].notna().sum()}")
say(f"localization merge -> " + ", ".join(
    f"{k}:{int((df['loc4'] == k).sum())}" for k in LOC_ORDER))

# =============================================================== f  copy-number
cn = df[(df['exp_ptn_cnt_2026'] > 0) & df['loc4'].isin(LOC_ORDER)]
fig, ax = plt.subplots(figsize=DIST, constrained_layout=True)
allv = np.log10(cn['exp_ptn_cnt_2026'].values)
lo, hi = allv.min(), allv.max()
xgrid = np.linspace(lo, hi, 300)
bw = (hi - lo) / 25.0
for loc in LOC_ORDER:
    v = np.log10(cn.loc[cn['loc4'] == loc, 'exp_ptn_cnt_2026'].values)
    if len(v) < 5:
        say(f"f) {loc}: n={len(v)} (<5, curve skipped)")
        continue
    n = len(v); med = float(np.median(v)); med_cn = 10 ** med
    kde = gaussian_kde(v)
    ycurve = kde(xgrid) * n * bw
    ax.plot(xgrid, ycurve, color=LOC_COLORS[loc], lw=1.2,
            label=f"{LOC_LABEL[loc]} ({med_cn:.0f}, n={n})")
    ax.fill_between(xgrid, ycurve, color=LOC_COLORS[loc], alpha=0.13)
    ax.vlines(med, 0, float(kde(med)[0]) * n * bw, color=LOC_COLORS[loc], lw=0.9, ls='--')
    say(f"f) {loc}: n={n} median={med_cn:.1f} copies")
ax.set_xlabel('Copies per cell ($\\log_{10}$)')
ax.set_ylabel('Proteins')
ax.set_xlim(lo, hi)
leg = ax.legend(frameon=False, handlelength=0.8, labelspacing=0.2, borderpad=0.2,
                fontsize=5, loc='upper left',
                title='(median, n=unique)', title_fontsize=5)
leg._legend_box.align = 'left'
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn3a')
fig.savefig(f'{OUT}/panel_f_copynumber_by_localization.pdf', dpi=300); plt.close(fig)

# =============================================================== g  TPM vs iPM
gd = df[(df['iPM_mean'] > 0) & (df['Illumina_sense_TPM'] > 0) & df['loc4'].isin(LOC_ORDER)].copy()
gd['log10_TPM'] = np.log10(gd['Illumina_sense_TPM'])
gd['log10_iPM'] = np.log10(gd['iPM_mean'])
fig, ax = plt.subplots(figsize=CORR, constrained_layout=True)
for loc in LOC_ORDER:
    s = gd[gd['loc4'] == loc]
    ax.scatter(s['log10_TPM'], s['log10_iPM'], s=11, alpha=0.8,           # open = syn3A
               facecolors='none', edgecolors=LOC_COLORS[loc], linewidths=0.5,
               label=f"{LOC_LABEL[loc]} (n={len(s)})")
r_all = pearsonr(gd['log10_TPM'], gd['log10_iPM'])[0]
cyto = gd[gd['loc4'] == 'cytoplasmic']
r_cyto = pearsonr(cyto['log10_TPM'], cyto['log10_iPM'])[0]
fit_line(ax, gd['log10_TPM'].values, gd['log10_iPM'].values, color='black')
ax.text(0.04, 0.96, f"all $r$ = {r_all:.2f} ($n$ = {len(gd)})\n"
                    f"cytoplasmic $r$ = {r_cyto:.2f} ($n$ = {len(cyto)})",
        transform=ax.transAxes, va='top', ha='left', fontsize=6)
ax.set_xlabel('mRNA Illumina TPM ($\\log_{10}$)')
ax.set_ylabel('Protein iPM ($\\log_{10}$)')
ax.legend(frameon=True, framealpha=0.9, edgecolor='0.7', handlelength=1.0,
          labelspacing=0.25, loc='lower right', fontsize=5)
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn3a')
fig.savefig(f'{OUT}/panel_g_TPM_vs_iPM.pdf', dpi=300); plt.close(fig)
say(f"g) all r={r_all:.3f} (n={len(gd)}); cytoplasmic r={r_cyto:.3f} (n={len(cyto)})")

# =============================================================== genome + CDS
STOP_CODONS = {"TAA", "TAG"}
GENETIC_CODE_4 = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"W","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def load_single_fasta(path):
    chunks = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                chunks.append(line.strip().upper())
    return "".join(chunks)

genome = load_single_fasta(GENOME_FA)
assert len(genome) == GENOME_LEN, f"genome length {len(genome)} != {GENOME_LEN}"
_COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s):
    return s.translate(_COMP)[::-1]

def _fetch_circular(g_start, g_end):
    g_start %= GENOME_LEN
    g_end = g_end % GENOME_LEN if g_end % GENOME_LEN != 0 else GENOME_LEN
    return genome[g_start:g_end] if g_start < g_end else genome[g_start:] + genome[:g_end]

def extract_cds(start0, end0, strand):
    seq = _fetch_circular(int(start0), int(end0))
    return revcomp(seq) if strand == "-" else seq

# protein-coding genes with valid in-frame CDS (coords from the TPM table)
cds_seqs = {}
n_bad = 0
for _, r in df.iterrows():
    if pd.isna(r['start0']):
        continue
    seq = extract_cds(r['start0'], r['end0'], r['strand'])
    if len(seq) < 6 or len(seq) % 3 != 0:
        n_bad += 1
        continue
    cds_seqs[r['locus_tag']] = seq
say(f"\nCDSs extracted: {len(cds_seqs)} (skipped non-multiple-of-3: {n_bad})")

# =============================================================== i  CAI (syn3A)
ipm_lookup = df.set_index('locus_tag')['iPM_mean'].to_dict()
prot_ranked = (df[df['locus_tag'].isin(cds_seqs) & (df['iPM_mean'] > 0)]
               .sort_values('iPM_mean', ascending=False))
n_ref = max(1, int(round(0.20 * len(prot_ranked))))
ref_loci = set(prot_ranked['locus_tag'].iloc[:n_ref])
say(f"CAI reference set (top 20% by iPM): {n_ref}/{len(prot_ranked)} genes")

SYNONYMOUS = defaultdict(list)
for codon, aa in GENETIC_CODE_4.items():
    if aa != "*":
        SYNONYMOUS[aa].append(codon)

ref_codon_counts = defaultdict(int)
for lt in ref_loci:
    codons = [cds_seqs[lt][i:i+3] for i in range(0, len(cds_seqs[lt]), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    for c in codons[1:]:
        if "N" in c or c in STOP_CODONS or GENETIC_CODE_4.get(c) is None:
            continue
        ref_codon_counts[c] += 1

w = {}
for aa, codons in SYNONYMOUS.items():
    counts = np.array([ref_codon_counts[c] for c in codons], dtype=float)
    if len(codons) == 1 or aa in ("M", "W") or counts.max() == 0:
        for c in codons:
            w[c] = 1.0
        continue
    rel = counts / counts.max()
    rel = np.where(rel == 0, 0.01, rel)
    for c, rv in zip(codons, rel):
        w[c] = float(rv)

def compute_cai(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP_CODONS:
        codons = codons[:-1]
    logs = []
    for c in codons[1:]:
        aa = GENETIC_CODE_4.get(c)
        if "N" in c or c in STOP_CODONS or aa is None or aa in ("M", "W"):
            continue
        logs.append(np.log(w[c]))
    return float(np.exp(np.mean(logs))) if logs else np.nan

cai_map = {lt: compute_cai(seq) for lt, seq in cds_seqs.items()}
df['CAI'] = df['locus_tag'].map(cai_map)

# =============================================================== h  TIR (OSTIR, gene-level, from scratch)
ANTI_SD, UTR_WINDOW, CDS_WINDOW = "ACCUCCUUU", 30, 30

def extract_initiation_seq(start0, end0, strand):
    """Fixed 30-nt window around the annotated start codon (circular)."""
    if strand == "+":
        seq = _fetch_circular(start0 - UTR_WINDOW, start0 + CDS_WINDOW)
    else:
        seq = revcomp(_fetch_circular(end0 - CDS_WINDOW, end0 + UTR_WINDOW))
    return seq, UTR_WINDOW   # start codon at 0-based index UTR_WINDOW

if os.path.exists(TIR_CACHE):
    tir_df = pd.read_csv(TIR_CACHE)
    say(f"\nTIR: loaded cache {TIR_CACHE} ({len(tir_df)} genes)")
else:
    from ostir import run_ostir
    say("\nTIR: running gene-level OSTIR on syn3A start windows ...")
    rows = []
    coord = df.dropna(subset=['start0']).set_index('locus_tag')
    for lt in cds_seqs:
        r = coord.loc[lt]
        seq, utr_len = extract_initiation_seq(int(r['start0']), int(r['end0']), r['strand'])
        start_1b = utr_len + 1
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = run_ostir(seq, aSD=ANTI_SD, threads=8)
        except Exception as e:
            say(f"  OSTIR error {lt}: {e}")
            continue
        if not out:
            continue
        best = next((h for h in out if h.get("start_position") == start_1b), None)
        if best is None:
            best = min(out, key=lambda x: abs(x.get("start_position", 0) - start_1b))
        rows.append({"locus_tag": lt, "TIR": float(best.get("expression", np.nan)),
                     "start_position": best.get("start_position", np.nan)})
    tir_df = pd.DataFrame(rows)
    tir_df.to_csv(TIR_CACHE, index=False)
    say(f"  TIR computed for {tir_df['TIR'].notna().sum()}/{len(cds_seqs)} genes -> {TIR_CACHE}")

df = df.merge(tir_df[['locus_tag', 'TIR']], on='locus_tag', how='left')

# =============================================================== residual helper
def residual_set(col, positive=True):
    """genes with iPM/TPM>0 and a valid `col`; return copy with residual of log10(iPM)~log10(TPM)."""
    s = df[(df['iPM_mean'] > 0) & (df['Illumina_sense_TPM'] > 0) & df[col].notna()].copy()
    if positive:
        s = s[s[col] > 0]
    s['log10_TPM'] = np.log10(s['Illumina_sense_TPM'])
    s['log10_iPM'] = np.log10(s['iPM_mean'])
    sl, ic = np.polyfit(s['log10_TPM'], s['log10_iPM'], 1)
    s['residual'] = s['log10_iPM'] - (sl * s['log10_TPM'] + ic)
    return s

# --- i  CAI vs residual
ci = residual_set('CAI', positive=False)
r_cai = pearsonr(ci['CAI'], ci['residual'])[0]
base_cai = r2_lstsq(ci[['log10_TPM']].values, ci['log10_iPM'].values)
full_cai = r2_lstsq(ci[['log10_TPM', 'CAI']].values, ci['log10_iPM'].values)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(ci['CAI'], ci['residual'], s=5, alpha=0.5, c=SYN3A_COL, edgecolors='none')
fit_line(ax, ci['CAI'].values, ci['residual'].values, color='black')
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_cai:.2f}\n$n$ = {len(ci)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('CAI')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn3a')
fig.savefig(f'{OUT}/panel_i_CAI_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"i) CAI vs residual r={r_cai:.3f} (n={len(ci)}); R2 {base_cai:.3f}->{full_cai:.3f} "
    f"(dR2={full_cai-base_cai:+.3f})")

# --- h  TIR vs residual
ti = residual_set('TIR', positive=True)
ti['log10_TIR'] = np.log10(ti['TIR'])
r_tir = pearsonr(ti['log10_TIR'], ti['residual'])[0]
base_tir = r2_lstsq(ti[['log10_TPM']].values, ti['log10_iPM'].values)
full_tir = r2_lstsq(ti[['log10_TPM', 'log10_TIR']].values, ti['log10_iPM'].values)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
ax.scatter(ti['log10_TIR'], ti['residual'], s=5, alpha=0.5, c=SYN3A_COL, edgecolors='none')
fit_line(ax, ti['log10_TIR'].values, ti['residual'].values, color='black')
ax.axhline(0, color='black', lw=0.5, ls=':')
ax.text(0.04, 0.96, f"$r$ = {r_tir:.2f}\n$n$ = {len(ti)}", transform=ax.transAxes, va='top', fontsize=6)
ax.set_xlabel('TIR ($\\log_{10}$)')
ax.set_ylabel('Proteome residual')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn3a')
fig.savefig(f'{OUT}/panel_h_TIR_vs_residual.pdf', dpi=300); plt.close(fig)
say(f"h) TIR vs residual r={r_tir:.3f} (n={len(ti)}); R2 {base_tir:.3f}->{full_tir:.3f} "
    f"(dR2={full_tir-base_tir:+.3f})")

# --- R with/without CAI, all + cytoplasmic (parallel to syn1 panel d; reported only)
def R_pair(d):
    y = d['log10_iPM'].values
    return np.sqrt(r2_lstsq(d[['log10_TPM']].values, y)), \
           np.sqrt(r2_lstsq(d[['log10_TPM', 'CAI']].values, y))
cyt = ci[ci['loc4'] == 'cytoplasmic']
Rb_all, Rc_all = R_pair(ci)
Rb_cyt, Rc_cyt = R_pair(cyt)
say(f"   R all {Rb_all:.3f}->{Rc_all:.3f} (n={len(ci)}); cytoplasmic {Rb_cyt:.3f}->{Rc_cyt:.3f} (n={len(cyt)})")

# --- l  Pearson R with/without CAI, all + cytoplasmic (mirrors syn1 panel d/f)
fig, ax = plt.subplots(figsize=(QUART, QUART), constrained_layout=True)
xpos = np.array([0, 1]); wbar = 0.38
b1 = ax.bar(xpos - wbar/2, [Rb_all, Rb_cyt], wbar, color='#bbbbbb', label='TPM only')
b2 = ax.bar(xpos + wbar/2, [Rc_all, Rc_cyt], wbar, color=SYN3A_COL, label='+ CAI')
for bars in (b1, b2):
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.012,
                f"{bar.get_height():.2f}", ha='center', va='bottom', fontsize=5)
ax.set_xticks(xpos); ax.set_xticklabels([f'All\n(n={len(ci)})', f'Cytosolic\n(n={len(cyt)})'])
ax.set_ylabel('$R$'); ax.set_ylim(0, 0.95)
ax.legend(frameon=False, handlelength=1.0, labelspacing=0.25, loc='upper left')
ax.spines[['top', 'right']].set_visible(False)
org_tag(ax, 'syn3a')
fig.savefig(f'{OUT}/panel_l_R_improvement.pdf', dpi=300); plt.close(fig)
say(f"l) R barplot: all {Rb_all:.3f}->{Rc_all:.3f}, cytosolic {Rb_cyt:.3f}->{Rc_cyt:.3f}")

# =============================================================== combined table
out = df.copy()
cols = ['locus_tag', 'gene_name', 'gene_product', 'loc4', 'Illumina_sense_TPM',
        'iPM_mean', 'exp_ptn_cnt_2026', 'CAI', 'TIR']
out = out[cols].rename(columns={'loc4': 'protein_localization',
                                'Illumina_sense_TPM': 'TPM_illumina',
                                'exp_ptn_cnt_2026': 'protein_copy_number'})
for c in ['TPM_illumina', 'iPM_mean', 'CAI', 'TIR']:
    out[c] = out[c].round(3)
out.to_csv('syn3A_genes_transcriptomics_proteomics.csv', index=False)
say(f"\nsyn3A_genes_transcriptomics_proteomics.csv: {len(out)} genes "
    f"(TPM {out['TPM_illumina'].notna().sum()}, iPM {out['iPM_mean'].notna().sum()}, "
    f"CAI {out['CAI'].notna().sum()}, TIR {out['TIR'].notna().sum()})")

with open(f'{OUT}/R3_syn3A.txt', 'w') as fh:
    fh.write("R3 syn3A PANELS (transcriptome-proteome correlation repeated on JCVI-syn3A)\n")
    fh.write("=" * 68 + "\n")
    fh.write("MAIN: g 7/3 x 7/3; f 7/3 x 7/6.  SI: h,i,l 7/4 x 7/4. Default fonts.\n\n")
    fh.write("\n".join(log) + "\n")
print(f"\nSaved 5 panels (2 main + 3 SI) + R3_syn3A.txt to {OUT}/")
