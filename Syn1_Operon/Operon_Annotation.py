#!/usr/bin/env python
# coding: utf-8

# ## Operon Annotation
#
# Analysis-only companion to the Syn1 operon map (operons.candidate_blocks.tsv).
# Per-operon visualization (the plot_one_operon driver) now lives in
# Operon_Visualization.ipynb; this script characterizes the operon set. Run from Syn1_Operon/.
#
# Sections and outputs (under annotation/):
#   1. Size of Operons -- operon-length + sense-genes-per-operon distributions
#      (R1 panels c/b: operon_length.pdf, genes_per_operon.pdf) + a QC scatter, and a
#      size / gene-count + antisense-coverage summary (Operon_Annotation.txt). Gene
#      counts are deduplicated from *_gene_loci at read time (the tsv over-counts genes
#      shared by merged operons; max 21, not the column's 23 -- see MANUSCRIPT.md TODO H).
#   2. Find Canonical Operons -- operons whose TSS and TTS both fall in intergenic regions
#      (trustworthy boundaries); 2x2 TSS/TTS intergenic-vs-intragenic matrix.
#   3. Spacing of genes inside canonical operons -- 5' and 3' UTR length distributions
#      (canonical/utr5_utr3.pdf) + long-UTR outlier lists (feeds R4 L4.2).
#   4. Promoters (leading and internal) -- -10 / -35 motif search at canonical TSS and
#      sequence logos (canonical/promoter_logo_{full,minus35,minus10,tss}.pdf), plus the
#      -10 classification (canonical/promoter_minus10_classification.{pdf,tsv}). (R1 L1.4, panel d.)
#   5. Terminators (internal and ending) -- intrinsic terminators near TTS and internal to
#      operons; hairpin sequence logos (canonical/tts_hairpins/all_hairpins.pdf) +
#      per-operon plots. (R1 L1.4, panel d.)
# In[40]:


import pandas as pd
import matplotlib
matplotlib.use("Agg")   # headless: plt.show() blocks on an X11 event loop if DISPLAY is set
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import Operon_Visualization as ov   # OperonCoord / plot_one_operon used in the terminator section

os.makedirs('annotation', exist_ok=True)

OPERON_TSV   = 'operons.candidate_blocks.tsv'                 # canonical 459-operon map (main folder)
GENE_COV_TSV = 'segmentation/gene_operon_coverage.tsv'       # moved to segmentation/ in the Operon_Segmentation output reorg

op = pd.read_csv(OPERON_TSV, sep='\t')
gc = pd.read_csv(GENE_COV_TSV, sep='\t')

print(f'Total operons: {len(op)}')
print(f'Total genes:   {len(gc)}')
print(f'Operon types:\n{op["segmentation_type"].value_counts().to_string()}')


# ## Size of Operons

# In[26]:


# ── Q1: Distribution of operon length and sense-gene count (R1 figure panels b, c) ──
# OUTPUT.md print spec: figsize (7/3, 7/3) in, Arial 5-7 pt, pdf.fonttype 42, no titles, vector.
# Panel b shows SENSE genes per operon only; antisense / total stats go to Operon_Annotation.txt.
#
# DEDUP AT CONSUMPTION: operons.candidate_blocks.tsv over-counts genes shared by merged operons
# (sense_gene_count inflated in 14 operons; max 23 vs 21 true -- see segmentation TO DO in
# MANUSCRIPT.md). The .tsv is left untouched; unique loci are recounted here so the figure and
# stats are correct.
import matplotlib as mpl
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

# Deduplicated per-operon gene counts (unique loci from the loci strings).
def _n_unique_loci(s):
    return len({x for x in str(s).split(',') if x and x != 'nan'})
sense_n_unique = op['sense_gene_loci'].apply(_n_unique_loci)
anti_n_unique  = op['antisense_gene_loci'].apply(_n_unique_loci)

# --- Panel c: operon length distribution ---
fig, ax = plt.subplots(figsize=(7/3, 7/3), constrained_layout=True)
ax.hist(op['length'], bins=50, color='steelblue', edgecolor='white', linewidth=0.3)
med_len = op['length'].median()
ax.axvline(med_len, color='crimson', lw=1.0, ls='--', label=f'median {med_len:.0f} bp')
ax.set_xlabel('Operon length (bp)')
ax.set_ylabel('Operon count')
ax.legend(loc='upper right')
fig.savefig('annotation/operon_length.pdf', dpi=300)
plt.show()

# --- Panel b: sense genes per operon (single series, log y, per-bar counts) ---
fig, ax = plt.subplots(figsize=(7/3, 7/3), constrained_layout=True)
sense_counts = sense_n_unique.value_counts().sort_index()
all_x = sorted(sense_counts.index)
x = np.arange(len(all_x))
vals = [int(sense_counts.get(k, 0)) for k in all_x]
FLOOR = 0.7  # finite bar baseline on the log axis: anchoring at 0 maps to log(-inf),
             # bloating the PDF path bbox to ~42,772 pt and breaking moves in Illustrator.
bars = ax.bar(x, [v - FLOOR for v in vals], 0.75, bottom=FLOOR, color='steelblue', edgecolor='white', lw=0.3)
ax.set_yscale('log')
ax.set_ylim(FLOOR, max(vals) * 4)
ax.set_xticks(x)
ax.set_xticklabels([str(v) for v in all_x])
ax.set_xlabel('Sense genes per operon')
ax.set_ylabel('Operon count')
for xi, v in zip(x, vals):
    if v > 0:
        ax.text(xi, v * 1.2, str(v), ha='center', va='bottom', rotation=90, fontsize=5)
fig.savefig('annotation/genes_per_operon.pdf', dpi=300)
plt.show()

# --- QC (not a manuscript panel): operon length vs sense-gene count ---
fig, ax = plt.subplots(figsize=(7/3, 7/3), constrained_layout=True)
jitter = np.random.default_rng(0).uniform(-0.2, 0.2, len(op))
ax.scatter(sense_n_unique + jitter, op['length'], alpha=0.4, s=6, color='steelblue', edgecolors='none')
ax.set_xlabel('Number of sense genes')
ax.set_ylabel('Operon length (bp)')
fig.savefig('annotation/length_vs_genecount.pdf', dpi=300)
plt.show()

# --- Size / gene-count summary + antisense coverage -> Operon_Annotation.txt (deduplicated) ---
summary_df = pd.DataFrame({
    'length': op['length'],
    'sense_gene_count': sense_n_unique,
    'antisense_gene_count': anti_n_unique,
    'gene_count': sense_n_unique + anti_n_unique,
})
antisense_counts = anti_n_unique.value_counts().sort_index()
n_anti = int((anti_n_unique >= 1).sum())
with open('Operon_Annotation.txt', 'w') as fh:
    fh.write("OPERON SIZE / GENE-COUNT SUMMARY\n")
    fh.write("=" * 60 + "\n\n")
    fh.write("(Gene counts deduplicated from *_gene_loci; the tsv's *_gene_count over-counts\n")
    fh.write(" genes shared by merged operons -- see segmentation TO DO.)\n\n")
    fh.write(f"n_operons = {len(op)}\n\n")
    fh.write("Per-operon summary (length / sense / antisense / total gene count):\n")
    fh.write(summary_df.describe().round(1).to_string() + "\n\n")
    fh.write("Antisense gene coverage per operon (removed from panel b):\n")
    fh.write(f"  operons with >=1 antisense gene: {n_anti} ({n_anti / len(op):.1%})\n")
    fh.write("  antisense_gene_count distribution (n_antisense_genes : n_operons):\n")
    fh.write(antisense_counts.to_string() + "\n")
print('Summary statistics:')
print(summary_df.describe().round(1).to_string())
print('Wrote size summary + antisense coverage stats to Operon_Annotation.txt')


# ## Find Canonical Operons
# 
# Canonical operons largely retain the information of transcrition upon the RNA processing.

# In[74]:


# ── Find Canonical Operons ──────────────────────────────────────────────────
#
# Canonical operon: segmentation_type == 'isoform_operon'  AND
#   TSS falls in an intergenic region (not inside any same-strand gene body) AND
#   TTS falls in an intergenic region (not inside any same-strand gene body)
#
# Intergenic = position is not contained within [gene.start0, gene.end0) on the
# same chrom and strand.
#
# For + strand: TSS = operon.start0, TTS = operon.end0
# For - strand: TSS = operon.end0,   TTS = operon.start0

def build_gene_intervals(gc: pd.DataFrame):
    intervals = {}
    for (chrom, strand), grp in gc.groupby(["chrom", "strand"]):
        intervals[(chrom, strand)] = list(
            zip(grp["start0"].astype(int), grp["end0"].astype(int))
        )
    return intervals

def is_intergenic(pos: int, chrom: str, strand: str, gene_intervals: dict) -> bool:
    """Return True if pos is NOT inside any gene body on (chrom, strand)."""
    for g0, g1 in gene_intervals.get((chrom, strand), []):
        if g0 < pos < g1:
            return False
    return True

gene_intervals = build_gene_intervals(gc)

# ── Add tss_intergenic / tts_intergenic columns to op ────────────────────────
tss_flags = []
tts_flags = []
for _, o in op.iterrows():
    chrom  = str(o["chrom"])
    strand = str(o["strand"])
    tss = int(o["start0"]) if strand == "+" else int(o["end0"])
    tts = int(o["end0"])   if strand == "+" else int(o["start0"])
    tss_flags.append(is_intergenic(tss, chrom, strand, gene_intervals))
    tts_flags.append(is_intergenic(tts, chrom, strand, gene_intervals))

op["tss_intergenic"] = tss_flags
op["tts_intergenic"] = tts_flags

# ── 2×2 contingency table for isoform_operon subset ──────────────────────────
iso_op = op[op["segmentation_type"] == "isoform_operon"].copy()

ct = pd.crosstab(
    iso_op["tss_intergenic"].map({True: "TSS intergenic", False: "TSS intragenic"}),
    iso_op["tts_intergenic"].map({True: "TTS intergenic", False: "TTS intragenic"}),
    rownames=[""], colnames=[""],
    margins=True, margins_name="Total",
)
# Reorder rows/cols: intergenic first
row_order = [r for r in ["TSS intergenic", "TSS intragenic", "Total"] if r in ct.index]
col_order = [c for c in ["TTS intergenic", "TTS intragenic", "Total"] if c in ct.columns]
ct = ct.loc[row_order, col_order]

print(f"isoform_operon total: {len(iso_op)}")
print()
print("2×2 matrix — TSS vs TTS boundary location (isoform_operon only):")
print(ct.to_string())
print()

# Canonical = isoform_operon + both boundaries intergenic
op_canonical = iso_op[iso_op["tss_intergenic"] & iso_op["tts_intergenic"]].copy()
print(f"Canonical operons (both intergenic): {len(op_canonical)}")
print("  Strand:", op_canonical["strand"].value_counts().to_dict())
print("  Sense gene count:", op_canonical["sense_gene_count"].value_counts().sort_index().to_dict())


# ## Spacing of genes inside canonical operons

# In[75]:


# ── 5'UTR and 3'UTR for canonical operons ────────────────────────────────────
#
# Canonical operons already have TSS and TTS guaranteed in intergenic regions.
# 5'UTR = distance from TSS to the 5' end of the first sense gene (tx direction)
# 3'UTR = distance from the 3' end of the last sense gene to TTS (tx direction)
#
# +strand: TSS = operon.start0, TTS = operon.end0
#   5'UTR = first_gene.start0 - operon.start0
#   3'UTR = operon.end0        - last_gene.end0
#
# -strand: TSS = operon.end0,   TTS = operon.start0
#   5'UTR = operon.end0        - last_gene.end0
#   3'UTR = first_gene.start0  - operon.start0

utr5_rows = []
utr3_rows = []

for _, o in op_canonical.iterrows():
    if pd.isna(o['sense_gene_loci']) or o['sense_gene_count'] < 1:
        continue
    loci = str(o['sense_gene_loci']).split(',')
    genes_in_op = gc[gc['locus_tag'].isin(loci)].sort_values('start0').copy()
    if genes_in_op.empty:
        continue

    strand     = o['strand']
    first_gene = genes_in_op.iloc[0]
    last_gene  = genes_in_op.iloc[-1]
    base = {'operon_id': o['operon_id'], 'strand': strand}

    if strand == '+':
        utr5 = int(first_gene['start0']) - int(o['start0'])
        utr3 = int(o['end0'])            - int(last_gene['end0'])
        utr5_rows.append({**base, 'locus_tag': first_gene['locus_tag'], 'gene_name': first_gene['gene_name'], 'utr5_bp': utr5})
        utr3_rows.append({**base, 'locus_tag': last_gene['locus_tag'],  'gene_name': last_gene['gene_name'],  'utr3_bp': utr3})
    else:
        utr5 = int(o['end0'])            - int(last_gene['end0'])
        utr3 = int(first_gene['start0']) - int(o['start0'])
        utr5_rows.append({**base, 'locus_tag': last_gene['locus_tag'],  'gene_name': last_gene['gene_name'],  'utr5_bp': utr5})
        utr3_rows.append({**base, 'locus_tag': first_gene['locus_tag'], 'gene_name': first_gene['gene_name'], 'utr3_bp': utr3})

utr5_df = pd.DataFrame(utr5_rows)
utr3_df = pd.DataFrame(utr3_rows)

# ── Export per-canonical-operon UTR table ────────────────────────────────────
# One row per canonical operon (isoform_operon, TSS+TTS both intergenic, >=1 sense
# gene) with both UTR lengths. Consumed by Syn1_Novel_ORF/R4_dist_panels.py
# (panel f) so its UTR distribution uses this exact canonical set rather than
# recomputing it (which had drifted from this definition).
utr_tbl = (utr5_df[['operon_id', 'strand', 'utr5_bp']]
           .merge(utr3_df[['operon_id', 'utr3_bp']], on='operon_id', how='outer')
           .sort_values('operon_id')
           .reset_index(drop=True))
os.makedirs('annotation/canonical', exist_ok=True)
utr_tbl.to_csv('annotation/canonical/operon_utr.tsv', sep='\t', index=False)
print(f"Saved canonical-operon UTR table: annotation/canonical/operon_utr.tsv ({len(utr_tbl)} operons)")

n = len(op_canonical)
print(f"Canonical operons with >=1 sense gene: {n}")
print()
print("5' UTR (TSS → first sense gene 5' end):")
print(utr5_df['utr5_bp'].describe().round(1).to_string())
print()
print("3' UTR (last sense gene 3' end → TTS):")
print(utr3_df['utr3_bp'].describe().round(1).to_string())


# In[76]:


# ── 5'UTR / 3'UTR plots — canonical operons ──────────────────────────────────
os.makedirs('annotation/canonical', exist_ok=True)

def hist_panel(ax, data, xlabel, title, color):
    data = data.dropna()
    p99 = data.quantile(1)
    ax.hist(data.clip(upper=p99), bins=40, color=color, edgecolor='white', linewidth=0.4)
    med = data.median()
    ax.axvline(med, color='black', lw=1.5, ls='--', label=f'Median {med:.0f} bp')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Operon Count')
    ax.set_title(title)
    ax.legend(fontsize=12)
    ax.text(0.9, 0.8, f'n_operon={len(data)}', transform=ax.transAxes,
            ha='right', va='top', fontsize=12)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
hist_panel(axes[0], utr5_df['utr5_bp'],
           "5' UTR length (bp)",
           "5' UTR\n(TSS → first sense gene 5' end)",
           '#2CA25F')
hist_panel(axes[1], utr3_df['utr3_bp'],
           "3' UTR length (bp)",
           "3' UTR\n(last sense gene 3' end → TTS)",
           '#E34A33')
# plt.suptitle("Canonical operons (isoform_operon, both boundaries intergenic)", fontsize=10)
plt.tight_layout()
plt.savefig("annotation/canonical/utr5_utr3.pdf", bbox_inches='tight')
plt.show()


# In[77]:


# ── Export operon IDs with extremely long 5' / 3' UTRs ───────────────────────
#
# "Extremely long" defined as > Q3 + 1.5 * IQR (standard Tukey upper fence).
# Operon plots are copied from operon_plots/ into annotation/canonical/long_5UTR
# and annotation/canonical/long_3UTR.

import shutil

os.makedirs('annotation/canonical/long_5UTR', exist_ok=True)
os.makedirs('annotation/canonical/long_3UTR', exist_ok=True)

def tukey_upper_fence(series: pd.Series) -> float:
    q1, q3 = series.quantile(0.25), series.quantile(0.75)
    return q3 + 1.5 * (q3 - q1)

fence5 = tukey_upper_fence(utr5_df['utr5_bp'])
fence3 = tukey_upper_fence(utr3_df['utr3_bp'])

long_5utr = utr5_df[utr5_df['utr5_bp'] > fence5].copy()
long_3utr = utr3_df[utr3_df['utr3_bp'] > fence3].copy()

print(f"5' UTR Tukey upper fence: {fence5:.0f} bp  →  {len(long_5utr)} outlier operons")
print(long_5utr[['operon_id', 'strand', 'locus_tag', 'gene_name', 'utr5_bp']]
      .sort_values('utr5_bp', ascending=False).to_string(index=False))
print()
print(f"3' UTR Tukey upper fence: {fence3:.0f} bp  →  {len(long_3utr)} outlier operons")
print(long_3utr[['operon_id', 'strand', 'locus_tag', 'gene_name', 'utr3_bp']]
      .sort_values('utr3_bp', ascending=False).to_string(index=False))

# Save operon ID lists
long_5utr[['operon_id', 'strand', 'locus_tag', 'gene_name', 'utr5_bp']] \
    .sort_values('utr5_bp', ascending=False) \
    .to_csv('annotation/canonical/long_5UTR/operon_ids.tsv', sep='\t', index=False)

long_3utr[['operon_id', 'strand', 'locus_tag', 'gene_name', 'utr3_bp']] \
    .sort_values('utr3_bp', ascending=False) \
    .to_csv('annotation/canonical/long_3UTR/operon_ids.tsv', sep='\t', index=False)

# Copy plots. operon_plots/ is written by Operon_Visualization.ipynb; if that has
# not been run the copies are skipped and only the operon_ids.tsv lists are produced.
def copy_operon_plots(operon_ids, dest):
    n = 0
    for opid in operon_ids:
        src = f'operon_plots/{opid}_wdepth.pdf'
        if os.path.exists(src):
            shutil.copy2(src, f'{dest}/{opid}_wdepth.pdf')
            n += 1
        else:
            print(f"  [warn] plot not found: {src}")
    return n

n5 = copy_operon_plots(long_5utr['operon_id'], 'annotation/canonical/long_5UTR')
n3 = copy_operon_plots(long_3utr['operon_id'], 'annotation/canonical/long_3UTR')

print(f"\nCopied {n5}/{len(long_5utr)} plots → annotation/canonical/long_5UTR/")
print(f"Copied {n3}/{len(long_3utr)} plots → annotation/canonical/long_3UTR/")


# ## Promoters (leading and internal)
# 
# - Signature of promoters
# 
# - Internal promoters.

# In[78]:


# ============================================================
# Promoter -10 motif analysis for canonical operon TSS
# ============================================================
#
# Algorithm mirrors End_Annotation/Peaks_Annotation.ipynb exactly:
#   1. TSS anchor: start0 (+ strand) or end0 (- strand) -- same geometry as summit_pos0
#   2. For each TSS, run best_shift_for_consensus independently for:
#        - 6-mer TANAAT  at [-12, -7]
#        - 9-mer TNNTANAAT at [-15, -7]
#      Selection: minimize IUPAC mismatches, tie-break by |shift| then sign
#      (identical to Peaks_Annotation best_shift_for_consensus)
#   3. Tier: strong_9mer > core_6mer > no_minus10  (np.select, same as Peaks_Annotation)
#   4. Build PFM from broad window [-40, +5] with no shift applied (unbiased logo)

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from Bio import SeqIO
from Bio.Seq import Seq

# -- Genome -------------------------------------------------------------------
GENOME_FASTA = "../Genomes_Input/syn1_genome.fasta"

genome = {}
for rec in SeqIO.parse(GENOME_FASTA, "fasta"):
    genome[rec.id] = str(rec.seq).upper()
print("Chromosomes loaded:", list(genome.keys()))

# -- Settings -----------------------------------------------------------------
SHIFT_RANGE = 2

CONS6 = "TANAAT";    REL6_START, REL6_END = -12, -7   # 6-mer -10 box
CONS9 = "TNNTANAAT"; REL9_START, REL9_END = -15, -7   # 9-mer -10 box

LOGO_START, LOGO_END = -40, +5   # broad window for sequence logo

# -- -10 scanner: single source in promoter_motif.py -------------------------
# These functions were extracted verbatim into promoter_motif.py so this canonical
# analysis and the R4 novel-transcription promoter scan share ONE implementation
# (verified: identical -10 classification, 0/127 mismatches). extract_tx_kmer reads
# promoter_motif.GENOME = the same syn1 FASTA loaded above as `genome`.
from promoter_motif import (
    IUPAC, consensus_to_regex, mismatches_iupac, circular_slice,
    extract_tx_kmer, best_shift_for_consensus, CONS6_RE, CONS9_RE,
)

# -- TSS coordinates ----------------------------------------------------------
def tss_pos(row):
    return int(row["start0"]) if row["strand"] == "+" else int(row["end0"])

# -- Run 6-mer and 9-mer independently (same as Peaks_Annotation) -------------
m6_best, m6_shift, m6_match, m6_mm = [], [], [], []
m9_best, m9_shift, m9_match, m9_mm = [], [], [], []

for _, o in op_canonical.iterrows():
    chrom = str(o["chrom"]); strand = str(o["strand"]); tss0 = tss_pos(o)
    k, sh, ma, mm = best_shift_for_consensus(tss0, chrom, strand, REL6_START, REL6_END, CONS6, CONS6_RE)
    m6_best.append(k); m6_shift.append(sh); m6_match.append(ma); m6_mm.append(mm)
    k, sh, ma, mm = best_shift_for_consensus(tss0, chrom, strand, REL9_START, REL9_END, CONS9, CONS9_RE)
    m9_best.append(k); m9_shift.append(sh); m9_match.append(ma); m9_mm.append(mm)

tss_class = op_canonical[["operon_id", "chrom", "strand"]].copy()
tss_class["tss0"]               = [tss_pos(o) for _, o in op_canonical.iterrows()]
tss_class["minus10_6mer_best"]  = m6_best
tss_class["minus10_6mer_shift"] = np.array(m6_shift, dtype=int)
tss_class["minus10_6mer_match"] = np.array(m6_match, dtype=bool)
tss_class["minus10_6mer_mm"]    = m6_mm
tss_class["minus10_9mer_best"]  = m9_best
tss_class["minus10_9mer_shift"] = np.array(m9_shift, dtype=int)
tss_class["minus10_9mer_match"] = np.array(m9_match, dtype=bool)
tss_class["minus10_9mer_mm"]    = m9_mm

# Tier (same np.select logic as Peaks_Annotation)
tss_class["motif_tier"] = np.select(
    [tss_class["minus10_9mer_match"], tss_class["minus10_6mer_match"]],
    ["strong_9mer", "core_6mer"],
    default="no_minus10"
)

print("-10 motif tier:")
print(tss_class["motif_tier"].value_counts().to_string())
print(f"\n9-mer match rate: {tss_class['minus10_9mer_match'].mean():.1%}")
print(f"6-mer match rate: {tss_class['minus10_6mer_match'].mean():.1%}")

# -- Broad-window PFM for logo (fixed window, no shift applied) ---------------
LOGO_LEN = LOGO_END - LOGO_START + 1
bases = ["A", "C", "G", "T"]
positions = list(range(LOGO_START, LOGO_END + 1))
pfm_counts = pd.DataFrame(0, index=bases, columns=positions, dtype=float)

valid_seqs = []
for _, row in tss_class.iterrows():
    seq = extract_tx_kmer(row["tss0"], row["chrom"], row["strand"], LOGO_START, LOGO_END)
    if len(seq) == LOGO_LEN and re.fullmatch(r"[ACGT]+", seq):
        valid_seqs.append(seq)
        for i, b in enumerate(seq):
            pfm_counts.iloc[pfm_counts.index.get_loc(b), i] += 1

pfm_freq = pfm_counts / pfm_counts.sum(axis=0)
print(f"\nValid sequences for PFM: {len(valid_seqs)} / {len(tss_class)}")


# In[31]:


# ============================================================
# Sequence logos (bits) -- -35 region, -10 region, and TSS
# ============================================================
#
# Two PFMs are built for different purposes:
#
#   pfm_freq      (TSS-aligned, no shift)
#     Built in the previous cell. All sequences anchored at the called TSS.
#     Use for: -35 logo and TSS context logo, where TSS is the fixed reference.
#
#   pfm_freq_10   (motif-aligned, 6-mer shift applied)
#     Each sequence is extracted at [REL6_START + shift, REL6_END + shift] so
#     all -10 boxes land at the same columns. This avoids the 1-2 bp smearing
#     that would otherwise wash out the logo when TSS calls vary slightly.
#     Use for: -10 region logo only.
#
# Y-axis is information content in bits (max 2 bits for DNA).
# IC formula: IC(pos) = 2 + sum_b [ f(b) * log2(f(b)) ]
# Pseudocount 0.5/4 per base to avoid log(0).

os.makedirs('annotation/canonical', exist_ok=True)

# -- Build motif-aligned PFM for the -10 region --------------------------------
# Window length is fixed at 14 positions (-17 to -4 relative to the motif center).
# We anchor by the best 6-mer shift, extending 5 bp on each side of [-12,-7].
MOTIF_FLANK = 5
M10_START = REL6_START - MOTIF_FLANK   # = -17 relative to shifted motif anchor
M10_END   = REL6_END   + MOTIF_FLANK   # = -2  relative to shifted motif anchor
M10_LEN   = M10_END - M10_START + 1

bases = ["A", "C", "G", "T"]
m10_positions = list(range(M10_START, M10_END + 1))
pfm10_counts = pd.DataFrame(0, index=bases, columns=m10_positions, dtype=float)

valid_m10 = []
for _, row in tss_class.iterrows():
    sh = int(row["minus10_6mer_shift"])
    seq = extract_tx_kmer(row["tss0"], row["chrom"], row["strand"],
                          M10_START + sh, M10_END + sh)
    if len(seq) == M10_LEN and re.fullmatch(r"[ACGT]+", seq):
        valid_m10.append(seq)
        for i, b in enumerate(seq):
            pfm10_counts.iloc[pfm10_counts.index.get_loc(b), i] += 1

pfm_freq_10 = pfm10_counts / pfm10_counts.sum(axis=0)
print(f"Motif-aligned sequences for -10 PFM: {len(valid_m10)} / {len(tss_class)}")

# -- Shared helpers ------------------------------------------------------------
def freq_to_bits(pfm_freq_df):
    """Convert frequency PFM to information content (bits) per position."""
    pseudo = 0.5 / 4
    pfm_ps = pfm_freq_df + pseudo
    pfm_ps = pfm_ps / pfm_ps.sum(axis=0)
    H = -(pfm_ps * np.log2(pfm_ps)).sum(axis=0)
    IC = 2 - H
    return pfm_ps.multiply(IC, axis=1)

def logo_bits(pfm_freq_df, col_slice, title, save_path, highlight_pos=None,
              figsize=None, hide_xticks=False):
    """
    Bits-scale sequence logo for columns col_slice (list of position ints).
    Y-axis in bits. When figsize is given the panel is born at that exact size
    (constrained_layout, no bbox-tight) for manuscript use; hide_xticks drops the
    per-position x tick labels.
    """
    sub = pfm_freq_df[col_slice].copy()
    bits = freq_to_bits(sub)
    logo_df = bits.T.copy()
    logo_df.columns = ["A", "C", "G", "T"]
    logo_df.index = logo_df.index.astype(int)
    logo_df = logo_df.sort_index()

    n_pos = len(logo_df)
    born_at_size = figsize is not None
    if figsize is None:
        figsize = (max(3, n_pos * 0.45), 2.5)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=born_at_size)
    logomaker.Logo(logo_df, ax=ax, color_scheme="classic",
                   shade_below=0.5, fade_below=0.5, stack_order="big_on_top")
    if highlight_pos is not None and highlight_pos in logo_df.index:
        ax.axvline(highlight_pos, color="black", lw=1, ls="-", alpha=0.35)
    ax.set_ylabel("bits", fontsize=6)
    ax.set_ylim(0, 1.0)
    if title:
        ax.set_title(title, fontsize=7)
    ax.set_xlabel("")
    if hide_xticks:
        ax.set_xticks([])
    else:
        ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if born_at_size:
        plt.savefig(save_path, dpi=300)            # exact figsize (e.g. 1x1 in)
    else:
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Saved: {save_path}")

# -- Sub-window definitions ----------------------------------------------------
W35  = list(range(-37, -31))    # TSS-aligned: -37 to -32  (broad -35 region)
WTSS = list(range(-1, 1))       # TSS-aligned: TSS context
W10  = list(range(REL6_START, REL6_END + 1))  # motif-aligned: full -10 window

# -35 and -10 promoter signatures: 1x1 in panels, no x-axis ticks (R1 panel d)
logo_bits(pfm_freq, W35, "", "annotation/canonical/promoter_logo_minus35.pdf",
          figsize=(1, 1), hide_xticks=True)
logo_bits(pfm_freq_10, W10, "", "annotation/canonical/promoter_logo_minus10.pdf",
          figsize=(1, 1), hide_xticks=True)
# TSS context logo (unchanged size)
logo_bits(pfm_freq, WTSS, "", "annotation/canonical/promoter_logo_tss.pdf", highlight_pos=0)

# -- Full-window overview logo (TSS-aligned) -----------------------------------
bits_all = freq_to_bits(pfm_freq)
logo_all = bits_all.T.copy()
logo_all.columns = ["A", "C", "G", "T"]
logo_all.index = logo_all.index.astype(int)
logo_all = logo_all.sort_index()

fig, ax = plt.subplots(figsize=(14, 2.5))
logomaker.Logo(logo_all, ax=ax, color_scheme="classic",
               shade_below=0.5, fade_below=0.5, stack_order="big_on_top")
for xv, ls, lab in [(-35, ":", "-35"), (-10, "--", "-10"), (0, "-", "TSS")]:
    ax.axvline(xv, color="gray", lw=1, ls=ls, alpha=0.7, label=lab)
ax.legend(fontsize=8, loc="upper left")
ax.set_ylabel("bits")
ax.set_ylim(0, 1.0)
ax.set_title(f"Promoter context -- canonical operon TSS (n={len(valid_seqs)})", fontsize=9)
ax.set_xlabel("")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("annotation/canonical/promoter_logo_full.pdf", dpi=300, bbox_inches="tight")
plt.show()
print("Saved: annotation/canonical/promoter_logo_full.pdf")

# -- Classification bar chart --------------------------------------------------
tiers = ["strong_9mer", "core_6mer", "no_minus10"]
color_map = {"strong_9mer": "#2CA25F", "core_6mer": "#FEC44F", "no_minus10": "#E34A33"}
counts = tss_class["motif_tier"].value_counts().reindex(tiers, fill_value=0)

fig, ax = plt.subplots(figsize=(4, 3))
bars = ax.bar(tiers, counts.values,
              color=[color_map[t] for t in tiers],
              edgecolor="white", linewidth=0.6)
for bar, v in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, v + 0.5, str(v),
            ha="center", va="bottom", fontsize=9)
ax.set_ylabel("Number of operons")
ax.set_title("-10 box match (canonical TSS)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("annotation/canonical/promoter_minus10_classification.pdf", dpi=300, bbox_inches="tight")
plt.show()
print("Saved: annotation/canonical/promoter_minus10_classification.pdf")

# -- Export table --------------------------------------------------------------
tss_class.to_csv("annotation/canonical/promoter_minus10_classification.tsv", sep="\t", index=False)
print("Saved: annotation/canonical/promoter_minus10_classification.tsv")


# ## Terminators (internal and ending)
# 
# For one-gene operon, the stop signal stronger?

# In[32]:


# ============================================================
# Parse TransTermHP output + annotate canonical operons
# ============================================================
#
# Questions:
#   1. Does every canonical operon have a predicted terminator near its TTS?
#   2. Are there internal terminators inside canonical operons?
#
# TransTermHP coordinate convention (1-based):
#   + strand: start1 < end1   (genomic left to right)
#   - strand: start1 > end1   (reported reversed; min=left, max=right in genome)
# We convert everything to 0-based half-open [start0, end0) with strand kept.
#
# Sequence line format (fixed-width blocks separated by 2+ spaces):
#   <5' tail>  <5' stem> <loop> <3' stem>  <3' tail>
# The middle block (stem5, loop, stem3) is single-space delimited.
# Gaps in stems are represented as '-'.

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

TRANSTERM_TXT = "syn1_TransTermHP.txt"
ISOFORMS_TSV  = "../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"

# -- Parse TransTermHP --------------------------------------------------------
TERM_RE = re.compile(
    r"^\s+TERM\s+(\d+)\s+(\d+)\s+-\s+(\d+)\s+([+-])\s+(\S+)\s+(\d+)\s+([\-\d\.]+)\s+([\-\d\.]+)"
)
SEQ_RE = re.compile(r"^\s+[ACGT]")   # sequence line: starts with whitespace then a nucleotide

def parse_seq_line(line):
    """
    Split TransTermHP sequence line into (tail5, stem5, loop, stem3, tail3).
    The line has three blocks separated by 2+ spaces:
      block0 = 5' tail
      block1 = '5stem loop 3stem'  (single-space separated)
      block2 = 3' tail
    All returned as RNA (T -> U).
    """
    blocks = re.split(r"  +", line.strip())
    if len(blocks) < 3:
        return None
    tail5 = blocks[0].replace("T", "U")
    tail3 = blocks[-1].replace("T", "U")
    mid   = blocks[1].strip().split(" ")
    if len(mid) < 3:
        return None
    stem5 = mid[0].replace("T", "U")
    loop  = mid[1].replace("T", "U")
    stem3 = mid[2].replace("T", "U")
    return tail5, stem5, loop, stem3, tail3


def rc_rna(s):
    """Reverse-complement an RNA/DNA block, preserving '-' gap characters."""
    return s.translate(str.maketrans("ACGUTacgut", "UGCAAugcaa"))[::-1]


terms = []
pending = None   # dict waiting for its sequence line

with open(TRANSTERM_TXT) as fh:
    chrom = None
    for line in fh:
        if line.startswith("SEQUENCE "):
            chrom = line.split()[1]
            pending = None
            continue

        m = TERM_RE.match(line)
        if m and chrom:
            term_id  = int(m.group(1))
            c1, c2   = int(m.group(2)), int(m.group(3))
            strand   = m.group(4)
            region   = m.group(5)
            conf     = int(m.group(6))
            hp_score = float(m.group(7))
            tail_sc  = float(m.group(8))
            g_min = min(c1, c2) - 1
            g_max = max(c1, c2)
            pending = {
                "term_id": term_id, "chrom": chrom,
                "strand": strand, "region": region,
                "start0": g_min, "end0": g_max,
                "conf": conf, "hp_score": hp_score, "tail_score": tail_sc,
                "tail5": "", "stem5": "", "loop": "", "stem3": "", "tail3": ""
            }
            continue

        # Sequence line immediately follows its TERM line
        if pending is not None and SEQ_RE.match(line):
            parsed = parse_seq_line(line)
            if parsed:
                t5, s5, lp, s3, t3 = parsed
                if pending["strand"] == "-":
                    # TransTermHP prints minus-strand terminators in FORWARD-genome
                    # orientation; reverse-complement to transcript 5'->3' so the poly-U
                    # becomes the 3' tail (matches the minus-strand transcript & panel d).
                    t5, s5, lp, s3, t3 = rc_rna(t3), rc_rna(s3), rc_rna(lp), rc_rna(s5), rc_rna(t5)
                pending["tail5"], pending["stem5"], pending["loop"], \
                    pending["stem3"], pending["tail3"] = t5, s5, lp, s3, t3
            terms.append(pending)
            pending = None

# Flush any trailing pending (shouldn't happen, but safe)
if pending is not None:
    terms.append(pending)

terms_df = pd.DataFrame(terms)
print(f"Total TransTermHP predictions: {len(terms_df)}")
print(f"Strand distribution:\n{terms_df['strand'].value_counts().to_string()}")
print(f"Confidence distribution:\n{terms_df['conf'].describe().round(1).to_string()}")
print(f"\nSample parsed sequences:")
print(terms_df[["term_id","strand","conf","tail5","stem5","loop","stem3","tail3"]].head(8).to_string(index=False))


# In[34]:


# -- Q1: Does every canonical operon have a terminator near its TTS? ----------
# "Near" = within TERM_WINDOW bp downstream of TTS (in transcription direction)
# We require same strand.
TERM_WINDOW = 50   # bp search window past TTS

def _term_3p(terms, strand):
    """Terminator 3' boundary in transcription direction:
    '+' strand -> end0 (larger coord); '-' strand -> start0 (smaller coord)."""
    return terms["end0"] if strand == "+" else terms["start0"]


def find_terminators_near_tts(op_row, terms_df, window=TERM_WINDOW):
    """
    Return same-strand terminators whose 3' boundary is within +/- window bp of
    the operon TTS.  The boundary must be strand-correct: for a '-' operon the
    terminator's 3' end is its SMALLER coordinate (start0), not end0.  The old
    rule tested end0 on both strands, which silently dropped almost every
    '-'-strand terminator (its end0 = the 5' stem-start sits ~15-30 bp past the
    TTS) -- e.g. OP_00099/0178 missed its conf-100 TransTermHP terminator.
    """
    chrom  = str(op_row["chrom"])
    strand = str(op_row["strand"])
    tts    = int(op_row["tts"])
    t = terms_df[(terms_df["chrom"] == chrom) & (terms_df["strand"] == strand)].copy()
    t["term_3p"] = _term_3p(t, strand)
    return t[(t["term_3p"] - tts).abs() <= window]

tts_hits = []
for _, o in op_canonical.iterrows():
    hit = find_terminators_near_tts(o, terms_df, window=TERM_WINDOW)
    rec = {
        "operon_id": o["operon_id"],
        "strand":    o["strand"],
        "n_sense_genes": o["sense_gene_count"],
        "operon_len": o["end0"] - o["start0"],
        "has_tts_term": len(hit) > 0,
        "n_tts_terms": len(hit),
        "best_conf":   hit["conf"].max() if len(hit) > 0 else np.nan,
        "best_hp":     hit["hp_score"].min() if len(hit) > 0 else np.nan,
    }
    # geometry of the highest-confidence TTS terminator -- the intrinsic-terminator
    # signature exported to the SI operon workbook (build_operon_xlsx.py)
    if len(hit) > 0:
        best = hit.loc[hit["conf"].idxmax()]
        rec.update({
            "term_id":        int(best["term_id"]),
            "term_stem_bp":   len(str(best["stem5"])),
            "term_loop_nt":   len(str(best["loop"])),
            "term_polyU_nt":  len(str(best["tail3"])),
            "term_stem5_seq": str(best["stem5"]),
            "term_loop_seq":  str(best["loop"]),
            "term_tail3_seq": str(best["tail3"]),
        })
    else:
        rec.update({"term_id": np.nan, "term_stem_bp": np.nan, "term_loop_nt": np.nan,
                    "term_polyU_nt": np.nan, "term_stem5_seq": "",
                    "term_loop_seq": "", "term_tail3_seq": ""})
    tts_hits.append(rec)

tts_df = pd.DataFrame(tts_hits)
# persist the per-operon terminator signature (single source of truth, read by
# build_operon_xlsx.py); symmetric with promoter_minus10_classification.tsv above
tts_df.to_csv("annotation/canonical/terminator_tts_classification.tsv", sep="\t", index=False)
print("Saved: annotation/canonical/terminator_tts_classification.tsv")
n_total = len(tts_df)
n_with  = tts_df["has_tts_term"].sum()
print(f"\nQ1 — Terminator within {TERM_WINDOW} bp of TTS:")
print(f"  Canonical operons: {n_total}")
print(f"  With predicted terminator: {n_with} ({n_with/n_total:.1%})")
print(f"  Without:                   {n_total - n_with} ({(n_total-n_with)/n_total:.1%})")
print()
print("By sense gene count:")
print(tts_df.groupby("n_sense_genes")["has_tts_term"]
      .agg(["sum","count"])
      .rename(columns={"sum":"with_term","count":"total"})
      .assign(frac=lambda d: (d["with_term"]/d["total"]).round(3))
      .to_string())


# ### Internal Promoters

# In[35]:


# -- Q2: Internal terminators inside canonical operons ------------------------
# "Internal" = terminator on the same strand, fully within [operon.start0, operon.end0),
#  AND NOT overlapping the TTS window (so it is upstream of the TTS).
def find_internal_terminators(op_row, terms_df, tts_window=TERM_WINDOW):
    chrom  = str(op_row["chrom"])
    strand = str(op_row["strand"])
    op_s   = int(op_row["start0"])
    op_e   = int(op_row["end0"])
    tts    = int(op_row["tts"])

    t = terms_df[(terms_df["chrom"] == chrom) & (terms_df["strand"] == strand)].copy()
    t["term_3p"] = _term_3p(t, strand)
    # 3' boundary inside the operon body, but not within the TTS window
    inside = t[t["term_3p"].between(op_s, op_e)]
    internal = inside[(inside["term_3p"] - tts).abs() > tts_window]
    return internal

int_hits = []
for _, o in op_canonical.iterrows():
    hit = find_internal_terminators(o, terms_df, tts_window=TERM_WINDOW)
    int_hits.append({
        "operon_id":      o["operon_id"],
        "strand":         o["strand"],
        "n_sense_genes":  o["sense_gene_count"],
        "operon_len":     o["end0"] - o["start0"],
        "tss":            o["tss"],
        "tts":            o["tts"],
        "has_internal":   len(hit) > 0,
        "n_internal":     len(hit),
        "internal_confs": ",".join(str(c) for c in hit["conf"].tolist()) if len(hit) > 0 else "",
        "internal_pos":   ",".join(str(s) for s in hit["start0"].tolist()) if len(hit) > 0 else "",
    })

int_df = pd.DataFrame(int_hits)
n_int = int_df["has_internal"].sum()
print(f"\nQ2 — Internal terminators (inside operon body, >={TERM_WINDOW} bp from TTS):")
print(f"  Operons with internal terminator: {n_int} / {len(int_df)} ({n_int/len(int_df):.1%})")
print()
print("By sense gene count:")
print(int_df.groupby("n_sense_genes")["has_internal"]
      .agg(["sum","count"])
      .rename(columns={"sum":"with_internal","count":"total"})
      .assign(frac=lambda d: (d["with_internal"]/d["total"]).round(3))
      .to_string())
print()
print("Operons with internal terminators:")
print(int_df[int_df["has_internal"]][
    ["operon_id","strand","n_sense_genes","operon_len","n_internal","internal_confs","internal_pos"]
].sort_values("n_internal", ascending=False).to_string(index=False))


# In[36]:


# ============================================================
# Visualize operons with internal terminators + their isoforms
# ============================================================
#
# Two-panel plot per operon:
#   Panel 1: gene arrow cartoon in genomic coords, with TSS/TTS and
#            internal/TTS terminator markers.
#   Panel 2: isoforms drawn via ov.draw_isoforms() in transcript coords
#            (5'->3' left->right for both strands), with terminator positions
#            overlaid as vertical lines in tx space.
#
# Isoform rendering (thickness ~ n_reads, row-packed, colored by TSS group)
# is inherited directly from Operon_Visualization.py.

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# -- Load isoforms ------------------------------------------------------------
isoforms_all = pd.read_csv(ISOFORMS_TSV, sep="\t")

MIN_ISOFORM_READS = 10

def get_operon_isoforms(op_row, iso_df, min_reads=MIN_ISOFORM_READS):
    chrom  = str(op_row["chrom"])
    strand = str(op_row["strand"])
    op_s   = int(op_row["start0"])
    op_e   = int(op_row["end0"])
    mask = (
        (iso_df["chrom"] == chrom) &
        (iso_df["strand"] == strand) &
        (iso_df["start0"] >= op_s - 50) &
        (iso_df["end0"]   <= op_e + 50) &
        (iso_df["n_reads"] >= min_reads)
    )
    return iso_df[mask].copy()

os.makedirs("annotation/canonical/internal_terminators", exist_ok=True)

# Operons with internal terminators
op_internal = op_canonical[op_canonical["operon_id"].isin(
    int_df[int_df["has_internal"]]["operon_id"]
)].copy()
op_internal = op_internal.merge(
    int_df[int_df["has_internal"]][["operon_id","n_internal","internal_confs","internal_pos"]],
    on="operon_id", how="left"
)

print(f"Plotting {len(op_internal)} operons with internal terminators...")

GENE_COLORS = plt.cm.tab10.colors

for _, op in op_internal.iterrows():
    opid   = op["operon_id"]
    chrom  = str(op["chrom"])
    strand = str(op["strand"])
    op_s   = int(op["start0"])
    op_e   = int(op["end0"])
    tss0   = int(op["tss"])
    tts0   = int(op["tts"])

    # OperonCoord for tx-space conversions (used by ov.draw_isoforms)
    oc = ov.OperonCoord(chrom=chrom, strand=strand, opid=opid, start0=op_s, end0=op_e)

    plot_s = op_s - 100
    plot_e = op_e + 100

    # Sense genes
    if pd.notna(op["sense_gene_loci"]) and op["sense_gene_count"] > 0:
        loci = str(op["sense_gene_loci"]).split(",")
        genes_in = gc[gc["locus_tag"].isin(loci)].sort_values("start0").copy()
    else:
        genes_in = pd.DataFrame()

    # Terminators
    int_terms = find_internal_terminators(op, terms_df, tts_window=TERM_WINDOW)
    tts_terms = find_terminators_near_tts(op, terms_df, window=TERM_WINDOW)

    # Isoforms
    iso = get_operon_isoforms(op, isoforms_all)

    # ── Layout ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(12, 5),
                              gridspec_kw={"height_ratios": [2, 3]},
                              sharex=False)
    ax_op  = axes[0]
    ax_iso = axes[1]

    # ── Panel 1: gene arrows in genomic coords ────────────────────────────────
    ax_op.hlines(0, op_s, op_e, color="black", lw=1.5, zorder=1)

    for gi, (_, g) in enumerate(genes_in.iterrows()):
        color = GENE_COLORS[gi % len(GENE_COLORS)]
        gname = g["gene_name"] if pd.notna(g["gene_name"]) else g["locus_tag"]
        rect = mpatches.FancyArrow(
            g["start0"], 0,
            g["end0"] - g["start0"], 0,
            width=0.25, head_width=0.35,
            head_length=max(20, (g["end0"] - g["start0"]) * 0.08),
            length_includes_head=True,
            color=color, zorder=2, alpha=0.85
        )
        ax_op.add_patch(rect)
        mid = (g["start0"] + g["end0"]) / 2
        ax_op.text(mid, 0.45, gname, ha="center", va="bottom", fontsize=7, color=color)

    ax_op.axvline(tss0, color="green",  lw=1.5, ls="--", alpha=0.7)
    ax_op.axvline(tts0, color="purple", lw=1.5, ls="--", alpha=0.7)

    for _, t in int_terms.iterrows():
        tmid = (t["start0"] + t["end0"]) / 2
        ax_op.annotate("", xy=(tmid, -0.5), xytext=(tmid, -0.15),
                       arrowprops=dict(arrowstyle="-|>", color="red", lw=1.5))
        ax_op.text(tmid, -0.62, f"T(c={t['conf']})", ha="center", va="top",
                   fontsize=6.5, color="red")

    for _, t in tts_terms.iterrows():
        tmid = (t["start0"] + t["end0"]) / 2
        ax_op.annotate("", xy=(tmid, -0.5), xytext=(tmid, -0.15),
                       arrowprops=dict(arrowstyle="-|>", color="gray", lw=1.2))
        ax_op.text(tmid, -0.62, f"T(c={t['conf']})", ha="center", va="top",
                   fontsize=6.5, color="gray")

    ax_op.set_xlim(plot_s, plot_e)
    ax_op.set_ylim(-1.0, 1.0)
    ax_op.set_yticks([])
    ax_op.set_title(
        f"{opid}  ({strand})  |  {op['sense_gene_count']} genes  |  "
        f"{op_e - op_s:,} bp  |  {len(int_terms)} internal terminator(s)",
        fontsize=9
    )
    handles = [
        Line2D([0],[0], color="green",  ls="--", label="TSS"),
        Line2D([0],[0], color="purple", ls="--", label="TTS"),
        mpatches.Patch(color="red",  label="Internal term"),
        mpatches.Patch(color="gray", label="TTS term"),
    ]
    ax_op.legend(handles=handles, fontsize=7, loc="upper right")
    ax_op.spines["top"].set_visible(False)
    ax_op.spines["left"].set_visible(False)
    ax_op.spines["right"].set_visible(False)

    # ── Panel 2: isoforms in transcript coords (via ov.draw_isoforms) ─────────
    ov.draw_isoforms(ax_iso, oc, iso, plot_s, plot_e)

    # Overlay terminator positions in tx coords
    # oc.tx_of_genome_pos0 maps genomic -> tx (0 = operon 5' end, op_e-op_s = 3' end)
    for _, t in int_terms.iterrows():
        tmid_g = (t["start0"] + t["end0"]) / 2
        tmid_tx = oc.tx_of_genome_pos0(int(tmid_g))
        ax_iso.axvline(tmid_tx, color="red",    lw=1.5, ls="--", alpha=0.8,
                       label=f"Internal T (c={t['conf']})")

    tts_tx = oc.tx_of_genome_pos0(tts0)
    ax_iso.axvline(tts_tx, color="purple", lw=1.5, ls="--", alpha=0.6, label="TTS")

    ax_iso.set_xlabel("Transcript position (bp, 5'→3')", fontsize=8)
    ax_iso.set_xlim(oc.tx_of_genome_pos0(plot_s), oc.tx_of_genome_pos0(plot_e))

    plt.tight_layout()
    save_path = f"annotation/canonical/internal_terminators/{opid}.pdf"
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()

print(f"Saved plots to annotation/canonical/internal_terminators/")
print()

# -- Summary: isoform termination at internal terminators ---------------------
summary_rows = []
for _, op in op_internal.iterrows():
    opid   = op["operon_id"]
    strand = str(op["strand"])
    int_terms = find_internal_terminators(op, terms_df, tts_window=TERM_WINDOW)
    iso = get_operon_isoforms(op, isoforms_all)
    if len(iso) == 0 or len(int_terms) == 0:
        continue
    for _, t in int_terms.iterrows():
        tmid = (t["start0"] + t["end0"]) / 2
        if strand == "+":
            ends3 = iso["end0"].values
        else:
            ends3 = iso["start0"].values
        n_at_term    = ((ends3 >= t["start0"] - 50) & (ends3 <= t["end0"] + 50)).sum()
        n_readthrough = (ends3 > t["end0"] + 50).sum() if strand == "+" else (ends3 < t["start0"] - 50).sum()
        summary_rows.append({
            "operon_id": opid, "strand": strand,
            "n_sense_genes": op["sense_gene_count"],
            "term_conf": t["conf"], "term_mid": int(tmid),
            "n_isoforms_total": len(iso),
            "n_at_internal_term": n_at_term,
            "n_readthrough": n_readthrough,
            "frac_at_term": n_at_term / max(1, len(iso)),
        })

if summary_rows:
    sum_df = pd.DataFrame(summary_rows)
    print("Isoform termination at internal terminators:")
    print(sum_df[["operon_id","strand","n_sense_genes","term_conf",
                  "n_isoforms_total","n_at_internal_term","n_readthrough","frac_at_term"]]
          .sort_values("frac_at_term", ascending=False)
          .to_string(index=False))


# In[37]:


# ============================================================
# Draw hairpin secondary structure for canonical operon TTS terminators
# ============================================================
#
# Classic hairpin/lollipop orientation:
#
#                    loop
#                   ┌────┐
#                   │    │
#       5' stem ────┘    └──── 3' stem
#                   (stem ascends bottom→top in the middle)
#              /                      \
#       5' tail (to lower-left)   3' tail (to lower-right, poly-U)
#
# Ts are already converted to Us by the parser.
# Base pairs are linked by short dashes drawn as text between the columns.

import os
import math
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("annotation/canonical/tts_hairpins", exist_ok=True)

# -- Collect TTS terminators for canonical operons + operon mapping -----------
term_to_operons = {}          # term_id -> list of operon_ids
tts_term_rows = []
for _, op in op_canonical.iterrows():
    hit = find_terminators_near_tts(op, terms_df, window=TERM_WINDOW)
    if len(hit) == 0:
        continue
    for _, h in hit.iterrows():
        term_to_operons.setdefault(int(h["term_id"]), []).append(str(op["operon_id"]))
    best = hit.loc[hit["conf"].idxmax()]
    tts_term_rows.append(best)

tts_terms_canonical = pd.DataFrame(tts_term_rows).drop_duplicates("term_id").reset_index(drop=True)
print(f"Unique TTS terminators for canonical operons: {len(tts_terms_canonical)}")

# -- Drawing helpers ----------------------------------------------------------
MONO = {"family": "monospace"}
# Exact logomaker 'classic' colours, so hairpin bases match the promoter logos.
BASE_COLOR = {
    "A": "#008000",   # green
    "C": "#0000ff",   # blue
    "G": "#ffa600",   # orange
    "U": "#ff0000",   # red
    "T": "#ff0000",
    "-": "#AAAAAA",
}

def base_col(ch):
    return BASE_COLOR.get(ch.upper(), "#333333")

def draw_terminator_hairpin(ax, term_row, title=None, operon_ids=None, minimal=False):
    """
    Draw a single terminator hairpin with:
      - stem rising vertically in the middle (base at bottom, loop on top)
      - base pairs linked by short connector lines between the two stem columns
      - 5' tail streaming down-and-left from the stem base
      - 3' tail (poly-U) streaming down-and-right from the stem base

    minimal=True  -> hairpin only: no title, no 5'/3' tail-end labels (for the
    standalone per-terminator panels used in the figure).
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    tail5 = term_row["tail5"]
    stem5 = term_row["stem5"]
    loop  = term_row["loop"]
    stem3 = term_row["stem3"]
    tail3 = term_row["tail3"]

    stem_len = max(len(stem5), len(stem3))
    n_loop   = max(1, len(loop))

    # ── Style / sizing ───────────────────────────────────────────────────────
    # The standalone per-terminator panels are born at 1x1 in (axes == figure),
    # so they need small fonts (Nature 5 pt floor) and tight, cw-scaled spacing;
    # the QC summary grid keeps the larger original style. cw = nucleotide
    # character width in axes units.
    if minimal:
        fs_nt, fs_loop = 5, 5
        cw       = 0.050
        row_cap  = 0.075
        tail_show = 6
    else:
        fs_nt, fs_loop = 10, 9
        cw       = 0.020
        row_cap  = 0.050
        tail_show = 10

    # Layout (axes [0,1] coords, y up). Vertically centered; row height adapts
    # so even a 19 bp stem fits inside the panel.
    cx        = 0.50
    stem_gap  = 0.9 * cw                 # horizontal gap between the two stem columns
    off       = 0.55 * cw                # base-letter offset from its column line
    Y0, Y1    = 0.07, 0.93
    row_h     = min(row_cap, (Y1 - Y0) / (stem_len + 2.0))
    loop_h    = 0.6 * row_h + 0.012 * n_loop
    block_h   = stem_len * row_h + loop_h
    stem_base = Y0 + ((Y1 - Y0) - block_h) / 2     # center the whole hairpin
    stem_top  = stem_base + stem_len * row_h

    x_L = cx - stem_gap / 2
    x_R = cx + stem_gap / 2

    # ── Stem bases: index 0 at bottom (nearest tails), increasing upward ─────
    # stem5 is 5'→3' entering the hairpin, so its first char is adjacent to tail5
    # (at stem_base) and its last char is adjacent to the loop (at stem_top).
    for i, ch in enumerate(stem5):
        y = stem_base + i * row_h
        ax.text(x_L - off, y, ch, ha="right", va="center",
                fontsize=fs_nt, color=base_col(ch), fontdict=MONO,
                transform=ax.transAxes)

    # stem3 is 5'→3' leaving the hairpin; its first char is adjacent to the
    # loop (top) and its last char is adjacent to tail3 (bottom). So drawn
    # top-to-bottom as-is on the right column.
    # Place so that stem3[0] sits at stem_top - row_h (paired with stem5[-1])
    for i, ch in enumerate(stem3):
        y = stem_top - row_h - i * row_h
        ax.text(x_R + off, y, ch, ha="left", va="center",
                fontsize=fs_nt, color=base_col(ch), fontdict=MONO,
                transform=ax.transAxes)

    # ── Base-pair dashes between columns ─────────────────────────────────────
    # Pair stem5[i] with stem3[stem_len-1-i] at row y_i
    for i in range(stem_len):
        ch5 = stem5[i] if i < len(stem5) else "-"
        j   = len(stem3) - 1 - i
        ch3 = stem3[j] if 0 <= j < len(stem3) else "-"
        if ch5 == "-" or ch3 == "-":
            continue
        y = stem_base + i * row_h
        pair = frozenset([ch5.upper(), ch3.upper()])
        # Connector line spans the gap between the two base columns; dashed for
        # a G-U wobble, solid for Watson-Crick. Scales with stem_gap.
        ls = (0, (1.5, 1.5)) if pair == frozenset(["G", "U"]) else "-"
        ax.plot([x_L - 0.35 * cw, x_R + 0.35 * cw], [y, y], color="#666666",
                lw=0.7, ls=ls, transform=ax.transAxes, zorder=0)

    # ── Loop bases on top of the stem (no arc drawn) ─────────────────────────
    # Loop ends sit roughly above the two stem columns (so it reads as connected)
    # and only widen slightly for large loops; shallow and anchored just above the
    # topmost stem base so the loop stays close to the stem on the y axis.
    loop_rx = (stem_gap / 2 + off) + 0.10 * cw * max(0, n_loop - 4)
    loop_ry = 0.022 + 0.5 * row_h
    y0_loop = stem_top - 0.6 * row_h
    for j, ch in enumerate(loop):
        t  = np.pi - (j + 0.5) / n_loop * np.pi  # left -> top -> right (5'->3')
        bx = cx + loop_rx * np.cos(t)
        by = y0_loop + loop_ry * np.sin(t)
        ax.text(bx, by, ch, ha="center", va="center",
                fontsize=fs_loop, color=base_col(ch), fontdict=MONO,
                transform=ax.transAxes)

    # ── 5' tail: stream from stem base down-and-left ─────────────────────────
    tail5_show = tail5[-tail_show:] if len(tail5) > tail_show else tail5
    tail5_label = ("…" if len(tail5) > tail_show else "") + tail5_show
    # Bases nearest the stem are the last chars of tail5_label (3' end of tail5)
    char_w = cw
    n5 = len(tail5_label)
    # Place characters so the LAST char sits just left of x_L at stem_base
    # and earlier chars extend toward the lower-left corner.
    for k, ch in enumerate(reversed(tail5_label)):
        # k=0 is closest to stem; start beyond the bottom stem letter so the
        # horizontal tail does not overlap the bottom row of the stem
        bx = x_L - (off + 1.7 * cw) - k * char_w
        by = stem_base
        if bx < 0.02 or by < 0.02:
            break
        col = base_col(ch) if ch != "…" else "#888888"
        ax.text(bx, by, ch, ha="center", va="center",
                fontsize=fs_nt, color=col, fontdict=MONO,
                transform=ax.transAxes)

    # 5' label near lower-left
    if not minimal:
        ax.text(0.015, max(0.02, stem_base - 0.012 * min(n5, 12) - 0.02),
                "5'", ha="left", va="center",
                fontsize=9, color="#555555", transform=ax.transAxes)

    # ── 3' tail: stream from stem base down-and-right (poly-U region) ────────
    # The 3' tail is the terminator poly-U at the mRNA 3' END, so NO trailing "…"
    # (the hidden part of tail3 is downstream genomic sequence past the transcript end).
    tail3_label = tail3[:tail_show] if len(tail3) > tail_show else tail3
    n3 = len(tail3_label)
    for k, ch in enumerate(tail3_label):
        bx = x_R + (off + 1.7 * cw) + k * char_w
        by = stem_base
        if bx > 0.98 or by < 0.02:
            break
        col = base_col(ch) if ch != "…" else "#888888"
        ax.text(bx, by, ch, ha="center", va="center",
                fontsize=fs_nt, color=col, fontdict=MONO,
                transform=ax.transAxes)

    if not minimal:
        ax.text(0.985, max(0.02, stem_base - 0.012 * min(n3, 12) - 0.02),
                "3'", ha="right", va="center",
                fontsize=9, color="#555555", transform=ax.transAxes)

    # ── Title with operon hit(s) ─────────────────────────────────────────────
    if not minimal:
        op_str = ""
        if operon_ids:
            op_str = "operon: " + ",".join(operon_ids)
        conf_str = (f"conf={term_row['conf']}  "
                    f"ΔG={term_row['hp_score']} kcal/mol  "
                    f"tail={term_row['tail_score']}")
        parts = [p for p in [title, op_str, conf_str] if p]
        ax.set_title("\n".join(parts), fontsize=7.5, pad=4)


# -- Individual plots (hairpin only: born at 1x1 in, no title/5'/3' labels) ---
# Clear stale TERM_*.pdf from earlier runs so the folder reflects the current set.
_hp_dir = "annotation/canonical/tts_hairpins"
for _f in os.listdir(_hp_dir):
    if _f.startswith("TERM_") and _f.endswith(".pdf"):
        os.remove(os.path.join(_hp_dir, _f))
# Name each file by its operon(s) (and term_id) so it is findable by operon, e.g.
# TERM_82_OP_00099.pdf for the 0178 operon.
for _, t in tts_terms_canonical.iterrows():
    fig, ax = plt.subplots(figsize=(1, 1))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)   # axes fills figure
    draw_terminator_hairpin(ax, t, minimal=True)
    _ops = term_to_operons.get(int(t["term_id"]), [])
    op_tag = "_".join(_ops) if _ops else "noop"
    save_path = f"{_hp_dir}/TERM_{int(t['term_id'])}_{op_tag}.pdf"
    plt.savefig(save_path)                                  # exact 1x1 in
    plt.close()

print(f"Saved {len(tts_terms_canonical)} individual hairpin plots.")

# -- Summary grid panel -------------------------------------------------------
n = len(tts_terms_canonical)
ncols = 5
nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 5))
axes_flat = axes.flatten() if n > 1 else [axes]

for i, (_, t) in enumerate(tts_terms_canonical.iterrows()):
    ops = term_to_operons.get(int(t["term_id"]), [])
    draw_terminator_hairpin(
        axes_flat[i], t,
        title=f"TERM {int(t['term_id'])} ({t['strand']})",
        operon_ids=ops,
    )

for j in range(i + 1, len(axes_flat)):
    axes_flat[j].axis("off")

fig.suptitle(f"Intrinsic terminators at canonical operon TTS  (n={n})",
             fontsize=11, y=1.01)
plt.tight_layout()
plt.savefig("annotation/canonical/tts_hairpins/all_hairpins.pdf", dpi=300, bbox_inches="tight")
plt.show()
print("Saved: annotation/canonical/tts_hairpins/all_hairpins.pdf")


# ============================================================
# Terminator statistics (R1 panel e): stem length, loop length, 3' poly-U tail
# ============================================================
#
# Three strip subplots, each born at (7/3, 7/9) in, to be assembled in
# Illustrator into one ~(7/3 x 7/3) panel. Computed over the unique intrinsic
# terminators at canonical operon TTSs (tts_terms_canonical). Base colours and
# the tail logo use the exact logomaker 'classic' scheme, matching the promoter
# logos in panel d.
#   1. term_stem_length.pdf  -- stem length (base pairs)
#   2. term_loop_length.pdf  -- loop length (nt)
#   3. term_tail3_logo.pdf   -- 3' tail composition logo (poly-U + A/C readthrough)

import logomaker

CLASSIC = {"A": "#008000", "C": "#0000ff", "G": "#ffa600", "U": "#ff0000"}

def _bp_count(s5, s3):
    """Number of base pairs = aligned non-gap positions (stem5[i] vs stem3[-1-i])."""
    return sum(1 for a, b in zip(str(s5), str(s3)[::-1]) if a != "-" and b != "-")

stem_bp  = tts_terms_canonical.apply(lambda r: _bp_count(r["stem5"], r["stem3"]), axis=1)
loop_len = tts_terms_canonical["loop"].astype(str).str.replace("-", "", regex=False).str.len()

# ── 1. Stem length (bp) ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7/3, 7/9), constrained_layout=True)
vc = stem_bp.value_counts().sort_index()
ax.bar(vc.index, vc.values, width=0.85, color="#377EB8", edgecolor="white", linewidth=0.3)
ax.axvline(stem_bp.median(), color="crimson", lw=0.8, ls="--")
ax.set_xlabel("Stem length (bp)", fontsize=6, labelpad=1.5)
ax.set_ylabel("Terminators", fontsize=6, labelpad=1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.savefig("annotation/canonical/term_stem_length.pdf", dpi=300)
plt.show()

# ── 2. Loop length (nt) ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7/3, 7/9), constrained_layout=True)
vc = loop_len.value_counts().sort_index()
ax.bar(vc.index, vc.values, width=0.85, color="#4DAF4A", edgecolor="white", linewidth=0.3)
ax.axvline(loop_len.median(), color="crimson", lw=0.8, ls="--")
ax.set_xlabel("Loop length (nt)", fontsize=6, labelpad=1.5)
ax.set_ylabel("Terminators", fontsize=6, labelpad=1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.savefig("annotation/canonical/term_loop_length.pdf", dpi=300)
plt.show()

# ── 3. 3' poly-U tail composition logo ───────────────────────────────────────
TAILLEN = 15   # TransTermHP reports a fixed 15-nt 3' tail window
tails = [str(t) for t in tts_terms_canonical["tail3"]
         if len(str(t)) == TAILLEN and set(str(t)) <= set("ACGU")]
tail_counts = pd.DataFrame(0.0, index=range(1, TAILLEN + 1), columns=["A", "C", "G", "U"])
for t in tails:
    for i, ch in enumerate(t):
        tail_counts.loc[i + 1, ch] += 1
tail_freq = tail_counts.div(tail_counts.sum(axis=1), axis=0)

# Poly-U tract length = leading run of consecutive U from the stem (stop at first
# non-U), the canonical intrinsic-terminator U-tail length.
def _polyU_run(t):
    n = 0
    for ch in t:
        if ch == "U":
            n += 1
        else:
            break
    return n
polyU_len = pd.Series([_polyU_run(t) for t in tails])
polyU_med = polyU_len.median()

fig, ax = plt.subplots(figsize=(7/3, 7/9), constrained_layout=True)
logomaker.Logo(tail_freq, ax=ax, color_scheme=CLASSIC, stack_order="big_on_top")
# median poly-U tract length (drawn at the position boundary it reaches)
ax.axvline(polyU_med + 0.5, color="crimson", lw=0.8, ls="--")
ax.set_ylim(0, 1)
ax.set_xticks([1, 5, 10, 15])
ax.set_xlabel("3' tail position (nt)", fontsize=6, labelpad=1.5)
ax.set_ylabel("Probability", fontsize=6, labelpad=1.5)
ax.tick_params(axis="both", labelsize=6, pad=1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.savefig("annotation/canonical/term_tail3_logo.pdf", dpi=300)
plt.show()

# ── 4. Stem G+C content (terminator stem vs the AT-rich genome) ──────────────
# Intrinsic-terminator stems are stabilized by G:C pairs; quantifying their G+C
# against the genome baseline tests whether the stems are genuinely GC-rich.
def _stem_gc(s5, s3):
    bases = [c for c in (str(s5) + str(s3)) if c in "ACGU"]
    return (sum(c in "GC" for c in bases) / len(bases)) if bases else np.nan
stem_gc   = tts_terms_canonical.apply(lambda r: _stem_gc(r["stem5"], r["stem3"]), axis=1)
genome_gc = (sum(genome[c].count("G") + genome[c].count("C") for c in genome)
             / sum(len(genome[c]) for c in genome))

fig, ax = plt.subplots(figsize=(7/3, 7/9), constrained_layout=True)
ax.hist(stem_gc * 100, bins=np.arange(0, 101, 10), color="#984ea3", edgecolor="white", linewidth=0.3)
ax.axvline(stem_gc.median() * 100, color="crimson", lw=0.8, ls="--")
ax.axvline(genome_gc * 100, color="black", lw=0.8, ls=":")   # genome baseline
ax.set_xlim(0, 100)
ax.set_xlabel("Stem G+C (%)", fontsize=6, labelpad=1.5)
ax.set_ylabel("Terminators", fontsize=6, labelpad=1.5)
ax.tick_params(axis="both", labelsize=6, pad=1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.savefig("annotation/canonical/term_stem_gc.pdf", dpi=300)
plt.show()

# ── 5. Distance from the mapped TTS to the 3' end of the poly-U tract ────────
# Hairpin spans [start0, end0); the poly-U tail begins at end0 (+ strand) or just
# below start0 (- strand), so the U-tract 3' end is end0 + polyU_run (+) or
# start0 - polyU_run (-). Signed distance in transcription direction:
#   d = strand_sign * (TTS - polyU_end);  d>0 = TTS lies 3' beyond the poly-U,
#   d=0 = TTS coincides with the poly-U end, d<0 = TTS falls short of it.
dist_rows = []
for _, o in op_canonical.iterrows():
    hit = find_terminators_near_tts(o, terms_df, window=TERM_WINDOW)
    if len(hit) == 0:
        continue
    best   = hit.loc[hit["conf"].idxmax()]
    run    = _polyU_run(str(best["tail3"]))
    s      = 1 if str(o["strand"]) == "+" else -1
    pu_end = (int(best["end0"]) + run) if s == 1 else (int(best["start0"]) - run)
    dist_rows.append(s * (int(o["tts"]) - pu_end))
tts_pu_dist = pd.Series(dist_rows, dtype=float)

lo = int(np.floor(tts_pu_dist.min() / 2.0) * 2)
hi = int(np.ceil(tts_pu_dist.max() / 2.0) * 2)
fig, ax = plt.subplots(figsize=(7/3, 7/9), constrained_layout=True)
ax.hist(tts_pu_dist, bins=np.arange(lo, hi + 2, 2), color="#FF7F00", edgecolor="white", linewidth=0.3)
ax.axvline(tts_pu_dist.median(), color="crimson", lw=0.8, ls="--")
ax.axvline(0, color="black", lw=0.6, ls=":")
ax.set_xlim(-15, 20)   # 96% lie within +/-10 nt; crop the few far outliers
ax.set_xlabel("TTS to poly-U end (nt)", fontsize=6, labelpad=1.5)
ax.set_ylabel("Terminators", fontsize=6, labelpad=1.5)
ax.tick_params(axis="both", labelsize=6, pad=1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.savefig("annotation/canonical/term_tts_polyU_distance.pdf", dpi=300)
plt.show()

# ── Stats -> append to Operon_Annotation.txt + console ───────────────────────
mean_U = tail_freq["U"].mean()
with open("Operon_Annotation.txt", "a") as fh:
    fh.write("\n\nTERMINATOR STATISTICS (canonical operon TTS intrinsic terminators)\n")
    fh.write("=" * 60 + "\n")
    fh.write(f"n unique terminators = {len(tts_terms_canonical)} "
             f"(tail logo built from {len(tails)} full 15-nt tails)\n\n")
    fh.write("Stem length (bp):\n" + stem_bp.describe().round(1).to_string() + "\n\n")
    fh.write("Loop length (nt):\n" + loop_len.describe().round(1).to_string() + "\n\n")
    fh.write("Poly-U tract length (leading consecutive U in the 3' tail):\n"
             + polyU_len.describe().round(1).to_string() + "\n")
    fh.write(f"  median poly-U tract length = {polyU_med:.0f} nt\n\n")
    fh.write(f"Mean U fraction across the 15-nt 3' tail: {mean_U:.2f}\n")
    fh.write("3' tail per-position U fraction:\n" + tail_freq["U"].round(2).to_string() + "\n\n")
    fh.write(f"Stem G+C content (%): median {stem_gc.median()*100:.0f}%, "
             f"mean {stem_gc.mean()*100:.0f}%  (genome G+C = {genome_gc*100:.0f}%)\n")
    fh.write((stem_gc * 100).describe().round(1).to_string() + "\n\n")
    within10 = (tts_pu_dist.abs() <= 10).mean()
    fh.write("Signed distance from mapped TTS to 3' end of poly-U tract "
             "(nt; + = TTS beyond poly-U, 0 = coincident):\n")
    fh.write(tts_pu_dist.describe().round(1).to_string() + "\n")
    fh.write(f"  median = {tts_pu_dist.median():.0f} nt; within +/-10 nt of the poly-U end: {within10:.0%}\n")
print(f"\nStem length (bp): median {stem_bp.median():.0f}, range {stem_bp.min()}-{stem_bp.max()}")
print(f"Loop length (nt): median {loop_len.median():.0f}, range {loop_len.min()}-{loop_len.max()}")
print(f"Poly-U tract length (nt): median {polyU_med:.0f}, range {polyU_len.min()}-{polyU_len.max()}")
print(f"Mean 3'-tail U fraction: {mean_U:.2f}")
print(f"Stem G+C: median {stem_gc.median()*100:.0f}% vs genome {genome_gc*100:.0f}%")
print(f"TTS to poly-U end (nt): median {tts_pu_dist.median():.0f}, "
      f"IQR {tts_pu_dist.quantile(.25):.0f} to {tts_pu_dist.quantile(.75):.0f}, "
      f"within +/-10 nt: {(tts_pu_dist.abs()<=10).mean():.0%}")
print("Saved: term_stem_length.pdf, term_loop_length.pdf, term_tail3_logo.pdf, "
      "term_stem_gc.pdf, term_tts_polyU_distance.pdf")


# In[ ]:




