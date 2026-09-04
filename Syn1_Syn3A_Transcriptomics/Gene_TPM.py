#!/usr/bin/env python3
"""Per-gene TPM for every RNA-seq library of JCVI-syn1.0 and JCVI-syn3A.

Background
----------
Both organisms were sequenced on several platforms: syn1 with Illumina (three
libraries), PacBio Iso-Seq, and two ONT direct-RNA runs; syn3A with Illumina and
one ONT run. This script quantifies all of them against their own reference
annotation and emits exactly one table per organism.

It replaces three earlier scripts that split the same work by organism and by
purpose (Syn1_Transcriptomics/Gene_TPM/Gene_Transcriptomics.py,
Syn3A_Transcriptomics/Gene_TPM/Syn3A_TPM.py, and compute_platform_TPM.py) and
emitted seven overlapping tables between them.

Algorithm
---------
Genes are parsed from the GFF3 `gene` and `pseudogene` feature types. syn1 has
911 `gene` and no `pseudogene` features; syn3A has 493 `gene` plus 3
`pseudogene` (0051, 0546, 0602), for 496 loci. syn3A carries `product=` on CDS
rows rather than gene rows, so products are mapped in from CDS.

Per library, per-base depth is read from the strand-split bedGraphs and turned
into a prefix sum, so a gene's mean depth is one subtraction. Sense is the plus
track for + genes and the minus track for - genes; antisense is the other.
Length normalisation is implicit in mean depth (read-bases / length), and TPM
uses a single sense-plus-antisense denominator per library:

    TPM = mean_depth / sum(all mean_depths, both strands) * 1e6

The three syn1 Illumina libraries are averaged in two steps -- technical
replicates within a biological sample first (SRR35996296 + SRR35996297 ->
sample_95; SRR35996298 -> sample_enr), then the two biological samples with
equal weight -> avg_sense_TPM. This keeps the two technical replicates from
double-weighting their biological sample.

All correlations use one convention: Pearson on log10(TPM), over loci with
TPM > 0 in both libraries being compared. No pseudo-count is added -- a constant
offset on a TPM scale inflates the dynamic range at the zero end and lets
unexpressed genes drive the fit. The cost is that the locus count varies by pair,
since a zero in either library drops that gene. Spearman is reported alongside for
the syn3A comparisons and is rank-based, so the transform does not affect it.
Every value is written to the log rather than left in memory.

Outputs
-------
- syn1_TPM.tsv   911 loci: annotation, per-library sense/antisense mean depth and
                 TPM for 3 Illumina libraries, PacBio, and ONT rep1/rep2/merged,
                 plus the two-step Illumina averages.
- syn3A_TPM.tsv  496 loci: annotation, sense/antisense mean depth and TPM for
                 Illumina and ONT.
- Gene_TPM.txt   loci counts and every correlation.

Run from Syn1_Syn3A_Transcriptomics/ in the Omics conda env.
"""
from __future__ import annotations

import os
from itertools import combinations
from typing import Dict

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

SYN1_GFF = os.path.join(ROOT, "Genomes_Input/syn1.genes.gff3")
SYN3A_GFF = os.path.join(ROOT, "Genomes_Input/syn3a_genome.gff3")

SYN1_ILLUMINA_DBG = os.path.join(ROOT, "Syn1_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph")
SYN1_PACBIO_DBG = os.path.join(ROOT, "Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph")
SYN1_ONT_DBG = os.path.join(ROOT, "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph")
SYN3A_ILLUMINA_DBG = os.path.join(ROOT, "Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph")
SYN3A_ONT_DBG = os.path.join(ROOT, "Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph")

PALSSON_CSV = os.path.join(HERE, "Processed_TPM_Palsson/GSM6204176_3A.csv")

OUT_SYN1 = os.path.join(HERE, "syn1_TPM.tsv")
OUT_SYN3A = os.path.join(HERE, "syn3A_TPM.tsv")
OUT_TXT = os.path.join(HERE, "Gene_TPM.txt")

# syn3A annotates 3 pseudogenes only under the `pseudogene` feature type; syn1
# has none, so one shared set gives 911 for syn1 and 496 for syn3A.
PRIMARY_FEATURES = {"gene", "pseudogene"}

ILLUMINA_SAMPLES = ["SRR35996296", "SRR35996297", "SRR35996298"]
SAMPLE_TO_BIO = {
    "SRR35996296": "sample_95",   # technical replicate pair
    "SRR35996297": "sample_95",
    "SRR35996298": "sample_enr",  # second biological sample
}
_log: list[str] = []


def say(msg: str = "") -> None:
    """Print and record for the .txt log."""
    print(msg)
    _log.append(msg)


# ---------------------------------------------------------------- annotation
def parse_gff_attributes(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in attr.split(";"):
        item = item.strip()
        if item and "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def read_genes_gff(path: str) -> pd.DataFrame:
    """Parse `gene` / `pseudogene` rows into a 0-based half-open frame."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in PRIMARY_FEATURES:
                continue
            chrom, _src, feature, start1, end1, _score, strand, _phase, attrs = parts
            a = parse_gff_attributes(attrs)
            rows.append((chrom, feature, int(start1) - 1, int(end1), strand,
                         a.get("locus_tag", ""), a.get("gene", ""),
                         a.get("rna_type", ""), a.get("product", "")))
    df = pd.DataFrame(rows, columns=["chrom", "feature", "start0", "end0", "strand",
                                     "locus_tag", "gene_name", "rna_type", "gene_product"])
    df = df.sort_values(["chrom", "start0", "end0"]).reset_index(drop=True)
    df["gene_len"] = df["end0"] - df["start0"]
    return df


def cds_product_map(path: str) -> Dict[str, str]:
    """locus_tag -> product taken from CDS rows (syn3A puts it there)."""
    out: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            a = parse_gff_attributes(parts[8])
            lt = a.get("locus_tag", "")
            if lt and lt not in out:
                out[lt] = a.get("product", "")
    return out


# --------------------------------------------------------------------- depth
def load_depth_prefix(path: str) -> Dict[str, np.ndarray]:
    """chrom -> prefix sum of per-base depth, so a gene mean is one subtraction."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start0", "end0", "depth"])
    prefix = {}
    for chrom, sub in df.groupby("chrom", sort=False):
        arr = np.zeros(int(sub["end0"].max()), dtype=np.float32)
        for s, e, d in zip(sub["start0"].to_numpy(dtype=int),
                           sub["end0"].to_numpy(dtype=int),
                           sub["depth"].to_numpy(dtype=np.float32)):
            arr[s:e] = d
        prefix[chrom] = np.concatenate([[0.0], np.cumsum(arr, dtype=np.float64)])
    return prefix


def add_strand_tpm(genes: pd.DataFrame, plus_bg: str, minus_bg: str,
                   prefix: str) -> pd.DataFrame:
    """Add {prefix}_{sense,antisense}_{avg_depth,TPM} columns in place."""
    plus_p = load_depth_prefix(plus_bg)
    minus_p = load_depth_prefix(minus_bg)

    sense = np.zeros(len(genes), dtype=np.float64)
    anti = np.zeros(len(genes), dtype=np.float64)
    for i, row in genes.iterrows():
        c, s0, e0 = row["chrom"], int(row["start0"]), int(row["end0"])
        glen = max(1, e0 - s0)
        pm = (plus_p[c][e0] - plus_p[c][s0]) / glen
        mm = (minus_p[c][e0] - minus_p[c][s0]) / glen
        sense[i], anti[i] = (pm, mm) if row["strand"] == "+" else (mm, pm)

    genes[f"{prefix}_sense_avg_depth"] = sense
    genes[f"{prefix}_antisense_avg_depth"] = anti
    denom = sense.sum() + anti.sum()
    genes[f"{prefix}_sense_TPM"] = sense / denom * 1e6 if denom > 0 else 0.0
    genes[f"{prefix}_antisense_TPM"] = anti / denom * 1e6 if denom > 0 else 0.0
    return genes


# -------------------------------------------------------------- correlations
def r_log10(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    """Pearson on log10(TPM) over loci expressed (TPM > 0) in both libraries."""
    m = (a > 0) & (b > 0)
    if m.sum() < 3:
        return float("nan"), int(m.sum())
    return float(pearsonr(np.log10(a[m]), np.log10(b[m]))[0]), int(m.sum())


def r_log10_sp(a: np.ndarray, b: np.ndarray) -> tuple[float, float, int]:
    """r_log10 plus Spearman (rank-based, so unaffected by the log transform)."""
    m = (a > 0) & (b > 0)
    return (float(pearsonr(np.log10(a[m]), np.log10(b[m]))[0]),
            float(spearmanr(a[m], b[m])[0]), int(m.sum()))


# ========================================================================
# syn1
# ========================================================================
say("=" * 74)
say("syn1 -- Illumina x3, PacBio Iso-Seq, ONT direct-RNA x2 (+ merged)")
say("=" * 74)

syn1 = read_genes_gff(SYN1_GFF)
say(f"loci parsed from syn1.genes.gff3 : {len(syn1)}")
say("  feature breakdown: " + ", ".join(
    f"{k}={v}" for k, v in syn1["feature"].value_counts().items()))

for s in ILLUMINA_SAMPLES:
    add_strand_tpm(syn1, f"{SYN1_ILLUMINA_DBG}/{s}.plus.bedGraph",
                   f"{SYN1_ILLUMINA_DBG}/{s}.minus.bedGraph", prefix=s)
    say(f"  Illumina {s} : TPM computed")

# Two-step averaging: technical replicates first, then biological samples.
bio_groups: Dict[str, list[str]] = {}
for s in ILLUMINA_SAMPLES:
    bio_groups.setdefault(SAMPLE_TO_BIO[s], []).append(s)
for bio, reps in bio_groups.items():
    for st in ("sense", "antisense"):
        syn1[f"{bio}_{st}_TPM"] = syn1[[f"{r}_{st}_TPM" for r in reps]].mean(axis=1)
for st in ("sense", "antisense"):
    syn1[f"avg_{st}_TPM"] = syn1[[f"{b}_{st}_TPM" for b in sorted(bio_groups)]].mean(axis=1)
say(f"  Illumina averaged: {' + '.join(sorted(bio_groups))} -> avg_sense_TPM")

add_strand_tpm(syn1, f"{SYN1_PACBIO_DBG}/syn1.PacBio.FLNC.HQ.plus.bedGraph",
               f"{SYN1_PACBIO_DBG}/syn1.PacBio.FLNC.HQ.minus.bedGraph", prefix="PacBio")
say("  PacBio       : TPM computed")

for rep in ("rep1", "rep2", "merged"):
    add_strand_tpm(syn1, f"{SYN1_ONT_DBG}/syn1.ONT.{rep}.plus.bedGraph",
                   f"{SYN1_ONT_DBG}/syn1.ONT.{rep}.minus.bedGraph", prefix=f"ONT_{rep}")
    say(f"  ONT {rep:<7}: TPM computed")

syn1_cols = (["chrom", "feature", "start0", "end0", "strand", "locus_tag",
              "gene_name", "rna_type", "gene_product", "gene_len"]
             + [f"{s}_{k}" for s in ILLUMINA_SAMPLES
                for k in ("sense_avg_depth", "antisense_avg_depth",
                          "sense_TPM", "antisense_TPM")]
             + [f"{b}_{st}_TPM" for b in sorted(bio_groups) for st in ("sense", "antisense")]
             + ["avg_sense_TPM", "avg_antisense_TPM"]
             + [f"{p}_{k}" for p in ("PacBio", "ONT_rep1", "ONT_rep2", "ONT_merged")
                for k in ("sense_avg_depth", "antisense_avg_depth",
                          "sense_TPM", "antisense_TPM")])
syn1[syn1_cols].to_csv(OUT_SYN1, sep="\t", index=False, float_format="%.4f")
say(f"\nwrote {os.path.basename(OUT_SYN1)}  ({len(syn1)} loci, {len(syn1_cols)} columns)")

# ---- syn1 correlations ----
say("\n-- Illumina replicate agreement (Pearson on log10 TPM, both > 0) --")
for a, b in combinations(ILLUMINA_SAMPLES, 2):
    kind = "technical" if SAMPLE_TO_BIO[a] == SAMPLE_TO_BIO[b] else "biological"
    r, n = r_log10(syn1[f"{a}_sense_TPM"].to_numpy(), syn1[f"{b}_sense_TPM"].to_numpy())
    say(f"  {a} vs {b}  r = {r:.4f}   (n = {n}, {kind})")

say("\n-- Cross-platform agreement (Pearson on log10 TPM, both > 0) --")
plats = [("Illumina", "avg_sense_TPM"), ("PacBio", "PacBio_sense_TPM"),
         ("ONT_rep1", "ONT_rep1_sense_TPM"), ("ONT_rep2", "ONT_rep2_sense_TPM"),
         ("ONT_merged", "ONT_merged_sense_TPM")]
for (na, ca), (nb, cb) in combinations(plats, 2):
    r, n = r_log10(syn1[ca].to_numpy(), syn1[cb].to_numpy())
    say(f"  {na:<11} vs {nb:<11} r = {r:.4f}   (n = {n})")

# ========================================================================
# syn3A
# ========================================================================
say("\n" + "=" * 74)
say("syn3A -- Illumina, ONT direct-RNA")
say("=" * 74)

syn3a = read_genes_gff(SYN3A_GFF)
syn3a["gene_product"] = (syn3a["locus_tag"].map(cds_product_map(SYN3A_GFF))
                         .fillna(syn3a["gene_product"]))
say(f"loci parsed from syn3a_genome.gff3 : {len(syn3a)}")
say("  feature breakdown: " + ", ".join(
    f"{k}={v}" for k, v in syn3a["feature"].value_counts().items()))

add_strand_tpm(syn3a, f"{SYN3A_ILLUMINA_DBG}/syn3A_rep1.plus.bedGraph",
               f"{SYN3A_ILLUMINA_DBG}/syn3A_rep1.minus.bedGraph", prefix="Illumina")
say("  Illumina : TPM computed")
add_strand_tpm(syn3a, f"{SYN3A_ONT_DBG}/syn3A.ONT.rep1.plus.bedGraph",
               f"{SYN3A_ONT_DBG}/syn3A.ONT.rep1.minus.bedGraph", prefix="ONT")
say("  ONT      : TPM computed")

syn3a_cols = (["locus_tag", "gene_name", "gene_product", "chrom", "start0", "end0",
               "strand", "gene_len"]
              + [f"{p}_{k}" for p in ("Illumina", "ONT")
                 for k in ("sense_avg_depth", "antisense_avg_depth",
                           "sense_TPM", "antisense_TPM")])
syn3a[syn3a_cols].to_csv(OUT_SYN3A, sep="\t", index=False, float_format="%.4f")
say(f"\nwrote {os.path.basename(OUT_SYN3A)}  ({len(syn3a)} loci, {len(syn3a_cols)} columns)")

say("\n-- syn3A agreement (Pearson on log10 TPM, both > 0) --")
pr, sr, n = r_log10_sp(syn3a["Illumina_sense_TPM"].to_numpy(),
                       syn3a["ONT_sense_TPM"].to_numpy())
say(f"  Illumina vs ONT      r = {pr:.3f}   rho = {sr:.3f}   (n = {n})")

palsson = (pd.read_csv(PALSSON_CSV)
           .rename(columns={"Geneid": "locus_tag", "Illumina_TPM": "Palsson_Illumina_TPM"}))
cmp_p = syn3a.merge(palsson[["locus_tag", "Palsson_Illumina_TPM"]],
                    on="locus_tag", how="inner")
pr, sr, n = r_log10_sp(cmp_p["Illumina_sense_TPM"].to_numpy(),
                       cmp_p["Palsson_Illumina_TPM"].to_numpy())
say(f"  Illumina vs Palsson  r = {pr:.3f}   rho = {sr:.3f}   (n = {n}, "
    f"overlap {len(cmp_p)} of {len(palsson)} reported loci)")

with open(OUT_TXT, "w") as fh:
    fh.write("\n".join(_log) + "\n")
print(f"\nwrote {os.path.basename(OUT_TXT)}")
