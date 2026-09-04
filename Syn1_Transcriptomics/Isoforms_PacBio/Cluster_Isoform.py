# ============================================================
# Cluster_Isoform_New.py
#
# Re-implementation of isoform clustering for PacBio FLNC reads.
# It takes ~10 minutes to run through.
#
# -------------------------------------------------------------
# Part 1 — WHY we need clustering (kept for the record)
#
# A "raw isoform" = unique (chrom, strand, pos5p0, pos3p0) tuple
# across primary FLNC alignments; pos5p0/pos3p0 are transcript-
# oriented aligned ends (soft-clips excluded, polyA already trimmed
# upstream in PacBio_Processing.ipynb).
#
# Empirical findings from 2,616,593 primary reads on Syn1:
#
#   * 724,453 unique raw tuples   (79.4% are singletons)
#   * Using a bare read-count threshold of n>=10 discards 988,591
#     reads (37.8% of all reads) — those reads are NOT noise in
#     aggregate; they live on tuples adjacent to high-count tuples
#     and represent biological ±wobble around one real TSS/TTS.
#   * Order-independent single-linkage clustering collapses the
#     17,291 n>=10 tuples to 4,002 clusters at eps=10 bp — nearly
#     identical to the old two-stage pipeline's 4,064 (<2% off),
#     confirming the greedy post-merge artifact is small but the
#     algorithm itself is sound.
#   * For the operon-segmentation set (n>=50), 3,736 raw tuples
#     collapse to 1,041 clusters at eps=10 bp — a 3.6x collapse
#     showing that most of the 3,736 are near-duplicates of the
#     same underlying transcript.
#
# Conclusion: clustering is required for (a) read-recovery for
# quantification, and (b) deduplication for downstream operon /
# gene-overlap analysis. Threshold alone is insufficient.
#
# -------------------------------------------------------------
# Part 2 — Choice of eps = 10 bp
#
# Tolerance sweep on a prototype run (single-linkage, Chebyshev):
#
#   tol   clusters (n>=50)   clusters (n>=10)
#    5    1,148              4,947
#   10    1,041              4,002
#   15      993              3,558
#   20      965              3,307
#   50      842              2,503
#
# The curve bends between 10 and 15 bp and flattens beyond 20 bp.
# 10 bp matches PacBio CCS end-position error, the operon-
# segmentation duplication tolerance already in use, and returns
# biologically sensible counts (~1 k distinct transcription units
# for ~911 genes organised into operons).
#
# -------------------------------------------------------------
# Part 3 — Why complete-linkage (not single-linkage)
#
# Single-linkage allows transitive chaining: if A-B and B-C are
# each within eps but A-C is not, {A,B,C} still forms one cluster.
# In dense highly-expressed operons (e.g. lpdA-pta-ackA), the dense
# carpet of degradation-product 5' ends creates unbroken chains of
# <=10-bp links that merge multiple real TSSs into one giant cluster
# (empirically: one single-linkage cluster spanning 4,233 bp on 5',
#  with 5' MAD = 839 bp).
#
# Complete-linkage under Chebyshev distance requires that every pair
# of members in a cluster be within eps on both axes. Equivalently,
# the cluster's bounding box has diameter <= eps on each axis. We
# enforce this during the sweep using a bounding-box-augmented
# union-find: a union is accepted only if the merged cluster's
# bounding box still fits in an eps-by-eps square. No chaining is
# possible, and the algorithm stays O(n log n) time, O(n) memory.
#
# -------------------------------------------------------------
# Part 4 — Algorithm implemented here
#
#   1. Stream the BAM once; for each primary read record:
#      (chrom, strand, pos5p0, pos3p0, g_start0, g_end0, read_id)
#   2. Aggregate to unique (chrom, strand, pos5p0, pos3p0) tuples
#      (still storing read_id / bounds lists for membership).
#   3. Per (chrom, strand): sweep-line + bounding-box union-find
#      (UnionFindBBox). For each tuple i, visit the next j with
#      |Δpos5p0| <= eps; if |Δpos3p0| <= eps, try_union(i, j, eps).
#      The union is rejected if it would blow the merged cluster's
#      bounding box past eps on either axis — this is exactly
#      complete-linkage with Chebyshev distance.
#   4. Aggregate each cluster:
#        n_reads           = sum over member tuples
#        pos5p0 / pos3p0   = weighted median (by read count)
#        start0 / end0     = weighted median of genomic bounds
#        MAD_5 / MAD_3     = weighted MAD of member positions
#      All clusters are retained (MIN_CLUSTER_READS=1) so the output
#      is complete — downstream code applies its own read-count
#      threshold (>=50 for operon segmentation, >=10 for plotting).
#   5. Annotate against GFF (sense/antisense/intergenic/mixed) —
#      reused verbatim from the original Cluster_Isoform.py.
#   6. Write isoform table (tsv + xlsx) and optional membership.
#   7. Emit a MAD diagnostic plot and a cluster-vs-raw read-count
#      distribution summary.
#
# Caching: if OUT_ISOFORMS_TSV already exists and SKIP_IF_EXISTS is
# True, clustering/annotation is skipped; the table is loaded from
# disk and only the downstream plots / summaries are regenerated.
# ============================================================

from __future__ import annotations
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pysam
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist


# ----------------------------
# Inputs / outputs
# ----------------------------
MOTHER_FOLDER = "../.."   # project root; this script lives at Syn1_Transcriptomics/Isoforms_PacBio/
BAM_PATH = MOTHER_FOLDER + "/Syn1_Transcriptomics/PacBio/PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam"
GFF_PATH = MOTHER_FOLDER + "/Genomes_Input/syn1.genes.gff3"

OUT_FOLDER = "."   # write outputs alongside this script (Syn1_Transcriptomics/Isoforms_PacBio/)
Path(OUT_FOLDER).mkdir(parents=True, exist_ok=True)

OUT_RAW_TSV        = OUT_FOLDER + "/raw_isoforms_precluster.tsv"
OUT_ISOFORMS_TSV   = OUT_FOLDER + "/isoform_clusters_annotated.tsv"
OUT_ISOFORMS_XLSX  = OUT_FOLDER + "/isoform_clusters_annotated.xlsx"
OUT_MEMBERSHIP_TSV = OUT_FOLDER + "/isoform_membership.tsv"   # None to disable

# ----------------------------
# Parameters
# ----------------------------
REQUIRE_PRIMARY    = True
EPS_BP             = 10         # single clustering tolerance (see Part 2)
MIN_CLUSTER_READS  = 1          # keep ALL clusters; filter downstream
RUN_ANALYSIS       = True       # emit the Part-1/Part-2 sweeps (can be disabled)
SKIP_IF_EXISTS     = False       # skip clustering+annotation if output CSV exists

# Annotation thresholds (reused from Cluster_Isoform.py)
GENE_FEATURE_TYPES = {"gene", "CDS"}
GENE_ID_KEYS       = ["locus_tag", "ID", "Name", "gene", "gene_id"]
GENE_NAME_KEYS     = ["gene", "Name", "product", "locus_tag", "ID"]
MIN_GENE_OVERLAP_BP = 50
PURE_FRAC           = 0.95


# ----------------------------
# Helpers
# ----------------------------
def pos_ends(read: pysam.AlignedSegment):
    r0 = read.reference_start
    r1 = read.reference_end
    if r0 is None or r1 is None or r1 <= r0:
        return None
    if read.is_reverse:
        return "-", int(r1), int(r0), int(r0), int(r1)   # strand, p5, p3, g_start0, g_end0
    return "+", int(r0), int(r1), int(r0), int(r1)


def overlap_len(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order].astype(float)
    cw = np.cumsum(w)
    cutoff = cw[-1] / 2.0
    return float(v[np.searchsorted(cw, cutoff)])


def weighted_mad(values: np.ndarray, weights: np.ndarray) -> float:
    if values.size <= 1:
        return 0.0
    med = weighted_median(values, weights)
    return float(weighted_median(np.abs(values - med), weights))


def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    if pd.isna(attr_str):
        return d
    for field in str(attr_str).split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        elif " " in field:
            k, v = field.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        d[k.strip()] = v.strip()
    return d


def first_attr(attrs: Dict[str, str], keys: List[str], default: str = "") -> str:
    for k in keys:
        if k in attrs and attrs[k] not in ("", ".", None):
            return str(attrs[k])
    return default


def load_genes_from_gff(gff_path: str, feature_types=GENE_FEATURE_TYPES) -> pd.DataFrame:
    rows = []
    with open(gff_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, source, feature, start1, end1, score, strand, phase, attrs = parts
            if feature not in feature_types or strand not in {"+", "-"}:
                continue
            try:
                start1 = int(start1); end1 = int(end1)
            except ValueError:
                continue
            if end1 < start1:
                continue
            a = parse_gff_attributes(attrs)
            gene_id = first_attr(a, GENE_ID_KEYS, default="")
            gene_name = first_attr(a, GENE_NAME_KEYS, default=gene_id)
            rows.append({
                "chrom": chrom, "start0": start1 - 1, "end0": end1, "strand": strand,
                "gene_id": gene_id or f"{feature}:{chrom}:{start1}-{end1}:{strand}",
                "gene_name": gene_name or gene_id,
                "feature_type": feature,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("No gene/CDS features parsed from GFF.")
    df = df[df["feature_type"] == ("gene" if (df["feature_type"] == "gene").any() else "CDS")].copy()
    return df.sort_values(["chrom", "start0", "end0", "strand"]).reset_index(drop=True)


def summarize_overlaps(chrom, strand, start0, end0, genes_chr) -> Dict[str, object]:
    start0 = int(start0); end0 = int(end0)
    iso_len = max(0, end0 - start0)
    if iso_len == 0:
        return {"class_main": "invalid", "class_detail": "invalid_interval",
                "isoform_len_bp": 0, "sense_overlap_bp": 0, "antisense_overlap_bp": 0,
                "intergenic_bp": 0, "frac_sense": 0.0, "frac_antisense": 0.0,
                "frac_intergenic": 0.0, "n_sense_genes": 0, "n_antisense_genes": 0,
                "sense_gene_ids": "", "sense_gene_names": "",
                "antisense_gene_ids": "", "antisense_gene_names": "",
                "nearest_left_gene_id": "", "nearest_left_gene_name": "",
                "nearest_left_dist_bp": np.nan,
                "nearest_right_gene_id": "", "nearest_right_gene_name": "",
                "nearest_right_dist_bp": np.nan}

    cand = genes_chr[(genes_chr["end0"] > start0) & (genes_chr["start0"] < end0)]
    sense_rows, anti_rows = [], []
    sense_bp = anti_bp = 0
    for _, g in cand.iterrows():
        ol = overlap_len(start0, end0, int(g["start0"]), int(g["end0"]))
        if ol < MIN_GENE_OVERLAP_BP:
            continue
        rec = {"gene_id": g["gene_id"], "gene_name": g["gene_name"], "overlap_bp": int(ol)}
        if g["strand"] == strand:
            sense_rows.append(rec); sense_bp += ol
        else:
            anti_rows.append(rec); anti_bp += ol
    sense_bp = min(sense_bp, iso_len); anti_bp = min(anti_bp, iso_len)

    merged = []
    for s, e in cand[["start0", "end0"]].sort_values(["start0", "end0"]).to_numpy():
        s = max(start0, int(s)); e = min(end0, int(e))
        if e <= s:
            continue
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    covered = sum(e - s for s, e in merged)
    intergenic_bp = max(0, iso_len - covered)

    frac_s = sense_bp / iso_len; frac_a = anti_bp / iso_len; frac_i = intergenic_bp / iso_len
    if frac_i >= PURE_FRAC:
        class_main = "intergenic"
    elif frac_s >= PURE_FRAC:
        class_main = "sense"
    elif frac_a >= PURE_FRAC:
        class_main = "antisense"
    else:
        class_main = "mixed"

    if class_main == "intergenic":
        detail = "pure_intergenic"
    elif class_main == "sense":
        detail = "sense_single_gene" if len(sense_rows) == 1 else "sense_multigene"
    elif class_main == "antisense":
        detail = "antisense_single_gene" if len(anti_rows) == 1 else "antisense_multigene"
    else:
        if frac_s > 0 and frac_i > 0 and frac_a == 0:
            detail = "mixed_sense_intergenic"
        elif frac_a > 0 and frac_i > 0 and frac_s == 0:
            detail = "mixed_antisense_intergenic"
        elif frac_s > 0 and frac_a > 0 and frac_i == 0:
            detail = "mixed_sense_antisense"
        else:
            detail = "mixed_complex"

    left = genes_chr[genes_chr["end0"] <= start0]
    right = genes_chr[genes_chr["start0"] >= end0]
    nl_id = nl_name = nr_id = nr_name = ""
    nl_dist = nr_dist = np.nan
    if not left.empty:
        lg = left.iloc[left["end0"].argmax()]
        nl_id = lg["gene_id"]; nl_name = lg["gene_name"]; nl_dist = int(start0 - lg["end0"])
    if not right.empty:
        rg = right.iloc[right["start0"].argmin()]
        nr_id = rg["gene_id"]; nr_name = rg["gene_name"]; nr_dist = int(rg["start0"] - end0)

    return {
        "isoform_len_bp": iso_len,
        "sense_overlap_bp": int(sense_bp), "antisense_overlap_bp": int(anti_bp),
        "intergenic_bp": int(intergenic_bp),
        "frac_sense": float(frac_s), "frac_antisense": float(frac_a), "frac_intergenic": float(frac_i),
        "n_sense_genes": len(sense_rows), "n_antisense_genes": len(anti_rows),
        "sense_gene_ids": ",".join(x["gene_id"] for x in sense_rows),
        "sense_gene_names": ",".join(x["gene_name"] for x in sense_rows),
        "antisense_gene_ids": ",".join(x["gene_id"] for x in anti_rows),
        "antisense_gene_names": ",".join(x["gene_name"] for x in anti_rows),
        "class_main": class_main, "class_detail": detail,
        "nearest_left_gene_id": nl_id, "nearest_left_gene_name": nl_name,
        "nearest_left_dist_bp": nl_dist,
        "nearest_right_gene_id": nr_id, "nearest_right_gene_name": nr_name,
        "nearest_right_dist_bp": nr_dist,
    }


# ============================================================
# Collect reads → tuple table (with per-tuple read-id / bounds)
# ============================================================
if not os.path.exists(BAM_PATH + ".bai"):
    pysam.index(BAM_PATH)

bam = pysam.AlignmentFile(BAM_PATH, "rb")

# Per-tuple: count + read_ids + genomic start/end lists
tuple_data: Dict[Tuple[str, str, int, int], Dict[str, list]] = defaultdict(
    lambda: {"read_ids": [], "gstart": [], "gend": []}
)
n_reads_used = 0
n_reads_skipped = 0

for chrom in bam.references:
    try:
        it = bam.fetch(chrom)
    except ValueError:
        continue
    for read in it:
        if read.is_unmapped:
            n_reads_skipped += 1; continue
        if REQUIRE_PRIMARY and (read.is_secondary or read.is_supplementary):
            n_reads_skipped += 1; continue
        info = pos_ends(read)
        if info is None:
            n_reads_skipped += 1; continue
        strand, p5, p3, gs, ge = info
        key = (chrom, strand, p5, p3)
        d = tuple_data[key]
        d["read_ids"].append(read.query_name)
        d["gstart"].append(gs)
        d["gend"].append(ge)
        n_reads_used += 1

bam.close()

tuple_rows = []
for (chrom, strand, p5, p3), d in tuple_data.items():
    tuple_rows.append({
        "chrom": chrom, "strand": strand, "pos5p0": p5, "pos3p0": p3,
        "n_reads": len(d["read_ids"]),
        "read_ids": d["read_ids"], "gstart": d["gstart"], "gend": d["gend"],
    })
df_tup = (
    pd.DataFrame(tuple_rows)
    .sort_values(["chrom", "strand", "pos5p0", "pos3p0"])
    .reset_index(drop=True)
)

n_unique = len(df_tup)
singletons = int((df_tup["n_reads"] == 1).sum())
print(f"Primary reads used       : {n_reads_used:,}")
print(f"Reads skipped            : {n_reads_skipped:,}")
print(f"Unique raw isoforms      : {n_unique:,}")
print(f"  singletons (n_reads=1) : {singletons:,}  ({singletons/n_unique:.1%})")
print(f"  multi-read (n_reads>=2): {n_unique-singletons:,}")

df_tup[["chrom", "strand", "pos5p0", "pos3p0", "n_reads"]].to_csv(
    OUT_RAW_TSV, sep="\t", index=False
)
print(f"\nWrote {OUT_RAW_TSV}  (rows={len(df_tup):,})")


# ============================================================
# Part-1 / Part-2 analysis sweeps (justification for clustering
# and for eps=10). Disable with RUN_ANALYSIS=False.
# ============================================================
if RUN_ANALYSIS:
    print("\nRead-count distribution across raw isoforms:")
    print(df_tup["n_reads"].describe(percentiles=[0.5, 0.9, 0.99]).to_string())

    print("\n--- Threshold sweep (no binning, no merging) ---")
    total_tuples = len(df_tup)
    print(f"{'min_reads':>9}  {'n_unique_isoforms':>18}  {'frac_of_raw':>12}")
    for t in [1, 2, 3, 5, 10, 20, 50, 100]:
        n = int((df_tup["n_reads"] >= t).sum())
        print(f"{t:>9}  {n:>18,}  {n/total_tuples:>12.4f}")

    total_reads = int(df_tup["n_reads"].sum())
    print("\n--- Reads kept vs discarded at each threshold ---")
    print(f"Total primary reads: {total_reads:,}")
    print(f"{'min_reads':>9}  {'n_isoforms':>12}  {'reads_kept':>14}  {'reads_discarded':>16}  {'frac_discarded':>14}")
    for t in [1, 2, 3, 5, 10, 20, 50, 100]:
        sel = df_tup["n_reads"] >= t
        kept = int(df_tup.loc[sel, "n_reads"].sum())
        disc = total_reads - kept
        print(f"{t:>9}  {int(sel.sum()):>12,}  {kept:>14,}  {disc:>16,}  {disc/total_reads:>14.4f}")

    def _count_clusters(d: pd.DataFrame, tol: int) -> int:
        if d.empty:
            return 0
        if tol <= 0:
            return d[["chrom", "strand", "pos5p0", "pos3p0"]].drop_duplicates().shape[0]
        total = 0
        for _, g in d.groupby(["chrom", "strand"], sort=False):
            pts = g[["pos5p0", "pos3p0"]].to_numpy()
            if len(pts) == 1:
                total += 1; continue
            Z = linkage(pdist(pts, metric="chebyshev"), method="single")
            total += int(fcluster(Z, t=tol, criterion="distance").max())
        return total

    for thr in (10, 50):
        sub = df_tup[df_tup["n_reads"] >= thr]
        print(f"\n--- Single-linkage clustering of n>={thr} tuples (n={len(sub):,}) ---")
        print(f"{'tol_bp':>6}  {'n_clusters':>11}")
        for tol in [0, 5, 10, 15, 20, 30, 50]:
            print(f"{tol:>6}  {_count_clusters(sub, tol):>11,}")


# ============================================================
# Part 3 — The clustering itself (single-linkage, eps=EPS_BP)
# Applied to ALL tuples; no read-count filter applied here so the
# output is complete — downstream users pick their own threshold.
#
# Implementation: sweep-line + union-find.
#   Sort tuples by pos5p0 within (chrom, strand). Two tuples i,j are
#   directly connected iff |Δpos5p0| ≤ eps AND |Δpos3p0| ≤ eps.
#   Sweep: for each i, union with every j > i whose pos5p0 is within
#   eps of i's (a bounded look-ahead window), checking the pos3p0
#   constraint. Equivalent to single-linkage with Chebyshev distance.
#
# If OUT_ISOFORMS_TSV already exists and SKIP_IF_EXISTS is True, we
# skip re-clustering / re-annotation and just load it for plotting.
# ============================================================
_need_compute = not (SKIP_IF_EXISTS and os.path.exists(OUT_ISOFORMS_TSV))
if not _need_compute:
    print(f"\n=== Loading existing {OUT_ISOFORMS_TSV} (skipping clustering) ===")
    isoforms_annot = pd.read_csv(OUT_ISOFORMS_TSV, sep="\t")
    print(f"Loaded {len(isoforms_annot):,} clusters")

class UnionFindBBox:
    """
    Union-find that tracks the bounding box (min5/max5/min3/max3) of each
    cluster at its root. Used for complete-linkage clustering under
    Chebyshev distance: two clusters merge only if the resulting bounding
    box has diameter <= eps on both axes — i.e., every pair of members
    (including across the two clusters) is within eps on both axes.
    """
    __slots__ = ("parent", "rank", "min5", "max5", "min3", "max3")
    def __init__(self, p5: np.ndarray, p3: np.ndarray):
        n = len(p5)
        self.parent = np.arange(n, dtype=np.int64)
        self.rank = np.zeros(n, dtype=np.int32)
        self.min5 = p5.copy().astype(np.int64)
        self.max5 = p5.copy().astype(np.int64)
        self.min3 = p3.copy().astype(np.int64)
        self.max3 = p3.copy().astype(np.int64)
    def find(self, x):
        p = self.parent
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    def try_union(self, a, b, eps):
        """Union a and b only if the merged cluster has diameter <= eps
        on both axes (= complete-linkage with Chebyshev metric)."""
        ra = self.find(a); rb = self.find(b)
        if ra == rb:
            return True
        new_min5 = min(self.min5[ra], self.min5[rb])
        new_max5 = max(self.max5[ra], self.max5[rb])
        new_min3 = min(self.min3[ra], self.min3[rb])
        new_max3 = max(self.max3[ra], self.max3[rb])
        if (new_max5 - new_min5) > eps or (new_max3 - new_min3) > eps:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        self.min5[ra] = new_min5; self.max5[ra] = new_max5
        self.min3[ra] = new_min3; self.max3[ra] = new_max3
        return True


if _need_compute:
    print(f"\n=== Complete-linkage clustering of {len(df_tup):,} tuples at eps={EPS_BP} bp ===")
    p5_all = df_tup["pos5p0"].to_numpy()
    p3_all = df_tup["pos3p0"].to_numpy()
    uf = UnionFindBBox(p5_all, p3_all)

    for (chrom, strand), g in df_tup.groupby(["chrom", "strand"], sort=False):
        idx = g.index.to_numpy()
        n = len(idx)
        if n == 1:
            continue
        p5 = p5_all[idx]
        p3 = p3_all[idx]
        for i in range(n):
            j = i + 1
            while j < n and p5[j] - p5[i] <= EPS_BP:
                if abs(p3[j] - p3[i]) <= EPS_BP:
                    uf.try_union(int(idx[i]), int(idx[j]), EPS_BP)
                j += 1

    roots = np.array([uf.find(i) for i in range(len(df_tup))], dtype=np.int64)
    _, cluster_labels = np.unique(roots, return_inverse=True)
    df_tup["cluster_id"] = cluster_labels.astype(np.int64)
    print(f"Total clusters: {int(cluster_labels.max()) + 1:,}")

    # --- Aggregate each cluster ---
    cluster_rows = []
    membership_rows = []
    for cid, g in df_tup.groupby("cluster_id", sort=False):
        chrom = g["chrom"].iloc[0]
        strand = g["strand"].iloc[0]
        p5 = g["pos5p0"].to_numpy()
        p3 = g["pos3p0"].to_numpy()
        w = g["n_reads"].to_numpy()

        gs_all = np.concatenate([np.asarray(x, dtype=int) for x in g["gstart"].tolist()])
        ge_all = np.concatenate([np.asarray(x, dtype=int) for x in g["gend"].tolist()])
        p5_per_read = np.concatenate([np.full(int(n), int(v), dtype=int) for v, n in zip(p5, w)])
        p3_per_read = np.concatenate([np.full(int(n), int(v), dtype=int) for v, n in zip(p3, w)])

        rep_p5 = int(np.median(p5_per_read))
        rep_p3 = int(np.median(p3_per_read))
        rep_gs = int(np.median(gs_all))
        rep_ge = int(np.median(ge_all))
        mad5 = float(np.median(np.abs(p5_per_read - rep_p5))) if len(p5_per_read) > 1 else 0.0
        mad3 = float(np.median(np.abs(p3_per_read - rep_p3))) if len(p3_per_read) > 1 else 0.0

        n_reads = int(w.sum())
        if n_reads < MIN_CLUSTER_READS or rep_ge <= rep_gs:
            continue

        cluster_rows.append({
            "cluster_id": int(cid),
            "chrom": chrom, "strand": strand,
            "pos5p0": rep_p5, "pos3p0": rep_p3,
            "start0": rep_gs, "end0": rep_ge,
            "n_reads": n_reads,
            "n_member_tuples": int(len(g)),
            "start_spread_mad_bp": mad5,
            "end_spread_mad_bp": mad3,
        })
        for _, row in g.iterrows():
            for rid, gs, ge in zip(row["read_ids"], row["gstart"], row["gend"]):
                membership_rows.append({
                    "cluster_id": int(cid),
                    "chrom": chrom, "strand": strand,
                    "read_id": rid,
                    "pos5p0": int(row["pos5p0"]), "pos3p0": int(row["pos3p0"]),
                    "start0": int(gs), "end0": int(ge),
                })

    clusters = pd.DataFrame(cluster_rows).reset_index(drop=True)
    clusters["isoform_id"] = [f"ISO_{i+1:06d}" for i in range(len(clusters))]
    print(f"Clusters retained (n_reads >= {MIN_CLUSTER_READS}): {len(clusters):,}")

    # --- Annotate against genes ---
    genes = load_genes_from_gff(GFF_PATH)
    print(f"Loaded gene annotations: n={len(genes)}")
    genes_by_chrom = {c: g.reset_index(drop=True) for c, g in genes.groupby("chrom", sort=False)}

    annot = []
    for _, r in clusters.iterrows():
        chrom = r["chrom"]
        if chrom not in genes_by_chrom:
            a = summarize_overlaps(chrom, r["strand"], int(r["start0"]), int(r["end0"]),
                                   pd.DataFrame(columns=["chrom", "start0", "end0", "strand",
                                                         "gene_id", "gene_name", "feature_type"]))
        else:
            a = summarize_overlaps(chrom, r["strand"], int(r["start0"]), int(r["end0"]),
                                   genes_by_chrom[chrom])
        rec = {
            "isoform_id": r["isoform_id"], "cluster_id": int(r["cluster_id"]),
            "chrom": chrom, "strand": r["strand"],
            "start0": int(r["start0"]), "end0": int(r["end0"]),
            "pos5p0": int(r["pos5p0"]), "pos3p0": int(r["pos3p0"]),
            "n_reads": int(r["n_reads"]),
            "n_member_tuples": int(r["n_member_tuples"]),
            "start_spread_mad_bp": float(r["start_spread_mad_bp"]),
            "end_spread_mad_bp": float(r["end_spread_mad_bp"]),
        }
        rec.update(a)
        annot.append(rec)

    isoforms_annot = pd.DataFrame(annot).sort_values(
        ["chrom", "start0", "end0", "strand", "n_reads"],
        ascending=[True, True, True, True, False]
    ).reset_index(drop=True)

    print("\nIsoform class counts:")
    print(isoforms_annot["class_main"].value_counts(dropna=False))

    # --- Output ---
    isoforms_annot.to_csv(OUT_ISOFORMS_TSV, sep="\t", index=False)
    with pd.ExcelWriter(OUT_ISOFORMS_XLSX, engine="openpyxl") as w:
        isoforms_annot.to_excel(w, sheet_name="isoforms", index=False)
        (isoforms_annot["class_main"].value_counts(dropna=False)
            .rename_axis("class_main").reset_index(name="n_isoforms")
            .to_excel(w, sheet_name="class_summary", index=False))
        (isoforms_annot["class_detail"].value_counts(dropna=False)
            .rename_axis("class_detail").reset_index(name="n_isoforms")
            .to_excel(w, sheet_name="detail_summary", index=False))
    print(f"\nwrote {OUT_ISOFORMS_TSV}")
    print(f"wrote {OUT_ISOFORMS_XLSX}")

    if OUT_MEMBERSHIP_TSV:
        cid_to_isoform = dict(zip(isoforms_annot["cluster_id"], isoforms_annot["isoform_id"]))
        mem_df = pd.DataFrame(membership_rows)
        mem_df = mem_df[mem_df["cluster_id"].isin(cid_to_isoform)].copy()
        mem_df["isoform_id"] = mem_df["cluster_id"].map(cid_to_isoform)
        mem_df = mem_df[["isoform_id", "cluster_id", "chrom", "strand",
                         "read_id", "pos5p0", "pos3p0", "start0", "end0"]]
        mem_df.to_csv(OUT_MEMBERSHIP_TSV, sep="\t", index=False)
        print(f"wrote {OUT_MEMBERSHIP_TSV}  (rows={len(mem_df):,})")

print(isoforms_annot.head(10))


# ============================================================
# Part 4 — MAD_5' and MAD_3' distributions across clusters
# (tells us how tight each cluster's ends are; a sanity check
# on EPS_BP and a diagnostic for residual 3'-end wobble, e.g.
# from imperfect polyA trimming or real TTS heterogeneity).
# ============================================================
import matplotlib.pyplot as plt

mad5 = isoforms_annot["start_spread_mad_bp"].to_numpy()
mad3 = isoforms_annot["end_spread_mad_bp"].to_numpy()

OUT_MAD_PDF = OUT_FOLDER + "/cluster_mad_distributions.pdf"

bin_max = int(np.ceil(max(np.nanmax(mad5), np.nanmax(mad3), EPS_BP))) + 1
bins = np.arange(0, bin_max + 1, 0.5)

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.hist(mad5, bins=bins, color="#2b6cb0", alpha=0.45,
        label=f"5' end MAD (n={len(mad5):,})")
ax.hist(mad3, bins=bins, color="#c0392b", alpha=0.45,
        label=f"3' end MAD (n={len(mad3):,})")
ax.axvline(EPS_BP, ls="--", color="grey", lw=0.8, label=f"eps = {EPS_BP} bp")
ax.set_xlabel("MAD (bp)")
ax.set_ylabel("# clusters (log)")
ax.set_yscale("log")
ax.set_title("Cluster end-position spread")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT_MAD_PDF)
plt.close(fig)
print(f"wrote {OUT_MAD_PDF}")

print("\n5' MAD summary (bp):")
print(pd.Series(mad5).describe(percentiles=[0.5, 0.9, 0.99]).to_string())
print("\n3' MAD summary (bp):")
print(pd.Series(mad3).describe(percentiles=[0.5, 0.9, 0.99]).to_string())


# ============================================================
# Part 5 — Read-count distribution: clusters vs. raw tuples
# Shows how clustering redistributes reads from a very long-tailed
# single-bp tuple distribution into a tighter cluster distribution.
# ============================================================
OUT_READCOUNT_PDF = OUT_FOLDER + "/readcount_distribution_cluster_vs_raw.pdf"
OUT_READCOUNT_CSV = OUT_FOLDER + "/readcount_distribution_cluster_vs_raw.csv"

raw_counts = df_tup["n_reads"].to_numpy()
cluster_counts = isoforms_annot["n_reads"].to_numpy()

# tabulate side-by-side at each min-reads threshold
thresholds = [1, 2, 3, 5, 10, 20, 50, 100, 500, 1000]
rc_rows = []
for t in thresholds:
    rc_rows.append({
        "min_reads": t,
        "raw_tuples_n":     int((raw_counts     >= t).sum()),
        "clusters_n":       int((cluster_counts >= t).sum()),
        "raw_reads_kept":   int(raw_counts[raw_counts >= t].sum()),
        "cluster_reads_kept": int(cluster_counts[cluster_counts >= t].sum()),
    })
rc_df = pd.DataFrame(rc_rows)
rc_df.to_csv(OUT_READCOUNT_CSV, index=False)
print(f"\nwrote {OUT_READCOUNT_CSV}")
print(rc_df.to_string(index=False))

print("\nRaw-tuple  n_reads summary:")
print(pd.Series(raw_counts).describe(percentiles=[0.5, 0.9, 0.99]).to_string())
print("\nCluster    n_reads summary:")
print(pd.Series(cluster_counts).describe(percentiles=[0.5, 0.9, 0.99]).to_string())

# Log-log histogram: raw tuples vs clusters, same axes
log_bins = np.logspace(0, np.log10(max(raw_counts.max(), cluster_counts.max())) + 0.05, 60)
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.hist(raw_counts, bins=log_bins, color="#888888", alpha=0.55,
        label=f"raw tuples (n={len(raw_counts):,})")
ax.hist(cluster_counts, bins=log_bins, color="#2b6cb0", alpha=0.65,
        label=f"clusters (n={len(cluster_counts):,})")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("n_reads per entity")
ax.set_ylabel("# entities")
ax.set_title("Read-count distribution: raw tuples vs clusters")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT_READCOUNT_PDF)
plt.close(fig)
print(f"\nwrote {OUT_READCOUNT_PDF}")
