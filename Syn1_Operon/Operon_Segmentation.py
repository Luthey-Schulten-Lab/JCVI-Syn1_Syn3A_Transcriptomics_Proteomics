#!/usr/bin/env python
# coding: utf-8

# ## Operon Segmentation — Isoform-based
# 
# Segment the *Syn1* minimal cell genome (~1.08 Mb, 911 genes) into operons directly from PacBio FLNC long-read RNA-seq.
# 
# **Biological rationale:** In bacteria, genes are co-transcribed in polycistronic operons from a shared promoter. A single full-length RNA isoform spanning multiple genes is therefore the natural evidence for operon membership. PacBio FLNC reads capture complete transcripts from TSS to TTS, making them ideal for operon discovery without relying on depth drop-offs or intergenic gap heuristics.
# 
# **Key insight:** Under containment-only clustering, each operon is defined by its **longest isoform** — the full-length transcript spanning the true TSS to TTS. All shorter co-clustered isoforms are nested subsets arising from RNase degradation or premature termination and do not change the boundary.
# 
# **Complications addressed:**
# - *RNase activity* — bacteria actively degrade mRNA from the 3′ end and internally; this fragments a single operon transcript into multiple shorter isoforms that can appear as separate clusters.
# - *Low-abundance transcripts* — essential genes under low expression may produce too few reads to pass the isoform threshold; these are rescued using BAM-level read evidence.
# - *rRNA operons* — rRNA transcripts are too abundant and structured for isoform clustering; these two operons (each encoding 5S–23S–16S) are annotated directly from the genome.
# 
# **Single tunable parameter:** `MIN_READS` — minimum read support per isoform cluster.
# 
# **Input:** `isoform_clusters_annotated.tsv` — 80 K clustered isoforms with genomic coordinates, strand, and read support.
#
# **Output:** `operons.candidate_blocks.tsv` — full operon table with gene annotations and
# segmentation type — is the ONLY file written to the main folder. Every other artifact
# (intermediate/decision TSVs, the merge-decision and conflict PDFs, IGV BED) is written to
# the `segmentation/` subfolder:
#   segmentation/operon_merge_decisions.tsv            — every merge candidate + verdict
#   segmentation/same_strand_conflicts.pdf             — OVERLAP merge-candidate panels
#   segmentation/operon_merge_decisions_gene_in_gap.pdf— GENE-IN-GAP merge-candidate panels
#   segmentation/same_strand_gene_conflicts.tsv        — shared-gene overlap pairs
#   segmentation/gene_operon_coverage.tsv              — per-gene coverage classification
#   segmentation/uncovered_genes.tsv                   — genes with no operon either strand
#   segmentation/operons.candidate_blocks.{plus,minus}.bed — IGV tracks
#
# **Step 5a merge rule (this revision):** overlapping OR gene-in-gap operon pairs are
# merged by direct co-transcription evidence — strand-specific PacBio bridging reads across
# the junction plus gap/flank read-depth continuity (params MERGE_*) — replacing the former
# TSS-location rule. Merging is pairwise with NO transitive chaining. Merged operons are then
# re-annotated from their spans (annotate_operons_with_genes), which lists any gene the merge
# now spans across a clustering coverage hole (e.g. rpsQ in the rPtn supercluster) and counts
# shared genes once (the final `dedup_operon_gene_lists` is then a no-op safety net). Every
# candidate's decision is written to `operon_merge_decisions.tsv`.

# ## Step 1 — Load isoforms
# 
# Filter the ~80 K clustered isoforms to those with at least `MIN_READS` supporting reads. Lower thresholds retain more isoforms but risk including sequencing noise and spurious clusters; higher thresholds miss lowly-expressed operons. At `MIN_READS = 50` the dominant transcription units are well-supported while RNase-generated fragments (which typically accumulate fewer reads than full-length transcripts) are partially suppressed.

# In[1]:


from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

# Headless execution: no-op the notebook display() and force a non-interactive backend.
try:
    from IPython.display import display  # noqa: F401
except Exception:
    def display(*args, **kwargs):
        return None
import matplotlib
matplotlib.use("Agg")

# ── Paths (organized server layout; run with Syn1_Operon/ as the working dir) ─
# OUT_FOLDER holds ONLY the canonical operon table (operons.candidate_blocks.tsv);
# every other artifact (intermediate TSVs, decision tables, QC/decision PDFs, BED)
# is written to SEG_FOLDER (./segmentation) to keep the main folder clean.
MOTHER_FOLDER = ".."
ISOFORMS_TSV  = MOTHER_FOLDER + "/Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
OUT_FOLDER    = "."
SEG_FOLDER    = OUT_FOLDER + "/segmentation"
Path(OUT_FOLDER).mkdir(parents=True, exist_ok=True)
Path(SEG_FOLDER).mkdir(parents=True, exist_ok=True)

# ── Single parameter ───────────────────────────────────────────────────────
MIN_READS = 50   # minimum read support per isoform — only tunable parameter

# ── Load ───────────────────────────────────────────────────────────────────
df = pd.read_csv(ISOFORMS_TSV, sep="\t")
df_iso = df[df["n_reads"] >= MIN_READS].copy()

print(f"Total isoforms loaded: {len(df)}")
print(f"Isoforms with n_reads >= {MIN_READS}: {len(df_iso)}")
print(f"\nStrand breakdown:")
print(df_iso["strand"].value_counts().to_string())


# ## Step 2 — Containment clustering → initial operons
# 
# **Biological basis:** Each operon produces one primary full-length transcript (longest isoform) and many degradation/truncation products. All products of the same operon are nested subsets of the full-length transcript. Two isoforms from *different* operons may overlap at their boundaries but one will not be fully contained within the other.
# 
# **Clustering rule:** Isoforms i and j join the same cluster iff one is fully contained within the other (within `BOUNDARY_TOL` bp tolerance to absorb soft-clipping). Partially overlapping isoforms — characteristic of adjacent operons whose 3′ and 5′ ends interleave — are **never** merged here.
# 
# **Mathematical guarantee:** The cluster boundary always equals the longest isoform. If any member extended beyond the longest isoform's boundaries, that member would itself be longer — a contradiction.
# 
# **Parameters:**
# - `MIN_READS` — isoform read-support threshold (Step 1)
# - `BOUNDARY_TOL` — boundary tolerance in bp; absorbs minor coordinate variation from soft-clipping at read ends
# 
# **Implementation:** Union-Find with path compression on containment pairs; sweep-line reduces candidate comparisons from O(n²) to O(n log n).

# In[2]:


BOUNDARY_TOL = 10   # bp — isoforms within this many bp of containment are treated as contained

def cluster_isoforms(isoforms: pd.DataFrame,
                     tol: int = BOUNDARY_TOL) -> list[dict]:
    """
    Cluster isoforms into operons by containment (with tolerance).

    Two isoforms i and j are in the same cluster iff one is within `tol` bp
    of being fully contained within the other:
        start_i >= start_j - tol  AND  end_i <= end_j + tol   (i inside j)
        OR vice-versa.

    The cluster boundary = the outermost (longest) isoform.
    Uses union-find + sweep-line for efficiency.
    """
    if isoforms.empty:
        return []

    iso    = isoforms.reset_index(drop=True)
    starts = iso["start0"].astype(int).tolist()
    ends   = iso["end0"].astype(int).tolist()
    n      = len(iso)

    # ── Union-Find with path compression ───────────────────────────────────
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)

    # ── Sweep: for each isoform i, check candidates j that start near i ────
    order = sorted(range(n), key=lambda i: starts[i])
    for idx, i in enumerate(order):
        si, ei = starts[i], ends[i]
        for j in order[idx + 1:]:
            sj = starts[j]
            if sj >= ei + tol:   # j starts well after i ends — no containment possible
                break
            ej = ends[j]
            # near-containment: i inside j (with tol), or j inside i (with tol)
            i_in_j = (si >= sj - tol) and (ei <= ej + tol)
            j_in_i = (sj >= si - tol) and (ej <= ei + tol)
            if i_in_j or j_in_i:
                union(i, j)

    # ── Collect connected components → one block each ──────────────────────
    components = defaultdict(list)
    for i in range(n):
        components[find(i)].append(i)

    chrom  = iso.loc[0, "chrom"]
    strand = iso.loc[0, "strand"]
    blocks = []
    for indices in components.values():
        members = [iso.loc[k, "isoform_id"] for k in indices]
        reads   = [int(iso.loc[k, "n_reads"])  for k in indices]
        s0 = min(starts[k] for k in indices)
        e0 = max(ends[k]   for k in indices)
        blocks.append({
            "chrom": chrom, "strand": strand,
            "start0": s0, "end0": e0,
            "members": members, "reads": reads,
        })
    return sorted(blocks, key=lambda b: b["start0"])


def blocks_to_df(blocks: list[dict]) -> pd.DataFrame:
    rows = [{
        "chrom":         b["chrom"],
        "strand":        b["strand"],
        "start0":        b["start0"],
        "end0":          b["end0"],
        "n_isoforms":    len(b["members"]),
        "n_reads_total": int(np.sum(b["reads"])),
        "member_ids":    ",".join(b["members"]),
    } for b in blocks]
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["chrom","strand","start0","end0","n_isoforms","n_reads_total","member_ids"])


# ── Cluster per strand ──────────────────────────────────────────────────────
operons_plus  = cluster_isoforms(df_iso[df_iso["strand"] == "+"])
operons_minus = cluster_isoforms(df_iso[df_iso["strand"] == "-"])

GENOME_LEN = 1_079_433   # CP002027.1 Syn1 genome length (bp)
def genome_coverage_pct(blocks, genome_len=GENOME_LEN):
    covered = sum(b["end0"] - b["start0"] for b in blocks)
    return covered, covered / genome_len * 100

cov_plus_bp,  cov_plus_pct  = genome_coverage_pct(operons_plus)
cov_minus_bp, cov_minus_pct = genome_coverage_pct(operons_minus)
print(f"BOUNDARY_TOL = {BOUNDARY_TOL} bp")
print(f"Operons — plus: {len(operons_plus)}, minus: {len(operons_minus)}, total: {len(operons_plus)+len(operons_minus)}")
print(f"Genome coverage — plus:  {cov_plus_bp:>8,} bp  ({cov_plus_pct:.1f}%)")
print(f"                  minus: {cov_minus_bp:>8,} bp  ({cov_minus_pct:.1f}%)")

# ── Build final DataFrame ───────────────────────────────────────────────────
all_blocks  = sorted(operons_plus + operons_minus, key=lambda b: (b["chrom"], b["start0"]))
operons_all = blocks_to_df(all_blocks)
operons_all.insert(0, "operon_id", [f"OP_{i+1:05d}" for i in range(len(operons_all))])
operons_all["length"] = operons_all["end0"] - operons_all["start0"]

# TSS/TTS: directional 5\'→3\' coordinates
# + strand: TSS = start0 (low), TTS = end0 (high)
# − strand: TSS = end0  (high), TTS = start0 (low)
operons_all["tss"] = np.where(operons_all["strand"] == "+",
                               operons_all["start0"], operons_all["end0"])
operons_all["tts"] = np.where(operons_all["strand"] == "+",
                               operons_all["end0"],  operons_all["start0"])

print(f"\nLength distribution (bp):")
print(operons_all["length"].describe().astype(int).to_string())

display(operons_all.drop(columns="member_ids").head(20))


# ## Step 3 — Gene annotation
# 
# Map the 911 annotated Syn1 genes (from GFF3) onto the isoform-derived operons. Each operon is annotated with:
# - **Sense genes** — annotated genes on the *same strand* that overlap the operon coordinates; these are the genes the operon is transcribing.
# - **Antisense genes** — annotated genes on the *opposite strand* that overlap the same genomic window; these reflect the dense gene packing of the minimal genome, where antisense transcription or convergent/divergent gene pairs are common.
# 
# Running gene annotation here — before the overlap and coverage analyses — makes gene context immediately available for diagnosing same-strand operon boundary conflicts in Step 4.

# In[3]:


import re

GFF3_FILE = MOTHER_FOLDER + "/Genomes_Input/syn1.genes.gff3"

# ── Parse GFF3 (1-based closed → 0-based half-open) ────────────────────────
gene_records = []
with open(GFF3_FILE) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        chrom  = parts[0]
        start1 = int(parts[3])
        end1   = int(parts[4])
        strand = parts[6]
        attrs  = parts[8]
        locus  = re.search(r"locus_tag=([^;]+)", attrs)
        name   = re.search(r"Name=([^;]+)",      attrs)
        prod   = re.search(r"product=([^;]+)",   attrs)
        gene_records.append({
            "chrom":     chrom,
            "strand":    strand,
            "start0":    start1 - 1,
            "end0":      end1,
            "locus_tag": locus.group(1) if locus else "",
            "gene_name": name.group(1)  if name  else "",
            "product":   prod.group(1)  if prod  else "",
        })

genes_df = pd.DataFrame(gene_records)
print(f"Genes loaded: {len(genes_df)}")
print(genes_df["strand"].value_counts().to_string())


# ── Annotate each operon with sense AND antisense overlapping genes ─────────
def annotate_operons_with_genes(operons: pd.DataFrame,
                                genes: pd.DataFrame) -> pd.DataFrame:
    """
    For each operon add:
      sense_gene_count / sense_gene_loci / sense_gene_names
          — genes on the SAME strand that overlap the operon
      antisense_gene_count / antisense_gene_loci / antisense_gene_names
          — genes on the OPPOSITE strand that overlap the operon
    """
    ops = operons.copy()
    opp = {"+": "-", "-": "+"}
    s_counts, s_loci, s_names = [], [], []
    a_counts, a_loci, a_names = [], [], []

    for _, op in ops.iterrows():
        overlaps  = (genes["start0"] < op["end0"]) & (genes["end0"] > op["start0"])
        sense     = genes[overlaps & (genes["strand"] == op["strand"])]
        antisense = genes[overlaps & (genes["strand"] == opp[op["strand"]])]
        s_counts.append(len(sense));  s_loci.append(",".join(sense["locus_tag"]));  s_names.append(",".join(sense["gene_name"]))
        a_counts.append(len(antisense)); a_loci.append(",".join(antisense["locus_tag"])); a_names.append(",".join(antisense["gene_name"]))

    ops["sense_gene_count"]     = s_counts
    ops["sense_gene_loci"]      = s_loci
    ops["sense_gene_names"]     = s_names
    ops["antisense_gene_count"] = a_counts
    ops["antisense_gene_loci"]  = a_loci
    ops["antisense_gene_names"] = a_names
    ops["gene_count"] = ops["sense_gene_count"] + ops["antisense_gene_count"]
    ops["gene_loci"]  = ops["sense_gene_loci"]
    ops["gene_names"] = ops["sense_gene_names"]
    return ops

operons_annotated = annotate_operons_with_genes(operons_all, genes_df)

print(f"\nOperons with sense genes:")
print(f"  0 sense genes:   {(operons_annotated['sense_gene_count'] == 0).sum()}")
print(f"  1 sense gene:    {(operons_annotated['sense_gene_count'] == 1).sum()}")
print(f"  2 sense genes:   {(operons_annotated['sense_gene_count'] == 2).sum()}")
print(f"  3+ sense genes:  {(operons_annotated['sense_gene_count'] >= 3).sum()}")
print(f"\nOperons with >=1 antisense gene: {(operons_annotated['antisense_gene_count'] > 0).sum()}")


# In[4]:


operons_annotated.head(20)


# In[5]:


operons_annotated.sort_values("gene_count", ascending=False).head(20)


# ## Step 4 — Operon overlap analysis
# 
# **Same-strand overlaps** arise when isoforms partially overlap without containment — typically RNase cleavage producing a 5′ fragment whose 3′ end falls inside the next operon. **Gene-in-gap** pairs arise when the clustering leaves a coverage hole over a gene (no isoform passes `MIN_READS`), so a same-strand gene sits in the gap between two consecutive operons (e.g. rpsQ inside the rPtn supercluster). Both are candidates resolved in Step 5a by a **co-transcription test** rather than by TSS location:
#
# - *Merge* if PacBio reads run continuously across the junction — `>= MERGE_MIN_BRIDGE` strand-specific bridging reads (spanning `[j-W, j+W]`) **and** gap/flank depth continuity `>= MERGE_MIN_CONT`.
# - *Separate* otherwise — a genuine terminator/promoter boundary, marked by few bridging reads and/or a depth valley at the junction.
#
# **Cross-strand (antisense) overlaps** are biologically expected in the compact Syn1 genome, where convergent and divergent gene pairs and naturally antisense transcription are common. These are reported for information but not resolved here.

# In[26]:


def find_overlapping_pairs(df_a: pd.DataFrame, df_b: pd.DataFrame,
                           same_strand: bool) -> pd.DataFrame:
    """
    Sweep-line overlap detection between two operon DataFrames.
    Returns one row per overlapping pair with overlap_bp and overlap fractions.
    If same_strand=True, df_a and df_b should be the same strand subset.
    """
    a = df_a.sort_values("start0").reset_index(drop=True)
    b = df_b.sort_values("start0").reset_index(drop=True)

    rows = []
    bi = 0
    for _, ra in a.iterrows():
        # advance b pointer past intervals that end before ra starts
        while bi < len(b) and b.loc[bi, "end0"] <= ra["start0"]:
            bi += 1
        j = bi
        while j < len(b) and b.loc[j, "start0"] < ra["end0"]:
            rb = b.loc[j]
            # skip self-comparison
            if same_strand and ra["operon_id"] == rb["operon_id"]:
                j += 1
                continue
            ovl = min(ra["end0"], rb["end0"]) - max(ra["start0"], rb["start0"])
            if ovl > 0:
                len_a = ra["end0"] - ra["start0"]
                len_b = rb["end0"] - rb["start0"]
                rows.append({
                    "operon_a":    ra["operon_id"],
                    "strand_a":    ra["strand"],
                    "start_a":     ra["start0"],
                    "end_a":       ra["end0"],
                    "operon_b":    rb["operon_id"],
                    "strand_b":    rb["strand"],
                    "start_b":     rb["start0"],
                    "end_b":       rb["end0"],
                    "overlap_bp":  ovl,
                    "frac_of_a":   round(ovl / len_a, 3),
                    "frac_of_b":   round(ovl / len_b, 3),
                    "relationship": (
                        "a_contains_b" if ra["start0"] <= rb["start0"] and ra["end0"] >= rb["end0"]
                        else "b_contains_a" if rb["start0"] <= ra["start0"] and rb["end0"] >= ra["end0"]
                        else "partial"
                    ),
                })
            j += 1

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "operon_a","strand_a","start_a","end_a",
        "operon_b","strand_b","start_b","end_b",
        "overlap_bp","frac_of_a","frac_of_b","relationship"])


# ══════════════════════════════════════════════════════════════════════════════
# Co-transcription merge primitives (PacBio bridging reads + depth continuity).
# These replace the TSS-location merge rule: a candidate junction is merged only
# if RNA runs continuously across it (direct read evidence), which is robust to
# where the clustering happens to place operon boundaries.
# ══════════════════════════════════════════════════════════════════════════════
MERGE_MAX_GAP    = 600    # bp: max gap for a gene-in-gap (coverage-hole) merge candidate
                          # (>= the 567 bp rPtn rpsQ gap; 500 would split rPtn)
MERGE_W          = 80     # bp: half-window for bridging + gap depth
MERGE_FLANK      = 300    # bp: flank window for continuity
MERGE_MIN_BRIDGE = 50     # min strand-specific reads spanning a junction to merge
MERGE_MIN_CONT   = 0.50   # min gap/flank depth ratio to merge

_MERGE_BAM_PATH = MOTHER_FOLDER + "/Syn1_Transcriptomics/PacBio/PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam"
_MERGE_DEPTH = {s: MOTHER_FOLDER + f"/Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.{n}.bedGraph"
                for s, n in (("+", "plus"), ("-", "minus"))}
import pysam
_MERGE_CHROM = str(genes_df["chrom"].iloc[0])
_merge_bam = pysam.AlignmentFile(_MERGE_BAM_PATH, "rb")
_merge_dep = {s: pd.read_csv(p, sep="\t", header=None, names=["c", "s", "e", "v"])
              for s, p in _MERGE_DEPTH.items()}


def _bridging_reads(j: int, strand: str, w: int = MERGE_W) -> int:
    """Strand-specific PacBio FLNC reads whose alignment spans [j-w, j+w]."""
    want_rev = (strand == "-")
    n = 0
    for r in _merge_bam.fetch(_MERGE_CHROM, max(0, j - w), j + w):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_reverse != want_rev:
            continue
        if r.reference_start <= j - w and r.reference_end >= j + w:
            n += 1
    return n


def _mean_depth(strand: str, s: int, e: int) -> float:
    d = _merge_dep[strand]
    sub = d[(d.c == _MERGE_CHROM) & (d.e > s) & (d.s < e)]
    if sub.empty or e <= s:
        return 0.0
    ov = np.minimum(sub.e.values, e) - np.maximum(sub.s.values, s)
    return float((sub.v.values * ov).sum()) / (e - s)


def junction_coexpressed(j: int, strand: str):
    """(merge_bool, n_bridge, continuity) for a junction at position j."""
    nbr = _bridging_reads(j, strand)
    gapd = _mean_depth(strand, j - MERGE_W, j + MERGE_W)
    flankd = 0.5 * (_mean_depth(strand, j - MERGE_FLANK, j - MERGE_W) +
                    _mean_depth(strand, j + MERGE_W, j + MERGE_FLANK))
    cont = gapd / flankd if flankd > 0 else 0.0
    return (nbr >= MERGE_MIN_BRIDGE and cont >= MERGE_MIN_CONT), nbr, cont


def find_merge_candidates(df_strand: pd.DataFrame, genes: pd.DataFrame,
                          max_gap: int = MERGE_MAX_GAP) -> pd.DataFrame:
    """Consecutive same-strand operon pairs eligible for a co-transcription merge:
      - OVERLAPPING (gap < 0, e.g. RNase 5'/3' fragments), or
      - separated by a gap (<= max_gap) that CONTAINS a same-strand gene -- a
        clustering coverage hole (e.g. rpsQ inside the rPtn supercluster).
    Clean intergenic gaps (no gene in the gap) are NOT candidates, so genuine
    operon boundaries are never merged. Schema mirrors find_overlapping_pairs
    plus gap_bp + junction_pos."""
    d = df_strand.sort_values("start0").reset_index(drop=True)
    cols = ["operon_a", "strand_a", "start_a", "end_a", "operon_b", "strand_b",
            "start_b", "end_b", "overlap_bp", "gap_bp", "frac_of_a", "frac_of_b",
            "relationship", "junction_pos"]
    if len(d) < 2:
        return pd.DataFrame(columns=cols)
    strand = d.iloc[0]["strand"]
    g = genes[genes["strand"] == strand]
    rows = []
    for i in range(len(d) - 1):
        ra, rb = d.loc[i], d.loc[i + 1]
        gap = int(rb["start0"]) - int(ra["end0"])         # < 0 => overlap
        if gap < 0:
            rel = "overlap"
        elif gap <= max_gap and bool(((g["end0"] > ra["end0"]) & (g["start0"] < rb["start0"])).any()):
            rel = "gene_in_gap"
        else:
            continue
        ovl = max(0, -gap)
        la, lb = int(ra["end0"]) - int(ra["start0"]), int(rb["end0"]) - int(rb["start0"])
        rows.append({
            "operon_a": ra["operon_id"], "strand_a": ra["strand"],
            "start_a": int(ra["start0"]), "end_a": int(ra["end0"]),
            "operon_b": rb["operon_id"], "strand_b": rb["strand"],
            "start_b": int(rb["start0"]), "end_b": int(rb["end0"]),
            "overlap_bp": ovl, "gap_bp": max(0, gap),
            "frac_of_a": round(ovl / la, 3) if la else 0.0,
            "frac_of_b": round(ovl / lb, 3) if lb else 0.0,
            "relationship": rel,
            "junction_pos": (int(ra["end0"]) + int(rb["start0"])) // 2,
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=cols)


# In[27]:


# ── 1. Same-strand overlap (should be zero by construction) ─────────────────
# operons_annotated carries sense_gene_loci so we can show gene context for
# each conflicting pair — useful for deciding which to merge in the RNase step.
plus_ops  = operons_annotated[operons_annotated["strand"] == "+"].copy()
minus_ops = operons_annotated[operons_annotated["strand"] == "-"].copy()

ss_plus  = find_overlapping_pairs(plus_ops,  plus_ops,  same_strand=True)
ss_minus = find_overlapping_pairs(minus_ops, minus_ops, same_strand=True)
ss_plus  = ss_plus[ss_plus["operon_a"]  < ss_plus["operon_b"]]
ss_minus = ss_minus[ss_minus["operon_a"] < ss_minus["operon_b"]]
same_strand_overlaps = pd.concat([ss_plus, ss_minus], ignore_index=True)

# Annotate each pair with sense gene loci from operons_annotated
op_gene_map = operons_annotated.set_index("operon_id")["sense_gene_loci"].to_dict()
if len(same_strand_overlaps):
    same_strand_overlaps["genes_a"] = same_strand_overlaps["operon_a"].map(op_gene_map)
    same_strand_overlaps["genes_b"] = same_strand_overlaps["operon_b"].map(op_gene_map)

print(f"Same-strand overlapping operon pairs: {len(same_strand_overlaps)}")
if len(same_strand_overlaps) > 0:
    print("  -> Caused by isoforms that partially overlap without containment.")
    print("  -> These pairs are candidates for merging in the RNase-rescue step.")
    # display(same_strand_overlaps[["operon_a","operon_b","strand_a",
    #                                "start_a","end_a","start_b","end_b",
    #                                "overlap_bp","frac_of_a","frac_of_b",
    #                                "relationship","genes_a","genes_b"]])


# In[28]:


# ── Pairs where genes_a and genes_b share at least one gene locus ────────────
# A shared gene means both operons claim the same gene on the same strand —
# a direct conflict that should be resolved by merging the two operons.

if len(same_strand_overlaps):
    def shared_loci(row):
        a = set(str(row["genes_a"]).split(",")) - {"", "nan"}
        b = set(str(row["genes_b"]).split(",")) - {"", "nan"}
        shared = a & b
        return ",".join(sorted(shared)) if shared else ""

    same_strand_overlaps["shared_genes"] = same_strand_overlaps.apply(shared_loci, axis=1)
    conflicts = same_strand_overlaps[same_strand_overlaps["shared_genes"] != ""].copy()

    print(f"Pairs with shared (conflicting) genes: {len(conflicts)} / {len(same_strand_overlaps)} total same-strand overlaps")

    if len(conflicts):
        display(conflicts[["operon_a","operon_b","strand_a",
                            "start_a","end_a","start_b","end_b",
                            "overlap_bp","frac_of_a","frac_of_b",
                            "relationship","genes_a","genes_b","shared_genes"]]
                .reset_index(drop=True))

        # Save to TSV for inspection
        OUT_SS = SEG_FOLDER + "/same_strand_gene_conflicts.tsv"
        conflicts.to_csv(OUT_SS, sep="\t", index=False)
        print(f"\nSaved: {OUT_SS}")
    else:
        print("No shared genes across conflicting pairs — overlaps are between adjacent distinct operons.")
        display(same_strand_overlaps[["operon_a","operon_b","strand_a",
                                       "start_a","end_a","start_b","end_b",
                                       "overlap_bp","relationship",
                                       "genes_a","genes_b"]].reset_index(drop=True))
else:
    print("No same-strand overlaps to analyse.")


# ## Step 5a — Merge possible RNase-Processed Operons

# In[29]:


# ══════════════════════════════════════════════════════════════════════════════
# Step 5a merge rule — CO-TRANSCRIPTION (PacBio bridging reads + depth continuity).
# (Replaces the former TSS-location rule. Parameters MERGE_* are defined above with
#  the merge primitives _bridging_reads / _mean_depth / junction_coexpressed.)
#
# Candidate pairs  (find_merge_candidates, per strand, consecutive operons):
#   - OVERLAP     : the two operons overlap (gap < 0; RNase 5'/3' fragments), or
#   - GENE_IN_GAP : gap <= MERGE_MAX_GAP that CONTAINS a same-strand gene (a
#                   clustering coverage hole, e.g. rpsQ inside the rPtn operon).
#   A clean intergenic gap (no gene in it) is NOT a candidate -> real operon
#   boundaries are never merged.
#
# Decision  (junction_coexpressed at the junction midpoint j):
#   bridging   = strand-specific PacBio FLNC primary reads spanning [j-W, j+W]
#                (reference_start <= j-W AND reference_end >= j+W)
#   continuity = mean sense-strand depth over [j-W, j+W]
#                / mean of (depth[j-FLANK, j-W], depth[j+W, j+FLANK])   (unbounded;
#                ~1 = uniform coverage, ~0 = a depth valley at a boundary)
#   MERGE iff  bridging >= MERGE_MIN_BRIDGE  AND  continuity >= MERGE_MIN_CONT.
#
# Apply: passing pairs are unioned with NO transitive chaining — strongest-bridging
#   first, each operon joins at most one merge, so every merged operon is a 2-operon
#   pair. Merged operons are re-annotated from their spans (so a gene the merge spans
#   across a coverage hole is listed, shared genes counted once). All candidate
#   decisions -> operon_merge_decisions.tsv.
# ══════════════════════════════════════════════════════════════════════════════
UTR_MAX_NT = 10   # (legacy, unused) nt upstream of a gene start; get_downstream_tss /
                  # tss_is_separate below are kept for reference only, not called.

def get_downstream_tss(row):
    """
    Return (downstream_operon_id, downstream_tss, upstream_operon_id).
    Downstream = further along the direction of transcription:
      + strand → larger start0 coordinate
      − strand → smaller end0 coordinate
    """
    strand = row["strand_a"]
    if strand == "+":
        if row["start_a"] >= row["start_b"]:
            return row["operon_a"], int(row["start_a"]), row["operon_b"]
        else:
            return row["operon_b"], int(row["start_b"]), row["operon_a"]
    else:
        if row["end_a"] <= row["end_b"]:
            return row["operon_a"], int(row["end_a"]), row["operon_b"]
        else:
            return row["operon_b"], int(row["end_b"]), row["operon_a"]


def tss_is_separate(tss: int, strand: str,
                    genes_on_strand: pd.DataFrame,
                    utr_max: int = UTR_MAX_NT):
    """
    Return (separate: bool, reason_str).
    Separate (genuine promoter) if EITHER:
      a) TSS is intergenic — not inside any same-strand gene body, OR
      b) TSS is within utr_max bp upstream of the nearest downstream gene start.
    Merge (RNase artifact) only when BOTH fail.
    """
    # Condition a: intergenic
    inside = genes_on_strand[
        (genes_on_strand["start0"] <= tss) &
        (genes_on_strand["end0"]   >  tss)
    ]
    if inside.empty:
        return True, f"intergenic (TSS {tss} not inside any gene)"

    # Condition b: within utr_max nt upstream of a downstream gene start
    # For + strand: downstream gene start > tss, distance = start - tss
    # For − strand: downstream gene end  < tss, distance = tss - end
    if strand == "+":
        nearby = genes_on_strand[
            (genes_on_strand["start0"] > tss) &
            (genes_on_strand["start0"] - tss <= utr_max)
        ]
    else:
        nearby = genes_on_strand[
            (genes_on_strand["end0"] < tss) &
            (tss - genes_on_strand["end0"] <= utr_max)
        ]
    if not nearby.empty:
        locus = nearby.iloc[0]["locus_tag"]
        dist  = (int(nearby.iloc[0]["start0"]) - tss
                 if strand == "+" else
                 tss - int(nearby.iloc[0]["end0"]))
        return True, (f"UTR window: TSS {tss} is {dist} nt upstream of "
                      f"{locus} (within {utr_max} nt)")

    loci = ",".join(inside["locus_tag"].tolist())
    return False, f"TSS {tss} inside gene(s) {loci} and >{utr_max} nt from any gene start"


# ── Evaluate each merge candidate by co-transcription (bridging + continuity) ─
# Candidate pairs: overlapping operons OR gapped operons with a same-strand gene
# in the gap (clustering coverage hole). Decision is read-evidence, not TSS
# location. (get_downstream_tss / tss_is_separate above are retained for
# reference but no longer drive the merge.)
op_gene_map = operons_annotated.set_index("operon_id")["sense_gene_loci"].to_dict()
merge_candidates = pd.concat([
    find_merge_candidates(operons_annotated[operons_annotated["strand"] == "+"], genes_df),
    find_merge_candidates(operons_annotated[operons_annotated["strand"] == "-"], genes_df),
], ignore_index=True)

merge_decisions = []
for _, row in merge_candidates.iterrows():
    strand = row["strand_a"]
    j = int(row["junction_pos"])
    do_merge, nbr, cont = junction_coexpressed(j, strand)
    ga = set(str(op_gene_map.get(row["operon_a"], "")).split(",")) - {"", "nan"}
    gb = set(str(op_gene_map.get(row["operon_b"], "")).split(",")) - {"", "nan"}
    merge_decisions.append({
        "operon_a":       row["operon_a"],
        "operon_b":       row["operon_b"],
        "strand":         strand,
        "upstream_op":    row["operon_a"],
        "downstream_op":  row["operon_b"],
        "downstream_tss": j,
        "relationship":   row["relationship"],
        "overlap_bp":     row["overlap_bp"],
        "gap_bp":         row["gap_bp"],
        "genes_a":        ",".join(sorted(ga)),
        "genes_b":        ",".join(sorted(gb)),
        "shared_genes":   ",".join(sorted(ga & gb)),
        "bridging_reads": nbr,
        "continuity":     round(cont, 3),
        "merge":          do_merge,
        "reason":         (f"co-transcribed (bridge={nbr}, continuity={cont:.2f})" if do_merge
                           else f"separate (bridge={nbr}, continuity={cont:.2f})"),
    })

decisions_df = pd.DataFrame(merge_decisions)
n_merge    = int(decisions_df["merge"].sum()) if len(decisions_df) else 0
n_separate = len(decisions_df) - n_merge

print(f"Co-transcription merge candidates (overlap or gene-in-gap): {len(decisions_df)}")
print(f"  Merge    (bridging >= {MERGE_MIN_BRIDGE} AND continuity >= {MERGE_MIN_CONT}): {n_merge}")
print(f"  Separate: {n_separate}")
if len(decisions_df):
    decisions_df.to_csv(SEG_FOLDER + "/operon_merge_decisions.tsv", sep="\t", index=False)
    display(decisions_df[["operon_a","operon_b","strand","relationship","downstream_tss",
                          "overlap_bp","gap_bp","bridging_reads","continuity",
                          "genes_a","genes_b","merge","reason"]])

# ── Apply merges to operons_annotated ────────────────────────────────────────
to_merge = decisions_df[decisions_df["merge"]]

if len(to_merge) == 0:
    print("\nNo merges to apply.")
    operons_merged = operons_annotated.copy()
else:
    from collections import defaultdict as _dd

    op_ids = list(operons_annotated["operon_id"])
    idx_of = {oid: i for i, oid in enumerate(op_ids)}
    parent = list(range(len(op_ids)))

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def _union(x, y):
        parent[_find(x)] = _find(y)

    # No transitive chaining: each operon may join at most ONE merge, so every
    # merged cluster is exactly a pair (2 operons). Process strongest first.
    skipped = []
    used = set()
    for _, row in to_merge.sort_values("bridging_reads", ascending=False).iterrows():
        a, b = row["operon_a"], row["operon_b"]
        a_ok, b_ok = a in idx_of, b in idx_of
        if not (a_ok and b_ok):
            skipped.append((a, b, "operon_a missing" if not a_ok else "operon_b missing"))
            continue
        if a in used or b in used:          # would chain -> skip
            continue
        _union(idx_of[a], idx_of[b])
        used.add(a); used.add(b)
    if skipped:
        print(f"WARNING: {len(skipped)} merge pair(s) skipped — operon ID not found in operons_annotated:")
        for a, b, reason in skipped:
            print(f"  {a} + {b}: {reason}")
        print("  These operons likely have IDs from a previous run. Re-run from Step 2.")

    # Build clusters and report sizes for full transparency
    clusters = _dd(list)
    for i, oid in enumerate(op_ids):
        clusters[_find(i)].append(oid)
    cluster_sizes = [len(v) for v in clusters.values()]
    n_merged_clusters  = sum(1 for s in cluster_sizes if s >= 2)
    n_singleton_ops    = sum(1 for s in cluster_sizes if s == 1)
    pairs_accounted    = sum(s - 1 for s in cluster_sizes if s >= 2)
    print(f"Merge cluster breakdown:")
    from collections import Counter as _Ctr
    for size, count in sorted(_Ctr(cluster_sizes).items()):
        if size == 1:
            print(f"  size 1 (unchanged):  {count} operons")
        else:
            print(f"  size {size} cluster:       {count} × ({size} operons → 1 merged) "
                  f"= {count*(size-1)} pairs consumed")
    print(f"  Total pairs consumed by Union-Find: {pairs_accounted}")
    print(f"  Merge pairs requested: {len(to_merge) - len(skipped)}")
    if pairs_accounted != len(to_merge) - len(skipped):
        print(f"  WARNING: mismatch — {len(to_merge)-len(skipped)} pairs requested "
              f"but {pairs_accounted} accounted for in clusters")

    ops = operons_annotated.set_index("operon_id")
    merged_rows = []
    for root, members in clusters.items():
        if len(members) == 1:
            merged_rows.append(ops.loc[members[0]].to_dict())
        else:
            group  = ops.loc[members]
            s0     = int(group["start0"].min())
            e0     = int(group["end0"].max())
            strand = ops.loc[members[0]]["strand"]
            merged_rows.append({
                **ops.loc[members[0]].to_dict(),
                "start0":               s0,
                "end0":                 e0,
                "length":               e0 - s0,
                "tss":                  s0 if strand == "+" else e0,
                "tts":                  e0 if strand == "+" else s0,
                "n_isoforms":           int(group["n_isoforms"].sum()),
                "n_reads_total":        int(group["n_reads_total"].sum()),
                "sense_gene_count":     int(group["sense_gene_count"].sum()),
                "sense_gene_loci":      ",".join(filter(None, group["sense_gene_loci"])),
                "sense_gene_names":     ",".join(filter(None, group["sense_gene_names"])),
                "antisense_gene_count": int(group["antisense_gene_count"].sum()),
                "antisense_gene_loci":  ",".join(filter(None, group["antisense_gene_loci"])),
                "antisense_gene_names": ",".join(filter(None, group["antisense_gene_names"])),
                "member_ids":           ",".join(filter(None, group["member_ids"])),
                "segmentation_type":          "isoform_operon_merged",
            })

    operons_merged = (pd.DataFrame(merged_rows)
                      .sort_values(["chrom", "start0"])
                      .reset_index(drop=True))

    if "operon_id" in operons_merged.columns:
        operons_merged = operons_merged.drop(columns=["operon_id"])
    operons_merged.insert(0, "operon_id",
        [f"OP_{i+1:05d}" for i in range(len(operons_merged))])

    # Re-annotate gene lists from the (possibly extended) spans so genes that a
    # merge now spans across a clustering gap (e.g. rpsQ in the rPtn supercluster)
    # are listed, and so shared genes are counted once (removes the merge
    # over-count). Span-based, idempotent for unchanged operons.
    operons_merged = annotate_operons_with_genes(operons_merged, genes_df)

    n_before = len(operons_annotated)
    n_after  = len(operons_merged)
    print(f"\nOperons before merge: {n_before}")
    print(f"Operons after merge:  {n_after}  (reduced by {n_before - n_after})")
    print(f"  = {n_singleton_ops} unchanged + {n_merged_clusters} merged clusters")
    print(f"  (isoform_operon_merged = number of merged OUTPUT operons, "
          f"not number of pairs consumed)")

    pm = operons_merged[operons_merged["strand"] == "+"]
    mm = operons_merged[operons_merged["strand"] == "-"]
    ss_check = pd.concat([
        find_overlapping_pairs(pm, pm, same_strand=True),
        find_overlapping_pairs(mm, mm, same_strand=True),
    ])
    ss_check = ss_check[ss_check["operon_a"] < ss_check["operon_b"]]
    print(f"Remaining same-strand overlaps after merge: {len(ss_check)}")


# In[30]:


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def plot_operon_conflict(row, genes_df, merged: bool = False, ax=None):
    """
    Draw a single same-strand operon conflict pair.

    Layout:
      Top track    — Operon A bar
      Middle track — gene bars
      Bottom track — Operon B bar
      Gold shading — overlap region

    merged=True  → frame and title in green  ("MERGED")
    merged=False → frame and title in orange ("KEPT SEPARATE")
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 2.2))

    strand   = row["strand_a"]
    s_a, e_a = int(row["start_a"]), int(row["end_a"])
    s_b, e_b = int(row["start_b"]), int(row["end_b"])
    ovl_s    = max(s_a, s_b)
    ovl_e    = min(e_a, e_b)

    shared   = set(str(row.get("shared_genes","")).split(",")) - {"","nan"}
    loci_a   = set(str(row.get("genes_a","")).split(","))      - {"","nan"}
    loci_b   = set(str(row.get("genes_b","")).split(","))      - {"","nan"}
    all_loci = loci_a | loci_b

    region_genes = genes_df[genes_df["locus_tag"].isin(all_loci)].copy()

    view_s = min(s_a, s_b)
    view_e = max(e_a, e_b)
    pad    = max((view_e - view_s) * 0.05, 50)
    view_s -= pad;  view_e += pad

    # Overlap shading
    ax.axvspan(ovl_s, ovl_e, color="gold", alpha=0.25, zorder=0)

    OPERON_Y = {"A": 1.0, "B": 0.0}
    OPERON_H = 0.25
    GENE_H   = 0.14
    COLORS   = {"A": "#4C72B0", "B": "#DD8452"}

    for label, (s, e) in [("A", (s_a, e_a)), ("B", (s_b, e_b))]:
        y = OPERON_Y[label]
        ax.broken_barh([(s, e - s)], (y - OPERON_H/2, OPERON_H),
                       facecolors=COLORS[label], alpha=0.85, zorder=2)
        mid = (s + e) / 2
        dx  = (e - s) * 0.08 * (1 if strand == "+" else -1)
        ax.annotate("", xy=(mid + dx, y), xytext=(mid, y),
                    arrowprops=dict(arrowstyle="->", color="white", lw=1.2), zorder=3)
        ax.text(s + (e - s) * 0.02, y + OPERON_H/2 + 0.04,
                f"Operon {label}: {row[f'operon_{label.lower()}']}  ({e-s:,} bp)",
                fontsize=7, va="bottom", color=COLORS[label], fontweight="bold")

    for _, gene in region_genes.iterrows():
        lt        = gene["locus_tag"]
        gs, ge    = int(gene["start0"]), int(gene["end0"])
        is_shared = lt in shared
        in_a, in_b = lt in loci_a, lt in loci_b
        color = "#C44E52" if is_shared else "#888888"
        ax.broken_barh([(gs, ge - gs)], (0.5 - GENE_H/2, GENE_H),
                       facecolors=color, alpha=0.75, zorder=2)
        gname = gene["gene_name"]
        claim = "A+B" if (in_a and in_b) else ("A" if in_a else "B")
        parts = [gname if gname else lt]
        if is_shared: parts.append("*")
        parts.append(f"[{claim}]")
        ax.text((gs + ge) / 2, 0.5 - GENE_H/2 - 0.06, " ".join(parts),
                fontsize=6, ha="center", va="top",
                color="#C44E52" if is_shared else "#444444")

    ax.set_xlim(view_s, view_e)
    ax.set_ylim(-0.45, 1.55)
    ax.set_xlabel("Genomic coordinate (bp)", fontsize=8)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_yticklabels(["Operon B", "Genes", "Operon A"], fontsize=7)
    ax.tick_params(axis="x", labelsize=7)

    # ── Merge / separate annotation ──────────────────────────────────────────
    decision_color = "#2ca02c" if merged else "#ff7f0e"   # green / orange
    decision_label = "MERGED" if merged else "KEPT SEPARATE"
    reason         = row.get("reason", "")

    shared_str = f"  shared: {', '.join(sorted(shared))}" if shared else ""
    ax.set_title(
        f"Strand {strand}  |  overlap {int(row['overlap_bp'])} bp  "
        f"({row['relationship']}){shared_str}\n"
        f"{reason}",
        fontsize=7, pad=3
    )

    # Coloured border to instantly signal decision
    for spine in ax.spines.values():
        spine.set_edgecolor(decision_color)
        spine.set_linewidth(2.0)

    # Decision badge (top-left corner)
    ax.text(0.01, 0.97, decision_label,
            transform=ax.transAxes, fontsize=8, fontweight="bold",
            color="white", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.2", facecolor=decision_color,
                      edgecolor="none", alpha=0.85))

    patches = [
        mpatches.Patch(color=COLORS["A"], label=f"Operon A ({row['operon_a']})"),
        mpatches.Patch(color=COLORS["B"], label=f"Operon B ({row['operon_b']})"),
        mpatches.Patch(color="#C44E52",   label="Shared gene"),
        mpatches.Patch(color="#888888",   label="Unique gene"),
        mpatches.Patch(color="gold",      alpha=0.5, label="Overlap region"),
    ]
    ax.legend(handles=patches, fontsize=6, loc="upper right",
              framealpha=0.7, ncol=2)
    return ax


# ── Plot all same-strand overlap pairs, annotated with merge decision ────────
if len(same_strand_overlaps) == 0:
    print("No same-strand overlaps to plot.")
else:
    # Join merge decision onto same_strand_overlaps via operon_a + operon_b key
    plot_df = same_strand_overlaps.merge(
        decisions_df[["operon_a","operon_b","merge","reason","downstream_tss",
                      "downstream_op","upstream_op"]],
        on=["operon_a","operon_b"], how="left"
    )
    plot_df["merge"] = plot_df["merge"].fillna(False)

    n   = len(plot_df)
    fig, axes = plt.subplots(n, 1, figsize=(12, 2.8 * n), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, (_, row) in zip(axes, plot_df.iterrows()):
        plot_operon_conflict(row, genes_df, merged=bool(row["merge"]), ax=ax)

    n_merged   = int(plot_df["merge"].sum())
    n_separate = n - n_merged
    fig.suptitle(
        f"Same-strand operon overlaps  "
        f"(MIN_READS={MIN_READS}, BOUNDARY_TOL={BOUNDARY_TOL})  |  "
        f"merged: {n_merged}  kept separate: {n_separate}",
        fontsize=10, y=1.01
    )

    OUT_FIG = SEG_FOLDER + "/same_strand_conflicts.pdf"
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    print(f"Saved: {OUT_FIG}")
    plt.show()


# ── Gene-in-gap merge decisions (the candidates same_strand_conflicts.pdf omits) ─
# same_strand_conflicts.pdf above shows only OVERLAP candidates. The other class of
# merge candidate — two same-strand operons separated by a gap that CONTAINS a
# same-strand gene (a clustering coverage hole, e.g. rpsQ inside the rPtn
# supercluster) — is visualized here. Each panel shows, around the junction:
# sense-strand PacBio FLNC depth (continuity), the same-strand genes (gap gene in
# purple), the junction (dashed), and the verdict (green=MERGE, red=SEPARATE) with
# bridging reads + continuity. Reuses the already-loaded depth (_merge_dep) + genes.
def _plot_gene_in_gap_decisions(decisions_df, genes, dep, chrom,
                                out_pdf, pad=1200):
    gig = (decisions_df[decisions_df["relationship"] == "gene_in_gap"]
           .sort_values("downstream_tss").reset_index(drop=True))
    print(f"\nGene-in-gap merge decisions: {len(gig)}  "
          f"(merge {int(gig['merge'].sum())}, "
          f"separate {int((~gig['merge'].astype(bool)).sum())})")
    if gig.empty:
        print("  Nothing to plot.")
        return

    def depth_array(dep_df, win_s, win_e):
        cov = np.zeros(win_e - win_s)
        sub = dep_df[(dep_df.c == chrom) & (dep_df.e > win_s) & (dep_df.s < win_e)]
        for _, r in sub.iterrows():
            a, b = max(int(r.s), win_s) - win_s, min(int(r.e), win_e) - win_s
            cov[a:b] = r.v
        return cov

    n    = len(gig)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 3.0, nrow * 1.7), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")

    for i, (_, r) in enumerate(gig.iterrows()):
        ax = axes[i // ncol][i % ncol]
        ax.axis("on")
        strand = r["strand"]
        j  = int(r["downstream_tss"])
        ws, we = j - pad, j + pad
        x   = np.arange(ws, we)
        cov = depth_array(dep[strand], ws, we)
        col = "#2ca02c" if bool(r["merge"]) else "#d62728"

        ax.fill_between(x, 0, cov, color="#9ecae1", lw=0)
        ax.plot(x, cov, color="#3182bd", lw=0.4)
        ymax = max(1.0, cov.max())
        ax.axvline(j, color=col, lw=1.0, ls="--")

        ab_loci = (set(str(r["genes_a"]).split(",")) |
                   set(str(r["genes_b"]).split(","))) - {"", "nan"}
        gw = genes[(genes.strand == strand) & (genes.end0 > ws) & (genes.start0 < we)]
        for _, g in gw.iterrows():
            gl, gr = max(g.start0, ws), min(g.end0, we)
            in_gap = g.locus_tag not in ab_loci
            gc = "#8856a7" if in_gap else "#999999"
            ax.add_patch(plt.Rectangle((gl, ymax * 1.02), gr - gl, ymax * 0.12,
                                       facecolor=gc, edgecolor="none", clip_on=False))
            if (gr - gl) > 0.15 * (we - ws):
                lbl = g.gene_name if g.gene_name else g.locus_tag
                ax.text((gl + gr) / 2, ymax * 1.20, lbl, ha="center", va="bottom",
                        fontsize=4.5, color=gc, clip_on=False)

        ax.set_xlim(ws, we); ax.set_ylim(0, ymax * 1.35)
        ax.set_yticks([0, int(ymax)])
        ax.set_xticks([ws, j, we])
        ax.set_xticklabels([f"{ws/1000:.1f}", "junc", f"{we/1000:.1f}"], fontsize=4.5)
        ax.tick_params(labelsize=4.5, length=2)
        ax.spines[["top", "right"]].set_visible(False)
        verdict = "MERGE" if bool(r["merge"]) else "SEPARATE"
        ax.set_title(f"{r['operon_a']}+{r['operon_b']} ({strand})  gap={int(r['gap_bp'])}\n"
                     f"bridge={int(r['bridging_reads'])}  cont={r['continuity']:.2f}  {verdict}",
                     fontsize=5, color=col)

    fig.suptitle("Gene-in-gap co-transcription merge decisions (sense-strand PacBio depth; "
                 "gap gene in purple; green=merge, red=separate)", fontsize=7, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    print(f"Saved: {out_pdf}  ({n} panels)")
    plt.close(fig)


if len(decisions_df):
    _plot_gene_in_gap_decisions(
        decisions_df, genes_df, _merge_dep, _MERGE_CHROM,
        SEG_FOLDER + "/operon_merge_decisions_gene_in_gap.pdf")


# ## Step 5b — Merge operons sharing the same sense gene
# 
# A separate failure mode occurs when RNase activity is so extensive that no isoform spans the full gene — the surviving isoform clusters land on either side of the cleavage site, each claiming the same gene. These appear as two same-strand operons with identical `sense_gene_loci`. Their coordinate union recovers the true operon span.

# In[31]:


# ── Find same-strand operon pairs with identical non-empty sense_gene_loci ──
from collections import defaultdict as _dd

def merge_duplicate_sense_loci(operons: pd.DataFrame) -> pd.DataFrame:
    """
    Merge same-strand operons that share identical sense_gene_loci into one
    operon labelled segmentation_type="isoform_gene_combined".
    Uses Union-Find so transitive chains (A=B, B=C → merge all three) are handled.
    """
    ops = operons.copy().reset_index(drop=True)
    n   = len(ops)

    parent = list(range(n))
    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def _union(x, y):
        parent[_find(x)] = _find(y)

    # Index by (strand, sense_gene_loci) — only non-empty loci
    loci_index = _dd(list)   # (strand, loci_str) → [row indices]
    for i, row in ops.iterrows():
        loci = str(row.get("sense_gene_loci", "") or "").strip()
        if loci:
            loci_index[(row["strand"], loci)].append(i)

    n_groups_merged = 0
    for (strand, loci), indices in loci_index.items():
        if len(indices) > 1:
            for idx in indices[1:]:
                _union(indices[0], idx)
            n_groups_merged += 1

    if n_groups_merged == 0:
        print("No duplicate sense_gene_loci found — nothing to merge.")
        return ops

    clusters = _dd(list)
    for i in range(n):
        clusters[_find(i)].append(i)

    merged_rows = []
    for root, members in clusters.items():
        if len(members) == 1:
            merged_rows.append(ops.iloc[members[0]].to_dict())
        else:
            group  = ops.iloc[members]
            s0     = int(group["start0"].min())
            e0     = int(group["end0"].max())
            strand = group.iloc[0]["strand"]
            base   = group.iloc[0].to_dict()
            base.update({
                "start0":               s0,
                "end0":                 e0,
                "length":               e0 - s0,
                "tss":                  s0 if strand == "+" else e0,
                "tts":                  e0 if strand == "+" else s0,
                "n_isoforms":           int(group["n_isoforms"].sum()),
                "n_reads_total":        int(group["n_reads_total"].sum()),
                "sense_gene_count":     int(group["sense_gene_count"].max()),
                "sense_gene_loci":      group.iloc[0]["sense_gene_loci"],
                "sense_gene_names":     group.iloc[0]["sense_gene_names"],
                "antisense_gene_count": int(group["antisense_gene_count"].max()),
                "antisense_gene_loci":  ",".join(filter(None,
                    group["antisense_gene_loci"].fillna("").tolist())),
                "antisense_gene_names": ",".join(filter(None,
                    group["antisense_gene_names"].fillna("").tolist())),
                "member_ids":           ",".join(filter(None,
                    group["member_ids"].fillna("").tolist())),
                "segmentation_type":    "isoform_gene_combined",
            })
            print(f"  Merged {len(members)} operons for gene(s) "
                  f"'{group.iloc[0]['sense_gene_loci']}' "
                  f"(strand {strand}, {s0}–{e0})")
            merged_rows.append(base)

    result = (pd.DataFrame(merged_rows)
              .sort_values(["chrom", "start0"])
              .reset_index(drop=True))
    if "operon_id" in result.columns:
        result = result.drop(columns=["operon_id"])
    result.insert(0, "operon_id",
        [f"OP_{i+1:05d}" for i in range(len(result))])
    return result


print("Merging operons with duplicate sense_gene_loci...")
operons_merged = merge_duplicate_sense_loci(operons_merged)

for stype in ["isoform_operon", "isoform_operon_merged", "isoform_gene_combined"]:
    n = (operons_merged["segmentation_type"] == stype).sum()
    print(f"  {stype:<30} {n}")


# In[32]:


# ── 2. Cross-strand overlap (antisense pairs) ───────────────────────────────
# Use operons_merged so cross-strand analysis reflects the post-merge state
plus_ops  = operons_merged[operons_merged["strand"] == "+"].copy()
minus_ops = operons_merged[operons_merged["strand"] == "-"].copy()
cross_strand = find_overlapping_pairs(plus_ops, minus_ops, same_strand=False)

print(f"\nCross-strand (antisense) overlapping operon pairs: {len(cross_strand)}")
if len(cross_strand) > 0:
    print(f"\nOverlap length distribution (bp):")
    print(cross_strand["overlap_bp"].describe().astype(int).to_string())
    print(f"\nRelationship breakdown:")
    print(cross_strand["relationship"].value_counts().to_string())
    print()
    display(cross_strand.sort_values("overlap_bp", ascending=False).head(20))


# ## Step 6 — Gene coverage
# 
# For each of the 911 annotated Syn1 genes, determine whether it is transcriptionally covered by the isoform-derived operons. Four coverage categories:
# 
# | Category | Meaning |
# |---|---|
# | `both` | Covered by operons on both strands (sense + antisense) |
# | `sense_only` | Covered only by a same-strand operon — normal transcription |
# | `antisense_only` | Only an antisense operon overlaps — the gene itself may be silent or very lowly expressed |
# | `uncovered` | No operon on either strand — candidate for Step 7 rescue |
# 
# Genes in the `uncovered` category are most likely lowly expressed (below `MIN_READS = 50`) rather than truly silent, since Syn1 is a minimal cell containing only essential genes.

# In[33]:


# ── Gene coverage: derived by inverting operons_annotated ──────────────────
def find_gene_coverage(genes: pd.DataFrame,
                       operons: pd.DataFrame) -> pd.DataFrame:
    """
    Invert sense/antisense loci columns in `operons` to map each gene back to
    its covering operons.  O(n_operons) — no per-gene operon scan needed.
    Returns genes_df augmented with:
      sense_operon_count / sense_covering_ops
      antisense_operon_count / antisense_covering_ops
      coverage_type: "both" | "sense_only" | "antisense_only" | "uncovered"
    """
    sense_map     = defaultdict(list)
    antisense_map = defaultdict(list)
    for _, op in operons.iterrows():
        oid = op["operon_id"]
        for locus in str(op.get("sense_gene_loci", "") or "").split(","):
            if locus:
                sense_map[locus].append(oid)
        for locus in str(op.get("antisense_gene_loci", "") or "").split(","):
            if locus:
                antisense_map[locus].append(oid)

    gdf = genes.copy()
    s_counts, s_ops, a_counts, a_ops, cov_types = [], [], [], [], []
    for _, gene in gdf.iterrows():
        lt = gene["locus_tag"]
        s  = sense_map.get(lt, [])
        a  = antisense_map.get(lt, [])
        s_counts.append(len(s));  s_ops.append(",".join(s))
        a_counts.append(len(a));  a_ops.append(",".join(a))
        if s and a:    cov_types.append("both")
        elif s:        cov_types.append("sense_only")
        elif a:        cov_types.append("antisense_only")
        else:          cov_types.append("uncovered")

    gdf["sense_operon_count"]     = s_counts
    gdf["sense_covering_ops"]     = s_ops
    gdf["antisense_operon_count"] = a_counts
    gdf["antisense_covering_ops"] = a_ops
    gdf["coverage_type"]          = cov_types
    return gdf

genes_covered = find_gene_coverage(genes_df, operons_merged)

n_both      = (genes_covered["coverage_type"] == "both").sum()
n_sense     = (genes_covered["coverage_type"] == "sense_only").sum()
n_anti_only = (genes_covered["coverage_type"] == "antisense_only").sum()
n_uncov     = (genes_covered["coverage_type"] == "uncovered").sum()
total       = len(genes_df)

print(f"Gene coverage summary (out of {total} genes):")
print(f"  Both sense + antisense:              {n_both:>4}  ({n_both/total*100:.1f}%)")
print(f"  Sense only:                          {n_sense:>4}  ({n_sense/total*100:.1f}%)")
print(f"  Antisense only:                      {n_anti_only:>4}  ({n_anti_only/total*100:.1f}%)")
print(f"  Uncovered (no operon either strand): {n_uncov:>4}  ({n_uncov/total*100:.1f}%)")

if n_anti_only > 0:
    print(f"\nGenes covered by antisense operon only:")
    display(genes_covered[genes_covered["coverage_type"] == "antisense_only"]
            [["locus_tag","gene_name","product","strand","start0","end0",
              "antisense_covering_ops"]].reset_index(drop=True))

uncovered_genes = genes_covered[genes_covered["coverage_type"] == "uncovered"].copy()
uncovered_genes["gene_length"] = uncovered_genes["end0"] - uncovered_genes["start0"]
print(f"\nUncovered genes by strand:")
print(uncovered_genes["strand"].value_counts().to_string())
print(f"\nUncovered gene length distribution (bp):")
print(uncovered_genes["gene_length"].describe().astype(int).to_string())
print(f"\nUncovered genes (first 30):")
display(uncovered_genes[["locus_tag","gene_name","product","chrom","strand",
                          "start0","end0","gene_length"]]
        .head(30).reset_index(drop=True))


# ## Step 7 — Rescue uncovered genes
# 
# Genes not covered by any isoform operon (Step 6 `uncovered`) are recovered through three complementary strategies:
# 
# **7a — BAM spanning read count:** Count PacBio FLNC reads that fully span each uncovered gene's coordinates in the raw BAM file. A read must map from ≤ gene start to ≥ gene end — this is stricter than counting any overlapping read and avoids contamination from adjacent high-expression operons. Even 1–2 spanning reads confirms basal transcription activity.
# 
# **7b — Consecutive gene rescue:** When multiple consecutive uncovered genes on the same strand are separated by gaps ≤ 500 bp, they likely form a single lowly-expressed operon. Any isoform (no read-depth filter) that spans the entire gene group is used as evidence; the best-supported spanning isoform defines the operon boundaries.
# 
# **7c — rRNA operons (hard-coded):** The two rRNA operons (5S–23S–16S on the minus strand, at ~92–97 kb and ~637–642 kb) are not captured by isoform clustering because rRNA is enormously abundant and heavily processed into mature subunits rather than read as full-length polycistronic transcripts. These operons are annotated directly from GFF3 coordinates.
# 
# All rescued operons are passed through `annotate_operons_with_genes` to populate `sense_gene_loci` / `antisense_gene_loci`, ensuring they are correctly counted in the final gene coverage summary.

# In[34]:


import pysam

BAM_FILE = MOTHER_FOLDER + "/Syn1_Transcriptomics/PacBio/PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam"

# ══════════════════════════════════════════════════════════════════════════════
# 4a — Count reads that span the ENTIRE gene region (reference_start ≤ gene_start
#      AND reference_end ≥ gene_end).  A simple bam.count() would include any
#      read that merely overlaps the region by even 1 bp — reads from adjacent
#      high-expression operons would inflate the count for low-expression neighbors.
#      Full-span counting is stricter and biologically more meaningful: it asks
#      "is there a transcript that covers this whole gene?"
# ══════════════════════════════════════════════════════════════════════════════
def count_spanning_reads(bam_path: str, chrom: str,
                         start0: int, end0: int) -> int:
    """Count reads whose mapped span fully covers [start0, end0)."""
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start0, end0):
            if (not read.is_unmapped and
                    read.reference_start <= start0 and
                    read.reference_end   >= end0):
                n += 1
    return n

spanning_counts = []
for _, gene in uncovered_genes.iterrows():
    spanning_counts.append(
        count_spanning_reads(BAM_FILE, gene["chrom"],
                             gene["start0"], gene["end0"])
    )

uncovered_genes = uncovered_genes.copy()
uncovered_genes["spanning_read_count"] = spanning_counts

print("Spanning-read distribution for uncovered genes:")
bins   = [0, 1, 5, 10, 20, 50, 100, 99999]
labels = ["0", "1-4", "5-9", "10-19", "20-49", "50-99", "≥100"]
for lo, hi, lbl in zip(bins, bins[1:], labels):
    n = ((uncovered_genes["spanning_read_count"] >= lo) &
         (uncovered_genes["spanning_read_count"] < hi)).sum()
    print(f"  {lbl:>6} reads: {n} genes")


# In[ ]:


uncovered_genes.sort_values("spanning_read_count", ascending=False).head(20)


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4b — Consecutive-gene rescue: check for spanning isoforms (any depth)
# ══════════════════════════════════════════════════════════════════════════════
# Sort uncovered genes by chrom, strand, start
ug = uncovered_genes.sort_values(["chrom", "strand", "start0"]).reset_index(drop=True)

# Group consecutive genes on same strand (gap ≤ MAX_GENE_GAP bp)
MAX_GENE_GAP = 500   # bp between adjacent gene ends — tune if needed

def group_consecutive(df):
    """Return list of gene-index groups that are consecutive on same strand."""
    groups = []
    current = [0]
    for i in range(1, len(df)):
        prev, cur = df.iloc[i - 1], df.iloc[i]
        if (cur["strand"] == prev["strand"] and
                cur["chrom"] == prev["chrom"] and
                cur["start0"] - prev["end0"] <= MAX_GENE_GAP):
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)
    return groups

cons_groups = group_consecutive(ug)
multi_groups = [g for g in cons_groups if len(g) > 1]
print(f"\nConsecutive uncovered gene groups (gap ≤ {MAX_GENE_GAP} bp): {len(multi_groups)}")

# For each multi-gene group, check if any isoform spans the full region
# (using ALL isoforms from df, regardless of read depth)
df_all_strand = {s: df[df["strand"] == s].reset_index(drop=True)
                 for s in ["+", "-"]}

spanning_rescues = []
for grp in multi_groups:
    genes_in_grp = ug.iloc[grp]
    strand  = genes_in_grp.iloc[0]["strand"]
    chrom   = genes_in_grp.iloc[0]["chrom"]
    reg_s   = genes_in_grp["start0"].min()
    reg_e   = genes_in_grp["end0"].max()
    loci    = ",".join(genes_in_grp["locus_tag"].tolist())

    iso_s = df_all_strand[strand]
    # isoform must span: start <= reg_s+BOUNDARY_TOL AND end >= reg_e-BOUNDARY_TOL
    spanning = iso_s[
        (iso_s["start0"] <= reg_s + BOUNDARY_TOL) &
        (iso_s["end0"]   >= reg_e - BOUNDARY_TOL)
    ]
    if len(spanning) > 0:
        best = spanning.loc[spanning["n_reads"].idxmax()]
        spanning_rescues.append({
            "chrom":      chrom,
            "strand":     strand,
            "start0":     int(reg_s),
            "end0":       int(reg_e),
            "n_genes":    len(grp),
            "gene_loci":  loci,
            "best_isoform": best["isoform_id"],
            "n_reads":    int(best["n_reads"]),
            "segmentation_type": "rescue_multiple",
        })

span_df = pd.DataFrame(spanning_rescues) if spanning_rescues else pd.DataFrame()
print(f"Consecutive groups rescued by spanning isoform: {len(span_df)}")
if len(span_df):
    display(span_df)


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4c — Build rescue operons
# ══════════════════════════════════════════════════════════════════════════════
rescued_loci = set()
if len(span_df):
    for loci_str in span_df["gene_loci"]:
        rescued_loci.update(loci_str.split(","))

# rRNA genes are handled separately in 4d — exclude them here
rrna_loci = {"MMSYN1_0067","MMSYN1_0068","MMSYN1_0069",
             "MMSYN1_0532","MMSYN1_0533","MMSYN1_0534"}
rescued_loci.update(rrna_loci)

remaining_ug = ug[~ug["locus_tag"].isin(rescued_loci)].copy()

single_gene_ops = []
for _, gene in remaining_ug.iterrows():
    single_gene_ops.append({
        "chrom":         gene["chrom"],
        "strand":        gene["strand"],
        "start0":        int(gene["start0"]),
        "end0":          int(gene["end0"]),
        "n_isoforms":    0,
        "n_reads_total": int(gene["spanning_read_count"]),
        "member_ids":    "",
        "gene_count":    1,
        "gene_loci":     gene["locus_tag"],
        "gene_names":    gene["gene_name"],
        "length":        int(gene["end0"] - gene["start0"]),
        "segmentation_type":   "rescue_single",
    })

print(f"Single-gene BAM-rescue operons: {len(single_gene_ops)}")
print(f"  (of which 0 spanning reads:  "
      f"{sum(1 for g in single_gene_ops if g['n_reads_total'] == 0)})")


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4d — rRNA operons (hard-coded)
# Both rRNA clusters are on the minus strand.
# Transcription goes 16S → 23S → 5S (high→low coordinate).
# So TSS = end0 (high coord), TTS = start0 (low coord).
# ══════════════════════════════════════════════════════════════════════════════
RRNA_OPERONS = [
    {   # cluster 1: MMSYN1_0067(5S)–0068(23S)–0069(16S), minus strand
        # GFF coords (1-based): rrf 92148-92256, rrl 92331-95225, rrs 95457-96980
        # 0-based: genomic min = 92147, genomic max = 96980
        "chrom": "CP002027.1", "strand": "-",
        "start0": 92147,   # genomic min (0-based) = TTS for minus strand
        "end0":   96980,   # genomic max           = TSS for minus strand
        "n_isoforms": 0, "n_reads_total": 0, "member_ids": "",
        "gene_count": 3,
        "gene_loci":  "MMSYN1_0067,MMSYN1_0068,MMSYN1_0069",
        "gene_names": "rrf,rrl,rrs",
        "segmentation_type": "rescue_multiple",
    },
    {   # cluster 2: MMSYN1_0532(5S)–0533(23S)–0534(16S), minus strand
        # GFF coords (1-based): rrf 637296-637404, rrl 637479-640373, rrs 640604-642127
        # 0-based: genomic min = 637295, genomic max = 642127
        "chrom": "CP002027.1", "strand": "-",
        "start0": 637295,
        "end0":   642127,
        "n_isoforms": 0, "n_reads_total": 0, "member_ids": "",
        "gene_count": 3,
        "gene_loci":  "MMSYN1_0532,MMSYN1_0533,MMSYN1_0534",
        "gene_names": "rrf,rrl,rrs",
        "segmentation_type": "rescue_multiple",
    },
]

# Count spanning reads across each rRNA operon
for rop in RRNA_OPERONS:
    rop["n_reads_total"] = count_spanning_reads(
        BAM_FILE, rop["chrom"], rop["start0"], rop["end0"])

print("rRNA operons:")
for rop in RRNA_OPERONS:
    print(f"  {rop['gene_loci']}  {rop['strand']}  "
          f"{rop['start0']}–{rop['end0']}  {rop['n_reads_total']} spanning reads")


# In[ ]:


# ══════════════════════════════════════════════════════════════════════════════
# 4e — Merge rescued operons + add TSS/TTS directional coordinates
# ══════════════════════════════════════════════════════════════════════════════
rescue_rows = []

if len(span_df):
    for _, r in span_df.iterrows():
        rescue_rows.append({
            "chrom":         r["chrom"],
            "strand":        r["strand"],
            "start0":        int(r["start0"]),
            "end0":          int(r["end0"]),
            "n_isoforms":    1,
            "n_reads_total": int(r["n_reads"]),
            "member_ids":    r["best_isoform"],
            "gene_count":    int(r["n_genes"]),
            "gene_loci":     r["gene_loci"],
            "gene_names":    "",
            "length":        int(r["end0"] - r["start0"]),
            "segmentation_type":   "rescue_multiple",
        })

rescue_rows.extend(single_gene_ops)
rescue_rows.extend(RRNA_OPERONS)

rescue_df = pd.DataFrame(rescue_rows)
rescue_df["length"] = rescue_df["end0"] - rescue_df["start0"]

# ── Populate sense/antisense gene loci for find_gene_coverage ──────────────
# rescue_df rows already have gene_loci set to exactly the intended locus/loci
# (one gene for single_gene_bam, N genes for consecutive_spanning/rRNA).
# We must NOT let annotate_operons_with_genes overwrite gene_loci/gene_count
# because it would pick up any gene that genomically overlaps the coordinates,
# potentially adding unintended neighbours to a single_gene_bam operon.
# Instead, copy gene_loci directly into sense_gene_loci (all rescue operons
# are sense by construction) and compute antisense_gene_loci independently.
opp = {"+": "-", "-": "+"}
s_loci, s_names, s_counts = [], [], []
a_loci, a_names, a_counts = [], [], []
for _, op in rescue_df.iterrows():
    overlaps  = (genes_df["start0"] < op["end0"]) & (genes_df["end0"] > op["start0"])
    antisense = genes_df[overlaps & (genes_df["strand"] == opp[op["strand"]])]
    # sense: use the already-correct gene_loci, not a fresh scan
    s_loci.append(str(op.get("gene_loci", "") or ""))
    s_names.append(str(op.get("gene_names", "") or ""))
    s_counts.append(len([l for l in str(op.get("gene_loci","") or "").split(",") if l]))
    a_loci.append(",".join(antisense["locus_tag"].tolist()))
    a_names.append(",".join(antisense["gene_name"].tolist()))
    a_counts.append(len(antisense))
rescue_df["sense_gene_loci"]      = s_loci
rescue_df["sense_gene_names"]      = s_names
rescue_df["sense_gene_count"]      = s_counts
rescue_df["antisense_gene_loci"]   = a_loci
rescue_df["antisense_gene_names"]  = a_names
rescue_df["antisense_gene_count"]  = a_counts

# Tag segmentation_type on operons_merged:
# rows that went through ss_merge already have "isoform_operon_merged";
# rows that were never in a conflict pair have NaN — fill those now.
operons_merged = operons_merged.copy()
operons_merged["segmentation_type"] = (
    operons_merged.get("segmentation_type", pd.Series(dtype=str))
    .fillna("isoform_operon")
)

# Combine all operons — sort by genomic position (start0 always = min coord)
all_operons = pd.concat(
    [operons_merged, rescue_df], ignore_index=True
).sort_values(["chrom", "start0"]).reset_index(drop=True)

# Drop any pre-existing operon_id before re-assigning
if "operon_id" in all_operons.columns:
    all_operons = all_operons.drop(columns=["operon_id"])
all_operons.insert(0, "operon_id",
    [f"OP_{i+1:05d}" for i in range(len(all_operons))])

# ── Add TSS / TTS directional coordinates ───────────────────────────────────
# For + strand: TSS = start0 (5′ end, low coord), TTS = end0 (3′ end, high coord)
# For − strand: TSS = end0  (5′ end, high coord), TTS = start0 (3′ end, low coord)
all_operons["tss"] = np.where(all_operons["strand"] == "+",
                               all_operons["start0"], all_operons["end0"])
all_operons["tts"] = np.where(all_operons["strand"] == "+",
                               all_operons["end0"],  all_operons["start0"])

print(f"Final operon count: {len(all_operons)}")
for stype in all_operons["segmentation_type"].unique():
    n = (all_operons["segmentation_type"] == stype).sum()
    print(f"  {stype:<30} {n}")
n_nan = all_operons["segmentation_type"].isna().sum()
if n_nan:
    print(f"  WARNING: {n_nan} operons with no segmentation_type (NaN)")

print("\nRescue operons (non-isoform):")
display(all_operons[all_operons["segmentation_type"] != "isoform_operon"]
        [["operon_id","strand","tss","tts","length","gene_loci",
          "n_reads_total","segmentation_type"]].reset_index(drop=True))


# ## Update isoform memberships

# ## Step 8 — Save outputs
# 
# Write the final operon table and gene coverage files. The operon table (`operons.candidate_blocks.tsv`) includes all segmentation types:
# 
# | `segmentation_type` | Origin |
# |---|---|
# | `isoform_operon` | Directly from containment clustering (Step 2) |
# | `isoform_operon_merged` | Two overlapping isoform operons merged by TSS rule (Step 5a) |
# | `isoform_gene_combined` | Two isoform operons sharing a gene merged by duplicate-loci rule (Step 5b) |
# | `rescue_multiple` | Consecutive low-expression genes or rRNA operon (Step 7b/7c) |
# | `rescue_single` | Single lowly-expressed gene from BAM evidence (Step 7a) |
# 
# BED6 files (plus and minus strand separately) are written for direct loading into IGV.

# In[ ]:


# ── Final deduplication of per-operon gene lists ────────────────────────────
# The Step-5a overlap merge concatenates the gene lists of the two merged operons
# and SUMS their counts, so a gene shared by both is listed and counted twice; the
# Step-5b merge only catches fully-identical-loci pairs, not partial overlaps. This
# final order-preserving pass guarantees unique gene membership in the written table.
def dedup_operon_gene_lists(df):
    df = df.copy()
    def _dedup(loci_str, name_str):
        loci  = [x for x in str(loci_str or "").split(",") if x]
        names = str(name_str or "").split(",")
        seen, ul, un = set(), [], []
        for i, lc in enumerate(loci):
            if lc in seen:
                continue
            seen.add(lc); ul.append(lc); un.append(names[i] if i < len(names) else "")
        return ",".join(ul), ",".join(un), len(ul)
    for kind in ("sense", "antisense"):
        out = [_dedup(r[f"{kind}_gene_loci"], r[f"{kind}_gene_names"]) for _, r in df.iterrows()]
        df[f"{kind}_gene_loci"]  = [o[0] for o in out]
        df[f"{kind}_gene_names"] = [o[1] for o in out]
        df[f"{kind}_gene_count"] = [o[2] for o in out]
    df["gene_count"] = df["sense_gene_count"] + df["antisense_gene_count"]
    if "gene_loci" in df.columns:
        df["gene_loci"]  = df["sense_gene_loci"]
        df["gene_names"] = df["sense_gene_names"]
    return df

n_dup_before = int(sum(len([x for x in str(r["sense_gene_loci"] or "").split(",") if x])
                       - len({x for x in str(r["sense_gene_loci"] or "").split(",") if x})
                       for _, r in all_operons.iterrows()))
all_operons = dedup_operon_gene_lists(all_operons)
print(f"Final dedup: removed {n_dup_before} duplicated sense-locus slots across operons.")

# ── Final gene coverage check ───────────────────────────────────────────────
genes_covered_final = find_gene_coverage(genes_df, all_operons)

n_both      = (genes_covered_final["coverage_type"] == "both").sum()
n_sense     = (genes_covered_final["coverage_type"] == "sense_only").sum()
n_anti_only = (genes_covered_final["coverage_type"] == "antisense_only").sum()
n_uncov     = (genes_covered_final["coverage_type"] == "uncovered").sum()
total       = len(genes_df)

print(f"Final gene coverage ({total} genes):")
print(f"  Both sense + antisense:              {n_both:>4}  ({n_both/total*100:.1f}%)")
print(f"  Sense only:                          {n_sense:>4}  ({n_sense/total*100:.1f}%)")
print(f"  Antisense only:                      {n_anti_only:>4}  ({n_anti_only/total*100:.1f}%)")
print(f"  Uncovered (no operon either strand): {n_uncov:>4}  ({n_uncov/total*100:.1f}%)")

if n_uncov > 0:
    display(genes_covered_final[genes_covered_final["coverage_type"] == "uncovered"]
            [["locus_tag","gene_name","product","strand","start0","end0"]]
            .reset_index(drop=True))

# ── Save all outputs ─────────────────────────────────────────────────────────
OUT_TSV = OUT_FOLDER + "/operons.candidate_blocks.tsv"
all_operons.to_csv(OUT_TSV, sep="\t", index=False)
print(f"\nSaved: {OUT_TSV}  ({len(all_operons)} operons)")
print(f"Columns: {list(all_operons.columns)}")

OUT_GCOV = SEG_FOLDER + "/gene_operon_coverage.tsv"
genes_covered_final.to_csv(OUT_GCOV, sep="\t", index=False)
print(f"Saved: {OUT_GCOV}")

OUT_UNCOV = SEG_FOLDER + "/uncovered_genes.tsv"
still_uncov = genes_covered_final[genes_covered_final["coverage_type"] == "uncovered"]
still_uncov.to_csv(OUT_UNCOV, sep="\t", index=False)
print(f"Saved: {OUT_UNCOV}  ({len(still_uncov)} genes)")

# ── BED6 for IGV ─────────────────────────────────────────────────────────────
def to_bed6(df, score_col="n_reads_total", name_col="operon_id"):
    bed = df[["chrom", "start0", "end0", name_col, score_col, "strand"]].copy()
    bed.columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    bed["score"] = bed["score"].clip(upper=1000).astype(int)
    return bed

to_bed6(all_operons[all_operons["strand"] == "+"]).to_csv(
    SEG_FOLDER + "/operons.candidate_blocks.plus.bed",  sep="\t", index=False, header=False)
to_bed6(all_operons[all_operons["strand"] == "-"]).to_csv(
    SEG_FOLDER + "/operons.candidate_blocks.minus.bed", sep="\t", index=False, header=False)
print("BED files written for IGV inspection.")


# ## Step 9 — Run summary (Operon_Segmentation.txt)
#
# Persist the stage-by-stage operon counts to a plain-text summary (mirrors
# Operon_Annotation.txt) so the Methods/Results text can quote them without
# re-parsing stdout. Written LAST: operons.candidate_blocks.tsv is already saved
# above, so any issue here cannot corrupt the canonical operon table.

# In[ ]:


def _seg_get(name, default=None):
    """Fetch a module-level variable that may be defined only on some branches."""
    return globals().get(name, default)

# ── Stage tallies (recomputed from the surviving stage DataFrames) ──────────
_init_total = len(operons_all)                                  # Step 2 initial operons
_init_plus  = len(operons_plus)
_init_minus = len(operons_minus)
_n_s0 = int((operons_annotated["sense_gene_count"] == 0).sum())  # Step 3 sense breakdown
_n_s1 = int((operons_annotated["sense_gene_count"] == 1).sum())
_n_s2 = int((operons_annotated["sense_gene_count"] == 2).sum())
_n_s3 = int((operons_annotated["sense_gene_count"] >= 3).sum())
_sense_ge1   = _n_s1 + _n_s2 + _n_s3
_n_anti_init = int((operons_annotated["antisense_gene_count"] > 0).sum())

if len(decisions_df):                                            # Step 5a candidate classes
    _cand_overlap = int((decisions_df["relationship"] == "overlap").sum())
    _cand_gig     = int((decisions_df["relationship"] == "gene_in_gap").sum())
    _n_pass       = int(decisions_df["merge"].sum())             # pairs passing co-transcription test
else:
    _cand_overlap = _cand_gig = _n_pass = 0
# NB: compute the separate count locally — the global `n_separate` is reassigned
# later in the overlap-plotting block, so it no longer holds the merge-decision split.
_n_sep = len(decisions_df) - _n_pass

_iso_types    = operons_merged["segmentation_type"].value_counts().to_dict()  # after 5a/5b
_n_iso_plain  = int(_iso_types.get("isoform_operon", 0))
_n_iso_merged = int(_iso_types.get("isoform_operon_merged", 0))
_n_iso_comb   = int(_iso_types.get("isoform_gene_combined", 0))
_iso_total    = len(operons_merged)

_pre_cov   = genes_covered["coverage_type"].value_counts().to_dict()          # Step 6 (pre-rescue)
_pre_uncov = int(_pre_cov.get("uncovered", 0))

_n_multi_groups = len(_seg_get("multi_groups", []))             # Step 7 rescue tallies
_n_span_rescue  = len(_seg_get("span_df", []))
_n_single       = len(_seg_get("single_gene_ops", []))
_n_rrna         = len(RRNA_OPERONS)

_final_types = all_operons["segmentation_type"].value_counts()  # Step 8 final

_lines = []
def _w(s=""):
    _lines.append(str(s))

_w("OPERON SEGMENTATION SUMMARY")
_w("=" * 60)
_w("")
_w("Stage-by-stage operon counts from Operon_Segmentation.py "
   "(isoform-based segmentation of the syn1 PacBio FLNC transcriptome).")
_w("")
_w(f"Parameters: MIN_READS={MIN_READS}, BOUNDARY_TOL={BOUNDARY_TOL} bp; "
   f"merge: MERGE_MAX_GAP={MERGE_MAX_GAP}, MERGE_W={MERGE_W}, MERGE_FLANK={MERGE_FLANK}, "
   f"MERGE_MIN_BRIDGE={MERGE_MIN_BRIDGE}, MERGE_MIN_CONT={MERGE_MIN_CONT}; "
   f"rescue: MAX_GENE_GAP={MAX_GENE_GAP} bp")
_w("")
_w("Step 1 -- Load isoforms")
_w(f"  isoform clusters loaded:        {len(df):>6}")
_w(f"  isoforms with n_reads >= {MIN_READS}:   {len(df_iso):>6}  "
   f"(plus {int((df_iso['strand']=='+').sum())}, minus {int((df_iso['strand']=='-').sum())})")
_w("")
_w("Step 2 -- Containment clustering into initial operons")
_w(f"  initial operons:                {_init_total:>6}  (plus {_init_plus}, minus {_init_minus})")
_w(f"  operon length (bp):  median {operons_all['length'].median():.0f}, mean {operons_all['length'].mean():.0f}")
_w("")
_w("Step 3 -- Gene annotation of initial operons")
_w(f"  annotated genes loaded:         {len(genes_df):>6}")
_w(f"  operons with >=1 sense gene:    {_sense_ge1:>6}  "
   f"(1 gene {_n_s1}, 2 genes {_n_s2}, 3+ genes {_n_s3}; 0 sense {_n_s0})")
_w(f"  operons with >=1 antisense gene:{_n_anti_init:>6}")
_w("")
_w("Step 4/5a -- Same-strand overlap + co-transcription merge")
_w(f"  same-strand overlap pairs:      {len(same_strand_overlaps):>6}  "
   f"(shared-gene pairs: {len(_seg_get('conflicts', []))})")
_w(f"  merge candidates:               {len(decisions_df):>6}  "
   f"(overlap {_cand_overlap}, gene_in_gap {_cand_gig})")
_w(f"    passed co-transcription test: {_n_pass:>6}")
_w(f"    kept separate:                {_n_sep:>6}")
_w(f"    -> merged output operons (pairwise, no chaining): {_n_iso_merged}")
_w("")
_w("Step 5b -- Merge operons sharing identical sense loci")
_w(f"  isoform-derived operons after merges: {_iso_total}")
_w(f"    isoform_operon         {_n_iso_plain:>5}")
_w(f"    isoform_operon_merged  {_n_iso_merged:>5}")
_w(f"    isoform_gene_combined  {_n_iso_comb:>5}")
_w("")
_w("Step 6/7 -- Coverage + rescue of uncovered genes")
_w(f"  uncovered genes (pre-rescue):   {_pre_uncov:>6} / {len(genes_df)} "
   f"({_pre_uncov/len(genes_df)*100:.1f}%)")
_w(f"  consecutive uncovered groups (gap <= {MAX_GENE_GAP} bp): {_n_multi_groups}")
_w(f"    rescued by spanning isoform:  {_n_span_rescue:>6}")
_w(f"  rRNA operons added:             {_n_rrna:>6}")
_w(f"  single-gene BAM rescues:        {_n_single:>6}")
_w("")
_w("Step 8 -- Final operon map")
_w(f"  total operons:                  {len(all_operons):>6}")
for _stype, _n in _final_types.items():
    _w(f"    {_stype:<24} {_n:>5}")
_w("")
_w(f"  Final gene coverage ({total} genes):")
_w(f"    both sense + antisense:       {n_both:>6}  ({n_both/total*100:.1f}%)")
_w(f"    sense only:                   {n_sense:>6}  ({n_sense/total*100:.1f}%)")
_w(f"    antisense only:               {n_anti_only:>6}  ({n_anti_only/total*100:.1f}%)")
_w(f"    uncovered:                    {n_uncov:>6}  ({n_uncov/total*100:.1f}%)")

OUT_SUMMARY = OUT_FOLDER + "/Operon_Segmentation.txt"
with open(OUT_SUMMARY, "w") as _fh:
    _fh.write("\n".join(_lines) + "\n")
print(f"Saved: {OUT_SUMMARY}")

