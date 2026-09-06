#!/usr/bin/env python3
"""
coexpression_common.py

Shared primitives for the syn3A co-transcription tests used by
06_single_operon_coexpression.py (operon-internal pairs) and
07_operon_pair_coexpression.py (cross-junction operon pairs).

The one shared idea: a (gene_a, gene_b) pair on the same strand is "co-
transcribed" in syn3A when an ONT direct-RNA read crosses the boundary between
them. Two stringencies:
  strict (spanning) : >=1 ONT read fully encloses CDS_A + intergenic + CDS_B
  loose  (bridging) : ONT bridging-read count >= max(1, min(CAP, 20% * min ONT
                      gene depth)) AND (gap < 10 bp OR mean Illumina gap depth
                      >= 20% * min Illumina gene depth)
All depths are matching-strand, mean per-base.
"""

from __future__ import annotations

import re
from pathlib import Path
import numpy as np
import pandas as pd
import pysam


# ---------- paths (shared inputs) ----------
GR   = Path(__file__).resolve().parent          # Genome_Reduction/
ROOT = GR.parent

SYN1_GFF  = ROOT / "Genomes_Input" / "syn1.genes.gff3"
SYN3A_GFF = ROOT / "Genomes_Input" / "syn3a_genome.gff3"
SYN1_OPERONS_TSV = ROOT / "Syn1_Operon" / "operons.candidate_blocks.tsv"

ONT_BAM      = ROOT / "Syn3A_Transcriptomics" / "ONT" / "ONT_Processing" / "syn3A.ONT.rep1.sorted.bam"
ONT_BG_PLUS  = ROOT / "Syn3A_Transcriptomics" / "ONT" / "ONT_Processing" / "depth_bedgraph" / "syn3A.ONT.rep1.plus.bedGraph"
ONT_BG_MINUS = ROOT / "Syn3A_Transcriptomics" / "ONT" / "ONT_Processing" / "depth_bedgraph" / "syn3A.ONT.rep1.minus.bedGraph"
ILL_BG_PLUS  = ROOT / "Syn3A_Transcriptomics" / "Illumina" / "Illumina_Processing" / "depth_bedgraph" / "syn3A_rep1.plus.bedGraph"
ILL_BG_MINUS = ROOT / "Syn3A_Transcriptomics" / "Illumina" / "Illumina_Processing" / "depth_bedgraph" / "syn3A_rep1.minus.bedGraph"

# 04 / 05 outputs consumed by 06 / 07
OPERON_CLASS_TSV = GR / "deletion_overlaid_operon" / "operon_deletion_classification.tsv"
JUNCTIONS_TSV    = GR / "deletion_junction" / "deletion_junctions.tsv"


# ---------- co-expression decision thresholds ----------
MIN_SPAN_FOR_PRESERVED            = 1     # strict: >= this many spanning ONT reads
BRIDGE_FRACTION_OF_MIN_GENE_DEPTH = 0.20  # loose: bridging reads vs min ONT gene depth
BRIDGE_ABS_CAP                    = 10    # cap so high-depth pairs don't need impractical counts
GAP_DEPTH_FRACTION_OF_MIN_GENE_DEPTH = 0.20  # Illumina gap depth vs min Illumina gene depth
MIN_GAP_LEN_FOR_DEPTH_CHECK          = 10    # bp; below this skip the Illumina check
BRIDGE_MIN_OVERLAP_BP                = 50    # min bp a bridging read must cover inside each CDS


# ---------- loaders ----------

_PAT_LOCUS = re.compile(r"locus_tag=([^;]+)")
_PAT_NAME  = re.compile(r"(?:Name|gene)=([^;]+)")


def load_gff_genes(gff: Path) -> pd.DataFrame:
    rows = []
    with gff.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = _PAT_LOCUS.search(f[8])
            if not m:
                continue
            n = _PAT_NAME.search(f[8])
            rows.append({
                "locus_tag": m.group(1),
                "name":      n.group(1) if n else "",
                "chrom":     f[0],
                "strand":    f[6],
                "start0":    int(f[3]) - 1,
                "end0":      int(f[4]),
            })
    return pd.DataFrame(rows).drop_duplicates(subset="locus_tag")


def locus_suffix(lt: str) -> str:
    return lt.rsplit("_", 1)[-1]


def load_bedgraph(path: Path) -> dict:
    """bedGraph -> {chrom: np.float32 per-base depth array} (absent = 0)."""
    chrom_max: dict = {}
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line or line[0] == "#" or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]; s = int(parts[1]); e = int(parts[2]); v = float(parts[3])
            rows.append((chrom, s, e, v))
            if e > chrom_max.get(chrom, 0):
                chrom_max[chrom] = e
    depth = {c: np.zeros(L, dtype=np.float32) for c, L in chrom_max.items()}
    for chrom, s, e, v in rows:
        depth[chrom][s:e] = v
    return depth


def mean_depth(arr_by_chrom: dict, chrom: str, start: int, end: int) -> float:
    arr = arr_by_chrom.get(chrom)
    if arr is None:
        return 0.0
    s = max(0, start); e = min(len(arr), end)
    if e <= s:
        return 0.0
    return float(arr[s:e].mean())


def pick_strand(plus_arr: dict, minus_arr: dict, strand: str) -> dict:
    return plus_arr if strand == "+" else minus_arr


def load_depths() -> dict:
    """Load all four bedGraphs once; returns dict with keys ont_plus/ont_minus/ill_plus/ill_minus."""
    return {
        "ont_plus":  load_bedgraph(ONT_BG_PLUS),
        "ont_minus": load_bedgraph(ONT_BG_MINUS),
        "ill_plus":  load_bedgraph(ILL_BG_PLUS),
        "ill_minus": load_bedgraph(ILL_BG_MINUS),
    }


# ---------- ONT read primitives ----------

def count_spanning_reads(bam: pysam.AlignmentFile, chrom: str,
                         span_s: int, span_e: int, strand: str) -> int:
    """STRICT: primary ONT alignments on `strand` whose span encloses [span_s, span_e)."""
    n = 0
    for read in bam.fetch(chrom, span_s, span_e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if ("-" if read.is_reverse else "+") != strand:
            continue
        if read.reference_start <= span_s and read.reference_end >= span_e:
            n += 1
    return n


def count_bridging_reads(bam: pysam.AlignmentFile, chrom: str,
                         a_s: int, a_e: int, b_s: int, b_e: int, strand: str,
                         min_overlap: int = BRIDGE_MIN_OVERLAP_BP) -> int:
    """LOOSE: primary ONT alignments on `strand` overlapping each CDS by
    >= min_overlap bp and fully covering the intergenic region between them."""
    if a_s <= b_s:
        L_s, L_e, R_s, R_e = a_s, a_e, b_s, b_e
    else:
        L_s, L_e, R_s, R_e = b_s, b_e, a_s, a_e
    n = 0
    for read in bam.fetch(chrom, min(L_s, R_s), max(L_e, R_e)):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if ("-" if read.is_reverse else "+") != strand:
            continue
        r_s, r_e = read.reference_start, read.reference_end
        if max(0, min(r_e, L_e) - max(r_s, L_s)) < min_overlap:
            continue
        if max(0, min(r_e, R_e) - max(r_s, R_s)) < min_overlap:
            continue
        if L_e <= R_s and (r_s > L_e or r_e < R_s):
            continue
        n += 1
    return n


# ---------- the shared pair test ----------

def test_pair(bam, depths: dict, chrom: str,
              a_s: int, a_e: int, b_s: int, b_e: int, strand: str) -> dict:
    """Co-transcription test for one same-strand gene pair. Returns metrics +
    pair_preserved_strict / pair_preserved_loose. Caller supplies syn3A coords."""
    span_s, span_e = min(a_s, b_s), max(a_e, b_e)
    n_span   = count_spanning_reads(bam, chrom, span_s, span_e, strand)
    n_bridge = count_bridging_reads(bam, chrom, a_s, a_e, b_s, b_e, strand)

    ont_arr = pick_strand(depths["ont_plus"], depths["ont_minus"], strand)
    ill_arr = pick_strand(depths["ill_plus"], depths["ill_minus"], strand)
    gap_s, gap_e = min(a_e, b_e), max(a_s, b_s)
    gap_len = max(0, gap_e - gap_s)

    ont_a = mean_depth(ont_arr, chrom, a_s, a_e)
    ont_b = mean_depth(ont_arr, chrom, b_s, b_e)
    ill_a = mean_depth(ill_arr, chrom, a_s, a_e)
    ill_b = mean_depth(ill_arr, chrom, b_s, b_e)
    ill_gap = mean_depth(ill_arr, chrom, gap_s, gap_e) if gap_len > 0 else float("nan")

    bridge_threshold = max(1.0, min(BRIDGE_ABS_CAP,
                                    BRIDGE_FRACTION_OF_MIN_GENE_DEPTH * min(ont_a, ont_b)))
    gap_threshold = GAP_DEPTH_FRACTION_OF_MIN_GENE_DEPTH * min(ill_a, ill_b)

    spanned = n_span >= MIN_SPAN_FOR_PRESERVED
    bridge_pass = n_bridge >= bridge_threshold
    gap_pass = True if gap_len < MIN_GAP_LEN_FOR_DEPTH_CHECK else (ill_gap >= gap_threshold)
    bridged = bridge_pass and gap_pass

    return {
        "syn3A_chrom":         chrom,
        "syn3A_intergenic_bp": gap_len,
        "n_spanning_reads":    n_span,
        "n_bridging_reads":    n_bridge,
        "ont_depth_a":         round(ont_a, 3),
        "ont_depth_b":         round(ont_b, 3),
        "ill_depth_a":         round(ill_a, 3),
        "ill_depth_b":         round(ill_b, 3),
        "ill_gap_depth":       (round(ill_gap, 3) if gap_len > 0 else float("nan")),
        "bridge_threshold":    round(bridge_threshold, 3),
        "gap_depth_threshold": round(gap_threshold, 3),
        "bridge_pass":         bridge_pass,
        "gap_depth_pass":      gap_pass,
        "pair_preserved_strict": spanned,
        "pair_preserved_loose":  bridged,
    }
