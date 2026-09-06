#!/usr/bin/env python
# coding: utf-8
"""
Novel_translation.py -- novel transcription/translation discovery from syn1 PacBio
long-read RNA-seq (script form of Novel_translation.ipynb).

Purpose: annotate potential novel peptides/proteins from full-length PacBio
isoforms and infer translation trap/stall mechanisms. The output FASTA is used in
mass-spectrometry-based proteogenomics search.

Inputs:
  - Annotated isoform-cluster table (built from the mapped, QC'd BAM): the CURRENT
    canonical file, post the Apr-22 clustering revision,
    ../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv
  - Annotated genome GFF3 + FASTA; canonical proteome FASTA (tryptic uniqueness).
Outputs:
  - candidate_orfs/  : all candidate ORFs (OSTIR-scored), scored set, top-100
                       strict set, and FASTA (novel ORFs / unique peptides).
  - trypsin_digest/  : in-silico tryptic peptides + MS-identifiable FASTA.

Pipeline (faithful conversion of notebook cells 2 -> 4 -> 5 -> 9 -> 12 -> 13 -> 15):
  1. Load the isoform table and build the abnormal ORF-input set.
  2. For every abnormal isoform, run OSTIR over all start codons, build the
     downstream ORF (genetic code 4, TGA = Trp), annotate overlap/initiation.
  3. Score by synthesis rate (reads x TIR), take the top 100, collapse to unique
     peptides, and run in-silico trypsin digestion for MS-detectability.

Run from Syn1_Novel_ORF/.
"""

from __future__ import annotations
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
from ostir import run_ostir

# ======================================================================== paths
ISOFORMS_TSV    = '../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv'
FASTA_PATH      = '../Genomes_Input/syn1_genome.fasta'
GFF_PATH        = '../Genomes_Input/syn1.genes.gff3'
CANONICAL_FASTA = '../Genomes_Input/syn1_proteins.faa'

OUT_FOLDER        = 'candidate_orfs'
OUT_DIGEST_FOLDER = 'trypsin_digest'
ISO_FILTERED_XLSX = 'isoform_abnormal_filtered.xlsx'

# ======================================================================== params
# OSTIR
ANTI_SD       = 'ACCTCCTTT'
OSTIR_THREADS = 8
# TIR thresholds (from canonical Syn1 gene distribution)
TIR_STRONG_THRESHOLD = 10000
TIR_MOD_THRESHOLD    = 1000
TIR_WEAK_THRESHOLD   = 100
LEADERLESS_MAX_DIST  = 5
DG_MRNA_RRNA_SD_THRESHOLD = -7.2
# ORF length thresholds
MIN_AA_LEN_MAIN  = 24
MIN_AA_LEN_MICRO = 15
OVERLAP_PURE_FRAC = 0.95
# trypsin
MAX_MISSED_CLEAVAGES = 1
MIN_PEPTIDE_LEN = 7
MAX_PEPTIDE_LEN = 25

STOP_CODONS = {"TAA", "TAG"}                       # genetic code 4: TGA = Trp
GENETIC_CODE_4 = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ======================================================================== helpers
def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def translate_code4(nt_seq: str) -> str:
    seq = nt_seq.upper()
    usable = (len(seq) // 3) * 3
    seq = seq[:usable]
    return "".join(GENETIC_CODE_4.get(seq[i:i + 3], "X") for i in range(0, len(seq), 3))


def load_genes_from_gff(gff_path: str) -> pd.DataFrame:
    def parse_attrs(attr_str):
        d = {}
        for field in str(attr_str).split(";"):
            field = field.strip()
            if "=" in field:
                k, v = field.split("=", 1)
                d[k.strip()] = v.strip()
        return d

    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, source, feature, start1, end1, score, strand, phase, attrs = parts
            if feature != "gene" or strand not in {"+", "-"}:
                continue
            attrd = parse_attrs(attrs)
            gene_id = attrd.get("locus_tag", attrd.get("ID", ""))
            gene_name = attrd.get("gene", attrd.get("Name", gene_id))
            rows.append({
                "chrom": chrom, "start0": int(start1) - 1, "end0": int(end1),
                "strand": strand, "gene_id": gene_id, "gene_name": gene_name,
            })
    return pd.DataFrame(rows).sort_values(["chrom", "start0"]).reset_index(drop=True)


def tx_local_to_genomic(iso_start0, iso_end0, strand, local_start0, local_end0):
    if strand == "+":
        return iso_start0 + local_start0, iso_start0 + local_end0
    return iso_end0 - local_end0, iso_end0 - local_start0


def merged_union_len(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged = []
    for s, e in intervals:
        if e <= s:
            continue
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return sum(e - s for s, e in merged)


def annotate_orf_overlap(chrom, strand, g0, g1, genes_chr):
    orf_len = int(g1) - int(g0)
    if orf_len <= 0:
        return {
            "orf_bp_in_sense": 0, "orf_bp_in_antisense": 0, "orf_bp_in_intergenic": 0,
            "orf_frac_in_sense": 0.0, "orf_frac_in_antisense": 0.0, "orf_frac_in_intergenic": 0.0,
            "orf_sense_gene_ids": "", "orf_sense_gene_names": "",
            "orf_antisense_gene_ids": "", "orf_antisense_gene_names": "",
            "orf_n_sense_genes": 0, "orf_n_antisense_genes": 0,
            "orf_novelty_class": "invalid", "orf_abnormal_frac": 0.0, "orf_abnormal_bp": 0,
        }
    cand = genes_chr[(genes_chr["end0"] > g0) & (genes_chr["start0"] < g1)]
    sense_iv, anti_iv, any_iv = [], [], []
    sense_ids, sense_names, anti_ids, anti_names = [], [], [], []
    for _, g in cand.iterrows():
        s, e = max(g0, int(g["start0"])), min(g1, int(g["end0"]))
        if e <= s:
            continue
        any_iv.append((s, e))
        if g["strand"] == strand:
            sense_iv.append((s, e)); sense_ids.append(str(g["gene_id"])); sense_names.append(str(g["gene_name"]))
        else:
            anti_iv.append((s, e)); anti_ids.append(str(g["gene_id"])); anti_names.append(str(g["gene_name"]))
    bp_s, bp_a, bp_any = merged_union_len(sense_iv), merged_union_len(anti_iv), merged_union_len(any_iv)
    bp_ig = max(0, orf_len - bp_any)
    abn_bp = bp_a + bp_ig
    abn_frac = abn_bp / orf_len
    if bp_ig / orf_len >= OVERLAP_PURE_FRAC:
        nov = "novel_intergenic_orf"
    elif bp_a / orf_len >= OVERLAP_PURE_FRAC:
        nov = "antisense_orf"
    elif bp_s / orf_len >= OVERLAP_PURE_FRAC:
        nov = "known_sense_orf"
    elif abn_frac >= 0.20 and bp_s > 0:
        nov = "mixed_boundary_orf"
    else:
        nov = "mixed_complex_orf"
    return {
        "orf_bp_in_sense": int(bp_s), "orf_bp_in_antisense": int(bp_a), "orf_bp_in_intergenic": int(bp_ig),
        "orf_frac_in_sense": bp_s / orf_len, "orf_frac_in_antisense": bp_a / orf_len,
        "orf_frac_in_intergenic": bp_ig / orf_len,
        "orf_sense_gene_ids": ",".join(sorted(set(sense_ids))), "orf_sense_gene_names": ",".join(sorted(set(sense_names))),
        "orf_antisense_gene_ids": ",".join(sorted(set(anti_ids))), "orf_antisense_gene_names": ",".join(sorted(set(anti_names))),
        "orf_n_sense_genes": len(set(sense_ids)), "orf_n_antisense_genes": len(set(anti_ids)),
        "orf_novelty_class": nov, "orf_abnormal_frac": float(abn_frac), "orf_abnormal_bp": int(abn_bp),
    }


def low_complexity_peptide(aa_seq):
    if not aa_seq:
        return True
    freqs = pd.Series(list(aa_seq.upper())).value_counts(normalize=True)
    return float(freqs.max()) >= 0.60


def classify_initiation_ostir(tir, dg_mrna_rrna, leaderless_dist):
    has_sd = (not math.isnan(dg_mrna_rrna)) and (dg_mrna_rrna < DG_MRNA_RRNA_SD_THRESHOLD)
    is_ll = leaderless_dist <= LEADERLESS_MAX_DIST
    if math.isnan(tir) or tir <= 0:
        if is_ll:
            return {"initiation_mode": "leaderless", "initiation_score": 1, "has_sd": False}
        return {"initiation_mode": "none", "initiation_score": 0, "has_sd": False}
    if is_ll and tir >= TIR_WEAK_THRESHOLD:
        sc = 3 if tir >= TIR_STRONG_THRESHOLD else (2 if tir >= TIR_MOD_THRESHOLD else 1)
        return {"initiation_mode": "leaderless", "initiation_score": sc, "has_sd": has_sd}
    if is_ll:
        return {"initiation_mode": "leaderless", "initiation_score": 1, "has_sd": has_sd}
    if tir >= TIR_STRONG_THRESHOLD and has_sd:
        return {"initiation_mode": "sd_strong", "initiation_score": 3, "has_sd": True}
    if tir >= TIR_MOD_THRESHOLD and has_sd:
        return {"initiation_mode": "sd_moderate", "initiation_score": 2, "has_sd": True}
    if tir >= TIR_WEAK_THRESHOLD and has_sd:
        return {"initiation_mode": "sd_weak", "initiation_score": 1, "has_sd": True}
    if tir >= TIR_MOD_THRESHOLD and not has_sd:
        return {"initiation_mode": "sd_independent", "initiation_score": 2, "has_sd": False}
    if tir >= TIR_WEAK_THRESHOLD and not has_sd:
        return {"initiation_mode": "sd_independent_weak", "initiation_score": 1, "has_sd": False}
    return {"initiation_mode": "none", "initiation_score": 0, "has_sd": has_sd}


# ======================================================================== cell 4/5
# Locate / filter abnormal transcripts. Each isoform was labelled base-by-base
# against the known gene model (upstream annotation), giving sense / antisense /
# intergenic coverage fractions.
#
# Source of anti-sense transcription (Serrano & Lluch-Senar): the high-AT genome
# spawns many spurious promoters -> pervasive sRNAs. The ORF-input set below keeps
# both spurious-promoter and read-through antisense plus intergenic transcripts as
# the substrate for novel-ORF enumeration.
def build_orf_input(df: pd.DataFrame) -> pd.DataFrame:
    """Abnormal ORF-input set (notebook cell 5): core_abnormal_strict + exploratory_mixed."""
    df = df.copy()
    df["abnormal_frac"] = df["frac_antisense"] + df["frac_intergenic"]
    df["abnormal_bp"] = df["antisense_overlap_bp"] + df["intergenic_bp"]
    for c in ["n_reads", "isoform_len_bp", "frac_sense", "frac_antisense",
              "frac_intergenic", "abnormal_frac"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    core_abnormal = df[
        (df["class_main"].isin(["intergenic", "antisense"])) |
        (df["class_detail"].isin(["mixed_sense_intergenic", "mixed_antisense_intergenic"]))
    ].copy()
    core_abnormal_strict = core_abnormal[
        (core_abnormal["n_reads"] >= 10) &
        (core_abnormal["abnormal_frac"] >= 0.20) &
        (core_abnormal["abnormal_bp"] >= 30)
    ].copy()
    exploratory_mixed = df[
        (df["class_main"] == "mixed") &
        (df["class_detail"].isin(["mixed_sense_antisense", "mixed_complex"])) &
        (df["n_reads"] >= 10) &
        (df["abnormal_frac"] >= 0.30) &
        (df["abnormal_bp"] >= 60)
    ].copy()
    orf_input = pd.concat([core_abnormal_strict, exploratory_mixed],
                          ignore_index=True).drop_duplicates(subset=["isoform_id"])
    print(f"core_abnormal: {len(core_abnormal)}  ->  strict: {len(core_abnormal_strict)}")
    print(f"exploratory_mixed: {len(exploratory_mixed)}")
    print(f"Final ORF-input isoforms: {len(orf_input)}")
    return orf_input


# ======================================================================== cell 9
# Enumerate & evaluate novel ORFs via OSTIR (unified pipeline). For each abnormal
# isoform, OSTIR scans the full transcript for all start codons and predicts a TIR
# at each; then for every hit we:
#   1. find the downstream in-frame stop (table 4: TAA, TAG; TGA = Trp),
#   2. translate (genetic code 4),
#   3. annotate overlap with known genes (sense / antisense / intergenic),
#   4. classify the initiation mode from TIR + dG_mRNA:rRNA.
# Anti-SD = ACCTCCTTT (Syn1 16S 3' tail). Replaces the old manual ORF enumeration
# + custom SD scan entirely. ~2 min for ~840 transcripts.
def find_stop_codon(tx_seq: str, start_local0: int) -> Tuple[Optional[int], str]:
    seq = tx_seq.upper()
    j = start_local0 + 3
    while j + 3 <= len(seq):
        codon = seq[j:j + 3]
        if codon in STOP_CODONS:
            return j, codon
        j += 3
    return None, ""


def process_one_isoform(iso: pd.Series, fasta_file, genes_by_chrom) -> List[Dict]:
    isoform_id = iso["isoform_id"]
    chrom = str(iso["chrom"]); strand = str(iso["strand"])
    start0 = int(iso["start0"]); end0 = int(iso["end0"])

    tx_seq = fasta_file.fetch(chrom, start0, end0).upper()
    if strand == "-":
        tx_seq = revcomp(tx_seq)
    tx_len = len(tx_seq)
    if tx_len < 10:
        return []

    try:
        ostir_hits = run_ostir(tx_seq, aSD=ANTI_SD, threads=OSTIR_THREADS)
    except Exception as e:
        print(f"  OSTIR error for {isoform_id}: {e}")
        return []
    if not ostir_hits:
        return []

    orf_rows = []
    for k, hit in enumerate(ostir_hits):
        local_start0 = hit.get("start_position", 0) - 1   # 1-based -> 0-based
        if local_start0 < 0 or local_start0 + 3 > tx_len:
            continue
        tir = hit.get("expression", np.nan)
        start_codon = tx_seq[local_start0:local_start0 + 3]
        stop_pos, stop_codon = find_stop_codon(tx_seq, local_start0)

        aa_seq_full = ""
        if stop_pos is not None:
            orf_end_local0 = stop_pos + 3
            nt_seq = tx_seq[local_start0:orf_end_local0]
            aa_seq_full = translate_code4(nt_seq)
            aa_seq = aa_seq_full.rstrip("*")
            if aa_seq:
                aa_seq = "M" + aa_seq[1:]
            termination = "complete_stop"
        else:
            usable_end = local_start0 + ((tx_len - local_start0) // 3) * 3
            if usable_end <= local_start0:
                continue
            orf_end_local0 = usable_end
            nt_seq = tx_seq[local_start0:orf_end_local0]
            aa_seq = translate_code4(nt_seq)
            if aa_seq:
                aa_seq = "M" + aa_seq[1:]
            termination = "no_stop_within_isoform"

        aa_len = len(aa_seq)
        if aa_len < MIN_AA_LEN_MICRO:
            continue
        length_tier = "main" if aa_len >= MIN_AA_LEN_MAIN else "micro"

        g0, g1 = tx_local_to_genomic(start0, end0, strand, local_start0, orf_end_local0)
        if chrom in genes_by_chrom:
            ov = annotate_orf_overlap(chrom, strand, g0, g1, genes_by_chrom[chrom])
        else:
            orf_len_bp = g1 - g0
            ov = {
                "orf_bp_in_sense": 0, "orf_bp_in_antisense": 0, "orf_bp_in_intergenic": int(orf_len_bp),
                "orf_frac_in_sense": 0.0, "orf_frac_in_antisense": 0.0, "orf_frac_in_intergenic": 1.0,
                "orf_sense_gene_ids": "", "orf_sense_gene_names": "",
                "orf_antisense_gene_ids": "", "orf_antisense_gene_names": "",
                "orf_n_sense_genes": 0, "orf_n_antisense_genes": 0,
                "orf_novelty_class": "novel_intergenic_orf",
                "orf_abnormal_frac": 1.0, "orf_abnormal_bp": int(orf_len_bp),
            }

        dg_mrna_rrna = hit.get("dG_rRNA:mRNA", np.nan)
        init = classify_initiation_ostir(tir, dg_mrna_rrna, leaderless_dist=local_start0)

        orf_rows.append({
            "orf_id": f"{isoform_id}__ORF_{k + 1:03d}",
            "isoform_id": isoform_id, "chrom": chrom, "strand": strand,
            "isoform_start0": start0, "isoform_end0": end0,
            "isoform_len_bp": end0 - start0, "isoform_n_reads": int(iso["n_reads"]),
            "class_main": iso["class_main"], "class_detail": iso["class_detail"],
            "isoform_frac_sense": float(iso["frac_sense"]),
            "isoform_frac_antisense": float(iso["frac_antisense"]),
            "isoform_frac_intergenic": float(iso["frac_intergenic"]),
            "isoform_sense_overlap_bp": int(iso["sense_overlap_bp"]),
            "isoform_antisense_overlap_bp": int(iso["antisense_overlap_bp"]),
            "isoform_intergenic_bp": int(iso["intergenic_bp"]),
            "tx_seq": tx_seq, "tx_len_bp": tx_len,
            "orf_local_start0": local_start0, "orf_local_end0": orf_end_local0,
            "orf_genomic_start0": int(g0), "orf_genomic_end0": int(g1),
            "start_codon": start_codon, "stop_codon": stop_codon if stop_codon else "",
            "termination_status": termination, "orf_nt_len": len(nt_seq),
            "orf_aa_len": aa_len, "length_tier": length_tier,
            "orf_nt_seq": nt_seq, "orf_aa_seq": aa_seq,
            "orf_aa_seq_with_stop": aa_seq_full if termination == "complete_stop" else aa_seq,
            "ostir_TIR": tir,
            "ostir_dG_total": hit.get("dG_total", np.nan),
            "ostir_dG_mRNA_rRNA": dg_mrna_rrna,
            "ostir_dG_mRNA": hit.get("dG_mRNA", np.nan),
            "ostir_dG_spacing": hit.get("dG_spacing", np.nan),
            "ostir_dG_standby": hit.get("dG_standby", np.nan),
            "ostir_dG_start_codon": hit.get("dG_start_codon", np.nan),
            "ostir_RBS_distance_bp": hit.get("RBS_distance_bp", np.nan),
            "ostir_start_codon": hit.get("start_codon", ""),
            **init, **ov,
            "leaderless": local_start0 <= LEADERLESS_MAX_DIST,
            "leaderless_dist_from_5p_nt": local_start0,
            "is_novel_orf": ov["orf_novelty_class"] in [
                "novel_intergenic_orf", "antisense_orf", "mixed_boundary_orf", "mixed_complex_orf"],
            "low_complexity_peptide": bool(low_complexity_peptide(aa_seq)),
        })
    return orf_rows


def enumerate_orfs(orf_input, fasta, genes_by_chrom) -> pd.DataFrame:
    all_orf_rows = []
    t0 = time.time(); n_iso = len(orf_input)
    for idx, (_, iso) in enumerate(orf_input.iterrows()):
        if idx % 100 == 0:
            print(f"  [{time.time() - t0:6.0f}s]  isoform {idx + 1:>4d}/{n_iso}  "
                  f"ORFs so far: {len(all_orf_rows)}")
        all_orf_rows.extend(process_one_isoform(iso, fasta, genes_by_chrom))
    orf_annot = pd.DataFrame(all_orf_rows)
    print(f"\nFinished {n_iso} isoforms in {time.time() - t0:.0f}s. Total candidate ORFs: {len(orf_annot)}")
    return orf_annot


# ======================================================================== cell 12/13
# Filtering, scoring and deduplication.
#   synthesis_rate = isoform RNA reads x OSTIR TIR
#   filters: aa > 15, abnormal_frac > 0.2, not known_sense_orf
#   -> top 100 by synthesis rate -> collapse to unique peptide sequences.
# Benchmark: a canonical protein has ~20 M synthesis rate (~4k read depth x ~5k TIR).
# WARNING: exact-sequence dedup still does not rule out near-identical peptides.
def score_and_rank(orf_annot: pd.DataFrame):
    scored = orf_annot.copy()
    scored["ostir_TIR"] = pd.to_numeric(scored["ostir_TIR"], errors="coerce").fillna(0)
    scored["isoform_n_reads"] = pd.to_numeric(scored["isoform_n_reads"], errors="coerce").fillna(0)
    scored["synthesis_rate"] = scored["isoform_n_reads"] * scored["ostir_TIR"]

    filtered = scored[
        (scored["orf_aa_len"] > 15) &
        (scored["orf_abnormal_frac"] > 0.2) &
        (~scored["orf_novelty_class"].isin(["known_sense_orf"]))
    ].copy()
    print(f"After filtering (aa>15, abnormal>0.2, not known_sense): {len(filtered)}")

    strict_set = (filtered.sort_values("synthesis_rate", ascending=False)
                  .head(100).reset_index(drop=True).copy())
    strict_set["strict_rank"] = range(1, len(strict_set) + 1)
    strict_set["orf_label"] = [f"NOVEL_ORF_{i + 1:03d}" for i in range(len(strict_set))]
    print(f"Strict set (top 100): {len(strict_set)}; unique peptides: {strict_set['orf_aa_seq'].nunique()}")
    return scored, strict_set


def collapse_unique(strict_set: pd.DataFrame) -> pd.DataFrame:
    def join_unique(vals):
        vals = [str(v) for v in vals if str(v) not in {"", "nan", "None"}]
        return ",".join(sorted(set(vals)))

    unique_peptides = (
        strict_set.groupby("orf_aa_seq", dropna=False)
        .agg(
            n_entries=("orf_id", "size"),
            total_synthesis_rate=("synthesis_rate", "sum"),
            all_isoform_ids=("isoform_id", join_unique),
            n_isoforms=("isoform_id", lambda x: len(set(x))),
            all_orf_ids=("orf_id", join_unique),
            chrom=("chrom", "first"), strand=("strand", "first"),
            orf_genomic_start0=("orf_genomic_start0", "first"),
            orf_genomic_end0=("orf_genomic_end0", "first"),
            start_codon=("start_codon", "first"), stop_codon=("stop_codon", "first"),
            termination_status=("termination_status", "first"),
            orf_aa_len=("orf_aa_len", "first"), orf_nt_seq=("orf_nt_seq", "first"),
            orf_aa_seq_with_stop=("orf_aa_seq_with_stop", "first"),
            novelty_classes=("orf_novelty_class", join_unique),
            initiation_modes=("initiation_mode", join_unique),
        )
        .reset_index()
        .sort_values("total_synthesis_rate", ascending=False)
        .reset_index(drop=True)
    )
    unique_peptides["peptide_id"] = [f"NOVEL_PEP_{i + 1:03d}" for i in range(len(unique_peptides))]
    print(f"Unique peptide sequences: {len(unique_peptides)}")
    return unique_peptides


# ======================================================================== cell 15
# In-silico trypsin digestion to find MS-detectable ORFs. A novel ORF is
# MS-identifiable if it yields >=1 tryptic peptide (7-25 aa, <=1 missed cleavage)
# absent from the canonical proteome. (A collaborator's MS re-search of the
# augmented database later confirmed matches to ORFs from this pipeline.)
def trypsin_cleavage_sites(seq: str) -> list:
    seq = seq.upper()
    cuts = [0]
    for i in range(len(seq) - 1):
        if seq[i] in {"K", "R"} and seq[i + 1] != "P":
            cuts.append(i + 1)
    cuts.append(len(seq))
    return cuts


def digest_trypsin(seq: str, max_mc: int = 1) -> list:
    seq = seq.upper()
    cuts = trypsin_cleavage_sites(seq)
    peptides = []
    for i in range(len(cuts) - 1):
        for mc in range(max_mc + 1):
            j = i + mc + 1
            if j >= len(cuts):
                continue
            pep = seq[cuts[i]:cuts[j]]
            if len(pep) == 0 or "X" in pep:
                continue
            peptides.append({"peptide_seq": pep, "pep_len": len(pep), "missed_cleavages": mc})
    return peptides


def read_fasta(path: str):
    header, seq_chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]; seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def trypsin_analysis(unique_peptides: pd.DataFrame):
    canonical_peptide_set = set()
    for header, seq in read_fasta(CANONICAL_FASTA):
        for p in digest_trypsin(seq, max_mc=MAX_MISSED_CLEAVAGES):
            canonical_peptide_set.add(p["peptide_seq"])
    print(f"Canonical tryptic peptides: {len(canonical_peptide_set)}")

    rows = []
    novel_peptide_to_orfs = defaultdict(set)
    for _, r in unique_peptides.iterrows():
        pid = r["peptide_id"]
        for p in digest_trypsin(r["orf_aa_seq"], max_mc=MAX_MISSED_CLEAVAGES):
            pep = p["peptide_seq"]
            novel_peptide_to_orfs[pep].add(pid)
            rows.append({
                "novel_peptide_id": pid, "tryptic_peptide": pep, "pep_len": p["pep_len"],
                "missed_cleavages": p["missed_cleavages"],
                "detectable_len": MIN_PEPTIDE_LEN <= p["pep_len"] <= MAX_PEPTIDE_LEN,
                "in_canonical": pep in canonical_peptide_set,
            })
    tryptic_df = pd.DataFrame(rows)
    tryptic_df["unique_vs_canonical"] = ~tryptic_df["in_canonical"]
    tryptic_df["n_novel_orfs_sharing"] = tryptic_df["tryptic_peptide"].map(
        lambda x: len(novel_peptide_to_orfs[x]))
    tryptic_df["globally_unique"] = tryptic_df["unique_vs_canonical"] & (tryptic_df["n_novel_orfs_sharing"] == 1)
    tryptic_df["usable_for_ms"] = tryptic_df["detectable_len"] & tryptic_df["unique_vs_canonical"]

    orf_digest_summary = tryptic_df.groupby("novel_peptide_id").agg(
        n_tryptic_total=("tryptic_peptide", "size"),
        n_detectable=("detectable_len", "sum"),
        n_unique_vs_canonical=("usable_for_ms", "sum"),
        n_globally_unique=("globally_unique", lambda x: x.sum()),
        example_unique_peptides=("tryptic_peptide", lambda x: ",".join(
            x[tryptic_df.loc[x.index, "usable_for_ms"]].drop_duplicates().head(5).tolist())),
    ).reset_index()
    orf_digest_summary = unique_peptides.merge(
        orf_digest_summary, left_on="peptide_id", right_on="novel_peptide_id", how="left")
    ms_identifiable = orf_digest_summary[orf_digest_summary["n_unique_vs_canonical"] >= 1].copy()
    print(f"MS-identifiable (>=1 unique tryptic peptide): {len(ms_identifiable)} / {len(orf_digest_summary)}")
    return tryptic_df, orf_digest_summary, ms_identifiable


# ======================================================================== fasta writers
def write_strict_fasta(df, path):
    with open(path, "w") as fh:
        for _, r in df.iterrows():
            fh.write(f">{r['orf_label']}|orf_id={r['orf_id']}|aa={int(r['orf_aa_len'])}"
                     f"|term={r['termination_status']}|synth_rate={float(r['synthesis_rate']):.0f}"
                     f"|TIR={float(r['ostir_TIR']):.0f}|reads={int(r['isoform_n_reads'])}"
                     f"|novelty={r['orf_novelty_class']}|isoform={r['isoform_id']}\n")
            fh.write(str(r["orf_aa_seq"]) + "\n")


def write_unique_fasta(df, path):
    with open(path, "w") as fh:
        for _, r in df.iterrows():
            fh.write(f">{r['peptide_id']}|aa={int(r['orf_aa_len'])}|term={r['termination_status']}"
                     f"|total_synth_rate={float(r['total_synthesis_rate']):.0f}"
                     f"|n_isoforms={int(r['n_isoforms'])}|novelty={r['novelty_classes']}"
                     f"|isoforms={r['all_isoform_ids']}\n")
            fh.write(str(r["orf_aa_seq"]) + "\n")


def write_ms_fasta(df, path):
    with open(path, "w") as fh:
        for _, r in df.iterrows():
            fh.write(f">{r['peptide_id']}|aa={int(r['orf_aa_len'])}|term={r['termination_status']}"
                     f"|total_synth_rate={float(r['total_synthesis_rate']):.0f}"
                     f"|n_isoforms={int(r['n_isoforms'])}|n_unique_tryptic={int(r['n_unique_vs_canonical'])}"
                     f"|novelty={r['novelty_classes']}|isoforms={r['all_isoform_ids']}\n")
            fh.write(str(r["orf_aa_seq"]) + "\n")


# ======================================================================== main
def main():
    Path(OUT_FOLDER).mkdir(parents=True, exist_ok=True)
    Path(OUT_DIGEST_FOLDER).mkdir(parents=True, exist_ok=True)

    print("=" * 60, "\nLoad isoforms + build abnormal ORF-input set")
    df = pd.read_csv(ISOFORMS_TSV, sep="\t")
    print(f"isoform clusters: {len(df):,}  (from {ISOFORMS_TSV})")
    orf_input = build_orf_input(df)
    with pd.ExcelWriter(ISO_FILTERED_XLSX, engine="openpyxl") as w:
        orf_input.to_excel(w, sheet_name="orf_input", index=False)

    print("=" * 60, "\nOSTIR scan + ORF construction")
    fasta = pysam.FastaFile(FASTA_PATH)
    genes = load_genes_from_gff(GFF_PATH)
    genes_by_chrom = {c: g.reset_index(drop=True) for c, g in genes.groupby("chrom", sort=False)}
    print(f"Loaded {len(genes)} genes")
    orf_annot = enumerate_orfs(orf_input, fasta, genes_by_chrom)
    orf_annot.to_csv(f"{OUT_FOLDER}/candidate_orfs_ostir.tsv", sep="\t", index=False)

    print("=" * 60, "\nScore + rank -> top 100 -> unique peptides")
    scored, strict_set = score_and_rank(orf_annot)
    unique_peptides = collapse_unique(strict_set)
    scored.to_csv(f"{OUT_FOLDER}/candidate_orfs_scored.tsv", sep="\t", index=False)
    strict_set.to_csv(f"{OUT_FOLDER}/candidate_orfs_strict_top100.tsv", sep="\t", index=False)
    write_strict_fasta(strict_set, f"{OUT_FOLDER}/novel_orfs_strict.fasta")
    write_unique_fasta(unique_peptides, f"{OUT_FOLDER}/novel_peptides_unique.fasta")

    print("=" * 60, "\nIn-silico trypsin digestion")
    tryptic_df, orf_digest_summary, ms_identifiable = trypsin_analysis(unique_peptides)
    orf_digest_summary.to_csv(f"{OUT_DIGEST_FOLDER}/novel_peptide_digest_summary.tsv", sep="\t", index=False)
    write_ms_fasta(ms_identifiable, f"{OUT_DIGEST_FOLDER}/novel_peptides_ms_identifiable.fasta")
    with pd.ExcelWriter(f"{OUT_DIGEST_FOLDER}/novel_peptide_trypsin_analysis.xlsx", engine="openpyxl") as w:
        orf_digest_summary.to_excel(w, sheet_name="digest_summary", index=False)
        ms_identifiable.to_excel(w, sheet_name="ms_identifiable", index=False)
        tryptic_df.to_excel(w, sheet_name="all_tryptic_peptides", index=False)

    print("=" * 60, f"\nDONE. Outputs in {OUT_FOLDER}/ and {OUT_DIGEST_FOLDER}/")


if __name__ == "__main__":
    main()
