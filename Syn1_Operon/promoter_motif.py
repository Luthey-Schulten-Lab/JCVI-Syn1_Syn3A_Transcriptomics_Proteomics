"""
Syn1 promoter -10 box motif scanning (single source of truth).

Extracted verbatim from Operon_Annotation.py so the canonical-operon promoter
analysis AND the R4 novel-transcription promoter quantification use the SAME
algorithm (no drift between two copies).

For a transcription start site (tss0, chrom, strand) the upstream sequence is read
on the transcribed strand (reverse-complemented for the minus strand, circular-safe)
and matched to the bacterial -10 box by IUPAC consensus:
  - 6-mer  TANAAT     over the window [-12, -7] relative to the TSS
  - 9-mer  TNNTANAAT  over the window [-15, -7]
each allowed a register shift of +/- SHIFT_RANGE bp, choosing the offset with the
fewest IUPAC mismatches (tie-break: smaller |shift|, then sign). The per-TSS tier is
strong_9mer > core_6mer > no_minus10.

Public API:
    extract_tx_kmer(tss0, chrom, strand, rel_start, rel_end) -> str
    best_shift_for_consensus(tss0, chrom, strand, rel_start, rel_end, consensus, cons_re)
        -> (best_kmer, best_shift, match_bool, n_mismatches)
    scan_minus10(tss0, chrom, strand) -> dict  (convenience: both boxes + tier)
plus the IUPAC helpers, constants (CONS6/CONS9, REL6_*/REL9_*, SHIFT_RANGE,
CONS6_RE/CONS9_RE) and the genome dict GENOME.
"""

import os
import re

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

# -- Genome (path relative to this file, so import works from any cwd) --------
_HERE = os.path.dirname(os.path.abspath(__file__))
GENOME_FASTA = os.path.join(_HERE, "..", "Genomes_Input", "syn1_genome.fasta")

GENOME = {}
for _rec in SeqIO.parse(GENOME_FASTA, "fasta"):
    GENOME[_rec.id] = str(_rec.seq).upper()

# -- Settings -----------------------------------------------------------------
SHIFT_RANGE = 2
CONS6 = "TANAAT";    REL6_START, REL6_END = -12, -7   # 6-mer -10 box
CONS9 = "TNNTANAAT"; REL9_START, REL9_END = -15, -7   # 9-mer -10 box

# -- IUPAC helpers ------------------------------------------------------------
IUPAC = {
    "A":"A","C":"C","G":"G","T":"T","N":"[ACGT]",
    "R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]","K":"[GT]","M":"[AC]",
    "B":"[CGT]","D":"[AGT]","H":"[ACT]","V":"[ACG]"
}

def consensus_to_regex(consensus):
    pat = "".join(IUPAC.get(ch.upper(), ch.upper()) for ch in consensus)
    return re.compile("^" + pat + "$")

def mismatches_iupac(seq, consensus):
    """Count positions where seq violates IUPAC consensus."""
    mm = 0
    for s, c in zip(seq.upper(), consensus.upper()):
        if c not in IUPAC:
            mm += (s != c)
        else:
            if not re.fullmatch(IUPAC[c], s):
                mm += 1
    return int(mm)

CONS6_RE = consensus_to_regex(CONS6)
CONS9_RE = consensus_to_regex(CONS9)

# -- Circular-safe fetch ------------------------------------------------------
def circular_slice(seq, s0, e0):
    L = len(seq)
    span = e0 - s0
    if L == 0 or span <= 0:
        return ""
    s = s0 % L
    if s + span <= L:
        return seq[s:s+span]
    first = seq[s:]
    remain = span - (L - s)
    fc, tail = divmod(remain, L)
    return first + (seq * fc) + (seq[:tail] if tail else "")

def extract_tx_kmer(tss0, chrom, strand, rel_start, rel_end):
    """
    Extract a window in transcription orientation.
      + strand: g = tss0 + rel
      - strand: g = tss0 - rel  (then reverse-complemented)
    """
    tss0 = int(tss0)
    if strand == "+":
        g0 = tss0 + int(rel_start)
        g1 = tss0 + int(rel_end) + 1
    else:
        g_a = tss0 - int(rel_start)
        g_b = tss0 - int(rel_end)
        g0, g1 = min(g_a, g_b), max(g_a, g_b) + 1
    seq = circular_slice(GENOME[chrom], g0, g1)
    if len(seq) != (g1 - g0):
        return ""
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq

def best_shift_for_consensus(tss0, chrom, strand, rel_start, rel_end, consensus, cons_re):
    """
    Return (best_kmer, best_shift, matches, n_mismatches).
    Selects by minimum IUPAC mismatches; tie-breaks by |shift| then sign.
    """
    K = rel_end - rel_start + 1
    shifts = sorted(range(-SHIFT_RANGE, SHIFT_RANGE+1), key=lambda x: (abs(x), x))
    best = None
    for sh in shifts:
        kmer = extract_tx_kmer(tss0, chrom, strand, rel_start+sh, rel_end+sh)
        if len(kmer) != K or not re.fullmatch(r"[ACGT]+", kmer):
            mm, match = 10**9, False
        else:
            mm = mismatches_iupac(kmer, consensus)
            match = bool(cons_re.match(kmer))
        cand = (mm, abs(sh), sh, kmer, match)
        if best is None or cand < best:
            best = cand
    mm, _, sh, kmer, match = best
    if mm >= 10**8:
        return ("", 0, False, np.nan)
    return (kmer, int(sh), bool(match), int(mm))


def scan_minus10(tss0, chrom, strand):
    """Convenience: scan both the 6-mer and 9-mer -10 boxes at one TSS and tier it.
    Returns a dict with both best k-mers + shift/match/mismatches and `motif_tier`
    (strong_9mer > core_6mer > no_minus10)."""
    k6, s6, m6, mm6 = best_shift_for_consensus(tss0, chrom, strand, REL6_START, REL6_END, CONS6, CONS6_RE)
    k9, s9, m9, mm9 = best_shift_for_consensus(tss0, chrom, strand, REL9_START, REL9_END, CONS9, CONS9_RE)
    tier = "strong_9mer" if m9 else ("core_6mer" if m6 else "no_minus10")
    return {
        "tss0": int(tss0), "chrom": chrom, "strand": strand,
        "minus10_6mer": k6, "shift6": s6, "match6": m6, "mm6": mm6,
        "minus10_9mer": k9, "shift9": s9, "match9": m9, "mm9": mm9,
        "motif_tier": tier,
    }
