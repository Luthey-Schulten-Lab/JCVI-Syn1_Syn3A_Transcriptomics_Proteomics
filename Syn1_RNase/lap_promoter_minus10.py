#!/usr/bin/env python
"""
Promoter -10 signatures at the two Fig. 2 example operons (lap/0154 and 0178).

PART 1 -- lap/0154 (Fig. 2 panel c): which 5' end is the real TSS?
PART 2 -- 0178/neopullulanase (Fig. 2 panel b): does its shared 5' end have a -10 box?

BACKGROUND
----------
Figure 2 panel c uses the leucyl-aminopeptidase operon OP_00078 (lap/MMSYN1_0154,
+ strand) as the example of 5'-biased RNA processing: many isoforms share a 3'
terminator and erode from the 5' end. Two candidate 5' ends compete for the label
"transcription start site":

  A. the OPERON START called by the segmentation (operons.candidate_blocks.tsv,
     `tss` column) -- the most upstream 5' end supported by full-length reads;
  B. the DEPTH-JUMP coordinate -- an internal position where PacBio coverage rises
     several-fold and where a large block of isoform 5' ends piles up.

A depth jump alone is ambiguous: it can mark a genuine internal promoter, or it can
mark the 5' boundary of a processed/degraded population (an endonucleolytic cut, or
5'->3' exonucleolytic stalling) sitting inside a transcript that actually starts
further upstream. The promoter -10 box discriminates between the two: a genuine
sigma-factor TSS carries a recognizable -10 hexamer ~10 bp upstream; a processing
boundary has no reason to.

ALGORITHM
---------
1. COORDINATES (both derived from data, nothing hard-coded):
   A = OP_00078 `tss` from the canonical 459-operon map.
   B = the largest positive step in the PacBio + strand depth bedGraph inside the
       operon span, cross-checked against the read-weighted isoform 5'-end
       histogram (the modal 5' end excluding A).
2. -10 SCAN: promoter_motif.scan_minus10() -- the SAME scanner used for the 127
   canonical operons in Syn1_Operon/Operon_Annotation.py (imported, not copied):
   6-mer TANAAT over [-12,-7] and 9-mer TNNTANAAT over [-15,-7], register shift
   +/-2 bp, best offset = fewest IUPAC mismatches. Tier strong_9mer > core_6mer >
   no_minus10.
3. CONTEXT: the -45..+10 upstream sequence of each candidate is printed with the
   matched box marked, so the call can be checked by eye.
4. NULL MODEL: because the Syn1 genome is AT-rich, TANAAT arises by chance. Two
   baselines are reported -- (i) N random + strand positions scanned identically,
   (ii) the 127 canonical operon TSSs (Operon_Annotation output) -- so "matches
   TANAAT" can be read against how often that happens by accident.
5. LOCAL SPECIFICITY: every position in a window around lap is scanned, to show
   how many positions in the neighbourhood would also pass.

REJECTED LINE OF ARGUMENT (recorded so it is not re-tried)
----------------------------------------------------------
B sits 21 nt upstream of the lap start codon, i.e. near the upstream edge of a
canonical 30S initiation-complex footprint, which invites the story "5'->3' trimming
stalled at a bound ribosome". Tested transcriptome-wide (all 20,885 isoforms with
>=10 reads, dominant 5' end per gene relative to its start codon) this does NOT
discriminate: 24.8% of the 509 testable genes have their dominant 5' end at -25..-15
but a statistically indistinguishable 22.2% have it at -15..-5, inside the footprint.
There is a broad concentration of 5' ends just upstream of start codons, with no
sharp peak at the -20 footprint edge, and short 5'UTRs from genuine promoters land in
exactly the same window. Position relative to the RBS therefore cannot separate
"processed end" from "TSS with a short 5'UTR" -- only the -10 box does.

PART 2 -- 0178 (OP_00099, - strand), the 3'-erosion example
-----------------------------------------------------------
Unlike lap there is no competing coordinate: 2,901 of 3,121 operon-contained reads
(93%) share ONE 5' base and coverage collapses 2861x -> 4x across it, so the only
question is whether that sharp end is a promoter. The same scanner is applied, then
the surrounding 66 bp is swept for the nearest perfect box and the read coverage at
that alternative position is checked.

OUTCOME (see the .txt)
----------------------
PART 1 (lap):
A = 197657 (operon start, 102 reads / 4 clusters): core_6mer -- a PERFECT TANAAT
    hexamer (TAAAAT, 0 mismatches) at register 0, i.e. canonical -10 spacing. The
    9-mer extension misses by 1 (GTGTAAAAT), so the tier is core_6mer, which is the
    tier of 87/127 canonical Syn1 operons.
B = 197721 (depth jump 169x -> 1143x, 6.8-fold; 1003 reads / 26 clusters, the modal
    5' end): no_minus10 -- best 6-mer TATTAT and best 9-mer TGTTATTAT both miss by 1
    at every register within +/-2 bp.
Background: the same scanner matches a 6-mer at 2.4% of random positions and at
68.5% of the 127 canonical operon TSSs (29x enrichment), and the 10 positive
positions among the 185 spanning A..B collapse to 2 box loci, one sitting exactly on
A -- B is in neither.
Conclusion: transcription initiates at 197657, and the 6.8-fold step at 197721 is a
post-transcriptional boundary inside that primary transcript, not a second promoter.

PART 2 (0178) -- no consensus match, but NOT evidence against a promoter:
The dominant 5' end (234511, 93% of reads, 44 nt 5'UTR) does not match the consensus:
best 6-mer TAATAT misses by 1 and best 9-mer AATTAATAT by 2, at every register within
+/-2 bp. This is a scanner non-assignment, not a demonstrated absence of a promoter --
a single deviation from TANAAT is ordinary in real bacterial -10 elements, and 40/127
(31%) of the canonical Syn1 operons are unassigned by the same scanner for the same
reason. The nearest perfect hexamer (TACAAT, genomic 234538-234543) would serve a TSS
~21 nt further upstream, but coverage there is 4x against 2861x at 234511, so it is
not the promoter for the observed transcripts.
Read as: 0178's 5' end is sharp, dominant and intergenic; the -10 scan simply adds
nothing either way. The affirmative -10 argument is available for lap (Part 1) and is
not claimed here.

MINUS-STRAND CONVENTION NOTE
----------------------------
Operon_Annotation.py's tss_pos() returns end0 for minus-strand operons, but end0 is
the half-open exclusive bound, so the true 5' base is end0-1 (a 1 bp offset). Checked
here explicitly: both conventions return the same k-mer and the same tier for 0178
(the +/-2 register search absorbs the offset), so the published classification is
unaffected -- but the 1 bp should not be propagated into any coordinate that is
reported as a TSS.

Run in base env (pandas + numpy + biopython):
  /home/enguang/anaconda3/bin/python Syn1_RNase/lap_promoter_minus10.py
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "Syn1_Operon"))

# the canonical -10 scanner (single source of truth, shared with Operon_Annotation.py)
from promoter_motif import (                                   # noqa: E402
    GENOME, scan_minus10, best_shift_for_consensus, circular_slice,
    CONS6, CONS9, CONS6_RE, CONS9_RE,
    REL6_START, REL6_END, REL9_START, REL9_END, SHIFT_RANGE,
)

ISOF   = f"{ROOT}/Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
OPS    = f"{ROOT}/Syn1_Operon/operons.candidate_blocks.tsv"
DEPTH  = f"{ROOT}/Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.{{}}.bedGraph"
CANON  = f"{ROOT}/Syn1_Operon/annotation/canonical/promoter_minus10_classification.tsv"
CHROM  = "CP002027.1"
OPERON = "OP_00078"
# GFF is 1-based inclusive (197743-199098) -> 0-based half-open start = 197742
GENE   = dict(locus="MMSYN1_0154", name="lap", gstart=197743, gend=199098, gstart_0based=197742)

TXT = f"{ROOT}/Syn1_RNase/lap_promoter_minus10.txt"

N_NULL = 20000
NULL_SEED = 0

_log = []
def say(msg=""):
    print(msg)
    _log.append(msg)


# ---------------------------------------------------------------- loaders
def load_depth(lo, hi, strand="+"):
    """Strand-specific PacBio depth over [lo, hi)."""
    d = np.zeros(hi - lo)
    for ln in open(DEPTH.format("plus" if strand == "+" else "minus")):
        f = ln.split()
        if f[0] != CHROM:
            continue
        a, b = int(f[1]), int(f[2])
        if b < lo or a > hi:
            continue
        d[max(a, lo) - lo:min(b, hi) - lo] = float(f[3])
    return d


def fmt_scan(sc):
    m6 = "match" if sc["match6"] else f"{sc['mm6']} mm"
    m9 = "match" if sc["match9"] else f"{sc['mm9']} mm"
    return (f"6-mer {sc['minus10_6mer']} (shift {sc['shift6']:+d}, {m6}) | "
            f"9-mer {sc['minus10_9mer']} (shift {sc['shift9']:+d}, {m9}) | {sc['motif_tier']}")


def context_string(tss0, strand, back=45, fwd=10):
    """Upstream..downstream sequence in transcription orientation, TSS at index `back`."""
    from Bio.Seq import Seq
    if strand == "+":
        s = circular_slice(GENOME[CHROM], tss0 - back, tss0 + fwd + 1)
    else:
        s = circular_slice(GENOME[CHROM], tss0 - fwd, tss0 + back + 1)
        s = str(Seq(s).reverse_complement())
    return s


def marker_line(seq_len, back, rel_start, rel_end, shift, char):
    """Ruler line marking a [rel_start+shift, rel_end+shift] window under the context."""
    line = [" "] * seq_len
    for rel in range(rel_start + shift, rel_end + shift + 1):
        i = back + rel
        if 0 <= i < seq_len:
            line[i] = char
    return "".join(line)


# ---------------------------------------------------------------- step 1: coordinates
ops = pd.read_csv(OPS, sep="\t")
op = ops[ops.operon_id == OPERON].iloc[0]
op_lo, op_hi, strand = int(op.start0), int(op.end0), str(op.strand)
tss_operon = int(op.tss)

iso = pd.read_csv(ISOF, sep="\t")
s = iso[(iso.strand == strand) & (iso.start0 >= op_lo) & (iso.end0 <= op_hi)].copy()
p5hist = s.groupby("start0").n_reads.agg(n_reads="sum", n_clusters="count").sort_values("n_reads", ascending=False)

PAD = 100
d = load_depth(op_lo - PAD, op_hi + PAD)
step = np.diff(d)
# largest positive step strictly inside the operon (exclude the operon start itself)
inside = np.array([op_lo < (op_lo - PAD) + i + 1 < op_hi for i in range(len(step))])
i_jump = int(np.argmax(np.where(inside, step, -np.inf)))
tss_jump = (op_lo - PAD) + i_jump + 1
d_before, d_after = d[i_jump], d[i_jump + 1]

say("=" * 78)
say("lap/0154 (Fig. 2 panel c) -- promoter -10 box at the two candidate 5' ends")
say("=" * 78)
say()
say(f"Operon {OPERON}  {CHROM}  {strand} strand  span {op_lo}-{op_hi} ({op_hi - op_lo} bp)")
say(f"Gene   {GENE['locus']} ({GENE['name']})  {GENE['gstart']}-{GENE['gend']}")
say(f"Isoforms contained in the operon span: {len(s)}  ({int(s.n_reads.sum())} reads)")
say()
say("-- Step 1: the two candidate coordinates ----------------------------------")
say()
say(f"A  operon start (segmentation TSS)   = {tss_operon}")
say(f"   PacBio 5' ends exactly here       = {int(p5hist.loc[tss_operon, 'n_reads'])} reads "
    f"in {int(p5hist.loc[tss_operon, 'n_clusters'])} isoform clusters")
say(f"   position relative to lap start    = {tss_operon - GENE['gstart']} bp (intergenic, upstream of the ORF)")
say()
say(f"B  depth-jump coordinate             = {tss_jump}")
say(f"   PacBio depth step                 = {d_before:.0f}x -> {d_after:.0f}x "
    f"(+{d_after - d_before:.0f}, {d_after / max(d_before, 1):.1f}-fold), the largest inside the operon")
say(f"   PacBio 5' ends exactly here       = {int(p5hist.loc[tss_jump, 'n_reads'])} reads "
    f"in {int(p5hist.loc[tss_jump, 'n_clusters'])} isoform clusters")
say(f"   position relative to lap start    = {tss_jump - GENE['gstart']} bp (still intergenic, "
    f"{tss_jump - tss_operon} bp downstream of A)")
say()
say("   Read-weighted isoform 5'-end histogram, top 8 positions:")
say(f"      {'pos':>8}  {'reads':>6}  {'clusters':>8}   note")
for pos, row in p5hist.head(8).iterrows():
    note = "<- A (operon start)" if pos == tss_operon else ("<- B (depth jump)" if pos == tss_jump else "")
    say(f"      {pos:>8}  {int(row.n_reads):>6}  {int(row.n_clusters):>8}   {note}")
say()
say(f"   The modal 5' end is B ({int(p5hist.iloc[0].n_reads)} reads), i.e. the depth jump and the")
say("   dominant isoform 5' end are the same coordinate -- the two derivations agree.")
say()

# ---------------------------------------------------------------- step 2: -10 scan
say("-- Step 2: -10 box at each candidate (promoter_motif.scan_minus10) --------")
say()
say(f"   scanner: 6-mer {CONS6} over [{REL6_START},{REL6_END}], "
    f"9-mer {CONS9} over [{REL9_START},{REL9_END}], register shift +/-{SHIFT_RANGE} bp")
say("   (imported from Syn1_Operon/promoter_motif.py -- identical to the 127-operon analysis)")
say()

cands = [("A", "operon start", tss_operon), ("B", "depth jump", tss_jump)]
for tag, label, pos in cands:
    sc = scan_minus10(pos, CHROM, strand)
    say(f"   {tag} {pos} ({label})")
    say(f"      {fmt_scan(sc)}")
say()

# sequence context with the boxes marked
say("   Sequence context (transcription orientation, TSS = position 0, marked '^'):")
say()
BACK, FWD = 45, 10
for tag, label, pos in cands:
    sc = scan_minus10(pos, CHROM, strand)
    ctx = context_string(pos, strand, BACK, FWD)
    ruler = [" "] * len(ctx); ruler[BACK] = "^"
    m6 = marker_line(len(ctx), BACK, REL6_START, REL6_END, sc["shift6"], "6")
    m9 = marker_line(len(ctx), BACK, REL9_START, REL9_END, sc["shift9"], "9")
    say(f"      {tag} {pos} ({label})   tier = {sc['motif_tier']}")
    say(f"         {ctx}")
    say(f"         {''.join(ruler)}")
    say(f"         {m9}   9-mer window -> {sc['minus10_9mer']} "
        f"({'MATCH' if sc['match9'] else str(sc['mm9']) + ' mismatch(es)'})")
    say(f"         {m6}   6-mer window -> {sc['minus10_6mer']} "
        f"({'MATCH' if sc['match6'] else str(sc['mm6']) + ' mismatch(es)'})")
    say()

# ---------------------------------------------------------------- step 3: null model
say("-- Step 3: how often does this happen by chance? -------------------------")
say()
rng = np.random.default_rng(NULL_SEED)
L = len(GENOME[CHROM])
null_pos = rng.integers(0, L, size=N_NULL)
n6 = n9 = 0
for p in null_pos:
    _, _, m6ok, _ = best_shift_for_consensus(int(p), CHROM, "+", REL6_START, REL6_END, CONS6, CONS6_RE)
    _, _, m9ok, _ = best_shift_for_consensus(int(p), CHROM, "+", REL9_START, REL9_END, CONS9, CONS9_RE)
    n6 += bool(m6ok); n9 += bool(m9ok)
say(f"   Random + strand positions (n={N_NULL}, same scanner, same +/-{SHIFT_RANGE} shift):")
say(f"      6-mer {CONS6}     matched {n6:>5} / {N_NULL}  = {100 * n6 / N_NULL:.1f}%")
say(f"      9-mer {CONS9}  matched {n9:>5} / {N_NULL}  = {100 * n9 / N_NULL:.1f}%")

can = pd.read_csv(CANON, sep="\t")
say()
say(f"   The 127 canonical operon TSSs (Operon_Annotation.py, same scanner):")
say(f"      6-mer matched {int(can.minus10_6mer_match.sum()):>5} / {len(can)}  "
    f"= {100 * can.minus10_6mer_match.mean():.1f}%")
say(f"      9-mer matched {int(can.minus10_9mer_match.sum()):>5} / {len(can)}  "
    f"= {100 * can.minus10_9mer_match.mean():.1f}%")
say()
say(f"   Enrichment at real operon TSSs over background: "
    f"6-mer {can.minus10_6mer_match.mean() / (n6 / N_NULL):.0f}x, "
    f"9-mer {can.minus10_9mer_match.mean() / (n9 / N_NULL):.0f}x.")
say("   A 6-mer hit is therefore informative, not automatic, in this AT-rich genome:")
say(f"   it is the tier candidate A reaches and it occurs at only {100 * n6 / N_NULL:.1f}% of positions by chance.")
say()

# ---------------------------------------------------------------- step 4: local specificity
say("-- Step 4: local specificity around lap ----------------------------------")
say()
WIN_LO, WIN_HI = tss_operon - 60, tss_jump + 60
hits6, hits9 = [], []
for p in range(WIN_LO, WIN_HI + 1):
    _, _, m6ok, _ = best_shift_for_consensus(p, CHROM, "+", REL6_START, REL6_END, CONS6, CONS6_RE)
    _, _, m9ok, _ = best_shift_for_consensus(p, CHROM, "+", REL9_START, REL9_END, CONS9, CONS9_RE)
    if m6ok: hits6.append(p)
    if m9ok: hits9.append(p)
def runs_of(positions):
    """Collapse consecutive positions into intervals (the +/-2 register tolerance makes
    one real box light up a short run of neighbouring TSS positions)."""
    if not positions:
        return []
    out, start = [], positions[0]
    for a, b in zip(positions, positions[1:] + [None]):
        if b is None or b != a + 1:
            out.append((start, a)); start = b
    return out

r6, r9 = runs_of(hits6), runs_of(hits9)
say(f"   Scanning every position in {WIN_LO}-{WIN_HI} ({WIN_HI - WIN_LO + 1} bp):")
say(f"      positions with a 6-mer -10 box: {len(hits6)} "
    f"({100 * len(hits6) / (WIN_HI - WIN_LO + 1):.0f}%), in {len(r6)} interval(s): " +
    ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in r6))
say(f"      positions with a 9-mer -10 box: {len(hits9)} "
    f"({100 * len(hits9) / (WIN_HI - WIN_LO + 1):.0f}%)" +
    ("" if not r9 else ", in intervals: " +
     ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in r9)))
say()
say(f"      A ({tss_operon}) carries a 6-mer box: {tss_operon in hits6}")
say(f"      B ({tss_jump}) carries a 6-mer box: {tss_jump in hits6}")
say()
say(f"      The {len(hits6)} positive positions collapse to {len(r6)} distinct box loci, one of which")
say("      sits exactly at A (the run spans A +/- the 2 bp register tolerance). B falls in")
say("      neither, and is 60+ bp from the nearest one.")
say()

# ---------------------------------------------------------------- verdict
say("-- Verdict ---------------------------------------------------------------")
say()
scA = scan_minus10(tss_operon, CHROM, strand)
scB = scan_minus10(tss_jump, CHROM, strand)
say(f"   A  {tss_operon}  operon start   {scA['motif_tier']:>12}   "
    f"6-mer {scA['minus10_6mer']} ({scA['mm6']} mm), 9-mer {scA['minus10_9mer']} ({scA['mm9']} mm)")
say(f"   B  {tss_jump}  depth jump     {scB['motif_tier']:>12}   "
    f"6-mer {scB['minus10_6mer']} ({scB['mm6']} mm), 9-mer {scB['minus10_9mer']} ({scB['mm9']} mm)")
say()
if scA["match6"] and not scB["match6"]:
    say(f"   The operon start carries a perfect {CONS6} hexamer ({scA['minus10_6mer']}) at register 0,")
    say("   i.e. canonical -10 spacing; the depth-jump coordinate carries no -10 box at any")
    say(f"   register within +/-{SHIFT_RANGE} bp, and is one of the 6-mer-negative positions in its own")
    say(f"   neighbourhood. The {d_after / max(d_before, 1):.1f}-fold coverage step at {tss_jump} therefore marks a")
    say(f"   post-transcriptional boundary inside a transcript initiated at {tss_operon}, not a")
    say("   second promoter. (What kind of boundary -- endonucleolytic cut vs 5'->3' trimming")
    say("   stall -- is NOT settled by these data; see the rejected RBS argument in the")
    say("   docstring.)")
    say()
    say(f"   Caveat: A reaches the core_6mer tier, not strong_9mer -- the 9-mer extension misses")
    say(f"   by 1 ({scA['minus10_9mer']}). That is the majority situation among real Syn1 operons")
    say("   (87/127 match the 6-mer, only 52/127 the 9-mer), so it does not weaken the call.")
else:
    say("   NOTE: the expected pattern did not reproduce -- inspect the scans above.")
say()

# ================================================================ PART 2: 0178
say()
say("=" * 78)
say("PART 2 -- 0178 / neopullulanase (Fig. 2 panel b, the 3'-erosion example)")
say("=" * 78)
say()

OPERON_B = "OP_00099"
GENE_B = dict(locus="MMSYN1_0178", name="0178", gstart=232672, gend=234468)

opb = ops[ops.operon_id == OPERON_B].iloc[0]
b_lo, b_hi, b_strand = int(opb.start0), int(opb.end0), str(opb.strand)
# minus strand: the annotation stores tss = end0, but the true 5' base is end0-1
tss_annot = int(opb.tss)
tss_true = b_hi - 1

sb = iso[(iso.strand == b_strand) & (iso.start0 >= b_lo) & (iso.end0 <= b_hi)].copy()
sb["p5"] = sb.end0 - 1
hb = sb.groupby("p5").n_reads.agg(n_reads="sum", n_clusters="count").sort_values("n_reads", ascending=False)
top5p = int(hb.index[0])
db = load_depth(b_lo - 150, b_hi + 150, b_strand)
def depth_at(p):
    return db[p - (b_lo - 150)]

say(f"Operon {OPERON_B}  {CHROM}  {b_strand} strand  span {b_lo}-{b_hi} ({b_hi - b_lo} bp)")
say(f"Gene   {GENE_B['locus']}  {GENE_B['gstart']}-{GENE_B['gend']} (neopullulanase)")
say(f"Isoforms contained in the operon span: {len(sb)}  ({int(sb.n_reads.sum())} reads)")
say()
say("-- Step 1: is there a competing 5' end here at all? -----------------------")
say()
say(f"   dominant 5' base                   = {top5p}  "
    f"({int(hb.iloc[0].n_reads)} reads in {int(hb.iloc[0].n_clusters)} clusters, "
    f"{100 * hb.iloc[0].n_reads / sb.n_reads.sum():.0f}% of operon reads)")
say(f"   coverage across it                 = {depth_at(top5p):.0f}x -> {depth_at(top5p + 1):.0f}x")
# GFF end is 1-based inclusive, so the last gene base is gend-1 in 0-based coords;
# the UTR runs gend .. top5p inclusive
say(f"   5'UTR to the ORF                   = {top5p - GENE_B['gend'] + 1} nt")
say()
say("   Next strongest 5' ends:")
for pos, row in hb.iloc[1:5].iterrows():
    say(f"      {pos}  {int(row.n_reads):>5} reads  ({int(row.n_clusters)} clusters)")
say()
say("   One sharp 5' end carries almost all the reads, so -- unlike lap -- there is no")
say("   second coordinate to arbitrate. The only question is whether it is a promoter.")
say()

say("-- Step 2: -10 box at the 0178 5' end -------------------------------------")
say()
say(f"   Minus-strand convention check (annotation tss = end0 = {tss_annot}, "
    f"true 5' base = end0-1 = {tss_true}):")
for p, lab in ((tss_annot, "end0 (Operon_Annotation convention)"), (tss_true, "end0-1 (true 5' base)")):
    sc = scan_minus10(p, CHROM, b_strand)
    say(f"      {p}  {lab}")
    say(f"         {fmt_scan(sc)}")
say(f"   -> same k-mer and same tier either way; the +/-{SHIFT_RANGE} register search absorbs the 1 bp.")
say()

scb = scan_minus10(tss_true, CHROM, b_strand)
ctx = context_string(tss_true, b_strand, BACK, FWD)
ruler = [" "] * len(ctx); ruler[BACK] = "^"
say(f"   Sequence context   tier = {scb['motif_tier']}")
say(f"      {ctx}")
say(f"      {''.join(ruler)}")
say(f"      {marker_line(len(ctx), BACK, REL9_START, REL9_END, scb['shift9'], '9')}   "
    f"9-mer -> {scb['minus10_9mer']} ({'MATCH' if scb['match9'] else str(scb['mm9']) + ' mismatch(es)'})")
say(f"      {marker_line(len(ctx), BACK, REL6_START, REL6_END, scb['shift6'], '6')}   "
    f"6-mer -> {scb['minus10_6mer']} ({'MATCH' if scb['match6'] else str(scb['mm6']) + ' mismatch(es)'})")
say()

say("-- Step 3: is there a usable -10 box anywhere nearby? ---------------------")
say()
WB_LO, WB_HI = tss_true - 30, tss_true + 36
hitsb = []
for p in range(WB_LO, WB_HI + 1):
    k, sh, m, _ = best_shift_for_consensus(p, CHROM, b_strand, REL6_START, REL6_END, CONS6, CONS6_RE)
    if m:
        hitsb.append((p, k))
say(f"   Sweeping {WB_LO}-{WB_HI} ({WB_HI - WB_LO + 1} bp) for a 6-mer box:")
if hitsb:
    rb = runs_of([p for p, _ in hitsb])
    kmer = hitsb[0][1]
    say(f"      {len(hitsb)} positive positions in {len(rb)} locus/loci: " +
        ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in rb) + f"  (box {kmer})")
    centre = (rb[0][0] + rb[0][1]) // 2
    say(f"      That box would serve a TSS near {centre}, i.e. "
        f"{abs(centre - tss_true)} nt further upstream than the observed 5' end.")
    say(f"      PacBio coverage there = {depth_at(centre):.0f}x, against "
        f"{depth_at(top5p):.0f}x at {top5p}.")
    say()
    say("      So that box is essentially unused: almost no transcript actually starts")
    say("      there, and it cannot be the promoter for the observed 0178 transcripts.")
else:
    say("      none")
say()

say("-- Verdict (0178) --------------------------------------------------------")
say()
say(f"   {tss_true}  0178 5' end   {scb['motif_tier']:>12}   "
    f"6-mer {scb['minus10_6mer']} ({scb['mm6']} mm), 9-mer {scb['minus10_9mer']} ({scb['mm9']} mm)")
say()
n_no10 = int((~can.minus10_6mer_match).sum())
say(f"   The consensus is not matched: TAATAT is 1 mismatch from {CONS6} at every register")
say(f"   within +/-{SHIFT_RANGE} bp. Treat this as a scanner non-assignment, NOT as evidence that the")
say("   5' end is not a promoter -- a single deviation from consensus is ordinary in real")
say(f"   bacterial -10 elements, and {n_no10}/{len(can)} ({100 * n_no10 / len(can):.0f}%) of the canonical Syn1 operons fall in the")
say("   same bin. The scan is simply uninformative here, in either direction.")
say()
say("   What does stand on its own: the 0178 5' end is sharp, dominant (93% of reads) and")
say("   intergenic. The affirmative -10 argument is available for lap (Part 1); it is not")
say("   claimed for 0178.")
say()

say(f"[written] {TXT}")
with open(TXT, "w") as fh:
    fh.write("\n".join(_log) + "\n")
