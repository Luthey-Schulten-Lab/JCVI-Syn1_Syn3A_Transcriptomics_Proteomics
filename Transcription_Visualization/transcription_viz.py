#!/usr/bin/env python
"""Transcription-activity browser: give it any gene, get every RNA-seq library over that region.

Generalisation of Genome_Reduction/dcw_operon.py. That script needed a hand-written config per
region (anchor, window, gene list, deletion bands, TSS/terminator marks); here `build_region()`
derives all of it from the pipeline tables, so `plot_region("ptsG")` is the whole interface.

Figure grammar (unchanged from dcw_operon.py)
--------------------------------------------
  gene arrows, Syn1 then Syn3A   ->  read stacks + depth, one pair per library
All panels share a 5'->3' RELATIVE axis: rel = (pos - anchor) on a plus-strand region,
(anchor - pos) on a minus-strand one, so the transcript always reads left to right. Each
organism is anchored at ITS OWN copy of the anchor gene - that is what makes the two
comparable across a deletion, and it keeps the Syn3A gene track CONSECUTIVE (you see what
Syn3A actually has next to the gene, including the other side of a scar) rather than
projecting Syn3A onto Syn1 coordinates and punching gaps in it.

Where a deletion falls inside the window, everything past the FIRST junction is greyed on the
Syn3A panels: the signal there is real Syn3A signal, but from a different locus, so it must
never be compared base-for-base with the Syn1 panels above it.

How a region is derived from one gene name
------------------------------------------
  window    the Syn1 operon block(s) covering the gene (operons.candidate_blocks.tsv), spanned
            end to end, + FLANK bp each side. Genes in no block fall back to gene +- NO_OP_PAD.
  anchor    5' end of the 5'-most sense gene inside that span (Syn1); the Syn3A anchor is set
            so the first window gene that Syn3A still has lands at the same rel.
  genes     every gene/pseudogene overlapping the window, per organism, from the two GFFs.
  deletions syn1_deleted_regions.bed intersected with the window -> red "absent from Syn3A"
            bands on the Syn1 gene row + the grey Syn3A mask.
  marks     operon TSS (green) and TransTermHP intrinsic terminators (red) inside the window.

Depth / normalisation conventions, inherited unchanged
------------------------------------------------------
  * SENSE depth only (the region's own strand); antisense is not drawn.
  * scale="abs" (the default): raw per-base depth, each library on its own 0..in-window-max
    axis, with that library's both-strand genome-mean depth drawn as the reference line. No
    normaliser, so a panel cannot understate how deep a library really is here.
    "linear"/"log" divide by that genome mean instead and share one axis.
  * The three Syn1 Illumina replicates are combined as Syn1_Syn3A_Transcriptomics/Gene_TPM.py does
    (each on its own library-size scale first, then weights .25/.25/.5). Under "abs" there is
    nothing to normalise, so they are pooled (summed) instead.
  * Illumina gets no read-stack panel: 51-nt paired-end fragments are not transcript-length
    observations.

Usage
-----
    import transcription_viz as tv
    fig, txt, reg = tv.plot_region("ptsG")        # or "0779", "MMSYN1_0779", "JCVISYN3A_0779"
    tv.check_inputs()                             # verify every external file this needs

Writes <gene>_transcription.{pdf,png} + _stats.txt under Transcription_Visualization/regions/.
Needs samtools on PATH. External inputs live in other folders but are all declared in the
PATHS block below - edit that block (or set MINCELL_ROOT) to move this anywhere.
"""
import os
import re
import functools
import subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrow, FancyArrowPatch

# ============================== paths (the only block that knows about other folders) ========
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("MINCELL_ROOT", os.path.abspath(os.path.join(HERE, "..")))
D = lambda p: os.path.join(ROOT, p)
OUTDIR = os.path.join(HERE, "regions")

GFF_S1 = "Genomes_Input/syn1.genes.gff3"
GFF_S3 = "Genomes_Input/syn3a_genome.gff3"
OPERONS = "Syn1_Operon/operons.candidate_blocks.tsv"
TRANSTERM = "Syn1_Operon/syn1_TransTermHP.txt"
PROMOTER_MOD = "Syn1_Operon"                       # holds promoter_motif.py
DELBED = "Genome_Reduction/aln/raw/syn1_deleted_regions.bed"
PAIRED = "Genome_Reduction/Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv"
PROT_S1 = "Syn1_Syn3A_Corr_RNA_Proteins/syn1_omics.xlsx"                    # sheet syn1_genes
PROT_S3 = "Syn1_Syn3A_Proteomics/syn3A_proteome.tsv"
ISO_TSV = "Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"

_PB = "Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ"
_ON1 = "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn1.ONT.rep1"
_ON2 = "Syn1_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn1.ONT.rep2"
_ILD = "Syn1_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph"
_S3ONT = "Syn3A_Transcriptomics/ONT/ONT_Processing/depth_bedgraph/syn3A.ONT.rep1"
_S3ILL = "Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph/syn3A_rep1"
BAM = {
    "Syn1 PacBio": "Syn1_Transcriptomics/PacBio/PacBio_Processing/syn1.PacBio.FLNC.sorted.HQ.bam",
    "Syn1 ONT 1": "Syn1_Transcriptomics/ONT/ONT_Processing/syn1.ONT.rep1.sorted.bam",
    "Syn1 ONT 2": "Syn1_Transcriptomics/ONT/ONT_Processing/syn1.ONT.rep2.sorted.bam",
    "Syn3A ONT": "Syn3A_Transcriptomics/ONT/ONT_Processing/syn3A.ONT.rep1.sorted.bam",
}
S1_CH, S1L = "CP002027.1", 1_078_809
S3_CH, S3L = "CP016816.2", 543_379
CHROM = {"s1": S1_CH, "s3": S3_CH}
GLEN = {"s1": S1L, "s3": S3L}

# ============================== tunables =====================================================
FLANK = 500          # bp added each side of the operon span
NO_OP_PAD = 2000     # window padding for a gene that sits in no operon block
N_STACK = 190        # a stack panel is ~39 pt tall, so ~190 rows is the most that stay
                     # individually resolvable; above that the panel reads as read density
SYN1_COL, SYN3A_COL = "#3182bd", "#c0392b"
C_GENE, C_BREAK, GREEN, RED = "#cfcfcf", "#111111", "#2ca02c", "#c0392b"

mpl.rcParams.update({
    "font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})


def show(R, width=680):
    """Display a saved figure inside a notebook, scaled to `width` PIXELS.

    The figure itself must stay 7 in wide (Nature 2-column) and is born at that size, so the
    exported PDF/PNG are untouched; this only shrinks the <img> the browser draws, which is
    what makes a ~7 x 7.5 in panel fit on one notebook screen. Pass R from plot_region().
    """
    from IPython.display import Image, display
    display(Image(filename=R["files"]["png"], width=width))


def check_inputs():
    """Report which external inputs are present. Call once after moving the folder."""
    need = [GFF_S1, GFF_S3, OPERONS, TRANSTERM, DELBED, PAIRED, ISO_TSV, PROT_S1, PROT_S3] \
        + list(BAM.values())
    need += [f"{p}.{s}.bedGraph" for p in (_PB, _ON1, _ON2, _S3ONT, _S3ILL,
                                           f"{_ILD}/SRR35996296", f"{_ILD}/SRR35996297",
                                           f"{_ILD}/SRR35996298") for s in ("plus", "minus")]
    miss = [p for p in need if not os.path.isfile(D(p))]
    print(f"ROOT = {ROOT}\n{len(need) - len(miss)}/{len(need)} inputs found")
    for p in miss:
        print("  MISSING:", p)
    if subprocess.run(["which", "samtools"], capture_output=True).returncode:
        print("  MISSING: samtools on PATH (needed for read stacks)")
    return not miss


# ============================== annotation loaders (cached) ==================================
@functools.lru_cache(maxsize=None)
def load_genes(org):
    """DataFrame of gene/pseudogene records: num, start0, end0, strand, name."""
    path, pref = (GFF_S1, "MMSYN1_") if org == "s1" else (GFF_S3, "JCVISYN3A_")
    rows = []
    for line in open(D(path)):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] not in ("gene", "pseudogene"):
            continue
        a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        lt = a.get("locus_tag", "")
        if not lt.startswith(pref):
            continue
        nm = a.get("Name") or a.get("gene") or ""
        if nm.startswith(("JCVISYN3A_", "MMSYN1_")):
            nm = ""
        rows.append((lt[-4:], int(p[3]) - 1, int(p[4]), p[6], nm, a.get("product", "")))
    return pd.DataFrame(rows, columns=["num", "start0", "end0", "strand", "name", "product"])


@functools.lru_cache(maxsize=None)
def load_operons():
    d = pd.read_csv(D(OPERONS), sep="\t")
    return d[d.chrom == S1_CH][["operon_id", "strand", "start0", "end0", "tss", "sense_gene_loci"]]


@functools.lru_cache(maxsize=None)
def load_deletions():
    d = pd.read_csv(D(DELBED), sep="\t", comment="#", header=None,
                    names=["chrom", "s", "e", "len"])          # file carries a #chrom header
    return d[d.chrom == S1_CH].astype({"s": int, "e": int})


@functools.lru_cache(maxsize=None)
def load_syn3a_proteome():
    """The curated Syn3A proteome table, indexed by 4-digit locus number."""
    # snake_case table -> the Title-Case names this module and the notebook use
    b = pd.read_csv(D(PROT_S3), sep="\t").rename(columns={
        "locus_tag": "Locus Tag", "gene_name": "Gene Name",
        "gene_product": "Gene Product", "essentiality": "Essentiality",
        "localization": "Localization",
        "primary_function": "Primary Function",
        "secondary_function": "Secondary Function",
        "tertiary_function": "Tertiary Function",
        "exp_ptn_cnt_2019": "Exp. Ptn Cnt 2019",
        "exp_ptn_cnt_2026": "Exp. Ptn Cnt 2026",
        "mRNA_copies_per_cell": "mRNA Copies/Cell"})
    return b.assign(num=b["Locus Tag"].astype(str).str[-4:]).set_index("num")


@functools.lru_cache(maxsize=None)
def load_protein_copies():
    """Measured protein copies per cell, keyed by 4-digit locus number, for both organisms."""
    a = pd.read_excel(D(PROT_S1), sheet_name="syn1_genes")
    a = a.assign(num=a.locusTag.astype(str).str[-4:]).set_index("num").protein_copy_number
    b = load_syn3a_proteome()
    return a, b["Exp. Ptn Cnt 2019"], b["Exp. Ptn Cnt 2026"]


def resolve_gene(q):
    """'0779' | 'MMSYN1_0779' | 'JCVISYN3A_0779' | 'ptsG'  ->  4-digit locus number."""
    q = str(q).strip()
    m = re.fullmatch(r"(?:MMSYN1_|JCVISYN3A_)?(\d{4})", q)
    if m:
        return m.group(1)
    for org in ("s1", "s3"):                       # fall back to a gene NAME, Syn1 first
        g = load_genes(org)
        hit = g[g.name.str.lower() == q.lower()]
        if len(hit):
            return hit.iloc[0].num
    raise KeyError(f"gene {q!r} not found in either annotation")


# ============================== region construction =========================================
def build_region(gene, flank=FLANK):
    """Everything the plotter needs, derived from one gene. See module docstring."""
    num = resolve_gene(gene)
    g1 = load_genes("s1")
    g3 = load_genes("s3")
    row = g1[g1.num == num]
    if not len(row):
        raise KeyError(f"{num} has no Syn1 gene record (Syn3A-only genes are not supported)")
    row = row.iloc[0]
    strand = row.strand
    sgn = 1 if strand == "+" else -1

    # -- core span: the operon block(s) covering the gene, spanned; else the gene itself --
    ops = load_operons()
    cov = ops[(ops.strand == strand) & (ops.start0 < row.end0) & (ops.end0 > row.start0)]
    if len(cov):
        span0, span1 = int(cov.start0.min()), int(cov.end0.max())
    else:
        span0, span1 = int(row.start0), int(row.end0)
    # An operon boundary often cuts through the first or last gene it carries (the block is
    # defined by transcript ends, not ORF ends). Grow the core to whole gene bodies so no gene
    # is ever drawn truncated, and so the flank stays true flank.
    ov = g1[(g1.strand == strand) & (g1.start0 < span1) & (g1.end0 > span0)]
    if len(ov):
        span0, span1 = min(span0, int(ov.start0.min())), max(span1, int(ov.end0.max()))
    flank = flank if len(cov) else NO_OP_PAD

    # -- anchor: 5' end of the 5'-most sense gene inside the span --
    ins = g1[(g1.strand == strand) & (g1.start0 < span1) & (g1.end0 > span0)]
    five = ins.start0 if strand == "+" else ins.end0
    inside = ins[(five >= span0) & (five <= span1)]        # a gene may only overhang the span;
    src = inside if len(inside) else ins                   # its 5' end would sit outside it
    five = src.start0 if strand == "+" else src.end0
    i5 = five.idxmin() if strand == "+" else five.idxmax()
    anchor1, anchor_num = int(five[i5]), src.loc[i5, "num"]

    rel = lambda pos, anc=None: sgn * (pos - (anchor1 if anc is None else anc))
    r_span = sorted([rel(span0), rel(span1)])
    rel_lo = int(min(r_span[0], 0) - flank)                # rel 0 (the anchor) is always in frame
    rel_hi = int(max(r_span[1], 0) + flank)

    # -- genes in the window, per organism (Syn3A read consecutively off its own coordinates) --
    win = sorted([anchor1 + sgn * rel_lo, anchor1 + sgn * rel_hi])
    w1 = g1[(g1.start0 < win[1]) & (g1.end0 > win[0])].copy()

    # -- Syn3A anchor: line the first window gene Syn3A still has up at the same rel --
    anchor3, s3_note = None, ""
    order = w1.assign(r=[sorted([rel(s), rel(e)])[0] for s, e in zip(w1.start0, w1.end0)])
    order = pd.concat([order[order.num == anchor_num], order.sort_values("r")])   # anchor first
    for _, g in order.iterrows():
        h = g3[(g3.num == g.num) & (g3.strand == g.strand)]
        if len(h):
            h = h.iloc[0]
            c1 = g.start0 if strand == "+" else g.end0
            c3 = h.start0 if strand == "+" else h.end0
            anchor3 = int(c3 - sgn * rel(c1))
            if g.num != anchor_num:
                s3_note = (f"Syn3A anchored on {g.num}: the anchor gene {anchor_num} is absent "
                           f"from Syn3A, so rel 0 is set from the nearest retained gene.")
            break
    s3nm = dict(zip(g3.num, g3.name))
    w1["label"] = [n or s3nm.get(k, "") for n, k in zip(w1.name, w1.num)]
    # Product descriptions: the curated Syn3A proteome wins over the GFF wherever the locus
    # matches - the GFF still carries "conserved hypothetical protein" for genes that have
    # since been characterised (0779 -> ptsG), and the workbook is the maintained source.
    try:
        p3 = load_syn3a_proteome()["Gene Product"]
        w1["product"] = [p3.get(k, q) if pd.notna(p3.get(k, np.nan)) else q
                         for k, q in zip(w1.num, w1["product"])]
    except Exception:                                                 # noqa - annotation only
        pass
    w3 = pd.DataFrame(columns=list(g3.columns) + ["label"])
    if anchor3 is not None:
        w3win = sorted([anchor3 + sgn * rel_lo, anchor3 + sgn * rel_hi])
        w3 = g3[(g3.start0 < w3win[1]) & (g3.end0 > w3win[0])].copy()
        w3["label"] = w3.name

    # -- deletions in the window -> red bands + the grey Syn3A mask at the first junction --
    bands = []
    for _, d in load_deletions().iterrows():
        if d.s < win[1] and d.e > win[0]:
            a, b = sorted([rel(max(int(d.s), win[0])), rel(min(int(d.e), win[1]))])
            bands.append((a, b))
    bands.sort()
    brk = bands[0][0] if bands else None

    # -- marks: operon TSSs (green) + TransTermHP terminators (red) --
    op_rel = sorted([rel(span0), rel(span1)])            # transcripts must overlap THIS, not the flank
    op_marks = [(*sorted([rel(int(o.start0)), rel(int(o.end0))]), o.operon_id)
                for _, o in cov.iterrows()]
    marks = [(rel(int(t)), GREEN, f"TSS {int(t):,}") for t in cov.tss.dropna()
             if win[0] <= int(t) <= win[1]]
    terms = []
    for ln in open(D(TRANSTERM)):
        if ln.strip().startswith("TERM"):
            f = ln.split()
            a, b = int(f[2]), int(f[4])
            if min(a, b) < win[1] and max(a, b) > win[0]:
                terms.append(ln.strip())
                marks.append((rel(min(a, b)), RED, f"TERM {f[1]}"))

    return dict(num=num, gene=row, strand=strand, sgn=sgn, anchor1=anchor1, anchor3=anchor3,
                anchor_num=anchor_num,
                rel_lo=rel_lo, rel_hi=rel_hi, win=win, w1=w1, w3=w3, bands=bands, brk=brk,
                core=(span0, span1), op_rel=op_rel, op_marks=op_marks,
                marks=marks, terms=terms, s3_note=s3_note,
                operons=list(cov.operon_id), tss=[int(t) for t in cov.tss.dropna()],
                label=w1.set_index("num").label.get(num, "") or "hyp",
                product=w1.set_index("num")["product"].get(num, row["product"]),
                title=f"{w1.set_index('num').label.get(num, '') or row['product'] or 'hyp'}"
                      f" / MMSYN1_{num}"
                      f"  [{', '.join(cov.operon_id) if len(cov) else 'no operon block'}], "
                      f"{strand} strand")


# ============================== depth / read loaders ========================================
def _libs(scale):
    ill = [(f"{_ILD}/SRR35996296", .25), (f"{_ILD}/SRR35996297", .25), (f"{_ILD}/SRR35996298", .5)]
    if scale == "abs":
        ill = [(p, 1.0) for p, _ in ill]           # nothing to normalise -> pool the replicates
    return [("Syn1 PacBio", "s1", [(_PB, 1.)]), ("Syn1 ONT 1", "s1", [(_ON1, 1.)]),
            ("Syn1 ONT 2", "s1", [(_ON2, 1.)]), ("Syn1 Illumina", "s1", ill),
            ("Syn3A ONT", "s3", [(_S3ONT, 1.)]), ("Syn3A Illumina", "s3", [(_S3ILL, 1.)])]


@functools.lru_cache(maxsize=None)
def _bg_array(pfx, st, ch, glen):
    """Whole-chromosome per-base depth for one bedGraph, as a float32 array.

    Cached, and worth it: the alternative (an awk pass per gene per library) re-reads ~200 MB
    of bedGraph for every region drawn, which dominates a 74-gene batch. Built with a
    difference array + cumsum, exact because bedGraph intervals never overlap. ~4 MB per
    array, ~60 MB for the full set of libraries and strands.
    """
    d = pd.read_csv(D(f"{pfx}.{st}.bedGraph"), sep="\t", header=None, comment="#",
                    names=["chrom", "s", "e", "v"])
    d = d[d.chrom == ch]
    acc = np.zeros(glen + 1, dtype=np.float64)
    np.add.at(acc, d.s.to_numpy(), d.v.to_numpy())
    np.add.at(acc, np.minimum(d.e.to_numpy(), glen), -d.v.to_numpy())
    return np.cumsum(acc)[:glen].astype(np.float32)


def _genome_mean(pfx, ch, glen):
    """Both-strand genome-mean depth of one library."""
    return float(sum(_bg_array(pfx, st, ch, glen).sum() for st in ("plus", "minus"))) / glen


def _depth(pfx, st, ch, lo, hi):
    return _bg_array(pfx, st, ch, GLEN["s1" if ch == S1_CH else "s3"])[lo:hi].astype(float)


def rel_depth(pwl, org, anchor, R, scale):
    """Sense depth on the ascending rel axis + the library's genome mean."""
    ch, glen, sgn = CHROM[org], GLEN[org], R["sgn"]
    lo, hi = sorted([anchor + sgn * R["rel_lo"], anchor + sgn * R["rel_hi"]])
    sense_bg = "plus" if R["strand"] == "+" else "minus"
    out, m = np.zeros(hi - lo), 0.0
    for pfx, wt in pwl:
        mi = _genome_mean(pfx, ch, glen)
        s = _depth(pfx, sense_bg, ch, lo, hi)
        out += wt * (s if scale == "abs" else s / mi)
        m += wt * mi
    r = sgn * (np.arange(lo, hi) - anchor)
    o = np.argsort(r)
    return r[o], out[o], m


def read_spans(bam, org, anchor, R):
    """(rel5', rel3') per primary sense-strand alignment in the window, sorted by 5' end.

    PacBio FLNC and ONT direct-RNA reads are both 5'->3' sense, so BAM strand == transcript
    strand. -F 0x904 drops unmapped/secondary/supplementary; span comes from the CIGAR with
    soft clips excluded, so one row is one of the alignments the depth track is built from.
    """
    ch, sgn = CHROM[org], R["sgn"]
    lo, hi = sorted([anchor + sgn * R["rel_lo"], anchor + sgn * R["rel_hi"]])
    flt = ["-F", "0x914"] if R["strand"] == "+" else ["-F", "0x904", "-f", "0x10"]
    out = subprocess.run(["samtools", "view"] + flt + [D(bam), f"{ch}:{lo + 1}-{hi}"],
                         capture_output=True, text=True).stdout
    sp = []
    for ln in out.splitlines():
        f = ln.split("\t")
        s0 = int(f[3]) - 1
        e0 = s0 + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]) if op in "MDN=X")
        a, b = sorted([sgn * (s0 - anchor), sgn * (e0 - anchor)])
        if b > R["op_rel"][0] and a < R["op_rel"][1]:     # must touch the operon, not just flank
            sp.append((a, b))
    return np.array(sorted(sp), dtype=float).reshape(-1, 2)


# ============================== plotting ====================================================
def _subsample(sp, rows):
    """Even subsample of a 5'-sorted transcript list to `rows` rows ('all' keeps every one)."""
    if rows == "all" or len(sp) <= int(rows):
        return sp
    return sp[np.unique(np.linspace(0, len(sp) - 1, int(rows)).round().astype(int))]


def _pack(intervals, gap=60):
    """Greedy row packing so overlapping isoform arrows never collide."""
    ends, out = [], [0] * len(intervals)
    for i in sorted(range(len(intervals)), key=lambda i: intervals[i][0]):
        a, b = intervals[i]
        for r, e in enumerate(ends):
            if a >= e + gap:
                ends[r], out[i] = b, r
                break
        else:
            ends.append(b); out[i] = len(ends) - 1
    return out


def _gene_row(ax, gd, anchor, R, org_label, col, bands=(), ops=()):
    sgn, lo, hi = R["sgn"], R["rel_lo"], R["rel_hi"]
    for i, (_, g) in enumerate(gd.sort_values("start0").iterrows()):
        a, b = sorted([sgn * (g.start0 - anchor), sgn * (g.end0 - anchor)])
        fwd = (g.strand == R["strand"])            # points right if co-oriented with the region
        x0, dx = (a, b - a) if fwd else (b, a - b)
        ax.add_patch(FancyArrow(x0, 0, dx, 0, width=.34, head_width=.34,
                                head_length=min(90, abs(dx)), length_includes_head=True,
                                color=C_GENE, ec="#8a8a8a", lw=.4))
        va_, vb_ = max(a, lo), min(b, hi)          # label the VISIBLE part of a clipped gene
        if vb_ <= va_:
            continue
        lab = f"{g['label']}/{g.num}" if g["label"] else g.num
        ax.text((va_ + vb_) / 2, .48 if i % 2 == 0 else -.48, lab, ha="center", va="center",
                fontsize=5, clip_on=True)
        if (vb_ - va_) > .015 * (hi - lo):         # absolute strand, inside the arrow
            ax.text((va_ + vb_) / 2, 0, "+" if g.strand == "+" else "−", ha="center",
                    va="center", fontsize=5, color="#4d4d4d", clip_on=True)
    for a, b in bands:
        ax.axvspan(a, b, facecolor="#e8736a", alpha=.20, lw=0, zorder=0)
    top = 1.15
    lv = _pack([(a, b) for a, b, _ in ops], gap=.06 * (hi - lo))   # overlapping blocks stack
    for (a, b, lab), k in zip(ops, lv):        # |-----| bracket = the operon block's extent
        va_, vb_ = max(a, lo), min(b, hi)
        if vb_ <= va_:
            continue
        y = .95 + .52 * k
        ax.plot([va_, vb_], [y, y], color="#333", lw=.7, clip_on=True)
        for x in (a, b):                       # end caps only where the block really ends
            if lo <= x <= hi:
                ax.plot([x, x], [y - .17, y + .17], color="#333", lw=.7, clip_on=True)
        ax.text((va_ + vb_) / 2, y + .27, lab, ha="center", va="center", fontsize=5, color="#333",
                zorder=7, bbox=dict(fc="white", ec="none", pad=.5))   # sits above the guide lines
        top = max(top, y + .5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-1.15, top)
    ax.axis("off")
    ax.text(lo - .012 * (hi - lo), 0, org_label, ha="right", va="center", fontsize=6, color=col)


def plot_region(gene, scale="abs", rows="all", mode="stack", flank=FLANK, flank_depth=.25,
                outdir=OUTDIR, save=True):
    """Draw every library over the region around `gene`. Returns (fig, stats_text, region).

    scale  "abs" raw per-base depth (default) | "linear"/"log" x that library's genome mean
    rows   "all" one row per transcript (default) | int N: even subsample to N rows. Only
           transcripts overlapping the OPERON are ever eligible - one lying entirely in the
           flank is not evidence about this operon, so it is dropped before subsampling.
    mode   "stack" read stacks + depth (default) | "pacbio" one PacBio isoform-cluster panel
    flank_depth  how the depth OUTSIDE the operon is rendered: an alpha (default 0.25 = faded)
           or None to drop it entirely. Either way the y-axis is scaled to the operon only -
           a transcript sitting in the flank contributes nothing inside the operon, so the
           only thing it can do to this panel is squash it (0776 next to ptsG). Flank depth is
           still drawn faded rather than blanked so it can never be misread as zero coverage.
    """
    if rows != "all" and (not str(rows).lstrip("-").isdigit() or int(rows) < 1):
        raise ValueError('rows must be "all" or a positive integer')
    R = build_region(gene, flank)
    if mode == "pacbio":
        scale = "log"
    lo, hi, sgn = R["rel_lo"], R["rel_hi"], R["sgn"]
    anc = {"s1": R["anchor1"], "s3": R["anchor3"]}
    libs = [l for l in _libs(scale) if anc[l[1]] is not None]

    dep = {n: rel_depth(p, o, anc[o], R, scale) for n, o, p in libs}
    spans, iso = {}, None
    if mode == "stack":
        spans = {n: read_spans(BAM[n], o, anc[o], R) for n, o, _ in libs if n in BAM}
    else:
        d = pd.read_csv(D(ISO_TSV), sep="\t")
        d = d[(d.chrom == S1_CH) & (d.strand == R["strand"]) & (d.n_reads >= 10) &
              (d.start0 < R["core"][1]) & (d.end0 > R["core"][0])].copy()   # operon, not flank
        rr = [sorted([sgn * (s - R["anchor1"]), sgn * (e - R["anchor1"])])
              for s, e in zip(d.start0, d.end0)]
        d["r5"], d["r3"] = [x[0] for x in rr], [x[1] for x in rr]
        d = d.sort_values("n_reads", ascending=False)
        iso = d if rows == "all" else d.head(int(rows))

    # ---- row plan: two gene rows, then per library [spacer] [stack] [depth] ----
    ROWS = [("gene", "s1"), ("gene", "s3")] if anc["s3"] is not None else [("gene", "s1")]
    HGT = [1.5] * len(ROWS)
    if mode == "pacbio":
        ROWS.append(("iso", None)); HGT.append(1.3)
    for i, (n, o, _p) in enumerate(libs):
        if i:
            ROWS.append(("gap", None)); HGT.append(.45)
        if mode == "stack" and n in BAM:
            ROWS.append(("stack", (n, o))); HGT.append(1.25)
        ROWS.append(("depth", (n, o))); HGT.append(1.0)
    fig, axes = plt.subplots(len(ROWS), 1, figsize=(7, 1.0 + .40 * sum(HGT)), sharex=True,
                             height_ratios=HGT, constrained_layout=True)
    fig.get_layout_engine().set(h_pad=.02, hspace=.005)
    for ax, (k, _t) in zip(axes, ROWS):
        if k == "gap":
            ax.axis("off")
        if k in ("gene", "iso", "stack", "gap"):
            ax.tick_params(axis="x", length=0)

    _gene_row(axes[0], R["w1"], R["anchor1"], R, "Syn1", SYN1_COL, R["bands"], R["op_marks"])
    if anc["s3"] is not None:
        _gene_row(axes[1], R["w3"], R["anchor3"], R, "Syn3A", SYN3A_COL)

    for ax, (k, t) in zip(axes, ROWS):
        if k == "iso":
            ints = list(zip(iso.r5, iso.r3))
            pr = _pack(ints)
            for (xa, xb), ri, nr in zip(ints, pr, iso.n_reads):
                lw = float(np.clip(.3 + .7 * np.log10(max(1, nr)), .5, 2.8))
                ax.add_patch(FancyArrowPatch((xa, ri), (xb, ri), arrowstyle="-|>",
                                             mutation_scale=5, lw=lw, color=SYN1_COL,
                                             alpha=.8, shrinkA=0, shrinkB=0))
            ax.set_ylim(-1, (max(pr) if pr else 0) + 1.5)
            ax.set_ylabel(f"PacBio isoforms\n(n={len(iso)})", fontsize=5.5, rotation=0,
                          ha="right", va="center", color=SYN1_COL)
            ax.set_yticks([]); ax.spines[["top", "right", "left"]].set_visible(False)
        elif k == "stack":
            n, o = t
            sp = spans[n]
            shown = _subsample(sp, rows)
            col = SYN1_COL if o == "s1" else SYN3A_COL
            lw_, al_ = (.08, 1.0) if rows == "all" else (.35, .75)
            ax.add_collection(LineCollection([[(a, i), (b, i)] for i, (a, b) in enumerate(shown)],
                                             colors=col, linewidths=lw_, alpha=al_, zorder=2))
            ax.set_ylim(len(shown) + 1, -2)
            lab = f"{n} reads\nn={len(sp):,}"
            if len(shown) < len(sp):
                lab += f" ({len(shown)} drawn)"
            ax.set_ylabel(lab, fontsize=5.5, rotation=0, ha="right", va="center", color=col)
            ax.set_yticks([]); ax.spines[:].set_visible(False)
        elif k == "depth":
            n, o = t
            xg, s, m = dep[n]
            col = SYN1_COL if o == "s1" else SYN3A_COL
            # Split the trace at the operon boundary. The axis is set from the OPERON part
            # only: a transcript lying in the flank adds nothing inside the operon, so all it
            # can do to this panel is squash it (0776 drawn next to ptsG). The flank part is
            # still drawn, faded, so a quiet flank is never confused with an undrawn one.
            inop = (xg >= R["op_rel"][0]) & (xg <= R["op_rel"][1])
            s_in = np.where(inop, s, np.nan)
            s_out = np.where(inop, np.nan, s)
            fill = lambda y, yy, a: ax.fill_between(xg, y, yy, color=col, lw=0, alpha=a, zorder=2)
            if scale == "abs":
                top = float(max(np.nanmax(s_in) if inop.any() else s.max(), 1.)) * 1.08
                fill(0, np.clip(s_in, 0, top), 1.0)
                if flank_depth:
                    fill(0, np.clip(s_out, 0, top), flank_depth)
                ax.set_ylim(0, top)
                ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=3, integer=True))
                if m <= top:
                    ax.axhline(m, color="#555", lw=.6, ls=(0, (3, 2)), zorder=4)
                    txt_ = f"genome avg {m:,.0f}"
                else:
                    txt_ = f"genome avg {m:,.0f} (above panel, max here {s.max():,.0f})"
                ax.annotate(txt_, xy=(hi, min(m, top)), xytext=(-2, 1.5), textcoords="offset points",
                            ha="right", va="bottom", fontsize=5, color="#555", zorder=6,
                            bbox=dict(fc="white", ec="none", pad=.5, alpha=.8))
            elif scale == "log":
                ax.set_yscale("log")
                top = 50
                fill(.05, np.clip(s_in, .05, None), 1.0)
                if flank_depth:
                    fill(.05, np.clip(s_out, .05, None), flank_depth)
                ax.set_ylim(.05, top)
                ax.set_yticks([.1, 1, 10]); ax.set_yticklabels(["0.1×", "1×", "10×"])
            else:
                top = 4
                fill(0, np.clip(s_in, 0, top), 1.0)
                if flank_depth:
                    fill(0, np.clip(s_out, 0, top), flank_depth)
                ax.set_ylim(0, 4.1)
                ax.set_yticks([0, 1, 2, 4]); ax.set_yticklabels(["0", "1×", "2×", "4×"])
                ax.axhline(1, color="#555", lw=.5, ls=(0, (3, 3)), zorder=4)
            # a flank peak taller than the operon axis is named, so clipping is never silent
            if inop.any() and np.nanmax(np.where(inop, -np.inf, s)) > top:
                i = int(np.nanargmax(np.where(inop, -np.inf, s)))
                ax.annotate(f"flank peak {s[i]:,.0f}", xy=(np.clip(xg[i], lo, hi), top),
                            xytext=(0, -7), textcoords="offset points", ha="center", va="top",
                            fontsize=5, color=col, zorder=7,
                            bbox=dict(fc="white", ec="none", pad=.4, alpha=.85))
            # bases with genuinely ZERO coverage, so "low" is never misread as "absent"
            y0, y1 = ax.get_ylim()
            ax.fill_between(xg, y0, y0 * 1.7 if scale == "log" else y0 + .055 * (y1 - y0),
                            where=(s <= 0), color="#111", lw=0, zorder=5)
            ax.set_ylabel(n if scale == "abs" else f"{n}\n1× = {m:,.0f}", fontsize=6,
                          rotation=0, ha="right", va="center", color=col)
            ax.set_xlim(lo, hi)
            ax.spines[["top", "right"]].set_visible(False)

    # ---- guide lines + the Syn3A non-syntenic mask ----
    for ax, (k, t) in zip(axes, ROWS):
        if k == "gap":
            continue
        for xr, cc, _ in R["marks"]:
            ax.axvline(xr, color=cc, lw=.6, ls=(0, (2, 2)), alpha=.85, zorder=3)
        if R["brk"] is not None:
            ax.axvline(R["brk"], color=C_BREAK, lw=.8, ls=(0, (4, 2)), alpha=.9, zorder=6)
            is3 = (k == "gene" and t == "s3") or (k in ("depth", "stack") and t[1] == "s3")
            if is3:
                ax.axvspan(R["brk"], hi, facecolor="#9e9e9e", alpha=.22, lw=0, zorder=1)
    if R["brk"] is not None:
        axes[0].annotate("absent from Syn3A", xy=((R["brk"] + hi) / 2, -.95), ha="center",
                         va="top", fontsize=5, color="#b03a2e")
        if anc["s3"] is not None:
            axes[1].annotate("different sequence (other side of the scar)",
                             xy=((R["brk"] + hi) / 2, 1.0), ha="center", va="top",
                             fontsize=5, color="#555")
    axes[0].set_title(R["title"], fontsize=7, loc="left")
    axes[-1].set_xlabel(f"Relative genome position from the 5′ end of {R['anchor_num']} (nt)",
                        fontsize=7)
    axes[-1].ticklabel_format(axis="x", style="plain")

    txt = _stats(R, dep, spans, iso, scale, rows, mode)
    if save:
        os.makedirs(outdir, exist_ok=True)
        nm = re.sub(r"\W+", "_", R["label"])[:24].strip("_") or "gene"
        # non-default settings get their own filename so a variant never clobbers the default
        tag = "" if (mode, scale, rows) == ("stack", "abs", "all") else f"_{mode}_{scale}_{rows}"
        base = os.path.join(outdir, f"{R['num']}_{nm}_transcription{tag}")
        fig.savefig(base + ".pdf", dpi=300)
        fig.savefig(base + ".png", dpi=300)
        open(base + "_stats.txt", "w").write(txt + "\n")
        R["files"] = {k: base + v for k, v in
                      (("pdf", ".pdf"), ("png", ".png"), ("txt", "_stats.txt"))}
    return fig, txt, R


# ============================== stats =======================================================
def _stats(R, dep, spans, iso, scale, rows, mode):
    sgn, a1, a3 = R["sgn"], R["anchor1"], R["anchor3"]
    L = [R["title"], "=" * len(R["title"]),
         f"anchor: Syn1 {a1:,}" + (f", Syn3A {a3:,}" if a3 else ", Syn3A n/a (region absent)"),
         f"rel = {'pos - anchor' if sgn > 0 else 'anchor - pos'};  window rel {R['rel_lo']} .. "
         f"{R['rel_hi']}  (Syn1 {R['win'][0]:,}-{R['win'][1]:,})",
         f"operon block(s): {', '.join(R['operons']) or 'none'}"
         + (f";  TSS {', '.join(f'{t:,}' for t in R['tss'])}" if R["tss"] else ""), ""]
    if R["s3_note"]:
        L += [R["s3_note"], ""]
    if R["brk"] is not None:
        L += [f"Syn3A leaves the syntenic block at rel {R['brk']:,} "
              f"({len(R['bands'])} deletion(s) in the window):",
              "  Syn3A panels are greyed past that point - the signal there is real Syn3A",
              "  signal from the OTHER side of the scar, not the same sequence as Syn1.", ""]
    L.append("Both-strand genome-mean depth per library "
             + ("(reference line):" if scale == "abs" else "(normaliser):"))
    for n, (_x, _s, m) in dep.items():
        L.append(f"  {n:15s} mean = {m:8.1f}")

    key = lambda n: ("syn1" if n.startswith("Syn1") else "syn3a") + "".join(n.split()[1:])
    L += ["", "Per-gene mean SENSE depth "
          + ("(RAW per-base depth):" if scale == "abs" else "(x genome mean):"),
          f"\n{'gene':16s}" + "".join(f"{key(n):>14s}" for n in dep)]
    g3 = load_genes("s3")
    for _, g in R["w1"][R["w1"].strand == R["strand"]].sort_values(
            "start0", ascending=(R["strand"] == "+")).iterrows():
        L.append(f"{(g['label'] or 'hyp') + '/' + g.num:16s}" + "".join(
            _gene_mean(dep[n], n, g, g3, R) for n in dep))

    if mode == "stack" and spans:
        L += ["", "Read stacks - primary sense-strand alignments in the window "
              f"(rows={rows}):", f"{'run':14s}{'reads':>8s}{'med_len':>9s}{'p90_len':>9s}"]
        for n, sp in spans.items():
            ln_ = sp[:, 1] - sp[:, 0] if len(sp) else np.array([0])
            L.append(f"{n:14s}{len(sp):8d}{np.median(ln_):9.0f}{np.percentile(ln_, 90):9.0f}")
    elif iso is not None and len(iso):
        L += ["", f"Top PacBio isoform clusters (n_reads >= 10, top {len(iso)}):",
              f"{'isoform':12s}{'rel5':>8s}{'rel3':>8s}{'n_reads':>9s}"]
        for _, r in iso.head(12).iterrows():
            L.append(f"{r.isoform_id:12s}{int(r.r5):8d}{int(r.r3):8d}{int(r.n_reads):9d}")

    if os.path.isfile(D(PAIRED)):
        p = pd.read_csv(D(PAIRED), sep="\t")
        p["num"] = p.locus_syn1.str[-4:]
        s = p[p.num.isin(set(R["w1"].num))]
        if len(s):
            L += ["", "Reduction paired table (mean-normalised rel units; TPM=Illumina):",
                  f"{'gene':8s}{'relTPM_syn1':>13s}{'relTPM_syn3a':>13s}{'TPM_FC':>9s}"
                  f"{'relIPM_syn1':>13s}{'relIPM_syn3a':>13s}{'iPM_FC':>9s}"]
            for _, r in s.sort_values("num").iterrows():
                L.append(f"{r.num:8s}" + "".join(
                    f"{r.get(c, np.nan):13.2f}" if "FC" not in c else f"{r.get(c, np.nan):9.2f}"
                    for c in ("relTPM_syn1", "relTPM_syn3a", "TPM_fold_change",
                              "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change")))

    try:
        s1c, s3a, s3b = load_protein_copies()
        L += ["", "Measured protein copies per cell "
              "(Syn1 syn1_omics.xlsx; Syn3A syn3A_proteome.tsv):",
              f"{'gene':16s}{'syn1':>12s}{'syn3a_2019':>12s}{'syn3a_2026':>12s}"]
        f = lambda v: f"{v:12,.0f}" if pd.notna(v) else f"{'-':>12s}"
        for _, g in R["w1"].sort_values("start0", ascending=(R["strand"] == "+")).iterrows():
            L.append(f"{(g['label'] or 'hyp') + '/' + g.num:16s}"
                     + f(s1c.get(g.num, np.nan)) + f(s3a.get(g.num, np.nan))
                     + f(s3b.get(g.num, np.nan)))
    except Exception as e:                                            # noqa
        L.append(f"  [protein copies skipped: {e}]")

    if R["tss"]:
        try:
            import sys
            sys.path.insert(0, D(PROMOTER_MOD))
            import promoter_motif as pm
            L += ["", "Promoter (-10) scan at the operon TSS(s):"]
            for t in R["tss"]:
                r = pm.scan_minus10(t, S1_CH, R["strand"])
                L.append(f"  TSS {t:,} (rel {sgn * (t - a1)}): "
                         f"-10 6-mer {r['minus10_6mer']} vs TANAAT (mm{r['mm6']}, "
                         f"shift {r['shift6']:+d}); 9-mer {r['minus10_9mer']} vs TNNTANAAT "
                         f"(mm{r['mm9']}, shift {r['shift9']:+d}); tier {r['motif_tier']}")
            L.append("  (best register within +-2 bp of -12..-7 / -15..-7; a tier above "
                     "no_minus10 needs ZERO IUPAC mismatches)")
        except Exception as e:                                    # noqa
            L.append(f"  [promoter scan skipped: {e}]")
    L += ["", "TransTermHP intrinsic terminators in the window (Syn1 coords):"]
    L += ["  " + t for t in R["terms"]] or ["  (none predicted)"]
    return "\n".join(L)


def _gene_mean(d, name, g, g3, R):
    """Mean sense depth over one gene for one library, on that library's own organism."""
    xg, s, _m = d
    if name.startswith("Syn3A"):
        h = g3[(g3.num == g.num) & (g3.strand == g.strand)]
        if not len(h) or R["anchor3"] is None:
            return f"{'-':>14s}"
        gg, anc = h.iloc[0], R["anchor3"]
    else:
        gg, anc = g, R["anchor1"]
    a, b = sorted([R["sgn"] * (gg.start0 - anc), R["sgn"] * (gg.end0 - anc)])
    msk = (xg >= a) & (xg < b)
    return f"{np.nanmean(s[msk]):14.2f}" if msk.any() else f"{'-':>14s}"
