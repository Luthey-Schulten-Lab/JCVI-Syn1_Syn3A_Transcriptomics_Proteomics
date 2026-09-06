"""
Operon plot with RNA isoforms (publication style).

Per-operon visualizer for the Syn1 PacBio long-read transcriptome. Three stacked
tracks, all in transcript space (5'->3' left->right; operon strand flips genomic
coordinates):

  1. Gene track   -- every gene overlapping the plotting window (window-overlap, so
                     flanking neighbours give context). Genes are gray polygon-arrows
                     pointing by gene strand; pseudogenes are purple. Genes ANTISENSE
                     to the operon (opposite strand) are drawn with a dotted outline
                     and a more transparent fill. syn3A-deleted regions are shaded by
                     default (band from Genome_Reduction/aln/raw/syn1_deleted_regions.bed).
  2. Isoform track -- PacBio FLNC isoforms as arrows, COLORED BY 5' end (TSS group, a
                     clean Okabe-Ito palette); line width ~ log10(reads) (clipped),
                     alpha ~ reads (global). Same drawer (draw_isoforms) is reused by
                     Operon_Annotation.py.
  3. Depth track  -- (optional) strand-correct PacBio depth, light-blue fill + line,
                     loaded on demand by an awk slice of the strand bedGraph.

Axes: bottom x = Transcript coordinate (nt, 0 = operon TSS); a secondary top axis on
the gene track shows Genome position (kb).

Born at final print size per OUTPUT.md (Arial, 5-7 pt, pdf.fonttype 42); see the
"Operon plot with RNA isoforms" section there for the full spec.

Public API (kept stable -- Operon_Annotation.py imports OperonCoord + draw_isoforms;
Operon_Visualization.ipynb drives plot_one_operon):
    OperonCoord, draw_isoforms, plot_one_operon
"""

from __future__ import annotations

import os
import subprocess
import urllib.parse
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── OUTPUT.md figure preamble ────────────────────────────────────────────────
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

# ── Paths (run with Syn1_Operon/ as the working dir) ─────────────────────────
MOTHER_FOLDER = ".."
CHROM         = "CP002027.1"
GFF3_FILE     = MOTHER_FOLDER + "/Genomes_Input/syn1.genes.gff3"
SYN3A_GFF     = MOTHER_FOLDER + "/Genomes_Input/syn3a_genome.gff3"   # gene names for labels
ISOFORMS_TSV  = MOTHER_FOLDER + "/Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
DEL_BED       = MOTHER_FOLDER + "/Genome_Reduction/aln/raw/syn1_deleted_regions.bed"
DEPTH_FILES   = {
    "+": MOTHER_FOLDER + "/Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.plus.bedGraph",
    "-": MOTHER_FOLDER + "/Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.minus.bedGraph",
}

# ── Style knobs ──────────────────────────────────────────────────────────────
FIG_W, FIG_H        = 7.0, 7.0 / 2          # default born-at-size canvas (in); tunable per call
PAD_BP              = 200                    # plot window = operon span +/- PAD_BP, then EXPANDED
                                            # so every gene touching [s0-PAD, e0+PAD] is shown with
                                            # its full body (no edge slivers). Operon span itself is
                                            # marked by the dashed boundary guide lines (tx 0 / len).
MARK_SYN3A_DELETION = True                  # shade syn3A-deleted regions on the gene track
ALT_LABELS          = False                 # stagger gene labels on two rows (dense operons)
MAX_ISOFORMS_TO_PLOT = 100                  # cap per operon (top by n_reads)

GENE_COLOR   = "#7a7a7a"                    # genes: gray for BOTH strands
PSEUDO_COLOR = "#b0a0c8"                    # pseudogenes: purple
SENSE_ALPHA, ANTI_ALPHA = 0.95, 0.38        # antisense genes: more transparent fill
DEL_COLOR    = "#e8736a"                    # syn3A-deletion shading
DEPTH_FILL, DEPTH_LINE = "#9ecae1", "#3182bd"
# Okabe-Ito (colour-blind safe), cycled across TSS groups
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#000000"]


# ── Load genes (GFF3, with rna_type for pseudogenes) ─────────────────────────
def _load_genes(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] != "gene":
                continue
            ad = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            rows.append({
                "chrom": p[0], "start0": int(p[3]) - 1, "end0": int(p[4]), "strand": p[6],
                "locus_tag": ad.get("locus_tag", ""),
                "gene_name": ad.get("Name") or ad.get("gene") or ad.get("locus_tag", ""),
                "rna_type": ad.get("rna_type", ""),
            })
    return pd.DataFrame(rows)


def _load_deletions(path: str) -> List[Tuple[int, int]]:
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            out.append((int(p[1]), int(p[2])))
    return out


def _load_syn3a_names(path: str) -> dict:
    """{numeric locus suffix -> syn3A gene Name}. Genes are matched to their syn3A
    ortholog by the preserved locus number (MMSYN1_NNNN <-> JCVISYN3A_NNNN). Covers
    every syn3A gene/pseudogene (incl. rRNA/tRNA); locus-name placeholders are dropped."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] not in ("gene", "pseudogene"):
                continue
            ad = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            num = ad.get("locus_tag", "").split("_")[-1]
            nm = urllib.parse.unquote((ad.get("Name") or ad.get("gene") or "").strip())
            for sep in (";", ","):            # compound names (e.g. "rnsD; uptD") -> primary
                nm = nm.split(sep)[0]
            nm = nm.strip()
            if num and nm and not nm.startswith("JCVISYN3A"):
                out[num] = nm
    return out


GENES       = _load_genes(GFF3_FILE)
SYN3A_NAMES = _load_syn3a_names(SYN3A_GFF)
GENES["syn3a_name"] = GENES["locus_tag"].apply(lambda lt: SYN3A_NAMES.get(str(lt).split("_")[-1], ""))
ISO         = pd.read_csv(ISOFORMS_TSV, sep="\t")
DELETIONS   = _load_deletions(DEL_BED)
print(f"[Operon_Visualization] genes {len(GENES)} (syn3A-named {int((GENES['syn3a_name'] != '').sum())}) "
      f"| isoforms {len(ISO)} | deletions {len(DELETIONS)}")


def load_depth_window(strand: str, win_s: int, win_e: int, chrom: str = CHROM) -> np.ndarray:
    """Per-base depth over [win_s, win_e) on `strand`, via an awk slice of the
    strand bedGraph (fast -- avoids loading the whole genome track)."""
    cov = np.zeros(max(0, win_e - win_s))
    path = DEPTH_FILES[strand]
    cmd = ["awk", "-F", "\t", f'$1=="{chrom}" && $3>{win_s} && $2<{win_e}', path]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        s, e, v = int(f[1]), int(f[2]), float(f[3])
        a, b = max(s, win_s) - win_s, min(e, win_e) - win_s
        if b > a:
            cov[a:b] = v
    return cov


# ── Coordinate transform ──────────────────────────────────────────────────────
@dataclass(frozen=True)
class OperonCoord:
    chrom: str
    strand: str
    opid: str
    start0: int
    end0: int

    def tx_of_genome_pos0(self, pos0: int) -> int:
        """Genomic pos0 -> operon transcript coordinate (5'->3' left->right)."""
        return pos0 - self.start0 if self.strand == "+" else self.end0 - pos0

    @property
    def length(self) -> int:
        return self.end0 - self.start0


def gene_label(r) -> str:
    """name/locusNum, preferring the **syn3A ortholog's gene name** (matched by locus
    number); falls back to the syn1 name, then pseudo/locusNum, then the bare number."""
    num = str(r["locus_tag"]).split("_")[-1]
    syn3a = str(r.get("syn3a_name", "") or "").strip()
    if syn3a and syn3a.lower() != "nan":
        return f"{syn3a}/{num}"
    nm = str(r["gene_name"]).strip()
    if nm and nm != str(r["locus_tag"]) and nm.lower() != "nan":
        return f"{nm}/{num}"
    if str(r["rna_type"]) == "pseudo":
        return f"pseudo/{num}"
    return num


def get_xticks(left: int, right: int) -> np.ndarray:
    """Nice tick positions (ascending) spanning [left, right] in transcript nt."""
    span = abs(right - left)
    step = 200 if span <= 2000 else 500 if span <= 5000 else 1000 if span <= 15000 else 2000
    lo, hi = min(left, right), max(left, right)
    start = (lo // step) * step
    if start > lo:
        start -= step
    end = ((hi + step - 1) // step) * step
    return np.arange(start, end + step, step, dtype=int)


# ── Gene track ────────────────────────────────────────────────────────────────
def draw_gene_arrows(ax, oc: OperonCoord, genes_df: pd.DataFrame, fig_w: float = None):
    """Genes overlapping the window as polygon-arrows (gray; pseudogenes purple).
    Sense genes (same strand as operon) point right with a solid outline; antisense
    genes point left with a dotted outline + transparent fill. syn3A-deleted regions
    are shaded behind the genes when MARK_SYN3A_DELETION is set.

    EVERY gene is labelled via gene_label(); if the full label would be wider than the
    gene's box (estimated cheaply, no renderer), it collapses to the bare locus number."""
    xlo, xhi = ax.get_xlim()
    win_lo, win_hi = min(xlo, xhi), max(xlo, xhi)
    xlim_span_tx = win_hi - win_lo
    fig_w_inches = fig_w if fig_w is not None else ax.figure.get_size_inches()[0]
    GENE_LABEL_FONTSIZE = 5
    char_pt = 0.58 * GENE_LABEL_FONTSIZE          # approx per-character width (pt)

    # syn3A deletion overlay (shaded band, behind everything)
    if MARK_SYN3A_DELETION and DELETIONS:
        drew = False; first_span = None
        for d0, d1 in DELETIONS:
            xa, xb = oc.tx_of_genome_pos0(d0), oc.tx_of_genome_pos0(d1)
            a, b = max(min(xa, xb), win_lo), min(max(xa, xb), win_hi)
            if b > a:
                ax.axvspan(a, b, facecolor=DEL_COLOR, alpha=0.17, lw=0, zorder=0)
                drew = True
                if first_span is None:
                    first_span = (a, b)
        if drew and ALT_LABELS:                 # below-arrow blocks fill the bottom-left; label the band itself
            ax.text(sum(first_span) / 2, 0.18, "syn3A\ndeletion", ha="center", va="top",
                    fontsize=5, color="#c0392b", linespacing=0.85)
        elif drew:
            ax.text(0.005, 0.04, "syn3A deletion (shaded)", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=5, color="#c0392b")

    Y, H, TRI, HEAD_FRAC, HEAD_MIN = 0.5, 0.30, 0.30, 0.18, 30
    ax.hlines(Y, win_lo, win_hi, color="black", lw=0.8, zorder=1)

    if not genes_df.empty:
        for gi, (_, r) in enumerate(genes_df.sort_values("start0").iterrows()):
            x0 = oc.tx_of_genome_pos0(int(r["start0"]))
            x1 = oc.tx_of_genome_pos0(int(r["end0"]))
            xl, xr = min(x0, x1), max(x0, x1)
            width = xr - xl
            head = min(max(HEAD_MIN, width * HEAD_FRAC), width)

            is_pseudo = str(r["rna_type"]) == "pseudo"
            is_sense  = (str(r["strand"]) == oc.strand)
            base = PSEUDO_COLOR if is_pseudo else GENE_COLOR
            alpha = SENSE_ALPHA if is_sense else ANTI_ALPHA
            ls    = "-" if is_sense else ":"
            elw   = 0.3 if is_sense else 0.7

            if is_sense:   # points right (operon transcription direction)
                tip, bse = xr, xr - head
                v = [(xl, Y - H/2), (bse, Y - H/2), (bse, Y - TRI/2), (tip, Y),
                     (bse, Y + TRI/2), (bse, Y + H/2), (xl, Y + H/2)]
            else:          # antisense -> points left
                tip, bse = xl, xl + head
                v = [(xr, Y - H/2), (bse, Y - H/2), (bse, Y - TRI/2), (tip, Y),
                     (bse, Y + TRI/2), (bse, Y + H/2), (xr, Y + H/2)]
            ax.add_patch(mpatches.Polygon(v, closed=True, facecolor=base, alpha=alpha,
                                          edgecolor="black", lw=elw, linestyle=ls, zorder=2))

            # label EVERY gene: full label if it fits the box, else the bare locus number.
            # ALT_LABELS staggers labels on two rows (+ a thin connector) so a dense operon
            # shows every name without collapsing to the locus number.
            full_label = gene_label(r)
            locus_num  = str(r["locus_tag"]).split("_")[-1]
            box_width_pt = (width / xlim_span_tx) * (fig_w_inches * 72) * 0.92
            cx = (xl + xr) / 2
            col = "#333" if is_sense else "#999"
            if ALT_LABELS:                            # 2-line block (name over locusNum), alternating
                blk = (full_label.split("/")[0] + "\n" + locus_num) if "/" in full_label else locus_num
                if gi % 2 == 0:                       # above the arrow
                    ax.text(cx, Y + TRI/2 + 0.05, blk, ha="center", va="bottom",
                            fontsize=GENE_LABEL_FONTSIZE, color=col, clip_on=True, linespacing=0.82)
                else:                                 # below the arrow
                    ax.text(cx, Y - TRI/2 - 0.05, blk, ha="center", va="top",
                            fontsize=GENE_LABEL_FONTSIZE, color=col, clip_on=True, linespacing=0.82)
            else:
                label = full_label if (len(full_label) * char_pt) <= box_width_pt else locus_num
                ax.text(cx, Y + TRI/2 + 0.06, label, ha="center", va="bottom",
                        fontsize=GENE_LABEL_FONTSIZE, color=col, clip_on=True)

    ax.set_ylim(*((-0.15, 1.15) if ALT_LABELS else (0, 1.15)))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)


# ── Isoform track ─────────────────────────────────────────────────────────────
def _pack_rows(intervals: List[Tuple[int, int]]) -> List[int]:
    """Greedy interval packing -> row index per interval (no overlap within a row)."""
    order = sorted(range(len(intervals)), key=lambda i: (intervals[i][0], intervals[i][1]))
    ends, rows = [], [0] * len(intervals)
    for i in order:
        a, b = intervals[i]
        for r, e in enumerate(ends):
            if e <= a:
                ends[r] = b; rows[i] = r; break
        else:
            ends.append(b); rows[i] = len(ends) - 1
    return rows


def layout_isoform_tracks(iso_df: pd.DataFrame, oc: OperonCoord,
                          plot_s0: int, plot_e0: int) -> pd.DataFrame:
    """Top-N isoforms by reads, clipped to the window and converted to tx space.
    Adds tx_left, tx_right, y (packed row), group_id (=tx_left, the 5' end -> colour),
    lw (~log10 reads, clipped 0.5-3.0) and alpha (~reads, global)."""
    if iso_df.empty:
        return iso_df
    iso = iso_df.sort_values("n_reads", ascending=False).head(MAX_ISOFORMS_TO_PLOT).copy()

    txl, txr = [], []
    for _, r in iso.iterrows():
        s = min(max(int(r["start0"]), plot_s0), plot_e0)
        e = min(max(int(r["end0"]),   plot_s0), plot_e0)
        x0, x1 = oc.tx_of_genome_pos0(s), oc.tx_of_genome_pos0(e)
        txl.append(min(x0, x1)); txr.append(max(x0, x1))
    iso["tx_left"], iso["tx_right"] = txl, txr
    iso = iso[iso["tx_right"] > iso["tx_left"]].copy()
    if iso.empty:
        return iso

    iso = iso.sort_values(["tx_left", "tx_right"]).reset_index(drop=True)
    iso["y"] = [r + 1 for r in _pack_rows(list(zip(iso["tx_left"].astype(int),
                                                    iso["tx_right"].astype(int))))]
    iso["group_id"] = iso["tx_left"].astype(int)

    nmax = float(iso["n_reads"].max())
    iso["lw"]    = [float(np.clip(0.3 + 0.7 * np.log10(max(1.0, n)), 0.5, 3.0))
                    for n in iso["n_reads"].astype(float)]
    iso["alpha"] = [float(min(1.0, 0.45 + 0.5 * (n / nmax)))
                    for n in iso["n_reads"].astype(float)]
    return iso


def draw_isoforms(ax, oc: OperonCoord, iso_df: pd.DataFrame, plot_s0: int, plot_e0: int):
    """Isoform arrows in tx space, coloured by 5' end (TSS group, Okabe-Ito cycled).
    Signature kept stable -- Operon_Annotation.py calls this directly."""
    if iso_df.empty:
        ax.text(0.01, 0.5, "No isoforms", transform=ax.transAxes, va="center", fontsize=5)
        ax.set_yticks([]); return
    iso = layout_isoform_tracks(iso_df, oc, plot_s0, plot_e0)
    if iso.empty:
        ax.text(0.01, 0.5, "No isoforms (after filter)", transform=ax.transAxes, va="center", fontsize=5)
        ax.set_yticks([]); return

    starts = sorted(iso["group_id"].unique().tolist())
    cmap = {s: OKABE_ITO[i % len(OKABE_ITO)] for i, s in enumerate(starts)}
    for _, r in iso.iterrows():
        left, right, y = float(r["tx_left"]), float(r["tx_right"]), float(r["y"])
        # transcript space already runs 5'->3' left->right for BOTH strands, so the
        # arrow always points right (head at tx_right)
        ax.add_patch(mpatches.FancyArrowPatch(
            (left, y), (right, y), arrowstyle="-|>", linewidth=float(r["lw"]),
            color=cmap[int(r["group_id"])], alpha=float(r["alpha"]),
            shrinkA=0, shrinkB=0, mutation_scale=5))
    ax.set_ylim(-1, float(iso["y"].max()) + 1)
    ax.set_ylabel("RNA isoforms", fontsize=6)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)


# ── Depth track ───────────────────────────────────────────────────────────────
def draw_depth(ax, oc: OperonCoord, cov: np.ndarray, plot_s0: int, plot_e0: int, strand: str):
    """Strand-correct depth fill + line in tx space. `cov` is per-base over
    [plot_s0, plot_e0) (genome order); mapped to tx via the operon coordinate."""
    gp = np.arange(plot_s0, plot_e0)
    xtx = (gp - oc.start0) if strand == "+" else (oc.end0 - gp)
    order = np.argsort(xtx)
    x, y = xtx[order], cov[order]
    ax.fill_between(x, 0, y, color=DEPTH_FILL, lw=0, zorder=1)
    ax.plot(x, y, color=DEPTH_LINE, lw=0.4, zorder=2)

    m = float(cov.max()) if cov.size else 0.0
    if m > 0:                                          # top tick = max rounded to 1 sig fig
        mag = 10 ** int(np.floor(np.log10(m)))
        T = int(round(m / mag) * mag)
        ax.set_yticks([0, T]); ax.set_ylim(0, max(m, T) * 1.06)
    else:
        ax.set_ylim(0, 1)
    ax.set_ylabel(f"depth ({strand})", fontsize=6)
    ax.tick_params(labelsize=5)
    ax.spines[["top", "right"]].set_visible(False)


# ── Main entry point ──────────────────────────────────────────────────────────
def plot_one_operon(op, save_path: str, dpi: int = 300, PLOT_DEPTH: bool = True,
                    isoform_reads_threshold: int = 10, fig_w: float = None, fig_h: float = None):
    """Render one operon (a row with chrom/strand/start0/end0/operon_id) as a
    gene | isoform | [depth] figure and save it. Genes are taken by window-overlap
    (not the operon's locus list), so flanking context is shown."""
    chrom  = str(op["chrom"])
    strand = str(op["strand"])
    s0, e0 = int(op["start0"]), int(op["end0"])
    opid   = str(op["operon_id"])
    oc = OperonCoord(chrom=chrom, strand=strand, opid=opid, start0=s0, end0=e0)

    # touched genes = any gene overlapping the padded operon span [s0-PAD, e0+PAD];
    # the plot window is then EXPANDED so every touched gene is shown with its full
    # body (no clipping at the window edge).
    pad = PAD_BP
    touch_lo, touch_hi = s0 - pad, e0 + pad
    genes_sub = GENES[(GENES["chrom"] == chrom) &
                      (GENES["end0"] > touch_lo) & (GENES["start0"] < touch_hi)].copy()
    if not genes_sub.empty:
        plot_s0 = min(touch_lo, int(genes_sub["start0"].min()))
        plot_e0 = max(touch_hi, int(genes_sub["end0"].max()))
    else:
        plot_s0, plot_e0 = touch_lo, touch_hi

    # isoform selection is UNCHANGED -- still only isoforms overlapping the operon span
    # [s0, e0] (not the expanded window).
    iso_sub = ISO[(ISO["chrom"].astype(str) == chrom) & (ISO["strand"].astype(str) == strand) &
                  (ISO["start0"] < e0) & (ISO["end0"] > s0) &
                  (ISO["n_reads"] >= isoform_reads_threshold)].copy()

    n_panels = 3 if PLOT_DEPTH else 2
    g_hr = 2.1 if ALT_LABELS else 1.0     # taller gene track to fit 2-line blocks above + below
    hr = [g_hr, 2.6, 1.1] if PLOT_DEPTH else [g_hr + 0.2, 2.8]
    fig, axes = plt.subplots(n_panels, 1, figsize=(fig_w or FIG_W, fig_h or FIG_H),
                             height_ratios=hr, sharex=True, constrained_layout=True)
    ax_genes, ax_iso = axes[0], axes[1]
    ax_depth = axes[2] if PLOT_DEPTH else None

    # x-limits in tx space cover the expanded window [plot_s0, plot_e0]
    xa, xb = oc.tx_of_genome_pos0(plot_s0), oc.tx_of_genome_pos0(plot_e0)
    left, right = min(xa, xb), max(xa, xb)
    ax_genes.set_xlim(left, right)

    draw_gene_arrows(ax_genes, oc, genes_sub, fig_w=(fig_w or FIG_W))
    draw_isoforms(ax_iso, oc, iso_sub, plot_s0, plot_e0)
    if PLOT_DEPTH:
        cov = load_depth_window(strand, plot_s0, plot_e0, chrom)
        draw_depth(ax_depth, oc, cov, plot_s0, plot_e0, strand)

    # secondary genome-kb axis on the gene track (top)
    def _tx2kb(x):  # tx -> genome kb
        return ((s0 + x) if strand == "+" else (e0 - x)) / 1000.0
    def _kb2tx(g):  # genome kb -> tx
        return (g * 1000.0 - s0) if strand == "+" else (e0 - g * 1000.0)
    secax = ax_genes.secondary_xaxis("top", functions=(_tx2kb, _kb2tx))
    secax.set_xlabel("Genome position (kb)", fontsize=6)
    secax.tick_params(labelsize=5)

    # bottom x-axis: transcript coordinate (nt). Only the bottom panel shows x ticks;
    # the upper panels get their tick MARKS removed too (sharex would otherwise bleed
    # the bottom ticks up between the gene and isoform panels).
    bottom = ax_depth if PLOT_DEPTH else ax_iso
    bottom.set_xticks(get_xticks(left, right))
    bottom.set_xlabel("Transcript coordinate (nt)", fontsize=7)
    bottom.tick_params(axis="x", labelsize=5, bottom=True, labelbottom=True)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)

    # mark the operon span with dashed boundary guide lines (tx 0 and tx len) so the
    # flanking context does not blur where the operon begins/ends
    for ax in axes:
        ax.axvline(0, color="#888", lw=0.5, ls=(0, (4, 3)), zorder=0.5)
        ax.axvline(oc.length, color="#888", lw=0.5, ls=(0, (4, 3)), zorder=0.5)

    fig.savefig(save_path, dpi=dpi)
    plt.close(fig)
