#!/usr/bin/env python
"""
R2 (RNase / RNA-processing) figure panels -- born-at-size per OUTPUT.md.

Output letters follow the Fig. 2 legend in MANUSCRIPT.md (the CLI keys a/b/c/f below are the
function names, NOT the figure letters):
  panel_a -> R2c_lap_5p_erosion.pdf      (figure panel c) : lap/0154 (OP_00078, + strand) 5' erosion
            ladder; isoforms sorted by 5' end, coloured by endpoint-context category; complete-5'
            isoforms at the bottom get EXTRA spacing so the candidate endo cut can be annotated.
  panel_b -> R2b_0178_3p_erosion.pdf     (figure panel b) : 0178/neopullulanase (OP_00099, - strand)
            3' erosion ladder, drawn 5'->3' (minus strand, high genome coord on the LEFT).
  panel_c -> R2e_truncation_categories.pdf (figure panel e) : isoform truncation categories.
  panel_f -> R2g_atp_synthase.pdf        (figure panel g) : ATP-synthase operon, split at atpA/alpha.
(panel d = the 0178 3' secondary structure -> fold_3prime_terminator.py; panel f = Illustrator schematic.)

Only isoforms CONTAINED within the operon span are shown. Categories / colours are shared
across the erosion + truncation panels and emitted once as a standalone legend (R2_legend.pdf):
  unprocessed (grey) | 5' eroded (blue) | 3' eroded (orange = Exo 3'->5' icon) | both (purple)
Per-panel category percentages are written to R2_panels/R2_panels.txt.

Run in base env (pandas + matplotlib):
  /home/enguang/anaconda3/bin/python Syn1_RNase/R2_figure_panels.py            # all
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.patheffects as pe

mpl.rcParams.update({"font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ISOF  = f"{ROOT}/Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.tsv"
OPS   = f"{ROOT}/Syn1_Operon/operons.candidate_blocks.tsv"
GFF   = f"{ROOT}/Genomes_Input/syn1.genes.gff3"
EPCTX = f"{ROOT}/Syn1_RNase/RNase/isoform_endpoint_context.tsv"
DEPTH = f"{ROOT}/Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.{{}}.bedGraph"
CHROM = "CP002027.1"
OUT   = f"{ROOT}/Syn1_RNase/R2_panels"
os.makedirs(OUT, exist_ok=True)

CATS = ["unprocessed", "5p_intragenic_only", "3p_intragenic_only", "both_intragenic"]
CLAB = ["Unprocessed", "5′ eroded", "3′ eroded", "Both eroded"]
CCOL = {"unprocessed": "#9e9e9e", "5p_intragenic_only": "#3b6db3",
        "3p_intragenic_only": "#f5901f", "both_intragenic": "#7a4fa3"}   # 3' eroded = Exo 3'->5' icon orange (PAC_O)
KB = FuncFormatter(lambda x, _: f"{x/1000:.1f}")

# example loci; op = operon span (0-based half-open), the display + containment window
A = dict(locus="MMSYN1_0154", name="0154/leucyl aminopeptidase", strand="+", gstart=197743, gend=199098,
         tss=197657, tts=199153, op=(197657, 199153), flip=False)
B = dict(locus="MMSYN1_0178", name="0178/neopullulanase", strand="-", gstart=232672, gend=234468,
         tss=234512, tts=232650, op=(232650, 234512), flip=True)

STATS = {}   # panel -> breakdown, filled by panel_a / panel_b


# ---------------------------------------------------------------- shared loaders
def load_intragenic_mask():
    rows, L = [], 0
    for ln in open(GFF):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or "rna_type=mRNA" not in f[8]:
            continue
        s0, e1 = int(f[3]) - 1, int(f[4])
        rows.append((s0, e1, f[6]))
        L = max(L, e1)
    m = {"+": np.zeros(L, bool), "-": np.zeros(L, bool)}
    for s0, e1, st in rows:
        if st in m:
            m[st][s0:e1] = True
    return m


def endpoints(strand, s0, e0):
    return (s0, e0 - 1) if strand == "+" else (e0 - 1, s0)


def classify(strand, p5, p3, mask):
    in5 = bool(0 <= p5 < len(mask[strand]) and mask[strand][p5])
    in3 = bool(0 <= p3 < len(mask[strand]) and mask[strand][p3])
    return ("both_intragenic" if in5 and in3 else
            "5p_intragenic_only" if in5 else
            "3p_intragenic_only" if in3 else "unprocessed")


def operon_isoforms(iso, cfg, mask, min_reads):
    """Isoforms CONTAINED within the operon span (omit reads spilling into neighbours)."""
    lo, hi = cfg["op"]
    s = iso[(iso.strand == cfg["strand"]) & (iso.start0 >= lo) & (iso.end0 <= hi) &
            (iso.n_reads >= min_reads)].copy()
    p5, p3, cat = [], [], []
    for r in s.itertuples():
        a, b = endpoints(cfg["strand"], r.start0, r.end0)
        p5.append(a); p3.append(b); cat.append(classify(cfg["strand"], a, b, mask))
    s["p5"], s["p3"], s["cat"] = p5, p3, cat
    return s


def load_depth(strand, lo, hi):
    d = np.zeros(hi - lo)
    for ln in open(DEPTH.format("plus" if strand == "+" else "minus")):
        f = ln.split()
        if f[0] != CHROM:
            continue
        a, b = int(f[1]), int(f[2])
        if b < lo or a > hi:
            continue
        d[max(a, lo) - lo:min(b, hi) - lo] = float(f[3])
    return np.arange(lo, hi), d


def breakdown(s):
    """{cat: (n_iso, pct_iso, n_reads, pct_reads)} over the displayed isoforms."""
    ni, nr = len(s), int(s["n_reads"].sum())
    out = {}
    for c in CATS:
        sub = s[s["cat"] == c]
        out[c] = (len(sub), 100 * len(sub) / ni if ni else 0,
                  int(sub["n_reads"].sum()), 100 * sub["n_reads"].sum() / nr if nr else 0)
    return out


# ---------------------------------------------------------------- drawing
def gene_arrow(ax, gstart, gend, strand, y=0.5, h=0.55, fc="#d9d9d9", ec="#888"):
    tip = min((gend - gstart) * 0.18, 90)
    if strand == "+":
        pts = [(gstart, y - h / 2), (gend - tip, y - h / 2), (gend, y),
               (gend - tip, y + h / 2), (gstart, y + h / 2)]
    else:
        pts = [(gend, y - h / 2), (gstart + tip, y - h / 2), (gstart, y),
               (gstart + tip, y + h / 2), (gend, y + h / 2)]
    ax.add_patch(plt.Polygon(pts, closed=True, facecolor=fc, edgecolor=ec, lw=0.6))


def erosion_panel(out, cfg, iso, mask, min_reads=2, sort="5p",
                  space_bottom=False, gap_small=1.0, gap_big=2.4, sep=3.0,
                  pad=90, figsize=(7 / 2, 7 / 4)):
    strand, flip = cfg["strand"], cfg["flip"]
    op_lo, op_hi = cfg["op"]
    lo, hi = op_lo - pad, op_hi + pad           # view window: PAD beyond the operon ends
    s = operon_isoforms(iso, cfg, mask, min_reads)
    STATS[os.path.basename(out)] = (cfg, breakdown(s))

    key = "p5" if sort == "5p" else "p3"
    ascending = (strand == "+") if sort == "5p" else (strand == "-")   # complete end -> bottom
    s = s.sort_values(key, ascending=ascending).reset_index(drop=True)
    if sort == "5p":
        bottom = ~s["cat"].isin(["5p_intragenic_only", "both_intragenic"])
    else:
        bottom = ~s["cat"].isin(["3p_intragenic_only", "both_intragenic"])

    ys, y, prev_b = [], 0.0, True
    for i in range(len(s)):
        is_b = bool(bottom.iloc[i])
        if space_bottom and (not is_b) and prev_b and i > 0:
            y += sep
        y += (gap_big if (space_bottom and is_b) else gap_small)
        ys.append(y); prev_b = is_b
    s["y"] = ys

    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.7, 5.0, 0.9], hspace=0.10)
    axg, axi, axd = (fig.add_subplot(gs[r, 0]) for r in range(3))   # gene / isoforms / depth

    # gene arrow (top) -- the reference ORF
    gene_arrow(axg, cfg["gstart"], cfg["gend"], strand)
    axg.text((cfg["gstart"] + cfg["gend"]) / 2, 0.5, cfg["name"], ha="center", va="center",
             fontsize=6, fontstyle="italic", color="#333")
    axg.set_ylim(0, 1); axg.axis("off")

    # isoform stack (middle) -- hangs from the gene
    for r in s.itertuples():
        lw = 0.35 + 0.55 * np.log10(r.n_reads + 1)
        x0, x1 = (r.p5, r.p3) if strand == "+" else (r.p3, r.p5)
        axi.plot([min(x0, x1), max(x0, x1)], [r.y, r.y], color=CCOL[r.cat],
                 lw=min(lw, 2.6), solid_capstyle="round")
    ytop = max(ys) if ys else 1
    axi.set_ylim(-(0.06 * ytop + 3), ytop + 3)            # leave space below the lowest isoform
    axi.set_yticks([]); axi.set_xticks([])
    axi.set_ylabel(f"isoforms (n={len(s)})", fontsize=6)
    for sp in ("top", "right", "bottom"):
        axi.spines[sp].set_visible(False)

    # depth (bottom) -- coverage footer carrying the shared genome axis
    xd, dd = load_depth(strand, lo, hi)
    axd.fill_between(xd, dd, color="#cfe2f3", lw=0)
    axd.plot(xd, dd, color="#3b6db3", lw=0.4)
    dmax = dd.max() if dd.max() else 1
    axd.set_ylim(0, dmax * 1.1)
    kr = max(1, round(dmax / 1000))                       # depth tick rounded to kilo
    axd.set_yticks([kr * 1000]); axd.set_yticklabels([f"{kr}k"])
    axd.set_ylabel("depth", fontsize=5, rotation=0, ha="right", va="center")
    axd.set_xlabel("Syn1 Genome Position (kb)", fontsize=6)
    axd.xaxis.set_major_formatter(KB)
    axd.tick_params(axis="both", labelsize=5, length=2, pad=1)
    for sp in ("top", "right"):
        axd.spines[sp].set_visible(False)

    # operon boundaries (dashed) + shared x-limits on every row
    for ax in (axg, axi, axd):
        for xb in (op_lo, op_hi):
            ax.axvline(xb, color="#999", lw=0.6, ls=(0, (3, 2)), zorder=0)
        ax.set_xlim((hi, lo) if flip else (lo, hi))       # minus strand -> 5'->3' (high coord left)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[{os.path.basename(out)}] {len(s)} isoforms (min_reads={min_reads}, operon-contained)")


def panel_a(iso, mask):
    erosion_panel(f"{OUT}/R2c_lap_5p_erosion.pdf", A, iso, mask,    # figure panel c (lap, 5' erosion)
                  min_reads=2, sort="5p", space_bottom=True)


def panel_b(iso, mask):
    erosion_panel(f"{OUT}/R2b_0178_3p_erosion.pdf", B, iso, mask,    # figure panel b (0178, 3' erosion)
                  min_reads=2, sort="3p", space_bottom=False)


def panel_c(*_):
    d = pd.read_csv(EPCTX, sep="\t")
    kinds = d["category"].value_counts()
    reads = d.groupby("category")["n_reads"].sum()
    rows = [("Reads", reads, reads.sum()), ("Isoforms", kinds, len(d))]

    fig, ax = plt.subplots(figsize=(7 / 4, 7 / 4), constrained_layout=True)
    for yi, (lab, ser, tot) in enumerate(rows):
        left = 0.0
        for c in CATS:
            w = 100.0 * ser.get(c, 0) / tot
            ax.barh(yi, w, left=left, color=CCOL[c], edgecolor="white", lw=0.6, height=0.62)
            if w >= 6:
                ax.text(left + w / 2, yi, f"{w:.0f}", ha="center", va="center",
                        fontsize=5.5, color="white", fontweight="bold")
            left += w
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Reads", "Isoforms"], fontsize=6.5)
    ax.set_xlim(0, 100); ax.set_ylim(-0.7, 1.75)
    ax.set_xlabel("% of isoform pool", fontsize=6.5)
    ax.tick_params(length=2, pad=1, labelsize=5.5)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    r5 = reads.get("5p_intragenic_only", 0) + reads.get("both_intragenic", 0)
    r3 = reads.get("3p_intragenic_only", 0) + reads.get("both_intragenic", 0)
    ratio = r5 / r3
    ax.text(50, 1.5, f"3′ erosion dominates  (5′/3′ eroded reads = {ratio:.2f})",
            ha="center", va="bottom", fontsize=5, color="#444", fontstyle="italic")
    fig.savefig(f"{OUT}/R2e_truncation_categories.pdf", dpi=300)    # figure panel e
    plt.close(fig)
    print(f"[R2e] reads 5'/3' eroded ratio = {ratio:.3f}")


def make_legend():
    """Standalone shared legend for a/b/c -- place it top-centre over panels a+b."""
    fig = plt.figure(figsize=(7 / 2, 7 / 24), constrained_layout=True)
    ax = fig.add_subplot(111); ax.axis("off")
    handles = [plt.Line2D([0], [0], color=CCOL[c], lw=4, label=CLAB[i]) for i, c in enumerate(CATS)]
    ax.legend(handles=handles, fontsize=6, frameon=False, ncol=4, loc="center",
              handlelength=1.2, columnspacing=1.3, handletextpad=0.4)
    fig.savefig(f"{OUT}/R2_legend.pdf", dpi=300)
    plt.close(fig)
    print("[R2_legend] standalone shared legend written")


def write_stats():
    L = ["R2 panels -- per-gene erosion-category composition (operon-contained isoforms, n_reads>=2)",
         "=" * 78, ""]
    for fn, (cfg, bd) in STATS.items():
        ni = sum(v[0] for v in bd.values())
        nr = sum(v[2] for v in bd.values())
        L.append(f"{fn}  --  {cfg['name']} ({cfg['locus']}, {cfg['strand']} strand, operon {cfg['op'][0]}-{cfg['op'][1]})")
        L.append(f"   isoforms shown {ni} ; reads {nr}")
        L.append(f"   {'category':<20}{'isoforms':>12}{'% iso':>9}{'reads':>12}{'% reads':>10}")
        for c in CATS:
            n, pi, r, pr = bd[c]
            L.append(f"   {c:<20}{n:>12}{pi:>8.1f}%{r:>12}{pr:>9.1f}%")
        # union 5'/3' ratio
        r5 = bd['5p_intragenic_only'][2] + bd['both_intragenic'][2]
        r3 = bd['3p_intragenic_only'][2] + bd['both_intragenic'][2]
        L.append(f"   5'/3' eroded reads ratio = {r5 / r3:.3f}" if r3 else "   (no 3' erosion)")
        L.append("")
    txt = "\n".join(L) + "\n"
    open(f"{OUT}/R2_panels.txt", "w").write(txt)
    print(txt)


# ---------------------------------------------------------------- panel f: ATP synthase operon
# The atp operon splits into a 5' block (a,c,b,delta) and a 3' block (gamma,beta,eps) at atpA
# (0792, alpha), assigned to RNase III by homology to B. subtilis atpA (BSU_36830).  The two
# B. subtilis RNase III sites were transferred onto Syn1 by protein homology (whole-gene fold +
# transcript-fraction projection; map_bsub_rnase_to_syn1.py) -> Syn1 mapped cut sites 932767 / 932881.
# These are the HOMOLOGY-MAPPED cleavage sites (NOT a structurally confirmed duplex): in Syn1 the two
# cuts fold into separate local stem-loops (cross-dist 104 nt), so the long B. subtilis stem is not
# conserved -- each mapped cut still sits at a stem.  Source: rnaseIII_syn1_predicted_cleavage_pairs.tsv.
# Local structure drawn by RNase_Site_Mapping/render_rnaseIII_structure.py ->
# output/rnaseIII/stems/R2_MMSYN1_0792_atpA_rnaseIII_stem.pdf (overlay as the panel-g inset).
ATP = dict(op5="OP_00395", op3="OP_00394", cuts=[932767, 932881], win=(929400, 936600),
           genes=["MMSYN1_0789", "MMSYN1_0790", "MMSYN1_0791", "MMSYN1_0792",
                  "MMSYN1_0793", "MMSYN1_0794", "MMSYN1_0795", "MMSYN1_0796", "MMSYN1_0797"])
def _rgb(*c):
    return tuple(v / 255 for v in c)

SUBUNIT = {"0796": "a", "0795": "c", "0794": "b", "0793": "δ", "0792": "α",
           "0791": "γ", "0790": "β", "0789": "ε", "0797": ""}
# per-subunit gene-arrow fills (RGB); c subunit grey, 0797 (no subunit) black
SUBUNIT_COLORS = {"0796": _rgb(61, 132, 181), "0795": _rgb(174, 174, 174),
                  "0794": _rgb(202, 112, 199), "0793": _rgb(234, 52, 38),
                  "0792": _rgb(219, 120, 66), "0791": _rgb(244, 193, 66),
                  "0790": _rgb(158, 214, 126), "0789": _rgb(76, 124, 49),
                  "0797": _rgb(0, 0, 0)}
F0_COL, F1_COL = _rgb(74, 124, 179), _rgb(0, 146, 69)   # isoform colour: 5'/F0 vs 3'/F1 block
GENOME_LEN = 1_078_809                                   # Syn1 (CP002027.1)
_MEAN_DEPTH = {}


def genome_mean_depth(strand):
    """Genome-wide mean per-base PacBio depth on one strand (for x-mean normalisation)."""
    if strand not in _MEAN_DEPTH:
        tot = 0.0
        for ln in open(DEPTH.format("plus" if strand == "+" else "minus")):
            c = ln.split()
            if not c or c[0] != CHROM:
                continue
            tot += (int(c[2]) - int(c[1])) * float(c[3])
        _MEAN_DEPTH[strand] = tot / GENOME_LEN
    return _MEAN_DEPTH[strand]


def _lum_textcolor(rgb01):
    r, g, b = (v * 255 for v in rgb01)
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 140 else "black"


def load_genes(loci):
    want = set(loci); g = {}
    for ln in open(GFF):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if attr.get("locus_tag") in want:
            g[attr["locus_tag"]] = (int(f[3]) - 1, int(f[4]), f[6])
    return g


def panel_f(iso, mask=None):
    ops = pd.read_csv(OPS, sep="\t").set_index("operon_id")
    genes = load_genes(ATP["genes"])
    lo, hi = ATP["win"]; cuts = ATP["cuts"]
    isoi = iso.set_index("isoform_id")

    rows = []
    for opid, blk in ((ATP["op5"], "F0"), (ATP["op3"], "F1")):
        for m in str(ops.loc[opid, "member_ids"]).split(","):
            if m in isoi.index:
                r = isoi.loc[m]
                rows.append((int(r.start0), int(r.end0), int(r.n_reads), blk))
    df = pd.DataFrame(rows, columns=["start0", "end0", "n_reads", "block"])
    df["p5"] = df["end0"] - 1                                  # minus strand: 5' = high coord
    df = df.sort_values("p5", ascending=False).reset_index(drop=True)

    fig = plt.figure(figsize=(7, 7 / 3), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.7, 3.9, 1.7], hspace=0.10)
    axg, axi, axd = (fig.add_subplot(gs[r, 0]) for r in range(3))

    # gene arrows (top): one colour per subunit (atpA/alpha keeps its own colour; cut shown below)
    axg.hlines(0.5, lo, hi, color="#555", lw=0.7, zorder=0)
    for lt in ATP["genes"]:
        s0, e0, st = genes[lt]
        num = lt.split("_")[-1]
        fc = SUBUNIT_COLORS.get(num, "#d9d9d9")
        gene_arrow(axg, s0, e0, st, fc=fc, ec="#555")
        sub = SUBUNIT.get(num, "")
        lab = f"{num}\n{sub}" if sub else num
        axg.text((s0 + e0) / 2, 0.5, lab, ha="center", va="center",
                 fontsize=4.0, color="white", linespacing=0.9, zorder=3,
                 path_effects=[pe.withStroke(linewidth=0.7, foreground="#333")])
    axg.set_ylim(0, 1); axg.axis("off")

    # isoforms (middle): one continuous stack ordered by 5' end (high coord first), 3' end as tiebreak;
    # coloured by block (F0 trans-membrane / F1 peripheral). line width ~ sqrt(n_reads), not log
    nmax = max(int(df["n_reads"].max()), 1)
    LW_MIN, LW_MAX = 0.4, 3.3
    order = df.sort_values(["p5", "start0"], ascending=[False, True]).reset_index(drop=True)
    for i, r in enumerate(order.itertuples(index=False)):
        lw = LW_MIN + (LW_MAX - LW_MIN) * (r.n_reads / nmax) ** 0.5
        col = F0_COL if r.block == "F0" else F1_COL
        axi.plot([r.start0, r.end0], [i, i], color=col, lw=lw, solid_capstyle="round")
    axi.set_ylim(-1.0, len(order)); axi.set_yticks([]); axi.set_xticks([])
    axi.set_ylabel("RNA isoforms", fontsize=6)
    for sp in ("top", "right", "bottom"):
        axi.spines[sp].set_visible(False)

    # depth (bottom): minus strand, normalised to the BOTH-STRAND genome-mean coverage
    # (plus+minus; = avg sense-gene depth), matching the Fig 4-6 depth panels + Table S1
    xd, dd = load_depth("-", lo, hi)
    dd = dd / (genome_mean_depth("+") + genome_mean_depth("-"))
    axd.fill_between(xd, dd, color="#dcdcdc", lw=0); axd.plot(xd, dd, color="#7a7a7a", lw=0.5)
    dmax = dd.max() if dd.max() else 1
    axd.set_ylim(0, dmax * 1.1)
    step = max(1, round(dmax / 2))
    ticks = list(range(step, int(dmax) + 1, step)) or [round(dmax, 1)]
    axd.set_yticks(ticks); axd.set_yticklabels([f"{t}×" for t in ticks])
    axd.set_ylabel("Depth\n(× mean)", fontsize=5, rotation=0, ha="right", va="center")
    axd.set_xlabel("Syn1 Genome Position (kb)", fontsize=6)
    axd.xaxis.set_major_formatter(KB); axd.tick_params(labelsize=5, length=2, pad=1)
    for sp in ("top", "right"):
        axd.spines[sp].set_visible(False)

    # minus-strand 5'->3' runs from high coord (left) to low coord (right)
    for ax in (axg, axi, axd):
        ax.set_xlim(hi, lo)
    fig.savefig(f"{OUT}/R2g_atp_synthase.pdf", dpi=300)
    plt.close(fig)
    n5 = int((df["block"] == "F0").sum())
    print(f"[R2g] ATP synthase: {len(df)} isoforms ({n5} 5'/F0 / {len(df)-n5} 3'/F1); "
          f"both-strand mean depth={genome_mean_depth('+')+genome_mean_depth('-'):.0f}x, peak={dmax:.1f}x mean")


PANELS = {"a": panel_a, "b": panel_b, "c": panel_c, "f": panel_f}


def main():
    which = [k for k in sys.argv[1:] if k in PANELS] or list(PANELS)
    iso = mask = None
    if any(k in ("a", "b", "f") for k in which):
        iso = pd.read_csv(ISOF, sep="\t",
                          usecols=["isoform_id", "strand", "start0", "end0", "n_reads"])
        mask = load_intragenic_mask()
    for k in which:
        PANELS[k](iso, mask)
    make_legend()
    if STATS:
        write_stats()


if __name__ == "__main__":
    main()
