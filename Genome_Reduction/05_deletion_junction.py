#!/usr/bin/env python3
"""
05_deletion_junction.py

Reframe each of the 95 syn1->syn3A deletions as a JUNCTION between two retained
operon fragments, and draw it so the strand relationship / operon structure /
deletion region can be eyeballed before trusting the junction-type call.

For each deleted block [s1_del_s, s1_del_e] (syn1 coords):
  operon_L = syn1 operon of the nearest syn3A-RETAINED gene to the LEFT of the
             deletion (genes absent from the syn3A annotation are skipped)
  operon_R = syn1 operon of the nearest syn3A-retained gene to the RIGHT
A gene that is only ANTISENSE-covered by an operon is still assigned to that
operon and treated as being on the OPERON's strand (transcription at that locus
runs on the operon's strand).

WHICH END OF EACH OPERON FACES THE JUNCTION
  Picture the deletion as a gap with operon_L on its left and operon_R on its
  right. Each operon presents the end nearest the gap to the junction:
    operon_L (on the left)  faces with its genomic-RIGHT end (end0):
        L = '+'  -> right end is the 3' end  = TERMINATOR  (txn runs left->right)
        L = '-'  -> right end is the 5' end  = PROMOTER    (txn runs right->left)
    operon_R (on the right) faces with its genomic-LEFT  end (start0):
        R = '+'  -> left end is the 5' end   = PROMOTER
        R = '-'  -> left end is the 3' end   = TERMINATOR
  A facing regulator is "lost" when that end's coordinate falls inside the
  deletion [s1_del_s, s1_del_e].

strand_relationship (relative transcriptional orientation of the two operons):
  intra_operon : operon_L == operon_R  -> deletion is internal to one operon.
  tandem       : operon_L.strand == operon_R.strand. The two operons are
                 co-oriented, so one transcript could in principle run across
                 the junction (upstream terminator meets downstream promoter):
                   '+'/'+': L is upstream  -> R is downstream
                   '-'/'-': R is upstream  -> L is downstream
                 This is the only fusion-capable configuration.
  convergent   : L='+', R='-'. Both operons transcribe TOWARD the junction;
                 their two TERMINATORS face each other across the gap.
  divergent    : L='-', R='+'. Both operons transcribe AWAY from the junction;
                 their two PROMOTERS face each other across the gap.

junction_type (tandem only)
  Only two regulators sit at a tandem junction and decide everything:
    - the UPSTREAM operon's TERMINATOR (its 3' end, facing the gap)
    - the DOWNSTREAM operon's PROMOTER  (its 5' end, facing the gap)
  A regulator is "lost" if its coordinate fell inside the deletion. In
  transcription order (upstream -> downstream) the original layout was:

    [==upstream==>]--TERM   ...operon(s) deleted...   PROM--[==downstream==>]

  The four classes are the 2x2 of {terminator kept/lost} x {promoter kept/lost}:

  clean_excision        terminator KEPT, promoter KEPT
        [==upstream==>]--TERM | PROM--[==downstream==>]
      Whole operon(s) in the middle excised; both flanks keep their regulators
      and work independently. No cis-regulatory change.

  decapitation          terminator KEPT, promoter LOST
        [==upstream==>]--TERM | (x)--[==downstream==>]
      Downstream lost its TSS; the intact upstream terminator blocks rescue by
      read-through -> downstream genes cannot initiate.  Predict: expression DROP.

  readthrough_extension terminator LOST, promoter KEPT
        [==upstream==>]--(x) | PROM--[==downstream==>]
      Upstream transcript no longer stops and runs into the downstream region,
      which still has its own promoter.  Predict: extended/additive coverage.

  fusion                terminator LOST, promoter LOST
        [==upstream==>]--(x) | (x)--[==downstream==>]  =>  [==upstream==>downstream==>]
      Both barriers gone -> a single chimeric transcript; downstream genes fall
      under the upstream promoter. The only class creating new co-transcription
      (ties to the cross-junction co-transcription test in 07_operon_pair_coexpression.py).

CONSISTENCY CHECK vs 04 (deletion_overlaid_operon)
  04 is operon-indexed (per operon: which boundary truncated, which genes
  deleted); 05 is deletion-indexed (per junction: did the flanking operon lose
  its facing regulator). Both are strand-aware, so each facing regulator maps
  to one truncation side: promoter -> 5'_truncation, terminator -> 3'_truncation.
  Every run appends a report (and the same checks gate against regressions):
    Check 1  regression guard: reg_lost is True but 04's overlap_class shows NO
             matching-side truncation -> a genuine 04<->05 disagreement. The
             reverse (04 shows a truncation this junction's deletion did not
             cause) is expected when another deletion hits the operon elsewhere
             and is NOT flagged. reg_lost uses the same strict boundary overlap
             as 04's classify_hit, so an abutting deletion (del_s == operon end0)
             is correctly NOT counted as removing the regulator.
    Check 2  invariant: no flank operon may be 04 'all_deleted'. (Currently 2
             expected exceptions: antisense-only flank genes whose covering
             operon is sense-all_deleted -- a documented consequence of the
             antisense->operon fallback; reported, not auto-fixed.)
    Check 3  UTR-only regulator loss: reg_lost AND gene_deletion_pattern=='intact'
             -- promoter/terminator lost without losing a CDS. Equals 04's
             5'/3'_truncation_UTR and 06_delete_gene's 'promoter_deleted'; this
             is the term that reconciles all three scripts.

Outputs (Genome_Reduction/deletion_junction/):
  deletion_junctions.tsv          one row per deletion
  deletion_junction_summary.txt   counts + cross-tabs + the consistency report
  plots/<strand_relationship>/<scar_id>.pdf   one gene-map per deletion
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch


# ---------- paths ----------
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

EVENTS_XLSX = HERE / "aln" / "analysis" / "genome_reduction_summary.xlsx"
SYN1_GFF    = ROOT / "Genomes_Input" / "syn1.genes.gff3"
SYN3A_GFF   = ROOT / "Genomes_Input" / "syn3a_genome.gff3"
SYN1_OPERONS_TSV    = ROOT / "Syn1_Operon" / "operons.candidate_blocks.tsv"
OPERON_CLASS_TSV    = HERE / "deletion_overlaid_operon" / "operon_deletion_classification.tsv"

OUT_DIR        = HERE / "deletion_junction"
OUT_JUNCTIONS  = OUT_DIR / "deletion_junctions.tsv"
OUT_SUMMARY    = OUT_DIR / "deletion_junction_summary.txt"
PLOT_DIR       = OUT_DIR / "plots"

PLOT_PAD_FRAC = 0.05
PLOT_PAD_MIN  = 300


# ---------- loaders ----------

_PAT_LOCUS = re.compile(r"locus_tag=([^;]+)")
_PAT_NAME  = re.compile(r"(?:Name|gene)=([^;]+)")


def locus_suffix(lt: str) -> str:
    return lt.rsplit("_", 1)[-1]


def build_locus2op(syn1_ops: pd.DataFrame) -> dict:
    """Map each syn1 locus_tag -> operon_id. Sense membership takes priority;
    genes that are only ANTISENSE-covered by an operon fall back to that operon
    (transcription at that locus happens on the operon's strand, so for the
    junction framing we treat the gene as belonging to that operon)."""
    m = {}
    for r in syn1_ops.itertuples():
        sense = r.sense_gene_loci if isinstance(r.sense_gene_loci, str) else ""
        for lt in (x.strip() for x in sense.split(",") if x.strip()):
            m.setdefault(lt, r.operon_id)
    for r in syn1_ops.itertuples():
        anti = r.antisense_gene_loci if isinstance(r.antisense_gene_loci, str) else ""
        for lt in (x.strip() for x in anti.split(",") if x.strip()):
            m.setdefault(lt, r.operon_id)
    return m


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
    return pd.DataFrame(rows).sort_values("start0").reset_index(drop=True)


# ---------- junction table ----------

def build_junctions(events: pd.DataFrame,
                    syn1_genes: pd.DataFrame,
                    syn1_ops: pd.DataFrame,
                    class05: pd.DataFrame,
                    syn3a_suffixes: set) -> pd.DataFrame:
    ev = events.sort_values("S1").reset_index(drop=True)

    op_by_id = {r.operon_id: r for r in syn1_ops.itertuples()}
    locus2op = build_locus2op(syn1_ops)
    cls_by_op = {r.operon_id: r for r in class05.itertuples()}

    # A syn1 gene is RETAINED iff its locus suffix is present in syn3A
    # (authoritative: a partially-deleted boundary gene that JCVI dropped from
    # the syn3A annotation is correctly treated as absent).
    g_all = syn1_genes.sort_values("start0").reset_index(drop=True)
    g = g_all[[locus_suffix(str(lt)) in syn3a_suffixes
               for lt in g_all.locus_tag]].reset_index(drop=True)

    rows = []
    idx = 0
    for i, row in ev.iterrows():
        if row["Change Case"] != "deleted":
            continue
        idx += 1
        s, e = int(row["S1"]), int(row["E1"])

        # syn3A junction coordinate = E2 of the preceding retained_ordered block
        prev_ret = None
        for j in range(i - 1, -1, -1):
            if ev.at[j, "Change Case"] == "retained_ordered":
                prev_ret = ev.iloc[j]; break
        syn3A_junction = int(prev_ret["E2"]) if prev_ret is not None and not pd.isna(prev_ret["E2"]) else -1

        # flanking syn1 genes — must be FULLY outside the deletion so that
        # deleted / straddling genes can never be chosen as a retained flank.
        left_cand  = g[g.end0   <= s]
        right_cand = g[g.start0 >= e]
        left_gene  = left_cand.iloc[left_cand.end0.values.argmax()]    if len(left_cand)  else None
        right_gene = right_cand.iloc[right_cand.start0.values.argmin()] if len(right_cand) else None

        opL_id = locus2op.get(str(left_gene.locus_tag))  if left_gene  is not None else None
        opR_id = locus2op.get(str(right_gene.locus_tag)) if right_gene is not None else None
        opL = op_by_id.get(opL_id); opR = op_by_id.get(opR_id)

        # deleted gene count (fully inside the deletion) — from the full gene set
        n_del = int(((g_all.start0 >= s) & (g_all.end0 <= e)).sum())

        base = {
            "scar_id":       f"DEL_{idx:03d}",
            "syn1_del_s":    s,
            "syn1_del_e":    e,
            "syn1_del_len":  int(row["LEN1"]) if not pd.isna(row["LEN1"]) else e - s,
            "syn3A_junction": syn3A_junction,
            "n_deleted_genes": n_del,
            "left_gene":     str(left_gene.locus_tag)  if left_gene  is not None else "",
            "right_gene":    str(right_gene.locus_tag) if right_gene is not None else "",
            "operon_L_id":   opL_id or "",
            "operon_R_id":   opR_id or "",
        }

        if opL is None or opR is None:
            base.update({"strand_relationship": "unknown", "junction_type": "unknown"})
            rows.append(base)
            continue

        Lstr, Rstr = str(opL.strand), str(opR.strand)
        same_operon = (opL_id == opR_id)
        L_reg = "terminator" if Lstr == "+" else "promoter"
        R_reg = "promoter"   if Rstr == "+" else "terminator"
        # Strict boundary overlap, matching 04's classify_hit (a deletion that
        # merely abuts a boundary, e.g. del_s == operon end0, does NOT remove it).
        #   operon_L faces its right boundary (end0): 04's `right` = del_s < end0 and del_e >= end0
        #   operon_R faces its left  boundary (start0): 04's `left`  = del_s <= start0 and del_e > start0
        L_reg_lost = (s < int(opL.end0))   and (e >= int(opL.end0))
        R_reg_lost = (s <= int(opR.start0)) and (e > int(opR.start0))

        if same_operon:
            strand_rel = "intra_operon"; jtype = "intra_operon"
        else:
            if Lstr == Rstr:
                strand_rel = "tandem"
            elif Lstr == "+" and Rstr == "-":
                strand_rel = "convergent"
            else:
                strand_rel = "divergent"

            if strand_rel == "tandem":
                if Lstr == "+":
                    up_lost_term, down_lost_prom = L_reg_lost, R_reg_lost
                else:
                    up_lost_term, down_lost_prom = R_reg_lost, L_reg_lost
                if up_lost_term and down_lost_prom:
                    jtype = "fusion"
                elif down_lost_prom:
                    jtype = "decapitation"
                elif up_lost_term:
                    jtype = "readthrough_extension"
                else:
                    jtype = "clean_excision"
            else:
                jtype = strand_rel

        clsL = cls_by_op.get(opL_id); clsR = cls_by_op.get(opR_id)
        base.update({
            "operon_L_strand": Lstr,
            "operon_L_gene_deletion_pattern": (clsL.gene_deletion_pattern if clsL is not None else ""),
            "operon_L_overlap_class":         (clsL.overlap_class if clsL is not None else ""),
            "operon_L_facing_reg": L_reg,
            "operon_L_reg_lost":   L_reg_lost,
            "operon_R_strand": Rstr,
            "operon_R_gene_deletion_pattern": (clsR.gene_deletion_pattern if clsR is not None else ""),
            "operon_R_overlap_class":         (clsR.overlap_class if clsR is not None else ""),
            "operon_R_facing_reg": R_reg,
            "operon_R_reg_lost":   R_reg_lost,
            "same_operon":         same_operon,
            "strand_relationship": strand_rel,
            "junction_type":       jtype,
        })
        rows.append(base)
    return pd.DataFrame(rows)


# ---------- plotting ----------

def _draw_gene(ax, gr, color: str, op_label: str = ""):
    y = 1.0 if gr.strand == "+" else -1.0
    s0, e0 = int(gr.start0), int(gr.end0)
    body = max(1, e0 - s0)
    head = min(body * 0.4, 180)
    if gr.strand == "+":
        tail_x, dx = s0, body
    else:
        tail_x, dx = e0, -body
    ax.add_patch(FancyArrow(tail_x, y, dx, 0, width=0.34,
                            head_width=0.58, head_length=head,
                            length_includes_head=True,
                            color=color, ec="black", lw=0.4, zorder=3))
    suf = locus_suffix(str(gr.locus_tag))
    # locus number nearest the arrow, operon id stacked further out
    if gr.strand == "+":
        label = f"{op_label}\n{suf}" if op_label else suf
        ax.text((s0 + e0) / 2, y + 0.5, label, ha="center", va="bottom",
                fontsize=5.5, zorder=4, linespacing=0.9)
    else:
        label = f"{suf}\n{op_label}" if op_label else suf
        ax.text((s0 + e0) / 2, y - 0.5, label, ha="center", va="top",
                fontsize=5.5, zorder=4, linespacing=0.9)


def plot_junction(jrow, syn1_genes: pd.DataFrame, op_by_id: dict,
                  locus2op: dict, save_path: Path) -> None:
    s, e = int(jrow["syn1_del_s"]), int(jrow["syn1_del_e"])
    opL_id = jrow.get("operon_L_id", ""); opR_id = jrow.get("operon_R_id", "")
    opL = op_by_id.get(opL_id); opR = op_by_id.get(opR_id)

    xs = [s, e]
    for op in (opL, opR):
        if op is not None:
            xs += [int(op.start0), int(op.end0)]
    x0, x1 = min(xs), max(xs)
    pad = int(PLOT_PAD_FRAC * (x1 - x0)) + PLOT_PAD_MIN
    x0 -= pad; x1 += pad

    genes = syn1_genes[(syn1_genes.start0 < x1) & (syn1_genes.end0 > x0)]

    # Assign a distinct color to each operon present in the window (by position).
    ops_in_window = []
    seen = set()
    for gr in genes.sort_values("start0").itertuples():
        opid = locus2op.get(str(gr.locus_tag))
        if opid and opid not in seen:
            seen.add(opid); ops_in_window.append(opid)
    cmap = plt.get_cmap("tab20")
    op_color = {op: cmap(i % 20) for i, op in enumerate(ops_in_window)}

    fig, ax = plt.subplots(figsize=(16, 4))
    # operon spans — colored per operon, shaded only on the operon's OWN strand
    # lane (top half for +, bottom half for -) so the strand is unambiguous.
    for opid in ops_in_window:
        op = op_by_id.get(opid)
        if op is None:
            continue
        ymin, ymax = (0.5, 1.0) if str(op.strand) == "+" else (0.0, 0.5)
        ax.axvspan(int(op.start0), int(op.end0), ymin=ymin, ymax=ymax,
                   color=op_color[opid], alpha=0.22, zorder=0)
    # deletion region — dashed boundaries (no fill, so operon colors show through)
    ax.axvline(s, color="red", ls="--", lw=1.4, zorder=2)
    ax.axvline(e, color="red", ls="--", lw=1.4, zorder=2)
    ax.axhline(0, color="gray", lw=0.6, zorder=1)

    for gr in genes.itertuples():
        opid = locus2op.get(str(gr.locus_tag))
        color = op_color.get(opid, "lightgray")
        _draw_gene(ax, gr, color, op_label=(opid or ""))

    ax.set_xlim(x0, x1)
    ax.set_ylim(-2.6, 2.6)
    ax.set_yticks([1, -1]); ax.set_yticklabels(["+ strand", "- strand"], fontsize=8)
    ax.set_xlabel("syn1 coordinate (bp)")
    title = (f"{jrow['scar_id']}   {jrow.get('strand_relationship','?')} / {jrow.get('junction_type','?')}\n"
             f"L={opL_id}({jrow.get('operon_L_strand','?')},{jrow.get('operon_L_gene_deletion_pattern','?')})"
             f"  |  DEL {s}-{e} ({jrow['n_deleted_genes']} genes, {jrow['syn1_del_len']} bp)  |  "
             f"R={opR_id}({jrow.get('operon_R_strand','?')},{jrow.get('operon_R_gene_deletion_pattern','?')})")
    ax.set_title(title, fontsize=10)
    # legend: operon -> color, marking L / R
    legend = []
    for opid in ops_in_window:
        lbl = opid
        tags = ("L" if opid == opL_id else "") + ("R" if opid == opR_id else "")
        if tags:
            lbl += f" ({tags})"
        legend.append(Patch(color=op_color[opid], label=lbl))
    legend.append(Patch(color="lightgray", label="no operon"))
    legend.append(Patch(facecolor="none", edgecolor="red", ls="--", label="deletion"))
    ax.legend(handles=legend, loc="upper right", fontsize=7,
              ncol=max(1, len(legend) // 2), framealpha=0.85)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


# ---------- consistency validation against 04 ----------

def validate_against_overlap_class(junc: pd.DataFrame) -> list[str]:
    """Cross-check the per-junction regulator-loss calls (this script, 05)
    against 04's per-operon overlap_class (already joined into `junc`).
    Both are strand-aware, so the facing regulator maps to one truncation side:
        facing promoter   -> expect a 5'_truncation_* in overlap_class
        facing terminator -> expect a 3'_truncation_*
    Reports three things:
      Check 1  regression guard: reg_lost is True but 04 sees NO matching-side
               truncation at all (a genuine 04<->05 disagreement). The reverse
               (04 shows a truncation this junction's deletion didn't cause) is
               expected when another deletion hits the operon elsewhere, so it
               is not flagged.
      Check 2  invariant : no flank operon may be 04 'all_deleted'.
      Check 3  UTR-only  : reg_lost AND gene_deletion_pattern=='intact'
               (promoter/terminator lost without losing a CDS; == 04's
               5'/3'_truncation_UTR == 06_delete_gene's promoter_deleted).
    """
    def _expected(facing: str) -> str:
        return "5'_truncation" if facing == "promoter" else "3'_truncation"

    mism, all_del, utr_only = [], [], []
    for r in junc.itertuples():
        for side in ("L", "R"):
            opid   = getattr(r, f"operon_{side}_id", "")
            facing = getattr(r, f"operon_{side}_facing_reg", "")
            lost   = getattr(r, f"operon_{side}_reg_lost", None)
            ovc    = getattr(r, f"operon_{side}_overlap_class", "")
            gdp    = getattr(r, f"operon_{side}_gene_deletion_pattern", "")
            if not isinstance(opid, str) or not opid:
                continue
            if not isinstance(facing, str) or facing not in ("promoter", "terminator"):
                continue
            ovc = ovc if isinstance(ovc, str) else ""
            gdp = gdp if isinstance(gdp, str) else ""
            if gdp == "all_deleted":
                all_del.append((r.scar_id, side, opid))
            # Only the "lost but 04 sees no truncation at all" direction is a
            # genuine disagreement. The reverse (04 shows a truncation this
            # junction's deletion didn't cause) is explained by another deletion
            # hitting the operon elsewhere, so it is NOT flagged.
            has_trunc = _expected(facing) in ovc
            if bool(lost) and not has_trunc:
                mism.append((r.scar_id, side, opid, facing, "reg_lost but overlap_class lacks "
                             + _expected(facing), ovc or "<empty>"))
            if bool(lost) and gdp == "intact":
                utr_only.append((r.scar_id, side, opid, facing, ovc or "<empty>"))

    lines = []
    lines.append("")
    lines.append("=" * 64)
    lines.append("CONSISTENCY CHECK vs 04 (operon_deletion_classification)")
    lines.append("=" * 64)
    lines.append(f"Check 1  facing-regulator-loss <-> overlap_class side : "
                 f"{len(mism)} mismatch(es)")
    for scar, side, opid, facing, why, ovc in mism:
        lines.append(f"    [{scar} {side}] {opid} ({facing}): {why}  [overlap_class={ovc}]")
    lines.append(f"Check 2  no flank operon is 04 'all_deleted'          : "
                 f"{'PASS' if not all_del else f'FAIL ({len(all_del)})'}")
    for scar, side, opid in all_del:
        lines.append(f"    [{scar} {side}] {opid} is all_deleted but used as a flank")
    lines.append(f"Check 3  UTR-only regulator loss (lost, genes intact) : {len(utr_only)}")
    for scar, side, opid, facing, ovc in utr_only:
        lines.append(f"    [{scar} {side}] {opid} ({facing} lost, no CDS deleted)  [overlap_class={ovc}]")
    return lines


# ---------- main ----------

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading events, GFFs, operons, 05 classification ...")
    events      = pd.read_excel(EVENTS_XLSX, sheet_name="events")
    syn1_genes  = load_gff_genes(SYN1_GFF)
    syn3a_genes = load_gff_genes(SYN3A_GFF)
    syn3a_suffixes = {locus_suffix(str(lt)) for lt in syn3a_genes.locus_tag}
    syn1_ops    = pd.read_csv(SYN1_OPERONS_TSV, sep="\t")
    class05     = pd.read_csv(OPERON_CLASS_TSV, sep="\t")

    op_by_id = {r.operon_id: r for r in syn1_ops.itertuples()}
    locus2op = build_locus2op(syn1_ops)

    print("Building junction table ...")
    junc = build_junctions(events, syn1_genes, syn1_ops, class05, syn3a_suffixes)
    junc.to_csv(OUT_JUNCTIONS, sep="\t", index=False)
    print(f"  junctions: {len(junc)}  -> {OUT_JUNCTIONS}")

    # ---- plots, organised by strand_relationship ----
    # Wipe the plot tree first so renamed/old files don't linger.
    if PLOT_DIR.exists():
        shutil.rmtree(PLOT_DIR)
    print("Rendering junction plots ...")
    n_plot = 0
    for r in junc.itertuples():
        sr = str(r.strand_relationship)
        sub = PLOT_DIR / sr
        sub.mkdir(parents=True, exist_ok=True)
        # under tandem/, prefix with junction_type so the four types are easy
        # to tell apart; other folders are already single-type.
        fname = (f"{r.junction_type}_{r.scar_id}.pdf" if sr == "tandem"
                 else f"{r.scar_id}.pdf")
        jrow = {k: getattr(r, k) for k in junc.columns}
        try:
            plot_junction(jrow, syn1_genes, op_by_id, locus2op, sub / fname)
            n_plot += 1
        except Exception as ex:
            print(f"  WARN {r.scar_id} plot failed: {ex}")
    print(f"  plotted {n_plot} junctions into {PLOT_DIR}/<strand_relationship>/")

    # ---- summary ----
    lines = []
    lines.append("=" * 64)
    lines.append("DELETION-JUNCTION TAXONOMY (syn1 -> syn3A)")
    lines.append("=" * 64)
    lines.append("")
    lines.append(f"deletions analysed : {len(junc)}")
    lines.append("")
    lines.append("strand_relationship:")
    sr = junc.strand_relationship.value_counts()
    for k in ("tandem", "convergent", "divergent", "intra_operon", "unknown"):
        if k in sr:
            lines.append(f"  {k:<14s} : {int(sr[k])}")
    lines.append("")
    lines.append("junction_type:")
    lines.append("  fusion                : upstream lost terminator AND downstream lost promoter")
    lines.append("  decapitation          : downstream operon lost its promoter only")
    lines.append("  readthrough_extension : upstream operon lost its terminator only")
    lines.append("  clean_excision        : both regulators retained (whole operon(s) removed between)")
    jt = junc.junction_type.value_counts()
    for k in ("fusion", "decapitation", "readthrough_extension", "clean_excision",
              "convergent", "divergent", "intra_operon", "unknown"):
        if k in jt:
            lines.append(f"  {k:<22s} : {int(jt[k])}")

    tand = junc[junc.strand_relationship == "tandem"]
    if len(tand):
        lines.append("")
        lines.append("Flanking-operon gene_deletion_pattern (from 05) at tandem junctions:")
        combo = (tand.groupby(["operon_L_gene_deletion_pattern",
                               "operon_R_gene_deletion_pattern"]).size()
                 .sort_values(ascending=False))
        for (lp, rp), n in combo.items():
            lines.append(f"  L={lp:<16s} R={rp:<16s} : {int(n)}")

    lines.append("")
    lines.append("Per-deletion detail (scar_id  rel/type  L | DEL | R):")
    lines.append("-" * 64)
    for r in junc.itertuples():
        lines.append(f"  {r.scar_id}  {str(r.strand_relationship):<12s} {str(r.junction_type):<22s}  "
                     f"L={r.operon_L_id}({getattr(r,'operon_L_strand','?')})  "
                     f"DEL {r.syn1_del_s}-{r.syn1_del_e} [{r.n_deleted_genes}g]  "
                     f"R={r.operon_R_id}({getattr(r,'operon_R_strand','?')})")

    # ---- consistency validation against 04 ----
    check_lines = validate_against_overlap_class(junc)
    lines.extend(check_lines)

    lines.append("")
    lines.append(f"Wrote: {OUT_JUNCTIONS}")
    lines.append(f"Wrote: {OUT_SUMMARY}")
    lines.append(f"Plots: {PLOT_DIR}/<strand_relationship>/")

    text = "\n".join(lines) + "\n"
    OUT_SUMMARY.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
