#!/usr/bin/env python3
"""
06_single_operon_coexpression.py

Single-operon internal co-transcription in syn3A. Unit = ONE operon; question =
"are its consecutive retained genes still co-transcribed?" Two operon classes,
both routed from 05's junction table:

  pristine     gene_deletion_pattern == 'intact' AND the operon is NOT a flank
               of any deletion junction (no deletion abuts it -> proximity
               unchanged). This is the CONTROL: it establishes the baseline
               internal-preservation rate when the reduction never touched the
               operon, so the junction-affected cases (07) can be read against it.

  intra_operon the operon is the (operon_L == operon_R) flank of an
               'intra_operon' junction in 05 -- a deletion sits strictly inside
               it, removing middle genes. Question: did excising the middle
               break the operon? Tested on the retained genes flanking the scar.

For each operon we walk consecutive RETAINED sense genes in transcription order
and test every adjacent pair with the shared co-transcription primitive
(coexpression_common.test_pair: strict spanning / loose bridging + Illumina gap
depth). Operon verdict (strict > loose > split): an operon is preserved only
when EVERY consecutive retained pair passes.

Read-based only (ONT + Illumina depth); expression fold-change tests live in
08_delete_gene.py. Operon pairs at inter-operon junctions live in
07_operon_pair_coexpression.py.

Outputs (Genome_Reduction/single_operon_coexpression/):
  single_operon_pairs.tsv     one row per consecutive retained pair
  single_operon_verdicts.tsv  one row per operon (class + verdict)
  single_operon_summary.txt   pristine baseline vs intra_operon comparison
"""

from __future__ import annotations

import sys
from pathlib import Path
import pandas as pd
import pysam

import coexpression_common as cc

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

OUT_DIR        = HERE / "single_operon_coexpression"
OUT_PAIRS      = OUT_DIR / "single_operon_pairs.tsv"
OUT_VERDICTS   = OUT_DIR / "single_operon_verdicts.tsv"
OUT_SUMMARY    = OUT_DIR / "single_operon_summary.txt"


def operon_verdict(n_pairs: int, n_strict: int, n_loose: int) -> str:
    if n_pairs == 0:
        return "untestable_no_pairs"
    if n_strict == n_pairs:
        return "preserved_strict"
    if n_loose == n_pairs:
        return "preserved_loose"
    return "split"


def analyse(operons_class: pd.DataFrame, operon_category: dict,
            syn3a_by_lt: dict, bam, depths) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Walk each categorised operon's consecutive retained pairs."""
    cls_by_id = {r.operon_id: r for r in operons_class.itertuples()}
    pair_rows, op_rows = [], []

    for op_id, category in operon_category.items():
        op = cls_by_id.get(op_id)
        if op is None:
            continue
        gdp = op.gene_deletion_pattern if isinstance(op.gene_deletion_pattern, str) else ""
        retained_str = op.retained_genes if isinstance(op.retained_genes, str) else ""
        sense_order  = op.sense_gene_locusNums if isinstance(op.sense_gene_locusNums, str) else ""
        retained_set = {x.strip() for x in retained_str.split(",") if x.strip()}
        ordered = [ln.strip() for ln in sense_order.split(",") if ln.strip()]
        ordered_retained = [ln for ln in ordered if ln in retained_set]

        if len(ordered_retained) < 2:
            op_rows.append({
                "operon_id": op_id, "category": category, "strand": op.strand,
                "gene_deletion_pattern": gdp,
                "n_retained_sense_genes": len(ordered_retained),
                "n_pairs": 0, "n_pairs_strict": 0, "n_pairs_loose": 0,
                "operon_verdict": "untestable_<2_retained",
            })
            continue

        n_pairs = n_strict = n_loose = 0
        for i in range(len(ordered_retained) - 1):
            ln_a, ln_b = ordered_retained[i], ordered_retained[i + 1]
            ga = syn3a_by_lt.get(f"JCVISYN3A_{ln_a}")
            gb = syn3a_by_lt.get(f"JCVISYN3A_{ln_b}")
            if ga is None or gb is None:
                continue
            res = cc.test_pair(bam, depths, str(ga.chrom),
                               int(ga.start0), int(ga.end0),
                               int(gb.start0), int(gb.end0), str(ga.strand))
            n_pairs += 1
            n_strict += int(res["pair_preserved_strict"])
            n_loose  += int(res["pair_preserved_loose"])
            pair_rows.append({
                "operon_id": op_id, "category": category, "strand": op.strand,
                "gene_a_locusNum": ln_a, "gene_b_locusNum": ln_b,
                **res,
            })

        op_rows.append({
            "operon_id": op_id, "category": category, "strand": op.strand,
            "gene_deletion_pattern": gdp,
            "n_retained_sense_genes": len(ordered_retained),
            "n_pairs": n_pairs, "n_pairs_strict": n_strict, "n_pairs_loose": n_loose,
            "operon_verdict": operon_verdict(n_pairs, n_strict, n_loose),
        })

    return pd.DataFrame(pair_rows), pd.DataFrame(op_rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading genes, operon classification (04), junctions (05), BAM, bedGraphs ...")
    syn3a_genes = cc.load_gff_genes(cc.SYN3A_GFF)
    syn3a_by_lt = {r.locus_tag: r for r in syn3a_genes.itertuples()}
    operons_class = pd.read_csv(cc.OPERON_CLASS_TSV, sep="\t")
    junc = pd.read_csv(cc.JUNCTIONS_TSV, sep="\t")
    bam = pysam.AlignmentFile(str(cc.ONT_BAM), "rb")
    depths = cc.load_depths()

    # Operon routing from 05's junctions.
    flank_ids = set(junc.operon_L_id.dropna()) | set(junc.operon_R_id.dropna())
    flank_ids.discard("")
    intra_ids = set(junc.loc[junc.strand_relationship == "intra_operon", "operon_L_id"].dropna())
    intra_ids.discard("")

    operon_category: dict = {}
    for r in operons_class.itertuples():
        gdp = r.gene_deletion_pattern if isinstance(r.gene_deletion_pattern, str) else ""
        if r.operon_id in intra_ids:
            operon_category[r.operon_id] = "intra_operon"
        elif gdp == "intact" and r.operon_id not in flank_ids:
            operon_category[r.operon_id] = "pristine"
    print(f"  pristine operons     : {sum(v=='pristine' for v in operon_category.values())}")
    print(f"  intra_operon operons : {sum(v=='intra_operon' for v in operon_category.values())}")

    pair_df, op_df = analyse(operons_class, operon_category, syn3a_by_lt, bam, depths)
    bam.close()
    pair_df.to_csv(OUT_PAIRS, sep="\t", index=False)
    op_df.to_csv(OUT_VERDICTS, sep="\t", index=False)

    # ---- summary ----
    lines = []
    lines.append("=" * 64)
    lines.append("SINGLE-OPERON INTERNAL CO-TRANSCRIPTION (syn3A)")
    lines.append("=" * 64)
    lines.append("Unit = one operon. Verdict = every consecutive retained pair must pass.")
    lines.append("  preserved_strict : every pair has a spanning ONT read")
    lines.append("  preserved_loose  : every pair passes loose-bridge (but not all strict)")
    lines.append("  split            : at least one pair fails both")
    lines.append("")
    VERDICTS = ("preserved_strict", "preserved_loose", "split",
                "untestable_no_pairs", "untestable_<2_retained")
    for cat in ("pristine", "intra_operon"):
        sub = op_df[op_df.category == cat]
        testable = sub[sub.n_pairs > 0]
        lines.append("-" * 64)
        lines.append(f"{cat}: {len(sub)} operons ({len(testable)} testable with >=2 retained genes)")
        for v in VERDICTS:
            n = int((sub.operon_verdict == v).sum())
            if not n:
                continue
            if v in ("preserved_strict", "preserved_loose", "split") and len(testable):
                lines.append(f"    {v:<24s} : {n}  ({n/len(testable):.0%} of testable)")
            else:
                lines.append(f"    {v:<24s} : {n}")
        # pair-level rate
        psub = pair_df[pair_df.category == cat]
        if len(psub):
            ns = int(psub.pair_preserved_strict.sum())
            nl = int(psub.pair_preserved_loose.sum())
            lines.append(f"    pair-level: {len(psub)} pairs  strict={ns} ({ns/len(psub):.0%})  "
                         f"loose={nl} ({nl/len(psub):.0%})")

    lines.append("")
    lines.append("Interpretation: the pristine rate is the baseline 'preserved' frequency")
    lines.append("when the reduction never touched the operon. Compare intra_operon against it")
    lines.append("to see whether excising interior genes breaks co-transcription beyond baseline.")
    lines.append("")
    lines.append(f"Wrote: {OUT_PAIRS}")
    lines.append(f"Wrote: {OUT_VERDICTS}")
    lines.append(f"Wrote: {OUT_SUMMARY}")

    text = "\n".join(lines) + "\n"
    OUT_SUMMARY.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
