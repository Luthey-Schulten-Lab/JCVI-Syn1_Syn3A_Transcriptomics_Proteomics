#!/usr/bin/env python3
"""
07_operon_pair_coexpression.py

Operon PAIRS at inter-operon deletion junctions. Unit = the two operons that the
reduction joined across a cut (operon_L != operon_R in 05's junction table:
tandem / convergent / divergent). Question = is the NEW cross-junction gene pair
co-transcribed in syn3A?

The cross-junction pair is the two retained flank genes recorded by 05
(left_gene, right_gene) -- the genes now adjacent across the scar. Each is mapped
to its syn3A ortholog by locus suffix and tested with the shared co-transcription
primitive (coexpression_common.test_pair: strict spanning / loose bridging +
Illumina gap depth).

Stratified by junction_type (read-based predictions):
  fusion                -> PREDICT co-transcribed (both regulators gone -> one unit)
  clean_excision        -> PREDICT NOT co-transcribed (negative control: both
                           operons kept their own terminator/promoter)
  readthrough_extension -> upstream terminator gone; may show co-transcription
  decapitation          -> downstream promoter gone; co-transcription possible via
                           upstream read-through (the EXPRESSION-drop prediction is
                           tested in 08_delete_gene.py, not here)
  convergent / divergent-> opposite strands; same-strand co-transcription is not
                           defined, so reported as 'opposite_strand' (structural
                           negative control).

Read-based only (ONT + Illumina depth). Expression fold-change tests: 08.

Outputs (Genome_Reduction/operon_pair_coexpression/):
  operon_pair_coexpression.tsv   one row per inter-operon junction
  operon_pair_summary.txt        co-transcription rate by junction_type
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

OUT_DIR     = HERE / "operon_pair_coexpression"
OUT_TSV     = OUT_DIR / "operon_pair_coexpression.tsv"
OUT_SUMMARY = OUT_DIR / "operon_pair_summary.txt"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading junctions (05), genes, BAM, bedGraphs ...")
    junc = pd.read_csv(cc.JUNCTIONS_TSV, sep="\t")
    syn3a_genes = cc.load_gff_genes(cc.SYN3A_GFF)
    syn3a_by_lt = {r.locus_tag: r for r in syn3a_genes.itertuples()}
    bam = pysam.AlignmentFile(str(cc.ONT_BAM), "rb")
    depths = cc.load_depths()

    # inter-operon junctions only (operon_L != operon_R): tandem/convergent/divergent
    inter = junc[junc.strand_relationship.isin(["tandem", "convergent", "divergent"])].copy()

    rows = []
    for r in inter.itertuples():
        lg, rg = str(r.left_gene), str(r.right_gene)
        ga = syn3a_by_lt.get(f"JCVISYN3A_{cc.locus_suffix(lg)}") if lg else None
        gb = syn3a_by_lt.get(f"JCVISYN3A_{cc.locus_suffix(rg)}") if rg else None
        base = {
            "scar_id":             r.scar_id,
            "strand_relationship": r.strand_relationship,
            "junction_type":       r.junction_type,
            "operon_L_id":         r.operon_L_id,
            "operon_R_id":         r.operon_R_id,
            "left_gene":           lg,
            "right_gene":          rg,
        }
        if ga is None or gb is None:
            rows.append({**base, "testable": False, "note": "missing syn3A ortholog"})
            continue
        if str(ga.strand) != str(gb.strand):
            # convergent / divergent: opposite strands -> same-strand co-transcription undefined
            rows.append({**base, "testable": False, "note": "opposite_strand",
                         "co_transcribed_strict": False, "co_transcribed_loose": False})
            continue
        res = cc.test_pair(bam, depths, str(ga.chrom),
                           int(ga.start0), int(ga.end0),
                           int(gb.start0), int(gb.end0), str(ga.strand))
        rows.append({
            **base, "testable": True, "note": "",
            **res,
            "co_transcribed_strict": res["pair_preserved_strict"],
            "co_transcribed_loose":  res["pair_preserved_loose"],
        })
    bam.close()

    out = pd.DataFrame(rows)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    # ---- summary ----
    lines = []
    lines.append("=" * 64)
    lines.append("OPERON-PAIR CROSS-JUNCTION CO-TRANSCRIPTION (syn3A)")
    lines.append("=" * 64)
    lines.append("Unit = the two operons joined at a deletion junction; tests the new")
    lines.append("cross-junction gene pair (left_gene -> right_gene) for co-transcription.")
    lines.append("")
    lines.append(f"inter-operon junctions: {len(out)}")
    lines.append("")
    lines.append("Co-transcription by junction_type (strict spanning / loose bridging):")
    lines.append("-" * 64)
    order = ["fusion", "readthrough_extension", "decapitation", "clean_excision",
             "convergent", "divergent"]
    for jt in order:
        sub = out[out.junction_type == jt]
        if not len(sub):
            continue
        testable = sub[sub.testable == True]
        if len(testable):
            ns = int(testable.co_transcribed_strict.fillna(False).sum())
            nl = int(testable.co_transcribed_loose.fillna(False).sum())
            lines.append(f"  {jt:<22s} n={len(sub):>3d}  testable={len(testable):>3d}  "
                         f"strict={ns} ({ns/len(testable):.0%})  loose={nl} ({nl/len(testable):.0%})")
        else:
            lines.append(f"  {jt:<22s} n={len(sub):>3d}  testable=0  (opposite-strand / no ortholog)")

    lines.append("")
    lines.append("Reading the controls:")
    lines.append("  fusion should co-transcribe at a HIGH rate (both barriers removed);")
    lines.append("  clean_excision is the NEGATIVE control (both regulators intact -> low rate);")
    lines.append("  convergent/divergent are opposite-strand (structurally not co-transcribed).")

    # list the fusion candidates explicitly (the key new operons)
    fus = out[(out.junction_type == "fusion") & (out.testable == True)]
    if len(fus):
        lines.append("")
        lines.append("Fusion junctions (cross-junction pair, co-transcription evidence):")
        for r in fus.itertuples():
            lines.append(f"  {r.scar_id}  {r.operon_L_id}->{r.operon_R_id}  "
                         f"{r.left_gene}->{r.right_gene}  "
                         f"n_span={r.n_spanning_reads} n_bridge={r.n_bridging_reads}  "
                         f"strict={r.co_transcribed_strict} loose={r.co_transcribed_loose}")

    lines.append("")
    lines.append(f"Wrote: {OUT_TSV}")
    lines.append(f"Wrote: {OUT_SUMMARY}")

    text = "\n".join(lines) + "\n"
    OUT_SUMMARY.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
