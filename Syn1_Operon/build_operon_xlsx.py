#!/usr/bin/env python
"""
build_operon_xlsx.py
====================
Assemble Supplementary Data S1 (operon.xlsx) for the manuscript.

Pure assembly: it JOINS the per-operon tables already produced by the operon
pipeline -- it does not recompute any signature, so there is no algorithm drift.

Inputs (all under Syn1_Operon/):
  operons.candidate_blocks.tsv                            canonical 459-operon map (Operon_Segmentation.py)
  annotation/canonical/promoter_minus10_classification.tsv  -10 promoter signature, canonical operons (Operon_Annotation.py)
  annotation/canonical/terminator_tts_classification.tsv    intrinsic-terminator signature, canonical operons (Operon_Annotation.py)
  annotation/canonical/operon_utr.tsv                       5'/3' UTR lengths (Operon_Annotation.py)
  protein_complexes.xlsx                                    macromolecular-complex -> member gene loci (curated)

Promoter and terminator columns are populated for CANONICAL operons only (TSS and
TTS both intergenic); non-canonical operons get blanks, with is_canonical = False.

Output:
  operon.xlsx   sheet 'Operons'  (one row per operon, merged signatures + complexes)
                sheet 'Protein_complexes'  (the curated complex table, verbatim)

    python Syn1_Operon/build_operon_xlsx.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ANN = HERE / "annotation" / "canonical"

OPERONS = HERE / "operons.candidate_blocks.tsv"
PROMOTER = ANN / "promoter_minus10_classification.tsv"
TERMINATOR = ANN / "terminator_tts_classification.tsv"
UTR = ANN / "operon_utr.tsv"
COMPLEXES = HERE / "protein_complexes.xlsx"
OUT = HERE / "operon.xlsx"

LOCUS = lambda s: str(s).split("_")[-1].zfill(4)   # MMSYN1_0685 / 0685 -> '0685'


def complex_map(pc):
    """suffix '0685' -> sorted list of complex names that contain it."""
    m = {}
    for _, r in pc.iterrows():
        name = str(r["Complex"]).strip()
        for g in str(r["Compositions - Gene Products"]).split(";"):
            g = g.strip()
            if g:
                m.setdefault(LOCUS(g), set()).add(name)
    return {k: sorted(v) for k, v in m.items()}


def operon_complexes(loci_field, cmap):
    """semicolon list of complexes represented among an operon's sense genes."""
    if pd.isna(loci_field) or not str(loci_field).strip():
        return ""
    found = set()
    for locus in str(loci_field).split(","):
        found.update(cmap.get(LOCUS(locus), []))
    return ";".join(sorted(found))


def main():
    op = pd.read_csv(OPERONS, sep="\t")
    prom = pd.read_csv(PROMOTER, sep="\t")
    utr = pd.read_csv(UTR, sep="\t")
    pc = pd.read_excel(COMPLEXES)
    term = (pd.read_csv(TERMINATOR, sep="\t") if TERMINATOR.is_file()
            else pd.DataFrame(columns=["operon_id"]))
    if not TERMINATOR.is_file():
        print(f"WARNING: {TERMINATOR.name} missing -- run Operon_Annotation.py first; "
              "terminator columns will be blank.")

    canonical = set(prom["operon_id"])      # promoter scan == canonical operons
    df = op.copy()
    df["is_canonical"] = df["operon_id"].isin(canonical)

    # promoter signature (canonical only)
    prom_cols = {"motif_tier": "promoter_minus10_tier",
                 "minus10_6mer_best": "promoter_minus10_6mer",
                 "minus10_9mer_best": "promoter_minus10_9mer"}
    df = df.merge(prom[["operon_id", *prom_cols]].rename(columns=prom_cols),
                  on="operon_id", how="left")

    # terminator signature (canonical only)
    term_cols = {"has_tts_term": "has_terminator", "best_conf": "terminator_conf",
                 "term_stem_bp": "terminator_stem_bp", "term_loop_nt": "terminator_loop_nt",
                 "term_polyU_nt": "terminator_polyU_nt", "term_tail3_seq": "terminator_polyU_seq"}
    keep = ["operon_id"] + [c for c in term_cols if c in term.columns]
    df = df.merge(term[keep].rename(columns=term_cols), on="operon_id", how="left")

    # UTR lengths
    df = df.merge(utr[["operon_id", "utr5_bp", "utr3_bp"]], on="operon_id", how="left")

    # macromolecular-complex annotation
    cmap = complex_map(pc)
    df["operon_complexes"] = df["sense_gene_loci"].apply(lambda x: operon_complexes(x, cmap))

    order = ["operon_id", "chrom", "strand", "start0", "end0", "length",
             "segmentation_type", "is_canonical", "tss", "tts",
             "n_isoforms", "n_reads_total",
             "sense_gene_count", "sense_gene_loci", "sense_gene_names",
             "antisense_gene_count", "antisense_gene_loci", "antisense_gene_names",
             "utr5_bp", "utr3_bp",
             "promoter_minus10_tier", "promoter_minus10_6mer", "promoter_minus10_9mer",
             "has_terminator", "terminator_conf", "terminator_stem_bp",
             "terminator_loop_nt", "terminator_polyU_nt", "terminator_polyU_seq",
             "operon_complexes"]
    order = [c for c in order if c in df.columns]
    df = df[order + [c for c in df.columns if c not in order]]

    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        df.to_excel(xw, sheet_name="Operons", index=False)
        pc.to_excel(xw, sheet_name="Protein_complexes", index=False)

    n_can = int(df["is_canonical"].sum())
    n_prom = int((df["promoter_minus10_tier"].fillna("no_minus10") != "no_minus10").sum())
    n_term = int((df["has_terminator"] == True).sum()) if "has_terminator" in df else 0
    n_cplx = int((df["operon_complexes"] != "").sum())
    print(f"wrote {OUT.name}: {len(df)} operons "
          f"({n_can} canonical; {n_prom} with -10 box; {n_term} with terminator; "
          f"{n_cplx} touching a complex) + Protein_complexes sheet ({len(pc)} complexes)")


if __name__ == "__main__":
    main()
