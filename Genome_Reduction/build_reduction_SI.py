#!/usr/bin/env python
"""
build_reduction_SI.py  --  assemble the Genome-Reduction Supplementary workbook.

BACKGROUND
----------
The genome-reduction pipeline (01->10) writes ~10 separate tables that the
Methods (M9, genome_reduction.tex) and Results (R5/R6) reference individually.
This script READS those already-computed tables (it recomputes nothing) and
consolidates them into a single multi-sheet Excel workbook for Supplementary
Information, mirroring the consolidated-script pattern used for R3's
`syn1_omics.xlsx`.

ALGORITHM (entity-level organisation)
-------------------------------------
The source tables sit at distinct granularities, so each becomes one sheet at
its own level; the per-gene structural (08) and omics (09) tables are MERGED
into one master table keyed by locus (all genes, deleted ones flagged):

  0. README                 sheet index + column glossary
  1. Deletions          (95)   05 deletion_junctions.tsv (+ per-deletion gene list from 02)
  2. Operon_classification (459) 04 operon_deletion_classification.tsv
  3. Coexpression      (pairs)  07 operon_pair_coexpression.tsv (cross-junction)
                                + 06 single_operon_pairs.tsv (operon-internal), unioned
  4. Gene_table        (all)   09 coding + non-coding RNA tables, LEFT-joined with
                                08 retained_gene_context.tsv (gene_impact_class, neighbours)
                                and curated Primary/Secondary/Tertiary function + essentiality;
                                a deleted_in_syn3A flag marks loci absent from Syn3A.
  5. Function_categories (15+44) 09 TPM_change_by_{secondary,tertiary}.tsv (stacked, `level` col)
  6. Summary_stats             09 deleted_gene_occupancy.txt + 10 macromolecule_complex_abundance.tsv

OUTPUT
------
  genome_reduction.xlsx   the SI workbook
  build_reduction_SI.txt  build log (per-sheet row/column counts)

Run from Genome_Reduction/ (paths are relative to it). Read-only on all inputs.
"""

import os
import re
import pandas as pd

# ----------------------------------------------------------------------------- paths
GR        = os.path.dirname(os.path.abspath(__file__))
CRP       = os.path.join(GR, "Compare_RNA_Protein")
OUT_XLSX  = os.path.join(GR, "genome_reduction.xlsx")
OUT_LOG   = os.path.join(GR, "build_reduction_SI.txt")

JUNCTIONS   = os.path.join(GR, "deletion_junction/deletion_junctions.tsv")
DEL_GENES   = os.path.join(GR, "aln/raw/syn1_deleted_genes.tsv")
OP_CLASS    = os.path.join(GR, "deletion_overlaid_operon/operon_deletion_classification.tsv")
PAIR_CROSS  = os.path.join(GR, "operon_pair_coexpression/operon_pair_coexpression.tsv")
PAIR_INTRA  = os.path.join(GR, "single_operon_coexpression/single_operon_pairs.tsv")
CODING      = os.path.join(CRP, "syn1_vs_syn3a_RNA_protein.tsv")
NONCODING   = os.path.join(CRP, "syn1_vs_syn3a_noncoding_RNA.tsv")
IMPACT_CTX  = os.path.join(GR, "delete_gene/retained_gene_context.tsv")
TPM_SEC     = os.path.join(CRP, "TPM_change_by_secondary.tsv")
TPM_TER     = os.path.join(CRP, "TPM_change_by_tertiary.tsv")
COMPLEX     = os.path.join(CRP, "macromolecule_complex_abundance.tsv")
OCCUPANCY   = os.path.join(CRP, "deleted_gene_occupancy.txt")
FUNC_TSV    = os.path.join(GR, "../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv")

log_lines = []
def log(msg=""):
    print(msg)
    log_lines.append(msg)

def locus_num(tag):
    """MMSYN1_0350 / JCVISYN3A_0350 -> '0350' (numeric suffix, preserved across organisms)."""
    if pd.isna(tag):
        return None
    m = re.search(r"_(\d+)\s*$", str(tag))
    return m.group(1) if m else None


# ============================================================ sheet 1: Deletions (95)
def build_deletions():
    j = pd.read_csv(JUNCTIONS, sep="\t")
    # per-deletion deleted-gene list, parsed from 02's per-(deletion x gene) table
    dg = pd.read_csv(DEL_GENES, sep="\t")
    dg["locus_tag"] = dg["gff_attributes"].str.extract(r"locus_tag=([A-Za-z0-9_]+)")
    glist = (dg.dropna(subset=["locus_tag"])
               .groupby(["del_start0", "del_end"])["locus_tag"]
               .apply(lambda s: ";".join(dict.fromkeys(s)))   # order-preserving unique
               .reset_index()
               .rename(columns={"locus_tag": "genes_in_deletion_interval"}))
    out = j.merge(glist, left_on=["syn1_del_s", "syn1_del_e"],
                  right_on=["del_start0", "del_end"], how="left").drop(
                  columns=["del_start0", "del_end"], errors="ignore")
    out["genes_in_deletion_interval"] = out["genes_in_deletion_interval"].fillna("")
    out["n_genes_in_interval"] = out["genes_in_deletion_interval"].apply(
        lambda s: 0 if not s else len(s.split(";")))
    cols = ["scar_id", "syn1_del_s", "syn1_del_e", "syn1_del_len",
            "n_deleted_genes", "n_genes_in_interval", "genes_in_deletion_interval",
            "syn3A_junction", "left_gene", "right_gene",
            "operon_L_id", "operon_L_strand", "operon_R_id", "operon_R_strand",
            "strand_relationship", "junction_type",
            "operon_L_facing_reg", "operon_L_reg_lost",
            "operon_R_facing_reg", "operon_R_reg_lost",
            "operon_L_gene_deletion_pattern", "operon_L_overlap_class",
            "operon_R_gene_deletion_pattern", "operon_R_overlap_class"]
    cols = [c for c in cols if c in out.columns]
    return out[cols]


# ================================================ sheet 2: Operon_classification (459)
def build_operon_classification():
    return pd.read_csv(OP_CLASS, sep="\t")


# ===================================================== sheet 3: Coexpression (pairs)
def build_coexpression():
    cross = pd.read_csv(PAIR_CROSS, sep="\t")
    intra = pd.read_csv(PAIR_INTRA, sep="\t")
    common = ["syn3A_intergenic_bp", "n_spanning_reads", "n_bridging_reads",
              "ont_depth_a", "ont_depth_b", "ill_depth_a", "ill_depth_b",
              "ill_gap_depth", "bridge_threshold", "gap_depth_threshold",
              "bridge_pass", "gap_depth_pass",
              "pair_preserved_strict", "pair_preserved_loose"]
    c = cross.rename(columns={"left_gene": "gene_a", "right_gene": "gene_b"})
    c["coexpression_level"] = "cross_junction"
    c["context"] = c["operon_L_id"].astype(str) + "->" + c["operon_R_id"].astype(str)
    i = intra.rename(columns={"gene_a_locusNum": "gene_a", "gene_b_locusNum": "gene_b"})
    i["coexpression_level"] = "operon_internal"
    i["context"] = i["operon_id"].astype(str)
    i["strand_relationship"] = "intra_operon"
    i["junction_type"] = i["category"]
    head = ["coexpression_level", "context", "gene_a", "gene_b",
            "strand_relationship", "junction_type"]
    keep = head + common
    c = c[[col for col in keep if col in c.columns]]
    i = i[[col for col in keep if col in i.columns]]
    out = pd.concat([c, i], ignore_index=True)
    return out[[col for col in keep if col in out.columns]]


# ================================================== sheet 4: Gene_table (master, all)
def build_gene_table():
    coding = pd.read_csv(CODING, sep="\t")
    noncod = pd.read_csv(NONCODING, sep="\t")
    genes  = pd.concat([coding, noncod], ignore_index=True, sort=False)
    genes["locus_num"] = genes["locus_syn1"].map(locus_num)

    # deleted set (authoritative): loci overlapped by a deletion in 02's table
    dg = pd.read_csv(DEL_GENES, sep="\t")
    dg["locus_tag"] = dg["gff_attributes"].str.extract(r"locus_tag=([A-Za-z0-9_]+)")
    deleted = set(dg["locus_tag"].dropna())
    genes["deleted_in_syn3A"] = genes["locus_syn1"].isin(deleted)

    # 08 structural context (retained genes only)
    ctx = pd.read_csv(IMPACT_CTX, sep="\t")
    ctx["locus_num"] = ctx["locus_tag"].map(locus_num)
    ctx_cols = ["locus_num", "gene_impact_class",
                "syn1_upstream_locus", "syn1_downstream_locus",
                "syn3A_upstream_locus", "syn3A_downstream_locus",
                "unaltered_cw_bps", "unaltered_ccw_bps",
                "cw_context_changed", "ccw_context_changed"]
    genes = genes.merge(ctx[[c for c in ctx_cols if c in ctx.columns]],
                        on="locus_num", how="left")

    # curated function + essentiality (Syn3A annotation; by preserved locus suffix)
    # snake_case table -> the Title-Case names used below
    fn = pd.read_csv(FUNC_TSV, sep="\t").rename(columns={
        "locus_tag": "Locus Tag", "essentiality": "Essentiality",
        "primary_function": "Primary Function",
        "secondary_function": "Secondary Function",
        "tertiary_function": "Tertiary Function"})
    fn["locus_num"] = fn["Locus Tag"].map(locus_num)
    fn_cols = ["locus_num", "Primary Function", "Secondary Function",
               "Tertiary Function", "Essentiality"]
    genes = genes.merge(fn[fn_cols].drop_duplicates("locus_num"),
                        on="locus_num", how="left")

    order = ["locus_syn1", "locus_syn3a", "gene_name", "gene_product", "rna_type",
             "deleted_in_syn3A",
             "Primary Function", "Secondary Function", "Tertiary Function", "Essentiality",
             "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change", "TPM_abs_change",
             "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change",
             "PTR_syn1", "PTR_syn3a", "PTR_fold_change",
             "gene_impact_class", "sense_covering_ops",
             "syn1_upstream_locus", "syn1_downstream_locus",
             "syn3A_upstream_locus", "syn3A_downstream_locus",
             "unaltered_cw_bps", "unaltered_ccw_bps",
             "cw_context_changed", "ccw_context_changed"]
    order = [c for c in order if c in genes.columns]
    genes = genes[order].sort_values("locus_syn1").reset_index(drop=True)
    return genes


# ============================================= sheet 5: Function_categories (15 + 44)
def build_function_categories():
    sec = pd.read_csv(TPM_SEC, sep="\t")
    ter = pd.read_csv(TPM_TER, sep="\t")
    sec.insert(0, "level", "secondary")
    ter.insert(0, "level", "tertiary")
    out = pd.concat([sec, ter], ignore_index=True, sort=False)
    front = ["level", "Primary Function", "Secondary Function", "category",
             "n_genes", "n_detected_FC",
             "syn1_pool_share_pct", "syn3a_pool_share_pct",
             "pool_share_change", "pool_share_FC",
             "median_TPM_FC_corr", "median_TPM_abs_corr", "mwu_p_vs_rest"]
    front = [c for c in front if c in out.columns]
    rest = [c for c in out.columns if c not in front]
    return out[front + rest]


# ================================================ sheet 6: Summary_stats (occupancy + complexes)
def build_summary_stats():
    rows = []
    with open(OCCUPANCY) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if not s or s.lstrip().startswith("#"):
                continue
            if ":" in s and not s.strip().startswith("MMSYN1"):
                k, _, v = s.partition(":")
                rows.append({"section": "occupancy", "metric": k.strip(), "value": v.strip()})
    occ = pd.DataFrame(rows)
    cx = pd.read_csv(COMPLEX, sep="\t")
    return occ, cx


# ====================================================================== README sheet
def build_readme(counts):
    idx = pd.DataFrame([
        ["1. Deletions", "One row per Syn1->Syn3A deletion: coordinates, flanking operons, "
         "strand relationship and junction type, regulators lost. n_deleted_genes counts the "
         "sense genes excised between the retained flanks; genes_in_deletion_interval "
         "(n_genes_in_interval) lists every annotated feature the raw interval overlaps, so it "
         "can be larger (it also counts watermarks/pseudogenes).",
         counts["Deletions"], "05_deletion_junction.py (+02_analyze.py)", "R5 L5.1, L5.3"],
        ["2. Operon_classification", "One row per Syn1 operon: how each deletion hit it "
         "(span-level overlap_class) and which sense genes were lost (gene_deletion_pattern).",
         counts["Operon_classification"], "04_deletion_overlaid_operon.py", "R5 L5.2"],
        ["3. Coexpression", "One row per gene pair tested for co-transcription in Syn3A reads "
         "(cross-junction pairs from 07 and operon-internal pairs from 06), with ONT "
         "spanning/bridging counts, Illumina gap depth, and pass flags.",
         counts["Coexpression"], "07_operon_pair_coexpression.py, 06_single_operon_coexpression.py",
         "R5 L5.3"],
        ["4. Gene_table", "Master per-gene table (all genes): relative transcript (relTPM) and "
         "protein (relIPM) abundance in each organism with fold/absolute changes and PTR, the "
         "structural gene_impact_class and neighbours, curated function and essentiality, and a "
         "deleted_in_syn3A flag.",
         counts["Gene_table"], "09_Compare_RNA_Protein.py + 08_delete_gene.py + Syn3A annotation",
         "R5 L5.4, R6 L6.1-L6.4"],
        ["5. Function_categories", "Per curated function category (Secondary and Tertiary): "
         "retained-pool mRNA share in each organism, share change, median TPM fold change "
         "(deletion-corrected) and Mann-Whitney p vs the rest.",
         counts["Function_categories"], "09_Compare_RNA_Protein.py", "R6 L6.2"],
        ["6. Summary_stats", "Deleted-gene occupancy of the Syn1 transcriptome/proteome and "
         "limiting-subunit abundance of RNA polymerase and the degradosome.",
         counts["Summary_stats"], "09_Compare_RNA_Protein.py, 10_Compare_Ptn.py", "R6 L6.1, L6.3"],
    ], columns=["Sheet", "Contents", "Rows", "Source script", "Supports"])

    glossary = pd.DataFrame([
        ["strand_relationship", "tandem (co-oriented, fusion-capable) | convergent (terminators face) "
         "| divergent (promoters face) | intra_operon (deletion internal to one operon)"],
        ["junction_type", "tandem junctions only: fusion (both facing regulators lost) | decapitation "
         "(downstream promoter lost) | readthrough_extension (upstream terminator lost) | "
         "clean_excision (both kept; whole operon(s) removed between intact neighbours)"],
        ["overlap_class", "span-level truncation of an operon: fully_deleted | 5'/3'_truncation_gene | "
         "5'/3'_truncation_UTR | intra_truncated | intact (multi:* = several hits)"],
        ["gene_deletion_pattern", "gene-level loss in an operon: all_deleted | leading_deleted | "
         "lagging_deleted | intra_deleted | intact"],
        ["gene_impact_class", "promoter-source change for a retained gene (precedence): promoter_lost > "
         "promoter_disconnected > new_promoter_fusion > readthrough_exposed > "
         "promoter_proximity_changed > context_only > unaffected"],
        ["relTPM / relIPM", "transcript (TPM, Illumina) / protein (iPM) abundance mean-normalised within "
         "organism+subset, so an average gene = 1 and an unchanged gene gives fold change ~1 despite the "
         "reduced gene count"],
        ["*_fold_change / *_abs_change", "Syn3A relative value / Syn1 relative value, and their difference"],
        ["PTR", "protein-to-transcript ratio relIPM/relTPM (steady-state translation-efficiency proxy, "
         "not Ribo-seq); PTR_fold_change = iPM_FC / TPM_FC"],
        ["deleted_in_syn3A", "TRUE if the Syn1 locus is overlapped by a deletion / absent from Syn3A"],
        ["locus correspondence", "MMSYN1_NNNN (Syn1) <-> JCVISYN3A_NNNN (Syn3A); numeric suffix preserved"],
    ], columns=["Term", "Meaning"])
    return idx, glossary


# ======================================================================= write workbook
def main():
    log("=" * 70)
    log("Genome-Reduction SI workbook build")
    log("=" * 70)

    deletions   = build_deletions()
    operon_cls  = build_operon_classification()
    coexpr      = build_coexpression()
    gene_table  = build_gene_table()
    func_cats   = build_function_categories()
    occ, cx     = build_summary_stats()

    counts = {
        "Deletions": len(deletions),
        "Operon_classification": len(operon_cls),
        "Coexpression": len(coexpr),
        "Gene_table": len(gene_table),
        "Function_categories": len(func_cats),
        "Summary_stats": len(occ) + len(cx),
    }
    idx, glossary = build_readme(counts)

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
        idx.to_excel(xw, sheet_name="README", index=False, startrow=0)
        glossary.to_excel(xw, sheet_name="README", index=False, startrow=len(idx) + 3)
        deletions.to_excel(xw, sheet_name="Deletions", index=False)
        operon_cls.to_excel(xw, sheet_name="Operon_classification", index=False)
        coexpr.to_excel(xw, sheet_name="Coexpression", index=False)
        gene_table.to_excel(xw, sheet_name="Gene_table", index=False)
        func_cats.to_excel(xw, sheet_name="Function_categories", index=False)
        occ.to_excel(xw, sheet_name="Summary_stats", index=False, startrow=0)
        cx.to_excel(xw, sheet_name="Summary_stats", index=False, startrow=len(occ) + 3)

    log(f"\nWrote {OUT_XLSX}")
    log("\nSheets:")
    log(f"  README                 index ({len(idx)} sheets) + glossary ({len(glossary)} terms)")
    log(f"  Deletions              {counts['Deletions']} rows x {deletions.shape[1]} cols")
    log(f"  Operon_classification  {counts['Operon_classification']} rows x {operon_cls.shape[1]} cols")
    log(f"  Coexpression           {counts['Coexpression']} rows x {coexpr.shape[1]} cols "
        f"({(coexpr['coexpression_level']=='cross_junction').sum()} cross-junction + "
        f"{(coexpr['coexpression_level']=='operon_internal').sum()} operon-internal)")
    log(f"  Gene_table             {counts['Gene_table']} rows x {gene_table.shape[1]} cols "
        f"({int(gene_table['deleted_in_syn3A'].sum())} deleted, "
        f"{int((~gene_table['deleted_in_syn3A']).sum())} retained)")
    log(f"  Function_categories    {counts['Function_categories']} rows x {func_cats.shape[1]} cols "
        f"({(func_cats['level']=='secondary').sum()} secondary + "
        f"{(func_cats['level']=='tertiary').sum()} tertiary)")
    log(f"  Summary_stats          {len(occ)} occupancy metrics + {len(cx)} complexes")

    with open(OUT_LOG, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    print(f"\nLog -> {OUT_LOG}")


if __name__ == "__main__":
    main()
