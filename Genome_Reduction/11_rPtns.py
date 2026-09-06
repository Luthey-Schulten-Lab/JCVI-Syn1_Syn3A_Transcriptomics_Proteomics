#!/usr/bin/env python3
"""
Ribosomal proteins across the syn1 -> syn3A reduction: does the ABSOLUTE
proteome agree with the relative mRNA-pool story?

Question
--------
09/10 showed that the 21-gene ribosomal-protein operon (MMSYN1_0652-0672)
takes over a third of the syn3A coding mRNA pool, up from ~12% in syn1. That
statement is made in *shares* of a relative unit (TPM), because the cellular
RNA dry mass of neither organism was measured, so TPM cannot be converted into
molecules per cell.

The proteome does not have this limitation: absolute per-gene copy numbers were
derived for both organisms from a cellular dry mass and a protein mass fraction
(Methods; syn1 12.8 fg / 58.2%, syn3A 10.2 fg / 54.727%), giving ~127,833
proteins per syn1 cell and ~100,341 per syn3A cell. So we can ask the same
question in molecules: are there MORE ribosomal proteins per minimal cell?

What this script does
---------------------
1. Defines the ribosomal-protein (rPtn) gene set as the union of loci called a
   ribosomal protein in either annotation (gene name ^rp[slm][A-Z]?$ or product
   starting with "(30S|50S)? ribosomal protein"). Two curation calls, both kept
   as flagged rows but excluded from the paired core set:
     - MMSYN1_0298 "ribosomal protein L7A family" -> syn3A re-annotates it as an
       uncharacterized L7Ae-family protein, not a ribosomal subunit.
     - JCVISYN3A_0930 / _0932 rpmG (L33): syn3A-only paralogs with no syn1
       counterpart (see CLAUDE.md).
2. Joins, per locus, Illumina TPM (both organisms) with iPM and ABSOLUTE protein
   copies per cell (syn1 `ptn_copy_number`, syn3A `exp_ptn_cnt_2026`).
3. Reports three complementary views of the change, because they answer
   different questions:
     - absolute copies per cell        : syn3A / syn1, no normalization at all.
     - share of the FULL pool          : matches how the mRNA claim is phrased.
     - share of the RETAINED pool      : syn1 renormalized to loci syn3A kept,
       removing the mechanical inflation caused by deleting ~420 genes
       (CLAUDE.md "Deletion-corrected cross-organism comparisons"). This is the
       fair like-for-like share comparison, and it is applied to the mRNA side
       too so RNA and protein are compared on identical footing.
4. Converts copies into an implied ribosome count per cell (copies / ribosomal
   stoichiometry; 4 for L7/L12, 1 otherwise) as an internal consistency check.

Headline outcome
----------------
  - rPtn copies per cell are essentially FLAT: ~37.1k (syn1) vs ~38.7k (syn3A).
  - Their share of the full proteome rises 29.0% -> 38.7% (x1.33), but almost
    all of that is deletion dilution: against the RETAINED syn1 proteome the
    share is 37.3% -> 38.7% (x1.03), i.e. no real reallocation.
  - The 21-gene operon (20 rPtns + secY) behaves the same way: 19.2k -> 20.0k copies
    (x1.04), while its mRNA share goes 12.1% -> 33.2% of the coding pool (x2.7).
  - The mRNA side moves far more than the protein side even after the same
    deletion correction, so the transcript surge is NOT matched one-for-one by
    protein. The per-gene ratio (protein share FC / mRNA share FC) quantifies
    this gene by gene.

Caveats (carried into the report verbatim)
------------------------------------------
  - Ribosomal-protein iPM is the least reliable part of the proteomics: small,
    highly basic, lysine/arginine-rich proteins are over-digested by trypsin and
    their iBAQ is compressed. 09/10 already exclude rPtns from the iPM outlier
    and PTR analyses for this reason.
  - Neither absolute scale is measured in syn1 or syn3A: the syn1 dry mass and
    protein mass fraction are borrowed from the parent M. mycoides, the syn3A
    ones from the Breuer 2019 whole-cell reconstruction. The two organisms
    therefore sit on independently assumed scales, and the absolute copy ratio
    inherits that assumption (a 1.27x difference in total proteins per cell
    comes straight from those constants).
  - Shares (% of pool) are free of the dry-mass assumption; only the absolute
    copies depend on it. Read the share rows first.

Inputs
------
  ../Syn1_Syn3A_Transcriptomics/syn1_TPM.tsv (avg_sense_TPM)
  ../Syn1_Syn3A_Transcriptomics/Processed_TPM_Palsson/GSM6204176_3A.csv (Illumina_TPM;
      same source 09 uses for the syn3A Illumina side)
  ../Syn1_Syn3A_Transcriptomics/syn3A_TPM.tsv (our own strand-aware
      Illumina_sense_TPM, carried as a second column for transparency; r ~ 0.998)
  ../Syn1_Syn3A_Proteomics/syn1_proteome.tsv
  ../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv
  ../Syn1_Syn3A_Transcriptomics/syn3A_rna_abundances.tsv (syn3A-only absolute mRNA cross-check)
  ../Genomes_Input/syn3a_genome.gff3 (retained-locus set)
  ../Syn1_Operon/operons.candidate_blocks.tsv (operon assignment)

Outputs (rProtein_Absolute/)
----------------------------
  rProtein_omics_syn1_syn3A.csv   : per-locus TPM + iPM + absolute copies + changes,
                                    for every ribosomal-protein gene of both organisms
  rPtn_operon_protein_change.csv  : the 21 genes of OP_00341 (0652-0672), same columns,
                                    ordered along the transcript (rpsJ/0672 -> secY/0652)
  rProtein_Absolute.txt           : full report (aggregates, per-gene tables, caveats)
"""

import os
import re
import numpy as np
import pandas as pd

# ── Inputs ───────────────────────────────────────────────────────────────────
SYN1_TPM_CSV   = "../Syn1_Syn3A_Transcriptomics/syn1_TPM.tsv"
SYN3A_ILL_CSV  = "../Syn1_Syn3A_Transcriptomics/Processed_TPM_Palsson/GSM6204176_3A.csv"
SYN3A_OWN_TSV  = "../Syn1_Syn3A_Transcriptomics/syn3A_TPM.tsv"
SYN1_PTN_CSV   = "../Syn1_Syn3A_Proteomics/syn1_proteome.tsv"
SYN3A_PTN_CSV  = "../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv"
SYN3A_RNA_ABS  = "../Syn1_Syn3A_Transcriptomics/syn3A_rna_abundances.tsv"
SYN3A_GFF      = "../Genomes_Input/syn3a_genome.gff3"
OPERON_TSV     = "../Syn1_Operon/operons.candidate_blocks.tsv"

SYN1_TPM_COL   = "avg_sense_TPM"      # Illumina, merged across biological samples
SYN3A_TPM_COL  = "Illumina_TPM"       # Illumina (Sandberg et al.), as used by 09

# ── Outputs ──────────────────────────────────────────────────────────────────
OUTDIR      = "rProtein_Absolute"
OUT_ALL     = f"{OUTDIR}/rProtein_omics_syn1_syn3A.csv"
OUT_OPERON  = f"{OUTDIR}/rPtn_operon_protein_change.csv"
OUT_REPORT  = f"{OUTDIR}/rProtein_Absolute.txt"
os.makedirs(OUTDIR, exist_ok=True)

CODING_RNA_TYPES = {"mRNA", "pseudo"}
RPTN_OPERON      = set(range(652, 673))   # OP_00341, 21 genes, minus strand
L7L12_LOCUS      = 806                    # rplL/rplG: 4 copies per ribosome
EXCLUDE_FROM_CORE = {298}                 # L7Ae-family, not a ribosomal subunit
SYN3A_ONLY_PARALOGS = {930, 932}          # rpmG (L33) paralogs, no syn1 counterpart

rep: list[str] = []


def say(line: str = "") -> None:
    print(line)
    rep.append(line)


def _num(s: pd.Series) -> pd.Series:
    return s.str.extract(r"(\d+)$").astype(int)


def _is_rprotein(gene_name: pd.Series, product: pd.Series) -> pd.Series:
    """Ribosomal-protein call. The product pattern is start-anchored so that
    'Maturation protease for ribosomal protein L27' (JCVISYN3A_0500) does not
    qualify."""
    gn = gene_name.fillna("")
    gp = product.fillna("")
    return gn.str.match(r"^rp[slm][A-Z]?$") | gp.str.match(r"(?i)^(30S |50S )?ribosomal protein")


def _subunit(gene_name: str, product: str) -> str:
    """30S / 50S from the gene name (rps -> 30S, rpl|rpm -> 50S), else the product."""
    gn, gp = str(gene_name or ""), str(product or "")
    if gn.startswith("rps"):
        return "30S"
    if gn.startswith(("rpl", "rpm")):
        return "50S"
    if re.search(r"\b30S\b", gp) or re.search(r"protein S\d", gp):
        return "30S"
    if re.search(r"\b50S\b", gp) or re.search(r"protein L\d", gp):
        return "50S"
    return ""


def _load_syn3a_loci(path: str) -> set:
    """locus_num set of syn3A gene + pseudogene records (the retained-gene set)."""
    nums = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] in ("gene", "pseudogene"):
                m = re.search(r"locus_tag=JCVISYN3A_(\d+)", f[8])
                if m:
                    nums.add(int(m.group(1)))
    return nums


# ─────────────────────────────────────────────────────────────────────────────
# Load
# ─────────────────────────────────────────────────────────────────────────────
syn1_tpm = pd.read_csv(SYN1_TPM_CSV, sep="\t")
syn3a_tpm = pd.read_csv(SYN3A_ILL_CSV).rename(columns={"Geneid": "locus_tag"})
syn3a_own = pd.read_csv(SYN3A_OWN_TSV, sep="\t")
syn1_ptn = pd.read_csv(SYN1_PTN_CSV, sep="\t")
syn3a_ptn = pd.read_csv(SYN3A_PTN_CSV, sep="\t")
# The three syn3A pseudogenes carry MS intensity but no curated protein (no sequence,
# length, MW or localization), so they are kept out of the absolute protein pool. They
# remain in the transcriptome tables, where all 496 loci count.
syn3a_ptn = syn3a_ptn[syn3a_ptn["primary_function"] != "Pseudogene"].copy()
for _d in (syn1_tpm, syn3a_tpm, syn3a_own, syn1_ptn, syn3a_ptn):
    _d["locus_num"] = _num(_d["locus_tag"])

retained = _load_syn3a_loci(SYN3A_GFF)          # syn1 locus_nums that syn3A kept
locus_to_rnatype = dict(zip(syn1_tpm["locus_num"], syn1_tpm["rna_type"]))

# locus_num -> syn1 operon id(s); a gene can sit under more than one operon block,
# so the covering operons are joined with ','.
op_map = pd.read_csv(OPERON_TSV, sep="\t")
_ops: dict[int, list[str]] = {}
for _, r in op_map.iterrows():
    if not isinstance(r["sense_gene_loci"], str):
        continue
    for lt in r["sense_gene_loci"].split(","):
        lt = lt.strip()
        if lt:
            _ops.setdefault(int(lt.split("_")[-1]), []).append(str(r["operon_id"]))
locus_to_op = {k: ",".join(v) for k, v in _ops.items()}

# ── Ribosomal-protein gene set (union of both annotations) ───────────────────
rp1 = set(syn1_tpm.loc[_is_rprotein(syn1_tpm["gene_name"], syn1_tpm["gene_product"]), "locus_num"])
rp3 = set(syn3a_ptn.loc[_is_rprotein(syn3a_ptn["gene_name"], syn3a_ptn["gene_product"]), "locus_num"])
# rows = every ribosomal-protein locus PLUS the 4 non-rPtn members of the 21-gene
# operon (secY etc.), so the operon export is complete at 21 genes.
rp_union = rp1 | rp3
rp_all = sorted(rp_union | RPTN_OPERON)
core = rp_union - EXCLUDE_FROM_CORE - SYN3A_ONLY_PARALOGS          # paired, ribosome-forming
paired_core = sorted(n for n in core if n in retained)

# ─────────────────────────────────────────────────────────────────────────────
# Pool denominators.
#   full     : every gene of that organism
#   retained : syn1 restricted to loci syn3A kept (deletion-corrected); for syn3A
#              the retained pool IS the full pool.
# ─────────────────────────────────────────────────────────────────────────────
syn1_coding = syn1_tpm[syn1_tpm["rna_type"].isin(CODING_RNA_TYPES)]
TPM1_FULL = syn1_coding[SYN1_TPM_COL].sum()
TPM1_RET  = syn1_coding.loc[syn1_coding["locus_num"].isin(retained), SYN1_TPM_COL].sum()
TPM3_FULL = syn3a_tpm[SYN3A_TPM_COL].sum()

IPM1_FULL = syn1_ptn["iPM_mean"].sum()
IPM1_RET  = syn1_ptn.loc[syn1_ptn["locus_num"].isin(retained), "iPM_mean"].sum()
IPM3_FULL = syn3a_ptn["iPM_mean"].sum()

CN1_FULL = syn1_ptn["ptn_copy_number"].sum()
CN1_RET  = syn1_ptn.loc[syn1_ptn["locus_num"].isin(retained), "ptn_copy_number"].sum()
CN3_FULL = syn3a_ptn["exp_ptn_cnt_2026"].sum()

# ─────────────────────────────────────────────────────────────────────────────
# Build the per-locus table
# ─────────────────────────────────────────────────────────────────────────────
t1 = syn1_tpm.set_index("locus_num")
t3 = syn3a_tpm.set_index("locus_num")
o3 = syn3a_own.set_index("locus_num")
q1 = syn1_ptn.set_index("locus_num")
q3 = syn3a_ptn.set_index("locus_num")


def _get(df: pd.DataFrame, num: int, col: str):
    if num in df.index and col in df.columns:
        v = df.at[num, col]
        return v if not isinstance(v, pd.Series) else v.iloc[0]
    return np.nan


rows = []
for n in rp_all:
    gname = _get(q3, n, "gene_name")
    if not isinstance(gname, str) or not gname:
        gname = _get(q1, n, "gene_name")
    if not isinstance(gname, str) or not gname:
        gname = _get(t1, n, "gene_name")
    prod1 = _get(t1, n, "gene_product")
    prod3 = _get(q3, n, "gene_product")

    if n in SYN3A_ONLY_PARALOGS:
        pairing = "syn3A_only"
    elif n not in retained:
        pairing = "deleted_in_syn3A"
    else:
        pairing = "both"

    tpm1 = _get(t1, n, SYN1_TPM_COL)
    tpm3 = _get(t3, n, SYN3A_TPM_COL)
    tpm3_own = _get(o3, n, "Illumina_sense_TPM")
    ipm1, ipm3 = _get(q1, n, "iPM_mean"), _get(q3, n, "iPM_mean")
    cn1, cn3 = _get(q1, n, "ptn_copy_number"), _get(q3, n, "exp_ptn_cnt_2026")
    # stoichiometry is defined only for genuine ribosome-forming subunits; secY and the
    # L7Ae-family locus get NaN so they never enter the implied-ribosome estimate.
    stoich = (4 if n == L7L12_LOCUS else 1) if n in core else np.nan

    sh_tpm1_full = 100 * tpm1 / TPM1_FULL if pd.notna(tpm1) else np.nan
    sh_tpm1_ret  = 100 * tpm1 / TPM1_RET  if pd.notna(tpm1) else np.nan
    sh_tpm3      = 100 * tpm3 / TPM3_FULL if pd.notna(tpm3) else np.nan
    sh_cn1_full  = 100 * cn1 / CN1_FULL if pd.notna(cn1) else np.nan
    sh_cn1_ret   = 100 * cn1 / CN1_RET  if pd.notna(cn1) else np.nan
    sh_cn3       = 100 * cn3 / CN3_FULL if pd.notna(cn3) else np.nan

    mrna_fc_ret = sh_tpm3 / sh_tpm1_ret if (pd.notna(sh_tpm3) and pd.notna(sh_tpm1_ret) and sh_tpm1_ret > 0) else np.nan
    ptn_fc_ret  = sh_cn3 / sh_cn1_ret  if (pd.notna(sh_cn3) and pd.notna(sh_cn1_ret) and sh_cn1_ret > 0) else np.nan

    rows.append({
        "locus_syn1": f"MMSYN1_{n:04d}" if n in set(syn1_tpm["locus_num"]) else "",
        "locus_syn3a": f"JCVISYN3A_{n:04d}" if n in retained else "",
        "gene_name": gname if isinstance(gname, str) else "",
        "subunit": _subunit(gname, prod3 if isinstance(prod3, str) else prod1),
        "gene_product_syn1": prod1 if isinstance(prod1, str) else "",
        "gene_product_syn3a": prod3 if isinstance(prod3, str) else "",
        "is_rprotein": n in rp_union,
        "pairing": pairing,
        "in_rPtn_operon": n in RPTN_OPERON,
        "operon_syn1": locus_to_op.get(n, ""),
        "stoich_per_ribosome": stoich,
        # --- RNA (Illumina) ---
        "TPM_syn1": tpm1,
        "TPM_syn3a": tpm3,
        "TPM_syn3a_own": tpm3_own,
        "mRNApool_pct_syn1_full": sh_tpm1_full,
        "mRNApool_pct_syn1_retained": sh_tpm1_ret,
        "mRNApool_pct_syn3a": sh_tpm3,
        "mRNApool_pct_FC_retained": mrna_fc_ret,
        # --- protein (relative) ---
        "iPM_syn1": ipm1,
        "iPM_syn3a": ipm3,
        # --- protein (absolute copies per cell) ---
        "copies_syn1": cn1,
        "copies_syn3a": cn3,
        "copies_FC": cn3 / cn1 if (pd.notna(cn1) and pd.notna(cn3) and cn1 > 0) else np.nan,
        "copies_change": cn3 - cn1 if (pd.notna(cn1) and pd.notna(cn3)) else np.nan,
        "proteome_pct_syn1_full": sh_cn1_full,
        "proteome_pct_syn1_retained": sh_cn1_ret,
        "proteome_pct_syn3a": sh_cn3,
        "proteome_pct_FC_retained": ptn_fc_ret,
        # --- do the two layers agree? ---
        "protein_over_mRNA_shareFC": ptn_fc_ret / mrna_fc_ret if (pd.notna(ptn_fc_ret) and pd.notna(mrna_fc_ret) and mrna_fc_ret > 0) else np.nan,
        "ribosomes_implied_syn1": cn1 / stoich if (pd.notna(cn1) and pd.notna(stoich)) else np.nan,
        "ribosomes_implied_syn3a": cn3 / stoich if (pd.notna(cn3) and pd.notna(stoich)) else np.nan,
        "locus_num": n,
    })

tbl = pd.DataFrame(rows).sort_values("locus_num").reset_index(drop=True)
tbl["note"] = ""
tbl.loc[tbl["locus_num"] == 298, "note"] = (
    "syn1 calls it 'ribosomal protein L7A family'; syn3A re-annotates as uncharacterized "
    "L7Ae-family protein -> excluded from ribosomal-protein aggregates")
tbl.loc[tbl["locus_num"].isin(SYN3A_ONLY_PARALOGS), "note"] = (
    "syn3A-only rpmG (L33) paralog, no syn1 counterpart -> excluded from paired aggregates")

tbl.to_csv(OUT_ALL, index=False, float_format="%.6g")

# 21-gene operon table, ordered along the transcript (minus strand: 0672 -> 0652)
op_tbl = tbl[tbl["in_rPtn_operon"]].sort_values("locus_num", ascending=False).reset_index(drop=True)
op_tbl.insert(0, "transcript_order", np.arange(1, len(op_tbl) + 1))
op_tbl.to_csv(OUT_OPERON, index=False, float_format="%.6g")


# ─────────────────────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────────────────────
def _agg(nums, label):
    """Aggregate a locus set across both organisms in all three views."""
    s = tbl[tbl["locus_num"].isin(nums)]
    d = {
        "label": label,
        "n": len(s),
        "n_ptn_syn1": int(s["copies_syn1"].notna().sum()),
        "n_ptn_syn3a": int(s["copies_syn3a"].notna().sum()),
        "cn1": s["copies_syn1"].sum(skipna=True),
        "cn3": s["copies_syn3a"].sum(skipna=True),
        "tpm1": s["TPM_syn1"].sum(skipna=True),
        "tpm3": s["TPM_syn3a"].sum(skipna=True),
    }
    d["ptn_full1"] = 100 * d["cn1"] / CN1_FULL
    d["ptn_ret1"]  = 100 * d["cn1"] / CN1_RET
    d["ptn_3"]     = 100 * d["cn3"] / CN3_FULL
    d["rna_full1"] = 100 * d["tpm1"] / TPM1_FULL
    d["rna_ret1"]  = 100 * d["tpm1"] / TPM1_RET
    d["rna_3"]     = 100 * d["tpm3"] / TPM3_FULL
    return d


def _block(d):
    say(f"  {d['label']}  (n = {d['n']} genes; protein quantified in "
        f"{d['n_ptn_syn1']} syn1 / {d['n_ptn_syn3a']} syn3A)")
    say("                                        syn1        syn3A     FC (syn3A/syn1)")
    say(f"    protein, absolute copies/cell  {d['cn1']:11,.0f} {d['cn3']:11,.0f}   {d['cn3']/d['cn1']:8.2f}")
    say(f"    protein, % of FULL proteome    {d['ptn_full1']:11.2f} {d['ptn_3']:11.2f}   {d['ptn_3']/d['ptn_full1']:8.2f}")
    say(f"    protein, % of RETAINED proteome{d['ptn_ret1']:11.2f} {d['ptn_3']:11.2f}   {d['ptn_3']/d['ptn_ret1']:8.2f}   <- deletion-corrected")
    say(f"    mRNA,    % of FULL coding pool {d['rna_full1']:11.2f} {d['rna_3']:11.2f}   {d['rna_3']/d['rna_full1']:8.2f}")
    say(f"    mRNA,    % of RETAINED pool    {d['rna_ret1']:11.2f} {d['rna_3']:11.2f}   {d['rna_3']/d['rna_ret1']:8.2f}   <- deletion-corrected")
    say(f"    ==> deletion-corrected share change, protein / mRNA = "
        f"{(d['ptn_3']/d['ptn_ret1']) / (d['rna_3']/d['rna_ret1']):.2f}"
        "   (1.00 = protein follows mRNA)")
    say()


say("=" * 88)
say("RIBOSOMAL PROTEINS ACROSS THE syn1 -> syn3A REDUCTION")
say("absolute protein copies per cell + Illumina TPM")
say("=" * 88)
say()
say("Q: the 21-gene ribosomal-protein operon takes a third of the syn3A coding mRNA pool,")
say("   but TPM is relative and the cellular RNA dry mass was never measured. The proteome")
say("   HAS an absolute scale. Are there more ribosomal proteins per minimal cell?")
say()

say("-" * 88)
say("1. POOL SIZES (denominators)")
say("-" * 88)
say(f"  syn1  total protein copies / cell     : {CN1_FULL:12,.0f}   ({int(syn1_ptn['ptn_copy_number'].notna().sum())} proteins quantified)")
say(f"  syn1  copies on RETAINED loci only    : {CN1_RET:12,.0f}   ({100*CN1_RET/CN1_FULL:.1f}% of the syn1 proteome)")
say(f"  syn3A total protein copies / cell     : {CN3_FULL:12,.0f}   ({int(syn3a_ptn['exp_ptn_cnt_2026'].notna().sum())} proteins quantified)")
say(f"  syn3A / syn1 total proteins per cell  : {CN3_FULL/CN1_FULL:12.3f}")
say("     (this ratio is fixed by the assumed dry masses / protein mass fractions,")
say("      12.8 fg x 58.2% for syn1 vs 10.2 fg x 54.727% for syn3A -- it is NOT a measurement)")
say()
say(f"  syn1  coding TPM pool (full)          : {TPM1_FULL:12,.0f}")
say(f"  syn1  coding TPM pool (retained loci) : {TPM1_RET:12,.0f}   ({100*TPM1_RET/TPM1_FULL:.1f}% of the syn1 coding pool)")
say(f"  syn3A coding TPM pool                 : {TPM3_FULL:12,.0f}")
say()

say("-" * 88)
say("1b. DENOMINATOR AUDIT (what exactly is in each mRNA pool, and what the choice costs)")
say("-" * 88)
say("  Every pool share below is over CODING RNA only. The two organisms' pools are:")
_s1rt = syn1_tpm[syn1_tpm["rna_type"].isin(CODING_RNA_TYPES)].groupby("rna_type")[SYN1_TPM_COL].agg(["count", "sum"])
for rt, r in _s1rt.iterrows():
    say(f"    syn1  {rt:<16} {int(r['count']):4d} genes   {r['sum']:11,.0f} TPM")
say(f"    syn1  {'TOTAL':<16} {int(_s1rt['count'].sum()):4d} genes   {TPM1_FULL:11,.0f} TPM")
syn3a_tpm["_syn1_rt"] = syn3a_tpm["locus_num"].map(locus_to_rnatype).fillna("syn3A-only CDS")
_s3rt = syn3a_tpm.groupby("_syn1_rt")[SYN3A_TPM_COL].agg(["count", "sum"])
for rt, r in _s3rt.iterrows():
    say(f"    syn3A {rt:<16} {int(r['count']):4d} genes   {r['sum']:11,.0f} TPM   (syn1 annotation of the same locus)")
say(f"    syn3A {'TOTAL':<16} {int(_s3rt['count'].sum()):4d} genes   {TPM3_FULL:11,.0f} TPM")
say("  The syn3A source table (GSM6204176_3A.csv) is CDS-only by construction, so its whole")
say("  content IS the coding pool -- including the 3 CDS with no syn1 counterpart")
say("  (JCVISYN3A_0930/0931/0932; two of them are the rpmG/L33 paralogs).")
say()
say("  Fig 6 does not use one denominator throughout, so here is every convention side by side:")


def _conv_share(nums):
    """(syn1%, syn3A%) under each denominator convention."""
    n1 = syn1_tpm.loc[syn1_tpm["locus_num"].isin(nums), SYN1_TPM_COL].sum(skipna=True)
    n3 = syn3a_tpm.loc[syn3a_tpm["locus_num"].isin(nums), SYN3A_TPM_COL].sum(skipna=True)
    m_ret = syn1_tpm["locus_num"].isin(retained)
    d1_m   = syn1_tpm.loc[syn1_tpm["rna_type"] == "mRNA", SYN1_TPM_COL].sum()
    d1_mp  = TPM1_FULL
    d1_mpr = TPM1_RET
    d3_m   = syn3a_tpm.loc[syn3a_tpm["_syn1_rt"] == "mRNA", SYN3A_TPM_COL].sum()
    d3_mp  = syn3a_tpm.loc[syn3a_tpm["_syn1_rt"].isin(CODING_RNA_TYPES), SYN3A_TPM_COL].sum()
    d3_all = TPM3_FULL
    n3_m   = syn3a_tpm.loc[syn3a_tpm["locus_num"].isin(nums) & syn3a_tpm["_syn1_rt"].isin(CODING_RNA_TYPES),
                           SYN3A_TPM_COL].sum(skipna=True)
    return [
        ("R6 panel d as coded : syn1 mRNA-only, FULL",      100 * n1 / d1_m,   100 * n3_m / d3_m),
        ("mRNA + pseudo, FULL pool, syn1-typed",            100 * n1 / d1_mp,  100 * n3_m / d3_mp),
        ("Fig 6 panel a style : mRNA+pseudo, syn1 RETAINED", 100 * n1 / d1_mpr, 100 * n3_m / d3_mp),
        ("THIS SCRIPT : mRNA+pseudo / all syn3A CDS, FULL", 100 * n1 / d1_mp,  100 * n3 / d3_all),
        ("THIS SCRIPT : same, syn1 RETAINED",               100 * n1 / d1_mpr, 100 * n3 / d3_all),
    ]


for _lab, _nums in (("21-gene rPtn operon (0652-0672)", RPTN_OPERON),
                    ("all ribosomal proteins (paired core)", set(paired_core))):
    say(f"    {_lab}")
    say(f"      {'convention':<52}{'syn1':>7}{'syn3A':>8}{'FC':>7}")
    for cl, a, b in _conv_share(_nums):
        say(f"      {cl:<52}{a:7.2f}{b:8.2f}{b / a:7.2f}")
    say()
say("  Read: the pseudo genes cost ~0.4 pp on the syn1 side and the 3 syn3A-only CDS ~0.7 pp")
say("  on the syn3A side -- neither changes any conclusion. What DOES change the story is the")
say("  syn1 baseline: FULL vs RETAINED moves the operon's mRNA-share fold change from ~2.8x")
say("  to ~2.2x, because ~22% of syn1's coding pool sits on genes syn3A deleted. Fig 6 panel d")
say("  uses the FULL syn1 pool while panel a uses the RETAINED one, so the two panels of the")
say("  same figure are on different baselines.")
say()

say("-" * 88)
say("2. AGGREGATES")
say("-" * 88)
say()
_block(_agg(paired_core, "ALL ribosomal proteins (paired core set)"))
_block(_agg(RPTN_OPERON, "21-gene rPtn operon OP_00341 (0652-0672; 20 rPtns + secY)"))
_block(_agg(set(paired_core) - RPTN_OPERON, "ribosomal proteins OUTSIDE the 21-gene operon"))

# syn3A-only paralogs, reported separately so they are not silently dropped
extra = tbl[tbl["locus_num"].isin(SYN3A_ONLY_PARALOGS)]
if len(extra):
    say(f"  syn3A-only rpmG (L33) paralogs, excluded above: "
        f"{extra['copies_syn3a'].sum(skipna=True):,.0f} copies/cell "
        f"({100*extra['copies_syn3a'].sum(skipna=True)/CN3_FULL:.2f}% of the syn3A proteome)")
    say()

say("-" * 88)
say("3. IMPLIED RIBOSOMES PER CELL")
say("-" * 88)
say("  Each ribosomal protein sits at 1 copy per ribosome (L7/L12 at 4), so copies/stoichiometry")
say("  estimates the ribosome count independently from each protein. Spread across proteins")
say("  reflects the trypsin-digestion bias on small basic proteins, not real stoichiometry.")
pair = tbl[tbl["locus_num"].isin(paired_core)]
both = pair[pair["ribosomes_implied_syn1"].notna() & pair["ribosomes_implied_syn3a"].notna()]
for org, col in (("syn1 ", "ribosomes_implied_syn1"), ("syn3A", "ribosomes_implied_syn3a")):
    v = both[col]
    say(f"  {org}: median {v.median():8,.0f}   IQR {v.quantile(.25):,.0f}-{v.quantile(.75):,.0f}"
        f"   mean {v.mean():,.0f}   (n = {len(v)} proteins)")
say(f"  ratio of medians (syn3A / syn1)       : {both['ribosomes_implied_syn3a'].median()/both['ribosomes_implied_syn1'].median():.3f}")
say()

# syn3A-only absolute mRNA cross-check
if os.path.exists(SYN3A_RNA_ABS):
    ra = pd.read_csv(SYN3A_RNA_ABS, sep="\t")
    ra["locus_num"] = _num(ra["locus_tag"])
    rp_mrna = ra[ra["locus_num"].isin(paired_core)]["copies_per_cell"].sum()
    op_mrna = ra[ra["locus_num"].isin(RPTN_OPERON)]["copies_per_cell"].sum()
    tot_mrna = ra[ra["RNA Type"] == "mRNA"]["copies_per_cell"].sum()
    say("-" * 88)
    say("4. syn3A-ONLY CROSS-CHECK: absolute mRNA (Syn1_Syn3A_Transcriptomics/Calc_Abundances_syn3A.py)")
    say("-" * 88)
    say("  syn3A mRNA WAS converted to copies/cell there, using an RNA mass fraction that is")
    say("  itself borrowed (Breuer 2019) -- same class of assumption as the protein dry mass.")
    say("  No syn1 equivalent exists, so this is a within-syn3A sanity check only.")
    say(f"  syn3A total mRNA copies / cell            : {tot_mrna:10.1f}")
    say(f"  ribosomal-protein mRNA copies / cell      : {rp_mrna:10.1f}  ({100*rp_mrna/tot_mrna:.1f}% of the mRNA pool)")
    say(f"  21-gene operon mRNA copies / cell         : {op_mrna:10.1f}  ({100*op_mrna/tot_mrna:.1f}% of the mRNA pool)")
    rp_cn3 = tbl.loc[tbl["locus_num"].isin(paired_core), "copies_syn3a"].sum(skipna=True)
    say(f"  protein copies per rPtn mRNA (syn3A)      : {rp_cn3/rp_mrna:10.0f}")
    say(f"  ... vs whole-proteome protein per mRNA    : {CN3_FULL/tot_mrna:10.0f}")
    say()

say("-" * 88)
say("5. WHAT CARRIES THE AGGREGATE (largest single rPtn contributors)")
say("-" * 88)
_pc = tbl[tbl["locus_num"].isin(paired_core)]
for org, col, den in (("syn1  mRNA (retained pool)", "mRNApool_pct_syn1_retained", None),
                      ("syn3A mRNA", "mRNApool_pct_syn3a", None),
                      ("syn1  protein (retained pool)", "proteome_pct_syn1_retained", None),
                      ("syn3A protein", "proteome_pct_syn3a", None)):
    top = _pc.nlargest(5, col)[["locus_syn1", "gene_name", col]]
    items = "  ".join(f"{r['gene_name']}/{r['locus_syn1'][-4:]} {r[col]:.2f}%" for _, r in top.iterrows())
    say(f"  {org:<30}: {items}")
say()
say("  NOTE rpmE (L31, MMSYN1_0137) is a genuine single-gene outlier on the syn1 RNA side:")
say(f"  {_pc.loc[_pc['locus_num'] == 137, 'TPM_syn1'].iat[0]:,.0f} TPM in syn1 vs "
    f"{_pc.loc[_pc['locus_num'] == 137, 'TPM_syn3a'].iat[0]:,.0f} TPM in syn3A (our own syn3A Illumina")
say("  quantification agrees, 355 TPM), on a 279 bp gene whose neighbours are ordinary. It")
say("  alone carries several percent of the syn1 ribosomal-protein mRNA share, so the")
say("  'ALL ribosomal proteins' RNA aggregate is sensitive to it; the 21-gene operon block")
say("  and the protein aggregates are not.")
say()

say("-" * 88)
say("6. PER-GENE TABLE -- 21-gene rPtn operon OP_00341, transcript order (rpsJ/0672 -> secY/0652)")
say("-" * 88)
cols = ["transcript_order", "locus_syn1", "gene_name", "subunit", "is_rprotein",
        "TPM_syn1", "TPM_syn3a", "mRNApool_pct_syn1_retained", "mRNApool_pct_syn3a",
        "copies_syn1", "copies_syn3a", "copies_FC",
        "proteome_pct_syn1_retained", "proteome_pct_syn3a", "proteome_pct_FC_retained",
        "protein_over_mRNA_shareFC"]
say(op_tbl[cols].to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
say()

say("-" * 88)
say("7. PER-GENE TABLE -- all ribosomal proteins, ascending locus")
say("-" * 88)
cols_all = ["locus_syn1", "locus_syn3a", "gene_name", "subunit", "is_rprotein", "pairing", "in_rPtn_operon",
            "TPM_syn1", "TPM_syn3a", "copies_syn1", "copies_syn3a", "copies_FC",
            "proteome_pct_FC_retained", "mRNApool_pct_FC_retained", "protein_over_mRNA_shareFC"]
say(tbl[cols_all].to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
say()

say("-" * 88)
say("8. CAVEATS")
say("-" * 88)
say("  (a) Ribosomal-protein iPM is the least reliable slice of the proteomics. rPtns are")
say("      small, highly basic and lysine/arginine-rich, so trypsin over-digests them and")
say("      iBAQ-derived abundance is compressed. Scripts 09/10 already exclude rPtns from")
say("      the iPM outlier and PTR analyses for exactly this reason. Treat the TPM side as")
say("      the better-measured layer and this table as the corroborating check.")
say("  (b) Neither absolute scale is measured in syn1 or syn3A. syn1's dry mass (12.8 fg)")
say("      and protein mass fraction (58.2%) come from the parent M. mycoides; syn3A's")
say("      (10.2 fg, 54.727%) from the Breuer 2019 whole-cell reconstruction. The 1.27x")
say("      difference in total proteins per cell is a consequence of those two independently")
say("      assumed constants, so 'copies_FC' inherits that assumption directly.")
say("  (c) Shares (% of pool) do NOT depend on the dry-mass assumption -- only on the")
say("      within-organism iBAQ ranking. Read the share rows first; the absolute copies")
say("      are the same numbers rescaled by the assumed totals.")
say("  (d) Full-pool shares are inflated on the syn3A side simply because ~420 genes were")
say("      deleted and the same relative pool is divided among fewer genes. The")
say("      RETAINED-pool rows renormalize syn1 to the loci syn3A kept and are the fair")
say("      comparison; both the RNA and the protein layer are corrected identically here.")
say()

# ─────────────────────────────────────────────────────────────────────────────
# Figure: the ribosomal-protein landscape
#
# One column per ribosomal-protein gene (ascending locus, so the 21-gene operon
# forms a contiguous block), three stacked tracks:
#   top    gene panel  -- locus number / subunit-protein name (S.. or L..)
#   middle mRNA        -- syn1 / syn3A, single-hue TEAL, opacity = abundance
#   bottom protein     -- syn1 / syn3A, single-hue ORANGE, opacity = copies per cell
# Opacity is a sqrt ramp anchored at 0, shared across BOTH organisms within a
# track, so the two rows of a track are directly comparable. Genes with no
# measurement are drawn as an empty dashed-outlined cell.
#
# Two versions are written, differing only in the mRNA track:
#   ..._TPM.pdf    raw Illumina TPM, exactly as requested.
# The mRNA track is plotted as % of each organism's CODING TPM pool, not as raw
# TPM, because the two TPM tables are normalized over different feature sets:
# syn1's is normalized across ALL features, so its coding genes hold only 644k of
# the 1e6 (ncRNA 15.6%, tRNA 10.2%, rRNA 6.5%, tmRNA 2.4% take the rest), whereas
# the syn3A table used here (Sandberg et al., GSM6204176_3A.csv) contains ONLY the
# 458 CDS renormalized to 1e6 among themselves. That is a file-scope mismatch, not
# a library difference -- syn3A's Illumina still carries ~8.5% rRNA when we
# quantify it ourselves over all features. Raw TPM would therefore put a uniform
# ~1.55x inflation on the syn3A row; the coding-pool share removes it and is the
# unit panel a already reports.
# ─────────────────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import PowerNorm, ListedColormap, to_rgb
from matplotlib.cm import ScalarMappable

matplotlib.rcParams.update({
    "font.size": 7,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans", "Helvetica", "DejaVu Sans"],
    "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

# gene -> ribosomal-subunit protein designation, for the gene-panel labels
RP_DESIGNATION = {
    "rpsB": "S2", "rpsC": "S3", "rpsD": "S4", "rpsE": "S5", "rpsF": "S6", "rpsG": "S7",
    "rpsH": "S8", "rpsI": "S9", "rpsJ": "S10", "rpsK": "S11", "rpsL": "S12", "rpsM": "S13",
    "rpsN": "S14", "rpsO": "S15", "rpsP": "S16", "rpsQ": "S17", "rpsR": "S18", "rpsS": "S19",
    "rpsT": "S20", "rpsU": "S21",
    "rplA": "L1", "rplB": "L2", "rplC": "L3", "rplD": "L4", "rplE": "L5", "rplF": "L6",
    "rplG": "L7/L12", "rplL": "L7/L12", "rplI": "L9", "rplJ": "L10", "rplK": "L11",
    "rplM": "L13", "rplN": "L14", "rplO": "L15", "rplP": "L16", "rplQ": "L17", "rplR": "L18",
    "rplS": "L19", "rplT": "L20", "rplU": "L21", "rplV": "L22", "rplW": "L23", "rplX": "L24",
    "rpmA": "L27", "rpmB": "L28", "rpmC": "L29", "rpmE": "L31", "rpmF": "L32", "rpmG": "L33",
    "rpmH": "L34", "rpmI": "L35", "rpmJ": "L36",
}

# Track hues. Deliberately NOT the organism tokens (#3182bd syn1 / #c0392b syn3A):
# here colour encodes the LAYER and the two rows within a track encode the organism,
# so reusing the organism colours would read as a contradiction.
RNA_HUE   = "#00696e"   # mRNA track (teal)
CN_HUE    = "#d95f02"   # protein track (orange)
GENE_BG   = "#d9d9d9"   # gene-panel boxes, uniform (the operon carries no mark)
MISSING_EC = "#8c8c8c"                      # not-detected cells: dashed outline, no fill
MISSING_LS = (0, (1.4, 1.2))

# Opacity ramp. sqrt (PowerNorm gamma=0.5) anchored at 0 rather than log10: ribosomal
# protein abundances sit in a narrow band, and a log ramp compressed almost every cell
# into the same dark tone so up/down between the two organism rows was unreadable.
NORM_GAMMA = 0.5


def _alpha_cmap(hue: str, lo: float = 0.07):
    """Single-hue colormap where only the opacity varies (lo -> 1)."""
    r, g, b = to_rgb(hue)
    return ListedColormap([(r, g, b, lo + (1 - lo) * t) for t in np.linspace(0, 1, 256)])


def _cells(ax, values_by_row, norm, cmap, ncol):
    """One row of blocks per organism; NaN -> empty grey cell."""
    pad = 0.07
    for i, vals in enumerate(values_by_row):
        for j, v in enumerate(vals):
            if pd.isna(v) or v <= 0:
                ax.add_patch(Rectangle((j + pad, i + pad), 1 - 2 * pad, 1 - 2 * pad,
                                       facecolor="none", edgecolor=MISSING_EC, lw=0.35,
                                       linestyle=MISSING_LS))
            else:
                ax.add_patch(Rectangle((j + pad, i + pad), 1 - 2 * pad, 1 - 2 * pad,
                                       facecolor=cmap(norm(v)), lw=0))
    ax.set_xlim(0, ncol)
    ax.set_ylim(len(values_by_row), 0)
    ax.set_xticks([])
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["Syn1.0", "Syn3A"], fontsize=5)
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)


def rptn_landscape(metric: str, out_pdf: str) -> dict:
    """metric: 'TPM' (raw Illumina TPM) or 'share' (% of the coding TPM pool)."""
    d = tbl[tbl["is_rprotein"] & ~tbl["locus_num"].isin(EXCLUDE_FROM_CORE)].sort_values("locus_num")
    ncol = len(d)

    if metric == "TPM":
        rna = [d["TPM_syn1"].to_numpy(float), d["TPM_syn3a"].to_numpy(float)]
        rna_label, rna_unit = "mRNA", "Illumina TPM"
    else:
        rna = [d["mRNApool_pct_syn1_full"].to_numpy(float), d["mRNApool_pct_syn3a"].to_numpy(float)]
        rna_label, rna_unit = "mRNA", "% of coding TPM pool"
    cn = [d["copies_syn1"].to_numpy(float), d["copies_syn3a"].to_numpy(float)]

    W, H = 7.0, 7 / 3
    LEFT, RIGHT = 0.95, 0.10
    PW = W - LEFT - RIGHT

    # vertical layout, inches from the bottom edge (the canvas is only 7/3 in tall,
    # so every band is placed explicitly rather than by constrained_layout)
    y_cb, h_cb = 0.26, 0.095
    y_cn, h_cn = 0.56, 0.42
    y_tpm, h_tpm = 1.10, 0.42
    y_gene, h_gene = 1.58, 0.13
    h_lab = 0.50                      # rotated gene labels sit on top of the gene panel

    fig = plt.figure(figsize=(W, H))

    def _ax(x, y, w, h):
        return fig.add_axes([x / W, y / H, w / W, h / H])

    # --- gene panel -----------------------------------------------------------
    ag = _ax(LEFT, y_gene, PW, h_gene)
    for j, (_, r) in enumerate(d.iterrows()):
        ag.add_patch(Rectangle((j + 0.07, 0.06), 0.86, 0.88,
                               facecolor=GENE_BG, edgecolor="#7f7f7f", lw=0.2))
        name = RP_DESIGNATION.get(str(r["gene_name"]), str(r["gene_name"]))
        lab = f"{int(r['locus_num']):04d}/{name}"
        ag.text(j + 0.5, 1.30, lab, rotation=90, ha="center", va="bottom", fontsize=5,
                color="#888888" if r["pairing"] == "syn3A_only" else "#1a1a1a", clip_on=False)
    ag.set_xlim(0, ncol); ag.set_ylim(0, 1)
    ag.set_xticks([]); ag.set_yticks([])
    for sp in ag.spines.values():
        sp.set_visible(False)

    # --- the two data tracks --------------------------------------------------
    rna_v = np.concatenate([v[~np.isnan(v)] for v in rna])
    cn_v = np.concatenate([v[~np.isnan(v)] for v in cn])
    # Anchored at 0 so opacity reads as "how much", shared by the two organism rows.
    n_rna = PowerNorm(gamma=NORM_GAMMA, vmin=0, vmax=rna_v.max())
    n_cn = PowerNorm(gamma=NORM_GAMMA, vmin=0, vmax=cn_v.max())
    c_rna, c_cn = _alpha_cmap(RNA_HUE), _alpha_cmap(CN_HUE)

    at = _ax(LEFT, y_tpm, PW, h_tpm)
    _cells(at, rna, n_rna, c_rna, ncol)
    ac = _ax(LEFT, y_cn, PW, h_cn)
    _cells(ac, cn, n_cn, c_cn, ncol)

    fig.text(0.16 / W, (y_tpm + h_tpm / 2) / H, rna_label, rotation=90,
             ha="left", va="center", fontsize=6, color=RNA_HUE, fontweight="bold")
    fig.text(0.16 / W, (y_cn + h_cn / 2) / H, "Protein", rotation=90,
             ha="left", va="center", fontsize=6, color=CN_HUE, fontweight="bold")

    # --- colorbars ------------------------------------------------------------
    for x0, w, norm, cmap, lab in (
            (LEFT, 1.85, n_rna, c_rna, rna_unit),
            (LEFT + 2.75, 1.85, n_cn, c_cn, "protein copies per cell")):
        cax = _ax(x0, y_cb, w, h_cb)
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
        cb.outline.set_linewidth(0.3)
        cax.tick_params(labelsize=5, length=1.5, pad=1)
        cb.set_label(lab, fontsize=5, labelpad=1.5)

    # not-detected legend: one dashed swatch + label, bottom-right corner
    axl = _ax(5.90, y_cb, 0.13, h_cb)
    axl.add_patch(Rectangle((0.04, 0.10), 0.92, 0.80, facecolor="none",
                            edgecolor=MISSING_EC, lw=0.35, linestyle=MISSING_LS))
    axl.set_xlim(0, 1); axl.set_ylim(0, 1)
    axl.set_xticks([]); axl.set_yticks([])
    for sp in axl.spines.values():
        sp.set_visible(False)
    fig.text(6.07 / W, (y_cb + h_cb / 2) / H, "not detected",
             ha="left", va="center", fontsize=5, color="#4d4d4d")

    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    return {"n_genes": ncol,
            "rna_range": (float(rna_v.min()), float(rna_v.max())),
            "cn_range": (float(cn_v.min()), float(cn_v.max()))}


FIG_LANDSCAPE = f"{OUTDIR}/rPtn_landscape.pdf"
_m = rptn_landscape("share", FIG_LANDSCAPE)   # pass "TPM" instead for the raw-TPM variant

say("-" * 88)
say("9. FIGURE -- ribosomal-protein landscape (7 x 7/3 in)")
say("-" * 88)
say(f"  {_m['n_genes']} ribosomal-protein genes, ascending locus (MMSYN1_0298 dropped: the")
say("  L7Ae-family locus syn3A re-annotates as non-ribosomal). Three tracks: gene panel")
say("  (locus/subunit-protein name), mRNA (opacity = share of the coding TPM pool),")
say("  protein (opacity = copies per cell). Opacity is a sqrt ramp")
say("  anchored at 0 and shared between the two organism rows of a track, so the rows are")
say("  directly comparable; genes with no measurement are drawn as a dashed empty cell.")
say(f"  {FIG_LANDSCAPE}")
say(f"    mRNA track    : {_m['rna_range'][0]:.3f} - {_m['rna_range'][1]:.3f} % of the coding TPM pool")
say(f"    protein track : {_m['cn_range'][0]:,.1f} - {_m['cn_range'][1]:,.0f} copies per cell")
say("  The mRNA track is a coding-pool SHARE, not raw TPM, because the two TPM tables are")
say(f"  normalized over different feature sets: syn1's covers all features (coding = {TPM1_FULL:,.0f}")
say(f"  of 1e6; ncRNA/tRNA/rRNA/tmRNA take the rest) while the syn3A table is CDS-only,")
say(f"  renormalized to {TPM3_FULL:,.0f} among the 458 CDS. That is a file-scope mismatch, not a")
say("  library difference -- syn3A's Illumina still carries ~8.5% rRNA when quantified over")
say("  all features. Share is also the unit panel a reports.")
say()

say("-" * 88)
say("OUTPUT FILES")
say("-" * 88)
say(f"  {OUT_ALL}")
say(f"  {OUT_OPERON}")
say(f"  {OUT_REPORT}")
say(f"  {FIG_LANDSCAPE}")

with open(OUT_REPORT, "w") as fh:
    fh.write("\n".join(rep) + "\n")
