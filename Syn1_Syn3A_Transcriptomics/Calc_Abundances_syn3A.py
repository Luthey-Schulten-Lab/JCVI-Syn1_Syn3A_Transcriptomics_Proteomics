#!/usr/bin/env python
"""
Calc_Abundances_syn3A.py
========================
Convert relative Syn3A transcriptomics (TPM) into ABSOLUTE RNA copies per cell,
for all RNA classes (mRNA, rRNA, tRNA, ncRNA).

Cleaned and re-homed from
  MinimalCell_Motif-Identification_RNAseq/Omic_Quantification/Calc_Abundances.ipynb
and rewired to this project's data. Moved here from Syn3A_Corr_RNA_Proteins/ on
2026-09-04: every input is transcriptomic (the syn3A TPM table sits in this very
folder), it correlates nothing with protein, and both of its consumers --
Syn1_Syn3A_Proteomics and Genome_Reduction/11_rPtns.py -- are downstream of this
stage, so it belongs at the end of the TPM pipeline rather than inside a later
one. See the reorganization note in REPRODUCE.md.

Data sources:
  - TPM from our own Syn3A reads (Syn1_Syn3A_Transcriptomics/syn3A_TPM.tsv),
    not the Sandberg export (they agree at r=0.998). mRNA uses Illumina sense TPM;
    rRNA/tRNA/ncRNA use ONT sense TPM, because Ribo-Zero Illumina sees no rRNA.
  - RNA type + gene sequences from Genomes_Input/syn3a.gb and syn3A_genome.fasta.

Method (mass balance, after Breuer et al. eLife 2019):
  cell dry weight (gDW) from radius/density/water-model
   -> RNA dry mass = gDW * RNA mass fraction
   -> partitioned by RNA type (rRNA/tRNA/mRNA/ncRNA), with Ribo-Zero rRNA depletion
   -> per type, total copies = type_mass / (TPM-weighted average molecular weight)
   -> per gene, copies/cell = total * TPM_gene / sum(TPM within type).

All scaling assumptions are named constants below; your planned Syn3A dry-mass
measurement can replace the literature values directly.

Outputs (this folder):
  syn3A_rna_abundances.tsv   per-gene: RNA type, Illumina/ONT TPM, MW, copies/cell, mRNA pool share %
  Calc_Abundances_syn3A.txt  summary log (gDW, per-type mass/total counts, top mRNAs)

Run in the Omics env (needs biopython, pandas).
"""
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight

PROJ = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
TPM_TSV = HERE / "syn3A_TPM.tsv"          # written by Gene_TPM.py, this folder
GENBANK = PROJ / "Genomes_Input/syn3a.gb"
FASTA = PROJ / "Genomes_Input/syn3A_genome.fasta"
OUT_TSV = HERE / "syn3A_rna_abundances.tsv"
OUT_TXT = HERE / "Calc_Abundances_syn3A.txt"

# ---------------- scaling assumptions (literature; replace with measured dry mass) ----------------
NA = 6.0221409e23
SYN3A_RADIUS_NM = 400 / 2          # 200 nm, Syn3A (Breuer et al. eLife 2019)
RHO = 1.1                          # g/ml, cell density (Bratbak 1984)
PTN_MASS_FRAC = 54.727             # % protein dry-mass fraction (Breuer et al. 2019), for the water model
WATER_VOL_PER_PROTEIN = 4.8        # microl/mg cellular protein (Leblanc 1979)
RNA_MASS_FRAC = 0.16274            # RNA dry-mass fraction (Breuer et al. 2019)
RNA_SPLIT = {"mRNA": 0.05, "ncRNA": 0.01, "rRNA": 0.80, "tRNA": 0.15}  # of RNA mass (rRNA fixed; others rescaled)
RRNA_DEPLETED = True               # Ribo-Zero library
DEPLETION_EFF = 0.95               # assumed rRNA depletion efficiency (C. Fields)
NONCODING_TPM = "ONT_sense_TPM"    # rRNA/tRNA/ncRNA per-gene distribution (Illumina depletes rRNA)
MRNA_TPM = "Illumina_sense_TPM"    # mRNA per-gene distribution (the quantitative track)

L = []
def p(s=""): L.append(s); print(s)


def calc_gDW(radius_nm):
    """Cell dry weight (g) from radius, after the Breuer 2019 water/density model."""
    water_vol_per_gDW = WATER_VOL_PER_PROTEIN * 1e3 * PTN_MASS_FRAC / 100.0   # microl/gDW
    water_mass_per_gDW = water_vol_per_gDW * 1e-6 * 1e3                       # g/g (rho_water = 1 kg/l)
    dry_fraction = 1.0 / (1.0 + water_mass_per_gDW)
    volume = 4.0 / 3.0 * np.pi * radius_nm ** 3 * 1e-27                       # m^3
    cell_mass = volume * RHO * 1e6                                           # g (1 m^3 = 1e6 ml)
    return cell_mass * dry_fraction


def rna_type_map():
    """locus_tag -> RNA type from the GenBank (CDS->mRNA; pseudo->mRNA; tmRNA->ncRNA)."""
    rec = next(SeqIO.parse(str(GENBANK), "genbank"))
    typemap = {"CDS": "mRNA", "rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "ncRNA", "ncRNA": "ncRNA"}
    rt = {}
    for f in rec.features:
        lt = f.qualifiers.get("locus_tag", [None])[0]
        if lt is None:
            continue
        if f.type in typemap:
            rt[lt] = typemap[f.type]
        elif f.type == "gene" and "pseudo" in f.qualifiers and lt not in rt:
            rt[lt] = "mRNA"   # pseudogene transcripts grouped with mRNA, as in the original
    return rt


def main():
    p("=" * 64); p("Syn3A absolute RNA abundances (copies per cell)"); p("=" * 64)

    df = pd.read_csv(TPM_TSV, sep="\t")
    rt = rna_type_map()
    df["RNA Type"] = df["locus_tag"].map(rt).fillna("mRNA")

    # per-gene transcript molecular weight from the genome sequence
    genome = str(next(SeqIO.parse(str(FASTA), "fasta")).seq).upper()
    def mw(r):
        sub = genome[int(r.start0):int(r.end0)]
        rna = (sub if r.strand == "+" else
               sub.translate(str.maketrans("ACGT", "TGCA"))[::-1]).replace("T", "U")
        return molecular_weight(rna, seq_type="RNA")
    df["MW_g_per_mol"] = df.apply(mw, axis=1)

    # --- mass budget ---
    gDW = calc_gDW(SYN3A_RADIUS_NM)
    rna_mass = RNA_MASS_FRAC * gDW
    rescale = (1 - RNA_SPLIT["rRNA"]) / (RNA_SPLIT["mRNA"] + RNA_SPLIT["ncRNA"] + RNA_SPLIT["tRNA"])
    mass = {"rRNA": rna_mass * RNA_SPLIT["rRNA"],
            "mRNA": rna_mass * RNA_SPLIT["mRNA"] * rescale,
            "tRNA": rna_mass * RNA_SPLIT["tRNA"] * rescale,
            "ncRNA": rna_mass * RNA_SPLIT["ncRNA"] * rescale}
    if RRNA_DEPLETED:
        mass["rRNA"] *= (1 - DEPLETION_EFF)
    p(f"cell radius {SYN3A_RADIUS_NM:.0f} nm, density {RHO} g/ml -> dry weight {gDW:.3e} g ({gDW*1e15:.2f} fg)")
    p(f"assumed RNA dry-mass fraction = {RNA_MASS_FRAC} of cell dry weight (Breuer et al. 2019)")
    p(f"assumed rRNA depletion fraction = {DEPLETION_EFF} (Ribo-Zero applied: {RRNA_DEPLETED})")
    split_eff = {"rRNA": RNA_SPLIT["rRNA"],
                 "tRNA": RNA_SPLIT["tRNA"] * rescale,
                 "mRNA": RNA_SPLIT["mRNA"] * rescale,
                 "ncRNA": RNA_SPLIT["ncRNA"] * rescale}
    p(f"RNA dry mass = {rna_mass*1e15:.3f} fg, partitioned among the four RNA classes:")
    for k in ["rRNA", "tRNA", "mRNA", "ncRNA"]:
        pre = rna_mass * split_eff[k] * 1e15   # fg, before any depletion
        if k == "rRNA" and RRNA_DEPLETED:
            p(f"  {k:5s} {split_eff[k]*100:6.2f}% of RNA -> {pre:8.4f} fg, "
              f"reduced to {mass[k]*1e15:.4f} fg after Ribo-Zero depletion (x{1-DEPLETION_EFF:.2f})")
        else:
            p(f"  {k:5s} {split_eff[k]*100:6.2f}% of RNA -> {mass[k]*1e15:8.4f} fg")

    # --- per type: total copies + per-gene distribution ---
    df["copies_per_cell"] = 0.0
    p("\nTotal copies per cell by RNA type:")
    for rtype in ["mRNA", "rRNA", "tRNA", "ncRNA"]:
        sub = df[df["RNA Type"] == rtype]
        tpmcol = MRNA_TPM if rtype == "mRNA" else NONCODING_TPM
        tpm = sub[tpmcol].clip(lower=0)
        if tpm.sum() <= 0:
            p(f"  {rtype:6s}: no {tpmcol} signal, skipped ({len(sub)} genes)")
            continue
        rel = tpm / tpm.sum()
        avg_mw_per_molecule = float((rel * sub["MW_g_per_mol"]).sum()) / NA
        total = mass[rtype] / avg_mw_per_molecule
        df.loc[sub.index, "copies_per_cell"] = total * rel
        p(f"  {rtype:6s}: {total:9.1f} total over {len(sub)} genes ({tpmcol})")

    # mRNA pool share (%) for the intuitive relative metric
    mrna = df["RNA Type"] == "mRNA"
    df["mRNA_pool_share_pct"] = np.nan
    df.loc[mrna, "mRNA_pool_share_pct"] = (
        100 * df.loc[mrna, MRNA_TPM] / df.loc[mrna, MRNA_TPM].sum())

    out = df[["locus_tag", "gene_name", "gene_product", "RNA Type",
              "Illumina_sense_TPM", "ONT_sense_TPM", "MW_g_per_mol",
              "copies_per_cell", "mRNA_pool_share_pct"]].copy()
    out.to_csv(OUT_TSV, sep="\t", index=False)

    p(f"\nSanity: total mRNA per cell = {df.loc[mrna,'copies_per_cell'].sum():.1f} "
      "(M. florum ~420 mRNA; Syn3A expected to be small)")
    p("\nTop 10 mRNAs by copies/cell:")
    top = df[mrna].nlargest(10, "copies_per_cell")
    for _, r in top.iterrows():
        p(f"  {r['locus_tag']} {str(r['gene_name'])[:8]:8s} {r['copies_per_cell']:7.2f} copies  "
          f"({r['mRNA_pool_share_pct']:.2f}% of mRNA pool)")
    OUT_TXT.write_text("\n".join(L) + "\n")
    print(f"\nwrote {OUT_TSV.name} ({len(out)} genes) and {OUT_TXT.name}")


if __name__ == "__main__":
    main()
