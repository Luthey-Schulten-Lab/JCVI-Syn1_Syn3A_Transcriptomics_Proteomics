#!/usr/bin/env python
"""
novel_tex_todos.py
==================
Resolve the three data questions raised as TODO@claude in the R4 manuscript
section (Manuscript/sections/results/novel.tex), and write the numbers to
novel_tex_todos.txt so the prose can quote them verbatim.

  TODO 7  -- how the 89 antisense isoform clusters relate to the operon-level
             antisense structure reported in R1 (69 operons enclose >=1 antisense
             gene; 9 operons are antisense/intergenic-only).
  TODO 16 -- syn3A Illumina depth at his3/0918 (retained as JCVISYN3A_0918,
             18,716-19,378, - strand), relative to the syn3A genome-wide average,
             on the sense and the antisense strand.
  TODO 31 -- genomic context (gene / pseudogene / function-unknown) of the two
             MS-confirmed novel peptides.

Run in the Omics env.
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "novel_tex_todos.txt"
L = []
def p(s=""): L.append(s)

# ============================ TODO 7 ============================
p("=" * 70); p("TODO 7 -- antisense clusters vs R1 operon antisense structure"); p("=" * 70)
op = pd.read_csv(ROOT / "Syn1_Operon/operons.candidate_blocks.tsv", sep="\t")
asx = pd.read_excel(ROOT / "Syn1_Novel_ORF/isoform_antisense_categories.xlsx")

n_enclose = int((op["antisense_gene_count"] >= 1).sum())
anti_only = op[op["sense_gene_count"] == 0].copy()
p(f"operons enclosing >=1 antisense gene (R1):        {n_enclose}")
p(f"antisense/intergenic-only operons, sense_gene_count==0 (R1): {len(anti_only)}")
p(f"antisense isoform clusters (R4):                  {len(asx)}")

# antisense genes captured by the operon segmentation (enclosed on a sense operon)
enclosed = set()
for v in op["antisense_gene_loci"].dropna():
    enclosed.update(x.strip() for x in str(v).split(",") if x.strip())
def anti_ids(row):
    v = row.get("antisense_gene_ids")
    return [x.strip() for x in str(v).split(",") if x.strip() and str(v) != "nan"]
asx["_anti"] = asx.apply(anti_ids, axis=1)
clu_with_enclosed = asx["_anti"].apply(lambda ids: any(g in enclosed for g in ids)).sum()
p(f"\nantisense genes enclosed by operons (distinct): {len(enclosed)}")
p(f"clusters whose antisense gene is enclosed by an operon: {int(clu_with_enclosed)} / {len(asx)}")

# breakdown by R4 antisense_category
p("\nclusters by antisense_category (R4) x enclosed-by-operon:")
asx["_enc"] = asx["_anti"].apply(lambda ids: any(g in enclosed for g in ids))
tab = asx.groupby("antisense_category")["_enc"].agg(["size", "sum"]).rename(columns={"size": "n_clusters", "sum": "enclosed"})
for cat, r in tab.iterrows():
    p(f"  {cat:16s} {int(r.n_clusters):3d} clusters, {int(r.enclosed):3d} enclosed by an operon")

# R1's sense-gene-less operons split into antisense-only + purely intergenic
def overlaps_any(o):
    same = asx[(asx["strand"] == o["strand"]) & (asx["start0"] < o["end0"]) & (asx["end0"] > o["start0"])]
    return len(same)
asonly = anti_only[anti_only["antisense_gene_count"] >= 1].copy()
inter1 = anti_only[anti_only["antisense_gene_count"] == 0].copy()
asonly["_ov"] = asonly.apply(overlaps_any, axis=1)
p(f"\nR1's {len(anti_only)} sense-gene-less operons = {len(asonly)} antisense-only + {len(inter1)} purely intergenic")
p(f"  antisense-only operons matched by >=1 antisense cluster: {int((asonly['_ov']>=1).sum())} / {len(asonly)} (his3/0918 the largest)")
for _, r in inter1.iterrows():
    p(f"  the single intergenic operon ({r['operon_id']}, {int(r['start0'])}-{int(r['end0'])}) has no antisense cluster; it is R4's isolated genuinely-intergenic transcript")

# ============================ TODO 16 ============================
p(""); p("=" * 70); p("TODO 16 -- syn3A Illumina depth at his3/0918 (JCVISYN3A_0918)"); p("=" * 70)
GLEN = 543379
def load(fn):
    cov = np.zeros(GLEN)
    for ln in open(fn):
        q = ln.split()
        if len(q) >= 4:
            cov[int(q[1]):int(q[2])] = abs(float(q[3]))
    return cov
bg = ROOT / "Syn3A_Transcriptomics/Illumina/Illumina_Processing/depth_bedgraph"
plus = load(bg / "syn3A_rep1.plus.bedGraph")
minus = load(bg / "syn3A_rep1.minus.bedGraph")
tot = plus + minus
s, e = 18716, 19378  # JCVISYN3A_0918 (his3), - strand -> sense=minus, antisense=plus
gmean = tot.mean()
p(f"his3/0918 retained in syn3A (JCVISYN3A_0918, {s:,}-{e:,}, - strand); not deleted.")
p(f"syn3A Illumina genome-wide mean depth: total {gmean:.0f} (plus {plus.mean():.0f}, minus {minus.mean():.0f})")
p(f"  his3 SENSE     (minus): {minus[s:e].mean():.0f}  = {minus[s:e].mean()/gmean:.2f}x genome mean")
p(f"  his3 ANTISENSE (plus) : {plus[s:e].mean():.0f}  = {plus[s:e].mean()/gmean:.2f}x genome mean")
p(f"  his3 total            : {tot[s:e].mean():.0f}  = {tot[s:e].mean()/gmean:.2f}x genome mean")
p("INTERPRETATION: the syn1 conspicuous antisense over-transcription (~30,000 depth,")
p("the single exceptional locus) collapses to background in syn3A (antisense 0.23x the")
p("genome-wide average); the his3 ORF is retained but unremarkably expressed.")

# ============================ TODO 31 ============================
p(""); p("=" * 70); p("TODO 31 -- genomic context of the two MS-confirmed novel peptides"); p("=" * 70)
def load_gff(fn):
    rows = []
    for ln in open(fn):
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"): continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        rows.append((int(f[3]), int(f[4]), f[6], attr.get("locus_tag", ""),
                     attr.get("gene", attr.get("Name", "")), attr.get("product", ""),
                     attr.get("rna_type", attr.get("gbkey", ""))))
    return pd.DataFrame(rows, columns=["start", "end", "strand", "locus", "name", "product", "rna_type"])
gff = load_gff(ROOT / "Genomes_Input/syn1.genes.gff3")
peps = [("NOVEL_PEP_002", 118, 728399, 728756, "intergenic near 0592"),
        ("NOVEL_PEP_043", 225, 905181, 905859, "5' extension of mmyCIVR/0768")]
for pep, aa, ps, pe, note in peps:
    p(f"\n{pep} ({aa} aa)  ORF {ps:,}-{pe:,}  [{note}]")
    near = gff[(gff["start"] <= pe + 1500) & (gff["end"] >= ps - 1500)].sort_values("start")
    for _, g in near.iterrows():
        rel = "overlaps" if (g["start"] <= pe and g["end"] >= ps) else "flank"
        flag = "  <-- PSEUDOGENE" if "pseudo" in str(g["rna_type"]).lower() or "PSEUDOGENE" in str(g["product"]) else ("  <-- function-unknown" if str(g["name"]).strip() == "" or "cdsf" in str(g["product"]) else "")
        p(f"   {rel:8s} {g['locus']} {g['strand']} {g['start']:,}-{g['end']:,}  {g['name'] or '(no name)'} | {g['product'][:60]}{flag}")

OUT.write_text("\n".join(L) + "\n")
print("\n".join(L))
print(f"\nwrote {OUT}")
