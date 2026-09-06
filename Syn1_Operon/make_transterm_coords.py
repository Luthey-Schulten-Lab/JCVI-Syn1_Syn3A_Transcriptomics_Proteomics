#!/usr/bin/env python
"""
make_transterm_coords.py
========================
Build the TransTermHP gene-coordinate file for the syn1 genome.

TransTermHP scores a hairpin by where it sits relative to the flanking genes
(H2T / T2T / H2H intergenic classes), so it needs the gene layout alongside the
FASTA. Its coordinate format is one tab-separated line per gene:

    <gene_name>  <start>  <end>  <chrom>

1-based inclusive, and for a minus-strand gene start > end (that is how
TransTermHP learns the strand). Gene names carry the "gene:" prefix of the GFF3
ID field, so the terminator table can be read against the annotation directly.

Input:   ../Genomes_Input/syn1.genes.gff3   (911 gene features)
Output:  syn1_genes.coords

Then predict terminators with:

    transterm -p "$CONDA_PREFIX/data/expterm.dat" \
        ../Genomes_Input/syn1_genome.fasta syn1_genes.coords > syn1_TransTermHP.txt

    python Syn1_Operon/make_transterm_coords.py
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
GFF3 = HERE / ".." / "Genomes_Input" / "syn1.genes.gff3"
OUT = HERE / "syn1_genes.coords"


def main():
    rows = []
    with open(GFF3) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] != "gene":
                continue
            ad = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            name = ad.get("ID") or "gene:" + ad.get("locus_tag", "")
            start, end = int(p[3]), int(p[4])
            if p[6] == "-":                      # minus strand: start > end
                start, end = end, start
            rows.append((name, start, end, p[0]))

    with open(OUT, "w") as fh:
        for name, start, end, chrom in rows:
            fh.write(f"{name}\t{start}\t{end}\t{chrom}\n")
    print(f"wrote {OUT.name}: {len(rows)} genes")


if __name__ == "__main__":
    main()
