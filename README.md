# Full-length transcriptomics and proteomics of synthetic minimal cells

Analysis code accompanying the manuscript:

> **Full-length transcriptomics and proteomics reveal how genome minimization reshapes gene expression in synthetic bacteria**

We combine full-length long-read RNA sequencing (PacBio Iso-Seq and Oxford Nanopore direct-RNA for **JCVI-syn1.0**, Oxford Nanopore direct-RNA for **JCVI-syn3A**), short-read Illumina quantification, and matched proteomics to build a genome-wide, operon-resolved map of transcription, RNA processing, and translation for both synthetic cells, and to compare how halving the genome reshapes that map.

This repository holds the complete pipeline, from raw-read retrieval to the figures and Supplementary Data of the paper.

---

## Repository structure

```
.
├── Genomes_Input/            reference genomes (FASTA + GFF) for both organisms
│
├── Syn1_Transcriptomics/     syn1 PacBio + ONT (2 runs) + Illumina: read processing, isoforms
├── Syn3A_Transcriptomics/    syn3A ONT + Illumina: read processing, isoforms
├── Syn1_Syn3A_Transcriptomics/  per-gene TPM for every library of both organisms,
│                             and the cross-platform comparison panels
│
├── Syn1_Operon/              syn1 operon segmentation, annotation, and visualization (459 operons)
│
├── Syn1_Syn3A_Proteomics/    per-protein abundance and curated function, both organisms
├── Syn1_Syn3A_Corr_RNA_Proteins/  transcriptome × proteome correlation and
│                             residual analysis, both organisms
├── Syn1_RNase/               RNA-processing / ribonuclease analysis (3' erosion)
├── Syn1_Novel_ORF/           novel-ORF / antisense / intergenic transcription discovery
│
├── Genome_Reduction/         syn1 → syn3A comparison (deletions, operon remodeling, expression)
│
├── Transcription_Visualization/  per-gene browser: every RNA-seq library over any gene, both organisms
│
└── env/                      conda environment specification
```

Most scripts carry their method, parameters and a result summary in a header docstring; output conventions are in [`OUTPUT.md`](OUTPUT.md).

---

## Environment

All analysis runs in a single conda environment (Python 3.10):

```bash
conda env create -f env/environment.yml
conda activate Omics
```

It provides the Python scientific stack and every command-line tool the pipeline calls:
minimap2 and bowtie2 for mapping, samtools, bedtools and seqkit for read handling,
MUMmer for the syn1↔syn3A genome alignment, BLAST+, ViennaRNA, OSTIR, SRA-Toolkit,
TransTermHP, and FastQC/MultiQC for read QC. Everything is a conda package — there is no
system or pip install step.

The exact versions used for the published analysis are recorded in
[`env/versions_used.tsv`](env/versions_used.tsv).

---

## Workflow

Data flows from raw reads (left/top) through per-organism processing and into the cross-organism comparison:

```mermaid
flowchart TD
    GEN["Genomes_Input<br/>reference FASTA + GFF"]

    subgraph SYN1["Syn1 transcriptome"]
        PB["PacBio HiFi"] --> PBP["PacBio_Processing"] --> ISO["Isoforms_PacBio"]
        ONT1["ONT direct-RNA<br/>2 runs"] --> ONTP1["ONT_Processing"]
        ILL1["Illumina"] --> ILP1["Illumina_Processing"]
    end

    subgraph SYN3["Syn3A transcriptome"]
        ONT["ONT direct-RNA"] --> ONTP["ONT_Processing"]
        ILL3["Illumina"] --> ILP3["Illumina_Processing"]
    end

    PROT["Syn1_Syn3A_Proteomics"]

    ISO --> OPN["Syn1_Operon<br/>459 operons"]
    ISO --> RNASE["Syn1_RNase"]
    ISO --> NOV["Syn1_Novel_ORF"]
    SYN1 --> TPM["Syn1_Syn3A_Transcriptomics<br/>per-gene TPM"]
    SYN3 --> TPM
    TPM --> PROT
    TPM --> CORR["Syn1_Syn3A_Corr_RNA_Proteins"]
    PROT --> CORR
    OPN --> GR["Genome_Reduction<br/>syn1 → syn3A"]
    PROT --> GR
    SYN1 --> GR
    SYN3 --> GR
    ONTP1 --> GR
    ONTP --> GR

    SYN1 --> VIZ["Transcription_Visualization<br/>per-gene browser"]
    SYN3 --> VIZ
    OPN --> VIZ
    PROT --> VIZ
    GR --> VIZ

    PBP --> RC["Syn1_Syn3A_Transcriptomics<br/>per-gene TPM + cross-platform QC"]
    ILP1 --> RC
    ONTP1 --> RC
    ILP3 --> RC
    ONTP --> RC

    GEN -.-> SYN1
    GEN -.-> SYN3
    GEN -.-> OPN
    GEN -.-> GR
```

The syn1 ONT libraries are two independent direct-RNA runs that are **kept separate throughout** — they differ in rRNA-depletion chemistry and agree only moderately (Pearson *r* = 0.82 on log₁₀ TPM), so they are never averaged. The merged BAM is retained only as a genome-browser track and as the isoform source for the syn1 ONT panels.

### Run order

Run the stages in this order; each folder's scripts read the outputs of the stages above it.

1. **Retrieve raw reads** — bash scripts in the `*_Raw/` folders download FASTQs from the NCBI SRA (`*_Transcriptomics/{PacBio,Illumina,ONT}/*_Raw/`).
2. **Process reads** — the `*_Processing/` folders map reads (minimap2 for long reads, bowtie2 for Illumina), sort/index with samtools, and emit per-strand depth bedGraphs. The syn1 ONT folder processes both runs separately (`syn1.ONT.rep{1,2}.sorted.bam`) plus a merged browser track; run 2 is sequenced 3′→5′ and is reoriented during mapping.
3. **Cluster isoforms** — one per organism, same algorithm and parameters: `Syn1_Transcriptomics/Isoforms_PacBio/` collapses the full-length syn1 PacBio reads, and `Syn3A_Transcriptomics/Isoforms_ONT/` the syn3A ONT reads, into isoform clusters.
4. **Per-gene TPM and cross-platform comparison** — `Syn1_Syn3A_Transcriptomics/Gene_TPM.py` computes sense and antisense TPM per gene for every library (syn1 Illumina ×3, PacBio, ONT run 1 / run 2 / merged; syn3A Illumina, ONT) from the depth tracks, writing one table per organism: `syn1_TPM.tsv` (911 loci) and `syn3A_TPM.tsv` (496 loci). The platform-agreement panels are built from those two tables. `Calc_Abundances_syn3A.py` then converts the syn3A table into absolute RNA copies per cell (`syn3A_rna_abundances.tsv`) by the Breuer *et al.* 2019 mass balance, which the proteomics stages consume.
5. **Proteome** — `Syn1_Syn3A_Proteomics/Proteome_Syn1_Syn3A.ipynb` converts mass-spectrometry iBAQ into absolute protein copies per cell for both organisms, adds syn1 localization (DeepTMHMM + SignalP) and the curated syn3A function hierarchy, and compares the two. It writes `syn1_proteome.tsv` and `syn3A_proteome.tsv` for the downstream stages, matching workbooks for distribution, and two interactive pages (see below).
6. **Operons** — `Syn1_Operon/` segments the syn1 PacBio isoforms into 459 operons (`Operon_Segmentation.py`), then characterizes them (`Operon_Annotation.py`): canonical boundaries, UTR lengths, −10 promoter boxes and intrinsic terminators. Terminators come from TransTermHP — `make_transterm_coords.py` builds its gene-coordinate file, and the command that produced `syn1_TransTermHP.txt` is recorded at the top of that file. `build_operon_xlsx.py` assembles Supplementary Data S1. Run with `Syn1_Operon/` as the working directory. To look at any single operon or gene, use [`Transcription_Visualization/`](Transcription_Visualization/) (step 9).
7. **Per-organism analyses** — `Syn1_Syn3A_Corr_RNA_Proteins/` (RNA↔protein correlation, one script per organism), `Syn1_RNase/` (RNA processing + ribonuclease cleavage-site mapping), `Syn1_Novel_ORF/` (antisense / intergenic / novel ORFs).
8. **Genome reduction** — `Genome_Reduction/` runs scripts `01`→`11` in numeric order to recast the syn1→syn3A deletions as operon junctions and quantify the transcriptome/proteome reallocation; `build_reduction_SI.py` then consolidates their tables into Supplementary Data S3. Run with `Genome_Reduction/` as the working directory.
9. **Browse any gene** — `Transcription_Visualization/` draws every library over a region of your choosing (see below). Nothing downstream depends on it, so run it whenever you want to look at a gene.

> Most Python scripts carry their full method, parameters, and a result summary in a header docstring, and write a companion `.txt` log next to their outputs.

---

## Browsing transcription at any gene

To look at one gene rather than a genome-wide table, **navigate to [`Transcription_Visualization/`](Transcription_Visualization/)** and open [`Transcription_Visualization.ipynb`](Transcription_Visualization/Transcription_Visualization.ipynb):

```python
import transcription_viz as tv

tv.check_inputs()                       # confirms every upstream file this needs is reachable
fig, txt, R = tv.plot_region("ptsG")    # or "0779", "MMSYN1_0779", "JCVISYN3A_0779"
tv.show(R)
```

One call stacks **every RNA-seq library over that region in a single figure** — syn1 PacBio, ONT run 1, ONT run 2 and Illumina; syn3A ONT and Illumina — each as a read stack over its own depth track, above gene arrows for both organisms. The window, anchor, gene lists, operon bracket, TSS and terminator marks, and the red *absent from Syn3A* deletion bands are all derived from the pipeline tables, so a gene name is the only input. Every call writes a PDF, a PNG and a `_stats.txt` (per-gene depth, read counts, protein copies per cell, promoter −10 scan, predicted terminators) into `Transcription_Visualization/regions/`.

It reads the depth bedGraphs and BAMs from steps 2–3, the proteome tables from step 5, the operon map from step 6, and the deletion map from step 8 — so run those first.

---

## Browsing the proteome

Step 5 writes two self-contained HTML pages — open either in a browser, no server needed:

- `Syn1_Syn3A_Proteomics/syn3A_proteome.html` — the syn3A proteome under its curated
  Primary › Secondary › Tertiary function hierarchy. Click any function to list its
  proteins.
- `Syn1_Syn3A_Proteomics/syn1_vs_syn3A_proteome.html` — the same table joined across
  both organisms, grouped by function and by what genome reduction did to each locus.

Both give every column a filter box (text matches anywhere; `>5`, `<=0.1`, `1..10`, `*`
for any value, `-` for blanks), sortable headers, and CSV/Excel download of whatever is
on screen.

---

## Key outputs (Supplementary Data)

| | File | Built by |
|---|---|---|
| **S1** | `operon.xlsx` — per-operon table (boundaries, signals, complexes) | `Syn1_Operon/build_operon_xlsx.py` |
| **S2** | `syn1_omics.xlsx` — paired transcriptome + proteome for 911 syn1 genes | `Syn1_Syn3A_Corr_RNA_Proteins/` |
| **S3** | `genome_reduction.xlsx` — deletions, junctions, per-gene expression change | `Genome_Reduction/` |
| **S4** | `Supplementary_Data_S4_QC.zip` — RNA-sample QC reports | assembled from the `*_Processing/qc/` reports of step 2 |

---

## Data availability

Raw sequencing data are deposited in the NCBI SRA under **BioProject [PRJNA1359397](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1359397)** (syn1 Illumina + PacBio + ONT, syn3A ONT). The syn3A Illumina data are from Sandberg *et al.* (2023). Reference genomes: **CP002027.1** (JCVI-syn1.0) and **CP016816.2** (JCVI-syn3A). Mass-spectrometry proteomics are on the MassIVE repository under **MSV000099558**.

---

## Citation

If you use this code or data, please cite the manuscript above. (Full citation details will be added upon publication.)
