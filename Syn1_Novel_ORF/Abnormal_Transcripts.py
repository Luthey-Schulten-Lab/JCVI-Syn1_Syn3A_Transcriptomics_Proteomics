"""
Abnormal_Transcripts.py -- antisense and intergenic transcription in syn1 PacBio
long-read RNA-seq.

Two independent analyses over the annotated isoform clusters:

1. Antisense classification. Clusters with any antisense coverage and >= 10 reads
   are collapsed to unique transcript species (10 bp end tolerance) and sorted by
   the order of sense and antisense segments along the transcript:
     - spurious promoter  antisense-first, a transcript starting on the wrong strand
     - read-through       sense running on into an antisense gene
     - embedded           sense -> antisense -> sense
   Length and read-support distributions are reported per class, unweighted and
   read-weighted.

2. Watermark expression. Per-strand mean PacBio depth over the four synthetic
   watermark spans, against the genome-wide mean.

Inputs
------
- ../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.xlsx
- ../Genomes_Input/syn1.genes.gff3, syn1_genome.fasta
- ../Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.{plus,minus}.bedGraph

Outputs
-------
- isoform_antisense_categories.xlsx   one row per unique antisense cluster + class
- watermark_expression.txt            per-watermark depth report
- frac_{antisense,intergenic}_{count,weighted}.pdf

Run from Syn1_Novel_ORF/.
"""

from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

ISOFORM_XLSX  = '../Syn1_Transcriptomics/Isoforms_PacBio/isoform_clusters_annotated.xlsx'
MIN_READS     = 10


PERCENTILES = [0, 25, 50, 75, 90, 95, 96, 97, 98, 99, 100]


def _weighted_median(vals: np.ndarray, w: np.ndarray) -> float:
    """Weighted median: smallest value v such that cumulative weight >= 0.5 * total."""
    order = np.argsort(vals)
    vals_s, w_s = vals[order], w[order]
    cumw = np.cumsum(w_s)
    return float(vals_s[np.searchsorted(cumw, 0.5 * cumw[-1])])


def _weighted_percentile(vals: np.ndarray, w: np.ndarray, pct: float) -> float:
    """Weighted percentile (0–100) via sorted cumulative weight."""
    order = np.argsort(vals)
    vals_s, w_s = vals[order], w[order]
    cumw = np.cumsum(w_s)
    return float(vals_s[np.searchsorted(cumw, (pct / 100) * cumw[-1])])


def _print_percentiles(vals: np.ndarray, w, field: str, subtitle: str):
    """Print a percentile table to stdout."""
    print(f'\n  {field}  [{subtitle}]')
    header = '  ' + ''.join(f'  p{p:<5}' for p in PERCENTILES)
    print(header)
    if w is None:
        row = '  ' + ''.join(f'  {np.percentile(vals, p):<7.4f}' for p in PERCENTILES)
    else:
        row = '  ' + ''.join(
            f'  {_weighted_percentile(vals, w, p):<7.4f}' for p in PERCENTILES)
    print(row)


def _plot_one_frac(vals, weights, field, color, ylabel, subtitle,
                   n_total, min_reads, out_pdf):
    """
    Draw and save a single figure for one (field, weighting) combination.

    vals    : Series of frac values for all isoforms (including zeros).
    weights : None for isoform-count histogram; Series of n_reads for weighted.

    Markers drawn on the plot:
      - Unweighted median (dashed black line, always shown).
      - Weighted median   (dashed gray  line, only when weights are provided).
    Percentile table is printed to stdout for every call.
    """
    bins = np.linspace(0, 1, 51)   # 50 equal bins; first bin [0, 0.02] captures zeros

    v = vals.values.astype(float)
    w = weights.loc[vals.index].astype(float).values if weights is not None else None

    med_unweighted = float(np.median(v))

    fig, ax = plt.subplots(figsize=(6, 3))
    # fig.suptitle(
    #     f'{field}  —  {subtitle}\n'
    #     f'(n_reads > {min_reads},  n = {n_total:,})',
    #     fontsize=10, fontweight='bold')

    ax.hist(v, bins=bins, weights=w, color=color, alpha=0.75,
            edgecolor='white', lw=0.3)
    ax.set_yscale('log')
    ax.set_xlim(-0.01, 1.01)
    ax.set_xlabel(field, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)

    # Always show the unweighted median
    ax.axvline(med_unweighted, color='black', lw=1.2, ls='--',
               label=f'median (Isoform kinds) = {med_unweighted:.3f}')

    # Weighted median only when weights are provided
    if w is not None:
        med_weighted = _weighted_median(v, w)
        ax.axvline(med_weighted, color='dimgray', lw=1.2, ls=':',
                   label=f'median (weighted by count)   = {med_weighted:.3f}')

    ax.legend(fontsize=8, framealpha=0.8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    # plt.show()
    print(f'Saved: {out_pdf}')

    # Percentile report
    _print_percentiles(v, w, field, subtitle)


def plot_frac_distributions(
        xlsx_path: str = ISOFORM_XLSX,
        min_reads: int = MIN_READS):
    """
    Produce four separate figures (each saved to its own PDF):
      1. frac_antisense   — isoform count
      2. frac_intergenic  — isoform count
      3. frac_antisense   — read-count weighted
      4. frac_intergenic  — read-count weighted

    Zeros are included in the first histogram bin [0, 0.02].
    The median (over all isoforms including zeros) is marked.
    """
    df  = pd.read_excel(xlsx_path)
    iso = df[df['n_reads'] > min_reads].copy()
    n_total = len(iso)
    
    print(f'Total isoforms with n_reads > {min_reads}: {n_total:,}')
    
    specs = [
        # (field,             color,     weights,          ylabel,                        subtitle,              out_pdf)
        ('frac_antisense',  '#d62728', None,              'Unique isoforms (log)',     'Isoform count',       'frac_antisense_count.pdf'),
        ('frac_intergenic', '#1f77b4', None,              'Unique isoforms (log)',     'Isoform count',       'frac_intergenic_count.pdf'),
        ('frac_antisense',  '#d62728', iso['n_reads'],    'Isoform read count (log)',      'Read-count weighted', 'frac_antisense_weighted.pdf'),
        ('frac_intergenic', '#1f77b4', iso['n_reads'],    'Isoform read count (log)',      'Read-count weighted', 'frac_intergenic_weighted.pdf'),
    ]

    for field, color, weights, ylabel, subtitle, out_pdf in specs:
        vals = iso[field].dropna()
        _plot_one_frac(vals, weights, field, color, ylabel,
                       subtitle, n_total, min_reads, out_pdf)


# if __name__ == '__main__':
#     plot_frac_distributions()


GFF3_FILE    = '../Genomes_Input/syn1.genes.gff3'
BOUNDARY_TOL = 10   # bp — isoforms within this many bp of containment are treated as identical


def cluster_isoforms(isoforms: pd.DataFrame,
                     tol: int = BOUNDARY_TOL) -> list:
    """
    Cluster isoforms by containment (with tolerance).  Copied verbatim from
    Operon_Segmentation.ipynb.

    Two isoforms i and j belong to the same cluster iff one is within `tol` bp
    of being fully contained within the other.  The cluster boundary is the
    outermost (longest) isoform span.  Uses union-find + sweep-line.

    Returns a list of dicts, each with keys:
      chrom, strand, start0, end0, members (list of isoform_id), reads (list of n_reads).
    """
    if isoforms.empty:
        return []

    iso    = isoforms.reset_index(drop=True)
    starts = iso['start0'].astype(int).tolist()
    ends   = iso['end0'].astype(int).tolist()
    n      = len(iso)

    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)

    order = sorted(range(n), key=lambda i: starts[i])
    for idx, i in enumerate(order):
        si, ei = starts[i], ends[i]
        for j in order[idx + 1:]:
            sj = starts[j]
            if sj >= ei + tol:
                break
            ej = ends[j]
            i_in_j = (si >= sj - tol) and (ei <= ej + tol)
            j_in_i = (sj >= si - tol) and (ej <= ei + tol)
            if i_in_j or j_in_i:
                union(i, j)

    components = defaultdict(list)
    for i in range(n):
        components[find(i)].append(i)

    chrom  = iso.loc[0, 'chrom']
    strand = iso.loc[0, 'strand']
    blocks = []
    for indices in components.values():
        members = [iso.loc[k, 'isoform_id'] for k in indices]
        reads   = [int(iso.loc[k, 'n_reads']) for k in indices]
        s0 = min(starts[k] for k in indices)
        e0 = max(ends[k]   for k in indices)
        blocks.append({
            'chrom': chrom, 'strand': strand,
            'start0': s0,   'end0': e0,
            'members': members, 'reads': reads,
        })
    return sorted(blocks, key=lambda b: b['start0'])


def read_genes_gff3(gff3_path: str) -> pd.DataFrame:
    """
    Parse a GFF3 file and return a DataFrame with columns:
      chrom, start0, end0, strand, locus_tag, gene_name

    GFF3 coordinates are 1-based inclusive; converted here to
    0-based half-open [start0, end0).
    """
    rows = []
    with open(gff3_path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            chrom, _, ftype, s1, e1, _, strand, _, attrs = parts
            if ftype != 'gene':
                continue

            # Parse locus_tag and gene name from attributes
            locus_tag = ''
            gene_name = ''
            for token in attrs.split(';'):
                if token.startswith('locus_tag='):
                    locus_tag = token[len('locus_tag='):]
                elif token.startswith('Name='):
                    gene_name = token[len('Name='):]

            rows.append({
                'chrom':     chrom,
                'start0':    int(s1) - 1,   # 1-based → 0-based
                'end0':      int(e1),        # inclusive → half-open
                'strand':    strand,
                'locus_tag': locus_tag,
                'gene_name': gene_name,
            })

    return pd.DataFrame(rows)


def _in_any(pos: int, gdf: pd.DataFrame) -> bool:
    """True if genomic position pos falls inside any gene in gdf [start0, end0)."""
    if gdf.empty:
        return False
    return bool(((gdf['start0'] <= pos) & (gdf['end0'] > pos)).any())


def _classify_one(p5: int, p3: int, strand: str,
                  chrom: str, genes: pd.DataFrame) -> str:
    """
    Classify a single isoform by the ordered sequence of overlapping genes
    in transcript 5'→3' direction.

    Each overlapping gene is labelled S (sense, same strand as isoform) or
    A (antisense, opposite strand).  The sequence of S/A labels determines
    the category:

      spurious_prom  — first gene in transcript order is antisense (A…)
      read_through   — starts with sense gene(s) then antisense (S…A),
                       with no sense gene after the last antisense gene
      embedded       — sense gene(s), then antisense gene(s), then sense
                       gene(s) again (S…A…S); antisense is sandwiched
      other          — no annotated genes overlap, or only sense genes
                       (frac_antisense > 0 but from sub-gene-resolution gaps)

    Sorting key: for + strand, genes are ordered by their midpoint ascending
    (left = 5'); for − strand, by midpoint descending (right = 5').
    """
    iso_s = min(p5, p3)
    iso_e = max(p5, p3)

    ov = genes[(genes['chrom'] == chrom) &
               (genes['start0'] < iso_e) &
               (genes['end0']   > iso_s)].copy()

    if ov.empty:
        return 'other'

    # Sort genes in transcript 5'→3' order using genomic midpoint
    ov['_mid'] = (ov['start0'] + ov['end0']) / 2.0
    if strand == '-':
        ov['_mid'] = -ov['_mid']        # negate: higher genomic coord = 5' on − strand
    ov = ov.sort_values('_mid')

    # Sequence of S/A labels in transcript order
    seq = ['S' if g == strand else 'A' for g in ov['strand']]

    # Classification based on sequence pattern
    if seq[0] == 'A':
        return 'spurious_prom'

    if 'A' not in seq:
        return 'other'                  # sense-only overlap

    last_anti  = max(i for i, s in enumerate(seq) if s == 'A')
    sense_after = any(s == 'S' for s in seq[last_anti + 1:])

    return 'embedded' if sense_after else 'read_through'


def classify_antisense_structure(
        xlsx_path: str = ISOFORM_XLSX,
        gff3_path: str = GFF3_FILE,
        min_reads: int = MIN_READS) -> pd.DataFrame:
    """
    For isoforms with n_reads > min_reads and frac_antisense > 0:

      1. Cluster isoforms per (chrom, strand) by containment using
         cluster_isoforms (BOUNDARY_TOL = 10 bp).  Isoforms that are
         near-identical variants of the same transcript collapse into one
         cluster.  The cluster representative is the member with the most
         reads.

      2. Classify each representative by the ordered gene sequence it
         overlaps in transcript 5'→3' direction (see _classify_one).

    Categories
    ----------
    read_through   sense gene(s) → antisense gene(s); no return to sense.
    embedded       sense → antisense → sense; antisense sandwiched in operon.
    spurious_prom  first overlapping gene is antisense; spurious promoter.
    other          only sense genes overlap (antisense from intergenic gap).

    Prints a two-level summary (raw isoforms → unique clusters) and returns
    the representative DataFrame with columns added:
      cluster_size, cluster_reads, antisense_category.
    """
    df    = pd.read_excel(xlsx_path)
    iso   = df[(df['n_reads'] > min_reads) & (df['frac_antisense'] > 0)].copy()
    genes = read_genes_gff3(gff3_path)

    n_raw = len(iso)
    n_all = len(df[df['n_reads'] > min_reads])
    print(f'\nGene models loaded from GFF3 : {len(genes):,}')
    print(f'Raw isoforms (frac_anti > 0) : {n_raw:,}  ({100*n_raw/n_all:.1f}% of all filtered)')

    # ── Step 1: cluster per (chrom, strand), pick representative ──────────
    reps = []
    for (chrom, strand), grp in iso.groupby(['chrom', 'strand']):
        for cl in cluster_isoforms(grp, tol=BOUNDARY_TOL):
            members = grp[grp['isoform_id'].isin(cl['members'])]
            rep = members.loc[members['n_reads'].idxmax()].copy()
            rep['cluster_size']  = len(cl['members'])
            rep['cluster_reads'] = int(sum(cl['reads']))
            reps.append(rep)

    reps_df = pd.DataFrame(reps).reset_index(drop=True)
    n_uniq  = len(reps_df)
    print(f'Unique clusters (repr.)      : {n_uniq:,}  '
          f'(collapsed from {n_raw} raw,  tol = {BOUNDARY_TOL} bp)')

    # ── Step 2: classify each representative ──────────────────────────────
    reps_df['antisense_category'] = [
        _classify_one(int(r['pos5p0']), int(r['pos3p0']),
                      r['strand'], r['chrom'], genes)
        for _, r in reps_df.iterrows()
    ]

    # ── Summary table ─────────────────────────────────────────────────────
    order = ['read_through', 'embedded', 'spurious_prom', 'other']
    labels = {
        'read_through':  'Read-through   (sense → anti)',
        'embedded':      '  └ Embedded   (sense → anti → sense)',
        'spurious_prom': 'Spurious prom  (anti → ...)',
        'other':         'Other          (sense-only overlap)',
    }

    print(f'\n  {"Category":<42}  {"n":>5}  {"% unique":>8}  {"% all":>6}')
    print(f'  {"-"*67}')
    for cat in order:
        cnt = (reps_df['antisense_category'] == cat).sum()
        print(f'  {labels[cat]:<42}  {cnt:>5,}  '
              f'{100*cnt/n_uniq:>7.1f}%  {100*cnt/n_all:>5.2f}%')
    print(f'  {"-"*67}')
    rt = reps_df['antisense_category'].isin(['read_through', 'embedded']).sum()
    print(f'  Read-through total (read_through + embedded) : {rt:,}  '
          f'({100*rt/n_uniq:.1f}% of unique antisense clusters)')

    # ── Length and read distributions per category ─────────────────────────
    pcts = [0, 25, 50, 75, 100]
    for metric, col, unit in [
        ('Isoform length (bp)',       'isoform_len_bp', 'bp'),
        ('Repr. reads (n_reads)',     'n_reads',        ''),
        ('Cluster reads (all members)', 'cluster_reads',''),
    ]:
        print(f'\n  {metric}')
        hdr = f'  {"Category":<38}  {"n":>4}  {"mean":>7}'
        hdr += ''.join(f'  {"p"+str(p):>7}' for p in pcts)
        hdr += f'  {"min":>7}  {"max":>7}'
        print(hdr)
        print(f'  {"-"*85}')
        for cat in order:
            sub = reps_df.loc[reps_df['antisense_category'] == cat, col].dropna()
            if sub.empty:
                print(f'  {labels[cat]:<38}  {"—":>4}')
                continue
            row = (f'  {labels[cat]:<38}  {len(sub):>4}  {sub.mean():>7.0f}')
            row += ''.join(f'  {sub.quantile(p/100):>7.0f}' for p in pcts)
            row += f'  {sub.min():>7.0f}  {sub.max():>7.0f}'
            print(row)

    return reps_df

iso = classify_antisense_structure()
iso.to_excel('isoform_antisense_categories.xlsx', index=False)

# Water marks

W1 = "TTAACTAGCTAAGTTCGAATATTTCTATAGCTGTACATATTGTAATGCTGATAACTAATACTGTGCGCTTGACTGTGATCCTGATAAATAACTTCTTCTGTAGGGTAGAGTTTTATTTAAGGCTACTCACTGGTTGCAAACCAATGCCGTACATTACTAGCTTGATCCTTGGTCGGTCATTGGGGGATATCTCTTACTAATAGAGCGGCCTATCGCGTATTCTCGCCGGACCCCCCTCTCCCACACCAGCGGTGTAGCATCACCAAGAAAATGAGGGGAACGGATGAGGAACGAGTGGGGGCTCATTGCTGATCATAATGACTGTTTATATACTAATGCCGTCAACTGTTTGCTGTGATACTGTGCTTTCGAGGGCGGGAGATTCGTTTTTGACATACATAAATATCATGACAAAACAGCCGGTCATGACAAAACAGCCGGTCATAATAGATTAGCCGGTGACTGTGAAACTAAAGCTACTAATGCCGTCAATAAATATGATAATAGCAACGGCACTGACTGTGAAACTAAAGCCGGCACTCATAATAGATTAGCCGGAGTCGTATTCATAGCCGGTAGATATCACTATAAGGCCCAGGATCATGATGAACACAGCACCACGTCGTCGTCCGAGTTTTTTTGCTGCGACGTCTATACCACGGAAGCTGATCATAAATAGTTTTTTTGCTGCGGCACTAGAGCCGGACAAGCACACTACGTTTGTAAATACATCGTTCCGAATTGTAAATAATTTAATTTCGTATTTAAATTATATGATCACTGGCTATAGTCTAGTGATAACTACAATAGCTAGCAATAAGTCATATATAACAATAGCTGAACCTGTGCTACATATCCGCTATACGGTAGATATCACTATAAGGCCCAGGACAATAGCTGAACTGACGTCAGCAACTACGTTTAGCTTGACTGTGGTCGGTTTTTTTGCTGCGACGTCTATACGGAAGCTCATAACTATAAGAGCGGCACTAGAGCCGGCACACAAGCCGGCACAGTCGTATTCATAGCCGGCACTCATGACAAAACAGCGGCGCGCCTTAACTAGCTAA"
W2 = "TTAACTAGCTAACAACTGGCAGCATAAAACATATAGAACTACCTGCTATAAGTGATACAACTGTTTTCATAGTAAAACATACAACGTTGCTGATAGTACTCCTAAGTGATAGCTTAGTGCGTTTAGCATATATTGTAGGCTTCATAATAAGTGATATTTTAGCTACGTAACTAAATAAACTAGCTATGACTGTACTCCTAAGTGATATTTTCATCCTTTGCAATACAATAACTACTACATCAATAGTGCGTGATATGCCTGTGCTAGATATAGAACACATAACTACGTTTGCTGTTTTCAGTGATATGCTAGTTTCATCTATAGATATAGGCTGCTTAGATTCCCTACTAGCTATTTCTGTAGGTGATATACGTCCATTGCATAAGTTAATGCATTTAACTAGCTGTGATACTATAGCATCCCCATTCCTAGTGCATATTTTCATCCTAGTGCTACGTGATATAATTGTACTAATGCCTGTAGATAATTTAATGCCTGGCTCGTTTGTAGGTGATAATTTAGTGCCTGTAAAACATATACCTGAGTGCTCGTTGCGTGATAGTTCGTTCATGCATATACAACTAGGCTGCTGTGATATGGTCACTGCCCTTACTGTGCTACATATTACTGCGAGGGGGATGACGTATAAACCTGTTGTAAGTGATATGACGTATATAACTACTAGTGATATGACGTATAGGCTAGAACAACGTGATATGACGTATATGACTACTGTCCCAAACATCAGTGATATGACGTATACTATAATTTCTATAATAGTGATAAATAAACCTGGGCTAAATACGTTCCTGAATACGTGGCATAAACCTGGGCTAACGAGGAATACCCATAGTTTAGCAATAAGCTATAGTTCGTCATTTTTAAGGCGCGCCTTAACTAGCTAA"
W3 = "TTAACTAGCTAATTTAACCATATTTAAATATCATCCTGATTTTCACTGGCTCGTTGCGTGATATAGATTCTACTGTAGTGCTAGATAGTTCTGTACTAGGTGATACTATAGATTTCATAGATAGCACTACTGGCTTCATGCTAGGCATCCCAATAGCTAGTGATAGTTTAGTGCATACAACGTCATGTGATACAACGTTGCTGGCTGTAGATACAACGTCGTATTCTGTAAGTGATACAATAGCTATTGCTGTGCATAGGCCTATAGTGGCTGTAACTAGTGATATCACGTAACAACCATATAAGTTAGATTTAATGCCCCTGACTGAACGCTCGTTGCGTGATAGTTTAGGCTCGTTGCATACAACTGTGATTTTCATAAAACAACGTGATAATTTAGTGCTAGATAAGTTCCGCTTAGCAAGTGATAGTTTCCGCTTGACTGTGCATAGTTCGTTCATGCGCTCGTTGCGTGATAAACTAGGCAGCTTCACAACTGATAATTTAATTGCTGATATTGCTGGCTGTCTAGTGCTAGTGATCATAGTGCGTGATAGTTTAAGCTGCTCTGTTTTAGATATCACGTGCTTGATAATGAAACTAACTAGTGATACTACGTAGTTAACTATGAATAGGCCTACTGTAAATTCAATAGTGCGTGATATTGAACTAGATTCTGCAACTGCTAATATGCCGTGCTGCACGTTTGGTGATAGTTTAGCATGCTTCACTATAATAAATATGGTAGTTGTAACTACTGCGAATAGGGGGAGCTTAATAAATATGATCACTGTGCTACGCTATATGCCGTTGAATATAGGCTATATGATCATAACATATATAGCTATAAGTGATAAGTTCCTGAATATAGGCTATATGATCATAACATATACAACTGTACTCATGAATAAGTTAACGAGGATTAACTAGCTAA"
W4 = "TTAACTAGCTAATTTCATTGCTGATCACTGTAGATATAGTGCATTCTATAAGTCGCTCCCACAGGCTAGTGCTGCGCACGTTTTTCAGTGATATTATCCTAGTGCTACATAACATCATAGTGCGTGATAAACCTGATACAATAGGTGATATCATAGCAACTGAACTGACGTTGCATAGCTCAACTGTGATCAGTGATATAGATTCTGATACTATAGCAACGTTGCGTGATATTTTCACTACTGGCTTGACTGTAGTGCATATGATAGTACGTCTAACTAGCATAACTAGTGATAGTTATATTTCTATAGCTGTACATATTGTAATGCTGATAACTAGTGATATAATCCAACTAGATAGTCCTGAACTGATCCCTATGCTAACTAGTGATAAACTAACTGATACATCGTTCCTGCTACGTGATAGCTTCACTGAGTTCCATACATCGTCGTGCTTAAACATCAGTGATAACACTATAGAGTTCATAGATACTGCATTAACTAGTGATATGACTGCAAATAGCTTGACGTTTTGCAGTCTAAAACAACGTGATAATTCTGTAGTGCTAGATACTATAGATTTCCTGCTAAGTGATAAGTCTACTGATTTACTAATGAATAGCTTGGTTTTGGCATACACTGTGCGCTGCACTGGTGATAGCTTTTCGTTGATGAATAATTTCCCTAGCACTGTGCGTGATATGCTAGATTCTGTAGATAGGCTAAATTCGTCTACGTTTGTAGGTGATAGTTTAGTTGCTGTAACTAATATTATCCCTGTGCCGTTGCTAAGCTGTGATATCATAGTGCTGCTAGATATGATAAGCAAACTAATAGAGTCGAGGGGGAGTCTCATAGTGAATACTGATATTTTAGTGCTGCCGTTGAATAAGTTCCCTGAACATTGTGATACTGATATTTTAGTGCTGCCGTTGAATATCCTGCATTTAACTAGCTTGATAGTGCATTCGAGGAATACCCATACTACTGTTTTCATAGCTAATTATAGGCTAACATTGCCAATAGTGCGGCGCGCCTTAACTAGCTAA"

WATERMARK_DICT = {
    "W1": W1,
    "W2": W2,
    "W3": W3,
    "W4": W4
}


from Bio import SeqIO
from Bio.Seq import Seq

SYN1_FASTA  = '../Genomes_Input/syn1_genome.fasta'
DEPTH_PLUS  = '../Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.plus.bedGraph'
DEPTH_MINUS = '../Syn1_Transcriptomics/PacBio/PacBio_Processing/depth_bedgraph/syn1.PacBio.FLNC.HQ.minus.bedGraph'
OUT_TXT     = 'watermark_expression.txt'

record = next(SeqIO.parse(SYN1_FASTA, 'fasta'))
CHROM = record.id
seq   = str(record.seq).upper()


def locate(sub: str):
    """(start0, end0, genome_strand) of an exact match of sub or its reverse complement."""
    i = seq.find(sub)
    if i != -1:
        return i, i + len(sub), '+'
    j = seq.find(str(Seq(sub).reverse_complement()))
    if j != -1:
        return j, j + len(sub), '-'
    return None


def covered_genes(genes_df, s: int, e: int):
    hit = genes_df[(genes_df['chrom'] == CHROM) & (genes_df['start0'] < e) & (genes_df['end0'] > s)]
    rows = []
    for _, g in hit.iterrows():
        ov = min(e, int(g['end0'])) - max(s, int(g['start0']))
        name = g['gene_name'] if str(g['gene_name']).strip() else g['locus_tag']
        rows.append(f"{g['locus_tag']} ({name}, {g['strand']} strand) overlap {ov} bp")
    return rows


def mean_depth(dep_df, s: int, e: int) -> float:
    """Mean per-base depth over [s, e); positions absent from the bedGraph count as 0."""
    L = e - s
    sub = dep_df[(dep_df['chrom'] == CHROM) & (dep_df['end'] > s) & (dep_df['start'] < e)]
    if sub.empty or L <= 0:
        return 0.0
    ov = np.minimum(sub['end'].values, e) - np.maximum(sub['start'].values, s)
    return float((sub['val'].values * ov).sum()) / L


genes_df  = read_genes_gff3(GFF3_FILE)
dep_plus  = pd.read_csv(DEPTH_PLUS,  sep='\t', header=None, names=['chrom', 'start', 'end', 'val'])
dep_minus = pd.read_csv(DEPTH_MINUS, sep='\t', header=None, names=['chrom', 'start', 'end', 'val'])

GENOME_LEN = len(seq)
def genome_mean(dep_df) -> float:
    return float((dep_df['val'].values * (dep_df['end'].values - dep_df['start'].values)).sum()) / GENOME_LEN
gmean_p, gmean_m = genome_mean(dep_plus), genome_mean(dep_minus)

out = []
def w(msg=''):
    print(msg); out.append(msg)

w('WATERMARK EXPRESSION ANALYSIS')
w('=' * 60)
w(f'Genome: {CHROM} ({len(seq):,} bp)')
w(f'Genome-wide mean PacBio depth (baseline): + strand {gmean_p:.0f}, - strand {gmean_m:.0f}.')
w('Depth = mean per-base PacBio FLNC coverage over the watermark span, per strand.')
w('')

sum_p = sum_m = tot_len = 0.0
for wname, wm in WATERMARK_DICT.items():
    loc = locate(wm)
    w(f'{wname}  (length {len(wm)} bp)')
    if loc is None:
        w('  not found in the syn1 genome (no exact / reverse-complement match)')
        w('')
        continue
    s, e, strand = loc
    w(f'  location : {s:,}-{e:,} (0-based half-open), genome strand {strand}')
    cg = covered_genes(genes_df, s, e)
    if cg:
        w(f'  covered genes ({len(cg)}):')
        for c in cg:
            w(f'    - {c}')
    else:
        w('  covered genes : none (no annotated gene overlaps)')
    mp, mm = mean_depth(dep_plus, s, e), mean_depth(dep_minus, s, e)
    sum_p += mp * (e - s); sum_m += mm * (e - s); tot_len += (e - s)
    sense, anti = (mp, mm) if strand == '+' else (mm, mp)
    w(f'  mean depth: + strand = {mp:.2f}, - strand = {mm:.2f}'
      f'   (sense {sense:.2f}, antisense {anti:.2f})')
    w('')

if tot_len:
    wp, wmn = sum_p / tot_len, sum_m / tot_len
    w('-' * 60)
    w(f'All 4 watermarks (length-weighted mean): + strand {wp:.1f}, - strand {wmn:.1f}')
    w(f'  vs genome-wide average: + {100*wp/gmean_p:.1f}%, - {100*wmn/gmean_m:.1f}% '
      f'({gmean_p/wp:.1f}x lower on +, {gmean_m/wmn:.1f}x lower on -)')
    w('')

with open(OUT_TXT, 'w') as fh:
    fh.write('\n'.join(out) + '\n')
print(f'Wrote {OUT_TXT}')
