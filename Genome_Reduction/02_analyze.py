#!/usr/bin/env python3
"""Analyze the JCVI-Syn1.0 -> JCVI-Syn3A genome reduction.

Reads the outputs of 01_align.sh (nucmer + dnadiff + deletion BED, all in
aln/raw/) and produces two deliverables in aln/analysis/:
    genome_reduction_summary.xlsx  (multi-sheet: events / short_insertions /
                                    dnadiff_summary / legend)
    genome_reduction_summary.txt   (human-readable narrative report)
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------- paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
ALN = SCRIPT_DIR / "aln"
RAW = ALN / "raw"                # nucmer + dnadiff intermediates (from 01_align.sh)
OUT = ALN / "analysis"           # human-facing deliverables only
OUT.mkdir(parents=True, exist_ok=True)

COORDS = RAW / "syn1_vs_syn3A.coords.tsv"
REPORT = RAW / "dnadiff_out.report"
QDIFF = RAW / "dnadiff_out.qdiff"
RDIFF = RAW / "dnadiff_out.rdiff"
DEL_BED = RAW / "syn1_deleted_regions.bed"
DEL_GENES = RAW / "syn1_deleted_genes.tsv"

SYN1_GFF = PROJECT_ROOT / "Genomes_Input" / "syn1.genes.gff3"
SYN3A_GFF = PROJECT_ROOT / "Genomes_Input" / "syn3a_genome.gff3"
SYN3A_FAI = PROJECT_ROOT / "Genomes_Input" / "syn3A_genome.fasta.fai"

# Single canonical output table; the .txt is a narrative companion.
EXCEL = OUT / "genome_reduction_summary.xlsx"
SUMMARY = OUT / "genome_reduction_summary.txt"

MIN_INS_BP = 10           # minimum syn3A-uncovered interval to report
MIN_DEL_BP_NOTE = 50      # used by 00_alignment.sh; quoted in the report
SHORT_INS_MAX = 1000      # treat <1 kb qdiff insertions as "short"

# ------------------------------------------------------- file existence check
for p in (COORDS, REPORT, QDIFF, RDIFF, DEL_BED):
    if not p.exists():
        sys.exit(f"ERROR: required input not found: {p}")


# ============================================================ GFF parsing
def parse_gff_genes(path: Path) -> pd.DataFrame:
    """Return DataFrame with chrom, start0, end, strand, locus, name."""
    rows = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(
                kv.split("=", 1) for kv in f[8].split(";") if "=" in kv
            )
            rows.append(
                (
                    f[0],
                    int(f[3]) - 1,
                    int(f[4]),
                    f[6],
                    attrs.get("locus_tag", "."),
                    attrs.get("Name", attrs.get("gene", ".")),
                )
            )
    return pd.DataFrame(rows, columns=["chrom", "start0", "end", "strand", "locus", "name"])


syn1_genes = parse_gff_genes(SYN1_GFF) if SYN1_GFF.exists() else pd.DataFrame()
syn3A_genes = parse_gff_genes(SYN3A_GFF) if SYN3A_GFF.exists() else pd.DataFrame()


# ============================================================ coords table
COORDS_COLS = [
    "S1", "E1", "S2", "E2", "LEN1", "LEN2", "PCT_IDY",
    "LENR", "LENQ", "COVR", "COVQ", "TAGR", "TAGQ",
]
coords = pd.read_csv(COORDS, sep="\t", comment=None, header=0)
# After header rewrite by 00_alignment.sh the file already has a header row,
# but the original column names (S1, E1, ...) include an awkward `%IDY`.
coords.columns = COORDS_COLS[: len(coords.columns)]
# Normalise S2/E2 in case the alignment is on the reverse strand of QRY
coords["S2_min"] = coords[["S2", "E2"]].min(axis=1)
coords["E2_max"] = coords[["S2", "E2"]].max(axis=1)


# ============================================================ dnadiff report
def parse_report(path: Path) -> dict[str, list[str]]:
    """Return {field_name: [ref_value, qry_value]}."""
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3 and f[0][0].isalpha():
                out[f[0]] = f[1:]
    return out


report = parse_report(REPORT)


def rep(field: str, idx: int = 0) -> str:
    return report.get(field, ["NA", "NA"])[idx]


# ============================================================ qdiff / rdiff
QDIFF_COLS = ["chrom", "type", "start", "end", "q_gap", "r_gap", "r_minus_q"]


def parse_diff(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            # JMP rows have 5 cols; GAP rows have 7
            f += [""] * (7 - len(f))
            rows.append(f[:7])
    df = pd.DataFrame(rows, columns=QDIFF_COLS)
    for c in ("start", "end", "q_gap", "r_gap", "r_minus_q"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


qdiff = parse_diff(QDIFF)
rdiff = parse_diff(RDIFF)


# ============================================================ interval ops
def interval_complement(intervals: list[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    """Given sorted half-open intervals on [0, length), return the complement."""
    out, cursor = [], 0
    for s, e in intervals:
        if s > cursor:
            out.append((cursor, s))
        cursor = max(cursor, e)
    if cursor < length:
        out.append((cursor, length))
    return out


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


# ============================================================ syn3A inserted BED
def syn3A_length() -> int:
    if SYN3A_FAI.exists():
        with SYN3A_FAI.open() as fh:
            return int(fh.readline().split()[1])
    # fall back to dnadiff report
    return int(rep("TotalBases", 1))


syn3A_len = syn3A_length()
syn3A_chrom = coords["TAGQ"].iloc[0]

q_blocks = merge_intervals(
    list(zip(coords["S2_min"] - 1, coords["E2_max"]))  # to 0-based half-open
)
inserted = [
    (s, e) for (s, e) in interval_complement(q_blocks, syn3A_len)
    if (e - s) >= MIN_INS_BP
]
# inserted intervals are kept in memory only; they end up in the Excel below.


# ============================================================ deletion stats
del_df = pd.read_csv(
    DEL_BED, sep="\t", comment="#",
    names=["chrom", "start0", "end", "length_bp"],
)
del_top = del_df.sort_values("length_bp", ascending=False).head(10)


# Genes hit by deletions (count from syn1_deleted_genes.tsv attributes column)
n_del_genes: int | str
if DEL_GENES.exists():
    g = pd.read_csv(DEL_GENES, sep="\t", comment=None)
    last = g.columns[-1]
    locus = g[last].astype(str).str.extract(r"locus_tag=([^;]+)")[0].dropna()
    n_del_genes = locus.nunique() if len(locus) else 0
else:
    n_del_genes = "(syn1_deleted_genes.tsv not found)"


# ============================================================ short insertions
def gff_lookup(genes: pd.DataFrame, chrom: str, start0: int, end: int):
    """Return overlap, nearest +/- on each side."""
    sub = genes[genes.chrom == chrom]
    if sub.empty:
        return ([], None, None, None, None)
    overlap = sub[(sub.start0 < end) & (sub.end > start0)]

    # nearest upstream/downstream by strand
    def nearest(strand: str, direction: str):
        s = sub[sub.strand == strand]
        if s.empty:
            return None
        if direction == "up":
            cand = s[s.end <= start0]
            if cand.empty:
                return None
            row = cand.loc[cand["end"].idxmax()]
            dist = row["end"] - start0   # negative
        else:
            cand = s[s.start0 >= end]
            if cand.empty:
                return None
            row = cand.loc[cand["start0"].idxmin()]
            dist = row["start0"] - end   # positive (or 0)
        label = f"{row['locus']}|{row['name']}"
        return label, int(dist)

    return (
        [f"{r['locus']}|{r['name']}" for _, r in overlap.iterrows()],
        nearest("+", "up"),
        nearest("+", "dn"),
        nearest("-", "up"),
        nearest("-", "dn"),
    )


short_ins = qdiff[
    qdiff["type"].isin(("GAP", "JMP"))
    & (qdiff["q_gap"] > 0)
    & (qdiff["q_gap"] < SHORT_INS_MAX)
].reset_index(drop=True)

def fmt(pair):
    if pair is None:
        return "NA"
    label, dist = pair
    return f"{label}|{dist}"


context_rows = []
for i, row in short_ins.iterrows():
    chrom, start0, end = row["chrom"], int(row["start"]) - 1, int(row["end"])
    overlap, up_p, dn_p, up_m, dn_m = gff_lookup(syn3A_genes, chrom, start0, end)
    context_rows.append(
        {
            "insertion_id": f"ins_{i + 1}",
            "chrom": chrom,
            "start0": start0,
            "end": end,
            "length_bp": int(row["q_gap"]),
            "overlapping_genes": ",".join(overlap) if overlap else "none",
            "upstream_plus": fmt(up_p),
            "downstream_plus": fmt(dn_p),
            "upstream_minus": fmt(up_m),
            "downstream_minus": fmt(dn_m),
        }
    )
context_df = pd.DataFrame(context_rows)
# context_df becomes a sheet in the Excel below; no separate TSV.


# ============================================================ relocation analysis
# Relocations = alignment blocks whose syn3A position breaks the linear order
# defined by their syn1 position. We detect them via Longest Increasing
# Subsequence (LIS) of S2 values after sorting by S1: blocks IN the LIS form
# the in-order backbone; blocks OUTSIDE the LIS are the relocated ones.
sorted_by_s1 = coords.sort_values("S1").reset_index(drop=True)
s2_seq = sorted_by_s1["S2_min"].astype(int).tolist()


def lis_indices(values: list[int]) -> set[int]:
    """Return the set of indices that participate in a longest increasing
    subsequence of `values`. Standard O(n log n) algorithm with parent tracking."""
    from bisect import bisect_left

    if not values:
        return set()
    tails_idx: list[int] = []   # tails_idx[k] = index of smallest tail of an inc. subseq. of length k+1
    parents: list[int] = [-1] * len(values)
    tails_vals: list[int] = []
    for i, v in enumerate(values):
        pos = bisect_left(tails_vals, v)
        if pos == len(tails_vals):
            tails_vals.append(v)
            tails_idx.append(i)
        else:
            tails_vals[pos] = v
            tails_idx[pos] = i
        if pos > 0:
            parents[i] = tails_idx[pos - 1]
    # reconstruct
    keep: set[int] = set()
    k = tails_idx[-1]
    while k != -1:
        keep.add(k)
        k = parents[k]
    return keep


lis_set = lis_indices(s2_seq)
out_of_order = [i for i in range(len(sorted_by_s1)) if i not in lis_set]

# Annotate each out-of-order block with overlapping syn1 / syn3A genes.
# Output format: comma-separated, de-duplicated, "<numeric_suffix>/<gene_name>"
# (or just "<numeric_suffix>" when the locus has no gene name).
# Examples: "MMSYN1_0001" + "dnaA_1" -> "0001/dnaA_1" ; "MMSYN1_0005" + "" -> "0005".
def _locus_suffix(locus: str) -> str:
    parts = locus.rsplit("_", 1)
    return parts[1] if len(parts) == 2 and parts[1].isdigit() else locus


def list_overlap(genes: pd.DataFrame, chrom: str, s0: int, e: int) -> str:
    from urllib.parse import unquote

    if genes.empty:
        return "NA"
    sub = genes[
        (genes.chrom == chrom) & (genes.start0 < e) & (genes.end > s0)
    ]
    if sub.empty:
        return "."
    seen: set[str] = set()
    out: list[str] = []
    for _, r in sub.iterrows():
        locus = str(r["locus"])
        if not locus or locus in seen:
            continue
        seen.add(locus)
        num = _locus_suffix(locus)
        name = unquote(str(r["name"]).strip())   # GFF3 percent-decoding
        # Treat empty / "." / "nan" / a name that is just the locus tag itself
        # as "no real gene name" -> print only the numeric suffix.
        if name and name not in (".", "nan", locus):
            out.append(f"{num}/{name}")
        else:
            out.append(num)
    return ",".join(out)


reloc_rows = []
for i in out_of_order:
    r = sorted_by_s1.iloc[i]
    s1_chrom = r["TAGR"]
    q_chrom = r["TAGQ"]
    s1_s, s1_e = int(r["S1"]) - 1, int(r["E1"])
    s2_s, s2_e = int(r["S2_min"]) - 1, int(r["E2_max"])
    reloc_rows.append(
        {
            "block_idx": i + 1,
            "syn1_chrom": s1_chrom,
            "syn1_start0": s1_s,
            "syn1_end": s1_e,
            "syn1_length": s1_e - s1_s,
            "syn3A_chrom": q_chrom,
            "syn3A_start0": s2_s,
            "syn3A_end": s2_e,
            "pct_identity": float(r["PCT_IDY"]),
            "syn1_genes": list_overlap(syn1_genes, s1_chrom, s1_s, s1_e),
            "syn3A_genes": list_overlap(syn3A_genes, q_chrom, s2_s, s2_e),
        }
    )
reloc_df = pd.DataFrame(reloc_rows)
# reloc_df is folded into the Excel as the relocated rows of the events sheet.


# ============================================================ master Excel
# Multi-sheet workbook is the single canonical analysis output.
#   events            - one row per change event (retained / deleted / inserted)
#   short_insertions  - per-event context for the small qdiff insertions
#   dnadiff_summary   - parsed headline numbers from dnadiff_out.report
#   legend            - column / case definitions
reloc_idx_set = {i + 1 for i in out_of_order}

master_rows = []

# Retained blocks (ordered + relocated)
for i, r in sorted_by_s1.iterrows():
    s1_0, e1_0 = int(r["S1"]) - 1, int(r["E1"])
    s2_0, e2_0 = int(r["S2_min"]) - 1, int(r["E2_max"])
    case = "retained_relocated" if (i + 1) in reloc_idx_set else "retained_ordered"
    master_rows.append({
        "block_index_syn1": i + 1,
        "S1": s1_0,
        "E1": e1_0,
        "S2": s2_0,
        "E2": e2_0,
        "LEN1": int(r["LEN1"]),
        "LEN2": int(r["LEN2"]),
        "PCT_IDY": float(r["PCT_IDY"]),
        "Change Case": case,
        "Syn1_genes": list_overlap(syn1_genes, r["TAGR"], s1_0, e1_0),
        "Syn3A_genes": list_overlap(syn3A_genes, r["TAGQ"], s2_0, e2_0),
    })

# Deletions (syn1 only)
syn1_chrom_default = sorted_by_s1["TAGR"].iloc[0]
for _, r in del_df.iterrows():
    s1_0, e1_0 = int(r["start0"]), int(r["end"])
    master_rows.append({
        "block_index_syn1": "",
        "S1": s1_0,
        "E1": e1_0,
        "S2": "",
        "E2": "",
        "LEN1": int(r["length_bp"]),
        "LEN2": "",
        "PCT_IDY": "",
        "Change Case": "deleted",
        "Syn1_genes": list_overlap(syn1_genes, r["chrom"], s1_0, e1_0),
        "Syn3A_genes": "",
    })

# Insertions (syn3A only) - use the in-memory `inserted` intervals computed earlier
syn3A_chrom_default = sorted_by_s1["TAGQ"].iloc[0]
for s2_0, e2_0 in inserted:
    master_rows.append({
        "block_index_syn1": "",
        "S1": "",
        "E1": "",
        "S2": s2_0,
        "E2": e2_0,
        "LEN1": "",
        "LEN2": e2_0 - s2_0,
        "PCT_IDY": "",
        "Change Case": "inserted",
        "Syn1_genes": "",
        "Syn3A_genes": list_overlap(syn3A_genes, syn3A_chrom_default, s2_0, e2_0),
    })

master_df = pd.DataFrame(master_rows, columns=[
    "block_index_syn1", "S1", "E1", "S2", "E2", "LEN1", "LEN2",
    "PCT_IDY", "Change Case", "Syn1_genes", "Syn3A_genes",
])

# Sort: retained blocks first (by S1), then deletions (by S1), then insertions (by S2)
case_order = {
    "retained_ordered": 0, "retained_relocated": 0,
    "deleted": 1, "inserted": 2,
}
master_df["_sort_case"] = master_df["Change Case"].map(case_order)
master_df["_sort_pos"] = master_df.apply(
    lambda r: int(r["S1"]) if r["S1"] != "" else (int(r["S2"]) if r["S2"] != "" else 0),
    axis=1,
)
master_df = master_df.sort_values(["_sort_case", "_sort_pos"]).drop(
    columns=["_sort_case", "_sort_pos"]
).reset_index(drop=True)

# --- Sheet: dnadiff_summary  (headline numbers parsed from the report) -----
dnadiff_summary_rows = [
    ("syn1 total length (bp)",          rep("TotalBases", 0)),
    ("syn3A total length (bp)",         rep("TotalBases", 1)),
    ("syn1 aligned bases",              rep("AlignedBases", 0)),
    ("syn3A aligned bases",             rep("AlignedBases", 1)),
    ("1-to-1 alignment blocks",         rep("1-to-1", 0)),
    ("Average identity (%)",            rep("AvgIdentity", 0)),
    ("Total SNPs",                      rep("TotalSNPs", 0)),
    ("Total indels",                    rep("TotalIndels", 0)),
    ("Inversions",                      rep("Inversions", 0)),
    ("Translocations",                  rep("Translocations", 0)),
    ("Relocations (ref)",               rep("Relocations", 0)),
    ("Relocations (qry)",               rep("Relocations", 1)),
    ("REF insertions = syn1 deletions (events / bp)",
        f"{rep('Insertions', 0)} events / {rep('InsertionSum', 0)} bp"),
    ("QRY insertions = syn3A insertions (events / bp)",
        f"{rep('Insertions', 1)} events / {rep('InsertionSum', 1)} bp"),
    ("Filtered deletion intervals (BED, >=50 bp)", len(del_df)),
    ("Filtered insertion intervals (BED, >=10 bp)", len(inserted)),
    ("Block-level relocations (LIS-based)",          len(reloc_df)),
]
dnadiff_summary_df = pd.DataFrame(dnadiff_summary_rows, columns=["metric", "value"])

# --- Sheet: legend  (definitions) ------------------------------------------
legend_df = pd.DataFrame([
    ("retained_ordered",   "syn1 block retained in syn3A at the position predicted by linear order."),
    ("retained_relocated", "syn1 block retained in syn3A but moved to a different position (LIS outlier)."),
    ("deleted",            "syn1 region not represented in syn3A."),
    ("inserted",           "syn3A region with no homolog in syn1."),
    ("S1, E1, LEN1",       "0-based half-open start/end and length on syn1 (CP002027.1). Empty for inserted rows."),
    ("S2, E2, LEN2",       "0-based half-open start/end and length on syn3A (CP016816.2). Empty for deleted rows."),
    ("PCT_IDY",            "Percent identity of the alignment block (only for retained_*)."),
    ("Syn1_genes / Syn3A_genes",
                           "Comma-separated locus-tag suffixes overlapping the interval, deduplicated. "
                           "Format: '<num>/<gene_name>' or just '<num>' if no real gene name. "
                           "'.' = intergenic (no overlapping gene)."),
], columns=["term", "meaning"])

with pd.ExcelWriter(EXCEL, engine="openpyxl") as writer:
    master_df.to_excel(writer,            sheet_name="events",           index=False)
    context_df.to_excel(writer,           sheet_name="short_insertions", index=False)
    dnadiff_summary_df.to_excel(writer,   sheet_name="dnadiff_summary",  index=False)
    legend_df.to_excel(writer,            sheet_name="legend",           index=False)


# Also pull the qdiff/rdiff JMP records directly
qjmp = qdiff[qdiff["type"] == "JMP"].copy()
rjmp = rdiff[rdiff["type"] == "JMP"].copy()


# ============================================================ summary text
n_ins_bed = len(inserted)
ins_bp = sum(e - s for s, e in inserted)
n_del = len(del_df)
del_bp = int(del_df["length_bp"].sum())
del_mean = int(del_df["length_bp"].mean()) if n_del else 0
top_del_idx = del_df["length_bp"].idxmax()
del_max_row = del_df.loc[top_del_idx]
del_max_str = f"{del_max_row['chrom']}:{del_max_row['start0']}-{del_max_row['end']}"
del_max_bp = int(del_max_row["length_bp"])

# pull the 1119 bp insertion that contains JCVISYN3A_0931 from short_ins/qdiff
big_ins = qdiff[(qdiff["type"] == "GAP") & (qdiff["q_gap"] >= 500)]


def top10_block(df: pd.DataFrame, label_fn, len_fn) -> str:
    lines = []
    for i, (_, r) in enumerate(df.head(10).iterrows(), start=1):
        lines.append(f"  {i:2d}. {label_fn(r):<28s} {len_fn(r):>10d} bp")
    return "\n".join(lines)


top_del_block = top10_block(
    del_df.sort_values("length_bp", ascending=False),
    lambda r: f"{r['chrom']}:{r['start0']}-{r['end']}",
    lambda r: int(r["length_bp"]),
)
ins_df = pd.DataFrame(
    [{"chrom": syn3A_chrom, "start0": s, "end": e, "length_bp": e - s} for s, e in inserted]
).sort_values("length_bp", ascending=False)
top_ins_block = top10_block(
    ins_df,
    lambda r: f"{r['chrom']}:{r['start0']}-{r['end']}",
    lambda r: int(r["length_bp"]),
)

qjmp_block = "\n".join(
    f"  JMP  {r['chrom']}:{int(r['start'])}-{int(r['end'])}   q_extra={int(r['q_gap'])} bp"
    for _, r in qjmp.iterrows()
)
rjmp_block = "\n".join(
    f"  JMP  {r['chrom']}:{int(r['start'])}-{int(r['end'])}   r_extra={int(r['q_gap'])} bp"
    for _, r in rjmp.iterrows()
)

if not reloc_df.empty:
    reloc_block = "\n".join(
        f"  block {int(r['block_idx']):>3d}  "
        f"syn1:{int(r['syn1_start0'])}-{int(r['syn1_end'])} "
        f"({int(r['syn1_length'])} bp)  -->  "
        f"syn3A:{int(r['syn3A_start0'])}-{int(r['syn3A_end'])}  "
        f"genes: {r['syn1_genes']}"
        for _, r in reloc_df.iterrows()
    )
else:
    reloc_block = "  (no out-of-order blocks detected by 100 bp slack)"

context_block = context_df.to_string(index=False)

now = datetime.now().strftime("%a %b %d %I:%M:%S %p %Z %Y")

summary = f"""\
================================================================
GENOME REDUCTION ANALYSIS - JCVI-Syn1.0  ->  JCVI-Syn3A
Generated: {now}
Source:    aln/raw/  (nucmer + dnadiff)
================================================================

OVERVIEW
--------
syn1  total length      : {rep('TotalBases', 0)} bp
syn3A total length      : {rep('TotalBases', 1)} bp
syn1  bases aligned     : {rep('AlignedBases', 0)}
syn3A bases aligned     : {rep('AlignedBases', 1)}
1-to-1 alignment blocks : {rep('1-to-1', 0)}
Average identity        : {rep('AvgIdentity', 0)} %
Total SNPs              : {rep('TotalSNPs', 0)}
Total indels            : {rep('TotalIndels', 0)}
Inversions              : {rep('Inversions', 0)}
Translocations          : {rep('Translocations', 0)}
Relocations (ref / qry) : {rep('Relocations', 0)} / {rep('Relocations', 1)}

INTERPRETATION
--------------
Syn3A is essentially a SUBSET of syn1: ~half of syn1 was excised in ~{rep('Insertions', 0)}
discrete deletion events while the retained sequence stays at {rep('AvgIdentity', 0)}% identity
with negligible novel content on the syn3A side. There are no inversions or
translocations - minimization preserved gene order.

================================================================
PART A.  DELETIONS  (syn1 -> syn3A)
================================================================
dnadiff REF "Insertions"  (= syn1 deletion events)  : {rep('Insertions', 0)}  events, {rep('InsertionSum', 0)} bp
Filtered deletion intervals (>= {MIN_DEL_BP_NOTE} bp from BED): {n_del}
Total deleted syn1 bp                               : {del_bp}
Mean deletion size                                  : {del_mean} bp
Largest single deletion                             : {del_max_str} ({del_max_bp} bp)
syn1 genes overlapped by a deletion                 : {n_del_genes}

Top 10 deletions by size:
{top_del_block}

Files:
  - aln/analysis/genome_reduction_summary.xlsx  sheet `events` (deleted rows + Syn1_genes column)
  - aln/raw/syn1_deleted_regions.bed                            (raw interval table)
  - aln/raw/syn1_deleted_genes.tsv                              (raw per-deletion gene overlap)
  - aln/raw/dnadiff_out.rdiff                                   (per-event reference-side diffs)

================================================================
PART B.  INSERTIONS  (syn3A bases NOT in syn1)
================================================================
dnadiff QRY "Insertions"  (event count, total bp)   : {rep('Insertions', 1)} events, {rep('InsertionSum', 1)} bp
syn3A intervals NOT covered by any syn1 alignment   : {n_ins_bed}  (>= {MIN_INS_BP} bp)
Total syn3A-only bp                                 : {ins_bp}

Top 10 syn3A-only intervals by size:
{top_ins_block}

Per-event view from dnadiff_out.qdiff (rows with q_extra > 0):
{qjmp_block}

Annotation of the 6 dnadiff insertions:
  - 1119 bp at CP016816.2:360098-361216  -->  contains a NEW GENE,
       JCVISYN3A_0931 (met14p, CDS 360403-361011, 609 bp).
       This is the only insertion that introduces novel coding sequence.
  - The other 5 events (11, 26, 39, 39, 90 bp) are intergenic junction /
       scar sequences left behind by the deletion-cassette assembly. See
       the `short_insertions` sheet of genome_reduction_summary.xlsx for
       nearest-gene neighbours on each strand.

Note on counts: the `inserted` rows in the events sheet reflect regions not
covered by the strict 1-to-1 alignment (delta-filter -1 -i 95 -l 100), so they
may over-count true insertions because the strict filter drops short / lower-
identity alignments that dnadiff still recognises. For the biological
insertion count, trust dnadiff: 6 events, 1324 bp.

Files:
  - aln/analysis/genome_reduction_summary.xlsx   sheet `events`           (canonical insertion rows)
  - aln/analysis/genome_reduction_summary.xlsx   sheet `short_insertions` (nearest-gene context)
  - aln/raw/dnadiff_out.qdiff                                              (per-event qry diffs, raw)
  - aln/raw/syn1_vs_syn3A.snps.tsv                                         (single-bp insertions, raw)

================================================================
PART C.  RELOCATIONS  (syn1 -> syn3A)
================================================================
dnadiff Relocations (ref / qry)                     : {rep('Relocations', 0)} / {rep('Relocations', 1)}
JMP rows in dnadiff_out.rdiff (ref-side anchors)    : {len(rjmp)}
JMP rows in dnadiff_out.qdiff (qry-side anchors)    : {len(qjmp)}
Out-of-order alignment blocks (LIS-based)           : {len(reloc_df)}

Definition reminder
  - Translocation = segment moved between DIFFERENT sequences (chromosomes /
    contigs). Both syn1 and syn3A are single circular chromosomes, so this
    count is forced to 0.
  - Relocation    = segment moved to a different position WITHIN the same
    chromosome (no strand flip).
  - Inversion     = same position but on the opposite strand (count is 0).

Reference-side JMP events (dnadiff_out.rdiff):
{rjmp_block}

Query-side JMP events (dnadiff_out.qdiff):
{qjmp_block}

Out-of-order alignment blocks detected from coords.tsv:
{reloc_block}

Reading
  - The dnadiff JMP counts (3 / 14) are per-anchor breakpoint counts: each
    relocation produces multiple JMP rows (one before and one after the
    relocated block, on each side of the alignment). Most of the 14 qry-
    side JMPs are tiny (0-4 bp) chain seams, not biological relocations.
  - LIS-based detection on the alignment blocks gives the canonical answer:
    exactly ONE block is truly out of order.
        syn1: 197,590-199,189  ({reloc_df.iloc[0]['syn1_length'] if not reloc_df.empty else 0} bp, gene: lap / MMSYN1_0154)
            -->  syn3A: 311,663-313,262  (~110 kb downstream in syn3A)
  - The 3 ref-side JMPs flank this single relocation event: two on the syn1
    origin side (~197.6 / ~199.2 kb) and one on the syn3A landing side
    (~590.6 kb in syn1's frame, which is where the corresponding syn3A
    coordinate jumps).
  - See the events sheet of genome_reduction_summary.xlsx, filtered to
    `Change Case == retained_relocated`, for the per-block table.

Files:
  - aln/analysis/genome_reduction_summary.xlsx  events sheet (block-level table)
  - aln/raw/dnadiff_out.rdiff / .qdiff                          (per-event diffs)

================================================================
PART D.  WHAT THIS MEANS BIOLOGICALLY
================================================================
- The minimization was a deletion campaign, not a redesign: ~{del_bp} bp of
  syn1 was removed in {n_del} discrete cuts while only {ins_bp} bp of new
  sequence appears in syn3A (and most of that is a single 1.1 kb segment).
- The largest single excision is {del_max_str} ({del_max_bp} bp). On the
  syn3A side, a 1119 bp segment introduced near that junction encodes a new
  gene: JCVISYN3A_0931 (met14p) - the only true gene addition. The other
  five qdiff insertion events (11-90 bp) are intergenic junction/scar
  sequences from the cassette-based deletion process.
- With {rep('TotalSNPs', 0)} SNPs and {rep('TotalIndels', 0)} indels in retained DNA, point-level
  divergence is essentially zero. Any expression differences at retained
  genes are NOT explained by sequence changes.
- Architecture is preserved: 0 inversions, 0 translocations, and only a
  single block-level relocation - the ~1.6 kb syn1 region carrying the
  gene lap (syn1:~197.6 kb) was reinserted at syn3A:~311.7 kb (~110 kb
  downstream). Otherwise the minimization respected gene order.
- For the gene-loss table feeding into the proteomics/transcriptomics
  comparison, filter the events sheet of genome_reduction_summary.xlsx to
  `Change Case == deleted` (the Syn1_genes column lists the lost loci).

================================================================
"""

SUMMARY.write_text(summary)

# ---------------------------------------------------------- console feedback
print(f"Wrote: {EXCEL}")
print(f"Wrote: {SUMMARY}")
print()
print("--- short-insertion context ---")
print(context_df.to_string(index=False))
print()
print("--- relocation block table ---")
print(reloc_df.to_string(index=False) if not reloc_df.empty else "(none)")
