#!/usr/bin/env python
# coding: utf-8
"""
PacBio cDNA HiFi preprocessing pipeline for JCVI-Syn1.0.

Input : demultiplexed cDNA HiFi reads (PacBio Iso-seq)
Output: strand-oriented, primer-trimmed, polyA-trimmed FLNC FASTQ (then mapped
        externally via 01_map_sort.sh and QC-filtered to HQ BAM).

Stages:
  00  orient reads to sense strand (via H1 / BCRC detection)
  01  trim H2 (3') and BC (5') primers with zero-mismatch tolerance
  02  trim trailing poly(A) tail (<= 30 bp)
  03  quick-look distribution of trailing-A length
  04  map FLNC -> genome (external 01_map_sort.sh)
  05  per-read QC table + HQ BAM filter (MAPQ/aln_frac/clip_frac + concatemer)
  06  post-filter QC report
  07  sequencing depth export (Run 02_sequencing_depth.sh afterwards)

Shared pigz/Progress helpers live in pacbio_io.py.

Author: Enguang Fu
"""

from __future__ import annotations

import io
import os
import random
import re
import shutil
import signal
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam

from pacbio_io import pigz_reader, pigz_writer, Progress, revcomp


# ============================================================================
# Paths
# ============================================================================
HOME_DIR = "../../../"   # project root; this script lives at Syn1_Transcriptomics/PacBio/PacBio_Processing/
WORK_DIR = "."           # outputs land alongside this script
INTERMEDIATE_DIR = os.path.join(WORK_DIR, "intermediate_files")
os.makedirs(INTERMEDIATE_DIR, exist_ok=True)

# Raw input: merged demultiplexed HiFi reads
MERGED_RAW_FASTQ = os.path.join(HOME_DIR, "Syn1_Transcriptomics/PacBio/PacBio_Raw/merged.hifi_reads.fastq.gz")

# Per-stage outputs
ORIENTED_FASTQ = os.path.join(INTERMEDIATE_DIR, "all.oriented.fastq.gz")
HB_FASTQ       = os.path.join(INTERMEDIATE_DIR, "all.hb.fastq.gz")       # H2 + BC trimmed
FLNC_FASTQ     = os.path.join(INTERMEDIATE_DIR, "all.FLNC.fastq.gz")     # polyA trimmed
TRAILING_A_TSV = os.path.join(INTERMEDIATE_DIR, "all.trailingA.tsv")

# Primer sequences
# After demultiplexing each raw read carries one of two primer configurations:
#   antisense cDNA: 5' H1 + TTTTTT... + RC(RNA) + BCRC  3'
#   sense    cDNA: 5' BC + RNA + AAAAAA... + H2         3'
H1   = "TAAGCAGTGGTATCAACGCAGAGTAC"   # 26 nt, oligo-dT handle, antisense 5'
H2   = "GTACTCTGCGTTGATACCACTGCTTA"   # sense 3'  (== revcomp(H1))
BC   = "TGCAATGAAGTCGCAGGGTTGGG"      # 23 nt barcode, sense 5'
BCRC = "CCCAACCCTGCGACTTCATTGCA"      # antisense 3'


# ============================================================================
# Stage 00  Orientation
# ============================================================================
# Flip antisense reads to sense, then enforce the polyA + H2 suffix as
# evidence of a full-length sense cDNA. Strict: suffix match is exact (no
# mismatches) — this is the step that drops 2.94M -> 2.69M reads and is
# intentionally tight.

def orient_fastq_stream(
    fin: io.BufferedReader,
    fout: io.BufferedWriter,
    *,
    H1: bytes,
    BCRC: bytes,
    W: int = 80,
    flush_every: int = 5000,
    progress: Optional[Progress] = None,
) -> Tuple[int, int, int, int]:
    """Orient to sense strand, then keep only reads ending in 10A + revcomp(H1).

    Returns (total, h1_flips, bc_flips, filtered_out).
    """
    required_suffix = (b"A" * 10) + revcomp(H1)

    total = h1_flips = bc_flips = filtered_out = written = 0
    readline = fin.readline
    write = fout.write
    join = b"".join
    buf: List[bytes] = []
    append = buf.append

    if progress is None:
        progress = Progress(enabled=False)
    progress.start(0)

    while True:
        header = readline()
        if not header:
            break
        seq = readline(); plus = readline(); qual = readline()
        if not qual:
            break

        total += 1
        header = header.rstrip(b"\r\n")
        seq    = seq.rstrip(b"\r\n")
        plus   = plus.rstrip(b"\r\n")
        qual   = qual.rstrip(b"\r\n")

        head = seq[:W] if len(seq) > W else seq
        if head.find(H1) != -1:
            s2, q2 = revcomp(seq), qual[::-1]
            h1_flips += 1
        elif seq.endswith(BCRC):
            s2, q2 = revcomp(seq), qual[::-1]
            bc_flips += 1
        else:
            s2, q2 = seq, qual

        if len(s2) != len(q2):
            m = min(len(s2), len(q2))
            s2, q2 = s2[:m], q2[:m]

        if not s2.endswith(required_suffix):
            filtered_out += 1
            progress.maybe_report(total, written=written, h1_flips=h1_flips,
                                  bc_flips=bc_flips, filtered=filtered_out)
            continue

        append(header + b"\n"); append(s2 + b"\n")
        append(plus + b"\n");   append(q2 + b"\n")
        written += 1

        if total % flush_every == 0 and buf:
            write(join(buf)); buf.clear()

        progress.maybe_report(total, written=written, h1_flips=h1_flips,
                              bc_flips=bc_flips, filtered=filtered_out)

    if buf:
        write(join(buf))
    progress.final(total, written=written, h1_flips=h1_flips,
                   bc_flips=bc_flips, filtered=filtered_out)
    return total, h1_flips, bc_flips, filtered_out


def orient_fastq_pigz(
    in_gz: str,
    out_gz: str,
    *,
    H1: str,
    BCRC: str,
    W: int = 80,
    threads_in: int = 16,
    threads_out: int = 16,
    compresslevel: int = 6,
    progress_every_seconds: float = 5.0,
) -> Tuple[int, int, int, int]:
    prog = Progress(enabled=True, every_seconds=progress_every_seconds, label="orient")
    with pigz_reader(in_gz, threads=threads_in) as fin, \
         pigz_writer(out_gz, threads=threads_out, compresslevel=compresslevel) as fout:
        return orient_fastq_stream(
            fin, fout,
            H1=H1.encode("ascii"),
            BCRC=BCRC.encode("ascii"),
            W=W,
            progress=prog,
        )


# ============================================================================
# Stage 01  Primer trimming (H2 at 3', BC at 5')
# ============================================================================
# Zero mismatches on purpose: guarantees retained reads have both primers
# intact. Reads missing either primer are dropped.

def count_mismatches(a: bytes, b: bytes) -> int:
    mism = 0
    for x, y in zip(a, b):
        if x != y:
            mism += 1
    return mism


def trim_front_primer_bytes(
    seq: bytes, qual: bytes, primer: bytes,
    mismatches: int, min_overlap: int, slack: int = 0,
) -> Tuple[bytes, bytes, bool]:
    """Trim primer from 5' end with up to `mismatches` mismatches."""
    if len(seq) < min_overlap:
        return seq, qual, False
    w = min(len(seq), len(primer) + slack)
    seq5 = seq[:w]
    max_i = min(len(primer), len(seq5))
    for i in range(max_i, min_overlap - 1, -1):
        if count_mismatches(seq5[:i], primer[:i]) <= mismatches:
            return seq[i:], qual[i:], True
    return seq, qual, False


def trim_tail_primer_bytes(
    seq: bytes, qual: bytes, primer: bytes,
    mismatches: int, min_overlap: int, slack: int = 0,
) -> Tuple[bytes, bytes, bool]:
    """Trim primer from 3' end with up to `mismatches` mismatches."""
    if len(seq) < min_overlap:
        return seq, qual, False
    w = min(len(seq), len(primer) + slack)
    seq3 = seq[-w:]
    max_i = min(len(primer), len(seq3))
    for i in range(max_i, min_overlap - 1, -1):
        if count_mismatches(seq3[-i:], primer[-i:]) <= mismatches:
            return seq[:-i], qual[:-i], True
    return seq, qual, False


def trim_primers_stream(
    fin: io.BufferedReader,
    fout: io.BufferedWriter,
    *,
    H2: bytes,
    BC: bytes,
    min_len: int = 10,
    min_overlap: int = 18,
    h2_mismatches: int = 0,
    bc_mismatches: int = 0,
    flush_every: int = 5000,
    progress: Optional[Progress] = None,
) -> Tuple[int, int, int, int, int]:
    """Trim H2 (3') then BC (5'); drop reads missing either primer.

    Returns (count_in, count_out, h2_trimmed, bc_trimmed, dropped).
    """
    count_in = count_out = 0
    h2_trimmed = bc_trimmed = dropped = 0

    readline = fin.readline
    write = fout.write
    join = b"".join
    buf: List[bytes] = []
    append = buf.append

    if progress is None:
        progress = Progress(enabled=False)
    progress.start(0)

    while True:
        header = readline()
        if not header:
            break
        seq = readline(); plus = readline(); qual = readline()
        if not qual:
            break

        count_in += 1
        header = header.rstrip(b"\r\n")
        seq    = seq.rstrip(b"\r\n")
        plus   = plus.rstrip(b"\r\n")
        qual   = qual.rstrip(b"\r\n")

        if len(seq) != len(qual):
            m = min(len(seq), len(qual))
            seq, qual = seq[:m], qual[:m]

        # 1) H2 at 3'
        seq2, qual2, did_h2 = trim_tail_primer_bytes(
            seq, qual, H2, mismatches=h2_mismatches, min_overlap=min_overlap,
        )
        if did_h2:
            h2_trimmed += 1
        else:
            dropped += 1
            progress.maybe_report(count_in, written=count_out, h2_trim=h2_trimmed,
                                  bc_trim=bc_trimmed, dropped=dropped)
            continue

        # 2) BC at 5'
        seq3, qual3, did_bc = trim_front_primer_bytes(
            seq2, qual2, BC, mismatches=bc_mismatches, min_overlap=min_overlap,
        )
        if did_bc:
            bc_trimmed += 1
        else:
            dropped += 1
            progress.maybe_report(count_in, written=count_out, h2_trim=h2_trimmed,
                                  bc_trim=bc_trimmed, dropped=dropped)
            continue

        # 3) min length
        if len(seq3) < min_len:
            dropped += 1
            progress.maybe_report(count_in, written=count_out, h2_trim=h2_trimmed,
                                  bc_trim=bc_trimmed, dropped=dropped)
            continue

        append(header + b"\n"); append(seq3 + b"\n")
        append(plus + b"\n");   append(qual3 + b"\n")
        count_out += 1

        if count_in % flush_every == 0 and buf:
            write(join(buf)); buf.clear()

        progress.maybe_report(count_in, written=count_out, h2_trim=h2_trimmed,
                              bc_trim=bc_trimmed, dropped=dropped)

    if buf:
        write(join(buf))
    progress.final(count_in, written=count_out, h2_trim=h2_trimmed,
                   bc_trim=bc_trimmed, dropped=dropped)
    return count_in, count_out, h2_trimmed, bc_trimmed, dropped


def trim_primers_pigz(
    input_fastq_gz: str,
    output_fastq_gz: str,
    *,
    H2: str,
    BC: str,
    min_len: int = 10,
    min_overlap: int = 18,
    h2_mismatches: int = 0,
    bc_mismatches: int = 0,
    threads_in: int = 16,
    threads_out: int = 16,
    compresslevel: int = 6,
    progress_every_seconds: float = 5.0,
) -> Tuple[int, int, int, int, int]:
    prog = Progress(enabled=True, every_seconds=progress_every_seconds, label="trim")
    with pigz_reader(input_fastq_gz, threads=threads_in) as fin, \
         pigz_writer(output_fastq_gz, threads=threads_out, compresslevel=compresslevel) as fout:
        return trim_primers_stream(
            fin, fout,
            H2=H2.encode("ascii"),
            BC=BC.encode("ascii"),
            min_len=min_len,
            min_overlap=min_overlap,
            h2_mismatches=h2_mismatches,
            bc_mismatches=bc_mismatches,
            progress=prog,
        )


# ============================================================================
# Stage 02  Trailing poly(A) trim
# ============================================================================
# After primer removal ~30 trailing As remain on most reads. Count and trim
# up to `trim_A_limit` bases per read; write per-read A-tail length to TSV.

def count_trailing_A_bytes(seq: bytes) -> int:
    n = 0
    for b in reversed(seq):
        if b == 65 or b == 97:  # 'A' or 'a'
            n += 1
        else:
            break
    return n


def count_trim_A_tail_in_fastq_pigz(
    fastq_gz: str,
    *,
    output_tsv: str,
    trim_A: bool = False,
    trim_A_limit: int = 35,
    max_A_limit: int = 70,
    output_clean_fastq_gz: str = "FLNC.fastq.gz",
    threads_in: int = 16,
    threads_out: int = 16,
    compresslevel: int = 6,
    progress_every_seconds: float = 5.0,
    bufsize: int = 1024 * 1024,
    flush_every: int = 5000,
    store_counts_in_memory: bool = True,
) -> List[int]:
    """Count and optionally trim the trailing A-run of each read.

    Reads whose trailing A-run exceeds `max_A_limit` are dropped entirely
    (likely concatemers or internal-priming artifacts; Syn1 has no native
    A-stretch that long). Retained reads have up to `trim_A_limit` trailing
    As removed.
    """
    prog = Progress(enabled=True, every_seconds=progress_every_seconds, label="trailingA")

    counts: List[int] = [] if store_counts_in_memory else []
    trimmed_written = 0
    dropped_long = 0
    n_reads = 0
    tsv_f = open(output_tsv, "wt", buffering=1024 * 1024)

    with pigz_reader(fastq_gz, threads=threads_in, bufsize=bufsize) as fin:
        readline = fin.readline

        out_ctx = (
            pigz_writer(output_clean_fastq_gz, threads=threads_out,
                        bufsize=bufsize, compresslevel=compresslevel)
            if trim_A else None
        )

        try:
            if out_ctx is None:
                prog.start(0)
                while True:
                    header = readline()
                    if not header:
                        break
                    seq = readline(); plus = readline(); qual = readline()
                    if not qual:
                        break
                    n_reads += 1
                    seq_b = seq.rstrip(b"\r\n")
                    nA = count_trailing_A_bytes(seq_b)
                    if store_counts_in_memory:
                        counts.append(nA)
                    tsv_f.write(f"{n_reads}\t{nA}\n")
                    prog.maybe_report(n_reads, trimmed_written=trimmed_written)
                prog.final(n_reads, trimmed_written=trimmed_written)

            else:
                with out_ctx as fout:
                    write = fout.write
                    join = b"".join
                    buf: List[bytes] = []
                    append = buf.append
                    prog.start(0)
                    while True:
                        header = readline()
                        if not header:
                            break
                        seq = readline(); plus = readline(); qual = readline()
                        if not qual:
                            break
                        n_reads += 1
                        header_b = header.rstrip(b"\r\n")
                        seq_b    = seq.rstrip(b"\r\n")
                        plus_b   = plus.rstrip(b"\r\n")
                        qual_b   = qual.rstrip(b"\r\n")

                        if len(seq_b) != len(qual_b):
                            m = min(len(seq_b), len(qual_b))
                            seq_b, qual_b = seq_b[:m], qual_b[:m]

                        nA = count_trailing_A_bytes(seq_b)
                        if store_counts_in_memory:
                            counts.append(nA)
                        tsv_f.write(f"{n_reads}\t{nA}\n")

                        if nA > max_A_limit:
                            dropped_long += 1
                            prog.maybe_report(
                                n_reads, written=trimmed_written, dropped_long=dropped_long
                            )
                            continue

                        trim_count = min(nA, trim_A_limit)
                        if trim_count:
                            seq_b = seq_b[:-trim_count]
                            qual_b = qual_b[:-trim_count]

                        append(header_b + b"\n"); append(seq_b + b"\n")
                        append(plus_b + b"\n");   append(qual_b + b"\n")
                        trimmed_written += 1

                        if n_reads % flush_every == 0 and buf:
                            write(join(buf)); buf.clear()

                        prog.maybe_report(
                            n_reads, written=trimmed_written, dropped_long=dropped_long
                        )
                    if buf:
                        write(join(buf))
                    prog.final(
                        n_reads, written=trimmed_written, dropped_long=dropped_long
                    )
        finally:
            tsv_f.close()

    return counts if store_counts_in_memory else []


# ============================================================================
# Stage 05  BAM QC helpers
# ============================================================================
CIGAR_TOKEN_RE = re.compile(r"(\d+)([MIDNSHP=XB])")


def cigar_stats(cigar: str) -> Tuple[int, int, int]:
    """Return (qlen, alnq, clip) from a CIGAR string."""
    if cigar in ("*", "", None):
        return (0, 0, 0)
    qlen = alnq = clip = 0
    for n_str, op in CIGAR_TOKEN_RE.findall(cigar):
        n = int(n_str)
        if op in ("M", "I", "S", "=", "X"):
            qlen += n
        if op in ("M", "I", "=", "X"):
            alnq += n
        if op == "S":
            clip += n
    return (qlen, alnq, clip)


def first_last_clip_from_cigar_tuples(cigartuples) -> Tuple[int, int]:
    if not cigartuples:
        return 0, 0
    left = cigartuples[0][1] if cigartuples[0][0] in (4, 5) else 0
    right = cigartuples[-1][1] if cigartuples[-1][0] in (4, 5) else 0
    return left, right


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def interval_len(intervals) -> int:
    return sum(e - s for s, e in intervals)


def _run(cmd, check=True) -> str:
    print("+", " ".join(map(str, cmd)))
    p = subprocess.run(list(map(str, cmd)), check=check,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.stderr.strip():
        print(p.stderr.strip())
    return p.stdout


def _flagstat_counts(bam_path: Path, threads: int = 8) -> dict:
    """Parse `samtools flagstat` once instead of 5 separate `view -c` calls."""
    out = _run(["samtools", "flagstat", "-@", str(threads), str(bam_path)])
    # Extract primary counts from flagstat lines; field 0 is the count.
    def _n(label_substr: str) -> int:
        for line in out.splitlines():
            if label_substr in line:
                return int(line.split()[0])
        return 0
    total          = _n("in total")
    primary_mapped = _n("primary mapped")
    if primary_mapped == 0:                 # fallback for old samtools
        primary_mapped = _n(" mapped (")
    mapped         = _n(" mapped (")
    secondary      = _n("secondary")
    supplementary  = _n("supplementary")
    unmapped       = total - mapped
    return {
        "total_alignments": total,
        "mapped_alignments": mapped,
        "unmapped_alignments": unmapped,
        "primary_mapped": primary_mapped,
        "secondary": secondary,
        "supplementary": supplementary,
    }


def _summarize(arr) -> dict:
    if len(arr) == 0:
        return {"n": 0}
    a = np.asarray(arr, dtype=float)
    return {
        "n": int(a.size),
        "min": float(np.min(a)),
        "p05": float(np.quantile(a, 0.05)),
        "p50": float(np.quantile(a, 0.50)),
        "p95": float(np.quantile(a, 0.95)),
        "max": float(np.max(a)),
        "mean": float(np.mean(a)),
    }


def qc_report_from_df(bam_path: Path, read_qc_df: pd.DataFrame, threads: int = 8) -> dict:
    """Zero-rescan QC summary: counts from flagstat, distributions from the
    already-built per-read QC DataFrame (primary alignments only)."""
    counts = _flagstat_counts(Path(bam_path), threads=threads)
    return {
        "bam": str(bam_path),
        "counts": counts,
        "sample_n_used": int(len(read_qc_df)),
        "mapq":              _summarize(read_qc_df["primary_mapq"]),
        "aligned_fraction":  _summarize(read_qc_df["primary_aln_frac"].dropna()),
        "softclip_fraction": _summarize(read_qc_df["primary_clip_frac"].dropna()),
    }


def qc_report(bam_path: Path, sample_n: int = 200000, seed: int = 1,
              threads: int = 8) -> dict:
    """Sampling QC report (used for the post-filter BAM).

    Uses `samtools flagstat` for global counts (one pass) and reservoir-
    samples primary alignments for distribution summaries.
    """
    bam_path = Path(bam_path)
    assert bam_path.exists(), f"Missing BAM: {bam_path}"
    random.seed(seed)

    counts = _flagstat_counts(bam_path, threads=threads)

    p = subprocess.Popen(
        ["samtools", "view", "-@", str(threads), "-F", "2308", str(bam_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )

    mapq_vals: List[int] = []
    alnfrac_vals: List[float] = []
    clipfrac_vals: List[float] = []
    seen = kept = 0

    for line in p.stdout:
        seen += 1
        take = kept < sample_n or (random.random() < (sample_n / float(seen)))
        if not take:
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        mapq = int(fields[4]); cigar = fields[5]
        qlen, alnq, clip = cigar_stats(cigar)
        if qlen <= 0:
            continue
        alnfrac = alnq / qlen
        clipfrac = clip / qlen
        if kept < sample_n:
            mapq_vals.append(mapq); alnfrac_vals.append(alnfrac); clipfrac_vals.append(clipfrac)
            kept += 1
        else:
            j = random.randrange(sample_n)
            mapq_vals[j] = mapq; alnfrac_vals[j] = alnfrac; clipfrac_vals[j] = clipfrac

    _, stderr = p.communicate()
    if p.returncode != 0:
        raise RuntimeError(f"samtools view failed on {bam_path}:\n{stderr}")

    return {
        "bam": str(bam_path),
        "counts": counts,
        "sample_n_used": len(mapq_vals),
        "mapq": _summarize(mapq_vals),
        "aligned_fraction": _summarize(alnfrac_vals),
        "softclip_fraction": _summarize(clipfrac_vals),
    }


def print_qc(report: dict, title: str) -> None:
    print("=" * 90); print(title); print("=" * 90)
    c = report["counts"]
    print(f"BAM: {report['bam']}")
    print(f"Total alignments    : {c['total_alignments']:,}")
    print(f"Mapped alignments   : {c['mapped_alignments']:,}")
    print(f"Unmapped alignments : {c['unmapped_alignments']:,}")
    print(f"Primary mapped      : {c['primary_mapped']:,}")
    print(f"Secondary           : {c['secondary']:,}")
    print(f"Supplementary       : {c['supplementary']:,}")
    print(f"Sample used         : {report['sample_n_used']:,} (from primary mapped)")
    for k in ("mapq", "aligned_fraction", "softclip_fraction"):
        s = report[k]
        if s.get("n", 0) == 0:
            print(f"{k:18s}: n=0"); continue
        print(
            f"{k:18s}: min={s['min']:.4g}  p05={s['p05']:.4g}  p50={s['p50']:.4g}  "
            f"p95={s['p95']:.4g}  max={s['max']:.4g}  mean={s['mean']:.4g}"
        )


def build_read_qc_table(
    in_bam: Path,
    *,
    min_mapq: int,
    min_aln_frac: float,
    max_clip_frac: float,
    mismatch_bp_min: int,
    mismatch_frac_min: float,
    high_query_coverage_frac: float,
    multi_seg_min: int,
    repeat_locus_ratio: float,
    distant_gap_bp: int,
    max_qlen_refspan_delta: int,
    threads: int = 8,
) -> pd.DataFrame:
    """One pass over the BAM: group by qname (primary + supplementary),
    compute per-read mapping QC and concatemer heuristics, return a DataFrame.

    `max_qlen_refspan_delta`: reads with |qlen - primary_ref_span| exceeding
    this are labeled FAIL_QLEN_REFSPAN_DELTA (large negative = chimeric
    split across distant loci; large positive = unexplained query residue).
    """
    bam = pysam.AlignmentFile(str(in_bam), "rb", threads=threads)
    read_groups = defaultdict(list)

    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary:
            continue
        qname = aln.query_name
        qlen = aln.infer_read_length() or aln.query_length or 0
        chrom = bam.get_reference_name(aln.reference_id)
        strand = "-" if aln.is_reverse else "+"
        ref_start0 = aln.reference_start
        ref_end0 = aln.reference_end if aln.reference_end is not None else aln.reference_start
        ref_span = max(0, ref_end0 - ref_start0)
        qstart = aln.query_alignment_start if aln.query_alignment_start is not None else 0
        qend = aln.query_alignment_end if aln.query_alignment_end is not None else qstart
        q_aln_span = max(0, qend - qstart)
        cigar = aln.cigarstring
        cigartuples = aln.cigartuples or ()
        qlen_cigar, alnq_cigar, clip_cigar = cigar_stats(cigar)
        left_clip, right_clip = first_last_clip_from_cigar_tuples(cigartuples)

        read_groups[qname].append({
            "qname": qname, "qlen": qlen, "chrom": chrom, "strand": strand,
            "ref_start0": ref_start0, "ref_end0": ref_end0, "ref_span": ref_span,
            "q_aln_start0": qstart, "q_aln_end0": qend, "q_aln_span": q_aln_span,
            "mapq": aln.mapping_quality,
            "is_primary": int(not aln.is_supplementary),
            "is_supplementary": int(aln.is_supplementary),
            "qlen_cigar": qlen_cigar, "alnq_cigar": alnq_cigar, "clip_cigar": clip_cigar,
            "left_clip": left_clip, "right_clip": right_clip,
        })
    bam.close()

    rows = []
    for qname, segs in read_groups.items():
        segs = sorted(segs, key=lambda x: (x["q_aln_start0"], x["q_aln_end0"]))
        qlen = max(s["qlen"] for s in segs) if segs else 0
        n_segments = len(segs)
        n_supplementary = sum(s["is_supplementary"] for s in segs)

        primary_segs = [s for s in segs if s["is_primary"] == 1]
        primary = (
            sorted(primary_segs, key=lambda x: (x["q_aln_span"], x["mapq"]), reverse=True)[0]
            if primary_segs else segs[0]
        )

        primary_mapq = primary["mapq"]
        pq = primary["qlen_cigar"]
        if pq > 0:
            primary_aln_frac  = primary["alnq_cigar"] / pq
            primary_clip_frac = primary["clip_cigar"] / pq
        else:
            primary_aln_frac = primary_clip_frac = np.nan

        fail_mapq      = int(primary_mapq < min_mapq)
        fail_aln_frac  = int(pd.isna(primary_aln_frac)  or primary_aln_frac  < min_aln_frac)
        fail_clip_frac = int(pd.isna(primary_clip_frac) or primary_clip_frac > max_clip_frac)
        pass_mapping_qc = int(fail_mapq == 0 and fail_aln_frac == 0 and fail_clip_frac == 0)

        q_intervals = [(s["q_aln_start0"], s["q_aln_end0"])
                       for s in segs if s["q_aln_end0"] > s["q_aln_start0"]]
        q_merged = merge_intervals(q_intervals)
        q_unique_aligned = interval_len(q_merged)
        q_aligned_frac_all = (q_unique_aligned / qlen) if qlen > 0 else np.nan

        total_ref_span = sum(s["ref_span"] for s in segs)
        chrom_strand_set = {(s["chrom"], s["strand"]) for s in segs}
        same_chr_strand_all = (len(chrom_strand_set) == 1)

        if same_chr_strand_all:
            ref_intervals = [(s["ref_start0"], s["ref_end0"]) for s in segs]
            ref_merged = merge_intervals(ref_intervals)
            merged_ref_span_same_locus = interval_len(ref_merged)
            xs = sorted(ref_intervals)
            gaps = [xs[i+1][0] - xs[i][1] for i in range(len(xs)-1)]
            max_ref_gap_same_locus = max(gaps) if gaps else 0
            repeat_locus_ratio_v = (
                total_ref_span / merged_ref_span_same_locus
                if merged_ref_span_same_locus > 0 else np.nan
            )
        else:
            merged_ref_span_same_locus = np.nan
            max_ref_gap_same_locus = np.nan
            repeat_locus_ratio_v = np.nan

        primary_ref_span = primary["ref_span"]
        mismatch_bp = qlen - primary_ref_span
        mismatch_frac = (mismatch_bp / qlen) if qlen > 0 else np.nan

        has_big_mismatch = int(
            (mismatch_bp >= mismatch_bp_min)
            and (pd.notna(mismatch_frac) and mismatch_frac >= mismatch_frac_min)
        )
        likely_terminal_clip = int(
            has_big_mismatch == 1 and n_segments == 1
            and (primary["left_clip"] + primary["right_clip"]) >= mismatch_bp_min
        )
        likely_multi_segment = int(
            has_big_mismatch == 1 and n_segments >= multi_seg_min
            and pd.notna(q_aligned_frac_all) and q_aligned_frac_all >= high_query_coverage_frac
        )
        likely_repeat_same_locus = int(
            likely_multi_segment == 1 and same_chr_strand_all
            and pd.notna(repeat_locus_ratio_v) and repeat_locus_ratio_v >= repeat_locus_ratio
        )
        likely_split_distant = int(
            likely_multi_segment == 1 and (
                (not same_chr_strand_all)
                or (pd.notna(max_ref_gap_same_locus) and max_ref_gap_same_locus >= distant_gap_bp)
            )
        )
        likely_concatemer_strong = int(likely_repeat_same_locus == 1 or likely_split_distant == 1)

        fail_qlen_refspan = int(abs(mismatch_bp) > max_qlen_refspan_delta)

        if pass_mapping_qc == 0:
            if fail_mapq == 1:
                final_label = "FAIL_MAPQ"
            elif fail_aln_frac == 1:
                final_label = "FAIL_ALN_FRAC"
            elif fail_clip_frac == 1:
                final_label = "FAIL_CLIP_FRAC"
            else:
                final_label = "FAIL_MAPPING_OTHER"
        elif fail_qlen_refspan == 1:
            final_label = "FAIL_QLEN_REFSPAN_DELTA"
        else:
            if likely_concatemer_strong == 1:
                final_label = "FAIL_CONCATEMER_STRONG"
            elif likely_terminal_clip == 1:
                final_label = "FAIL_TERMINAL_CLIP_STRUCTURE"
            elif has_big_mismatch == 1:
                final_label = "FAIL_AMBIGUOUS_STRUCTURE"
            else:
                final_label = "PASS_STRICT"

        rows.append({
            "qname": qname, "qlen": qlen,
            "n_segments": n_segments, "n_supplementary": n_supplementary,
            "primary_chrom": primary["chrom"], "primary_strand": primary["strand"],
            "primary_ref_start_1b": primary["ref_start0"] + 1,
            "primary_ref_end_1b": primary["ref_end0"],
            "primary_ref_span": primary_ref_span,
            "primary_mapq": primary_mapq,
            "primary_qlen_cigar": primary["qlen_cigar"],
            "primary_alnq_cigar": primary["alnq_cigar"],
            "primary_clip_cigar": primary["clip_cigar"],
            "primary_aln_frac": primary_aln_frac,
            "primary_clip_frac": primary_clip_frac,
            "fail_mapq": fail_mapq,
            "fail_aln_frac": fail_aln_frac,
            "fail_clip_frac": fail_clip_frac,
            "pass_mapping_qc": pass_mapping_qc,
            "q_unique_aligned": q_unique_aligned,
            "q_aligned_frac_all": q_aligned_frac_all,
            "total_ref_span": total_ref_span,
            "merged_ref_span_same_locus": merged_ref_span_same_locus,
            "same_chr_strand_all": int(same_chr_strand_all),
            "repeat_locus_ratio": repeat_locus_ratio_v,
            "max_ref_gap_same_locus": max_ref_gap_same_locus,
            "mismatch_bp_vs_primary_ref_span": mismatch_bp,
            "mismatch_frac_vs_primary_ref_span": mismatch_frac,
            "fail_qlen_refspan": fail_qlen_refspan,
            "has_big_mismatch": has_big_mismatch,
            "likely_terminal_clip": likely_terminal_clip,
            "likely_multi_segment": likely_multi_segment,
            "likely_repeat_same_locus": likely_repeat_same_locus,
            "likely_split_distant": likely_split_distant,
            "likely_concatemer_strong": likely_concatemer_strong,
            "final_label": final_label,
        })
    return pd.DataFrame(rows)


def filter_bam_hq(
    in_bam: Path, out_bam: Path, exclude_qname_txt: Path,
    *, min_mapq: int, min_aln_frac: float, max_clip_frac: float,
    threads: int, tmp_prefix: str,
) -> None:
    """Stream `samtools view | awk | samtools sort` to produce HQ BAM."""
    awk_prog = r'''
BEGIN {
  OFS="\t"
  while ((getline line < EXCLUDE_FILE) > 0) { bad[line] = 1 }
  close(EXCLUDE_FILE)
}
/^@/ { print; next }
{
  qname = $1
  if (qname in bad) next
  mapq = $5 + 0
  if (mapq < MIN_MAPQ) next
  cigar = $6
  if (cigar == "*" || cigar == "") next
  qlen = 0; alnq = 0; clip = 0
  while (match(cigar, /^([0-9]+)([MIDNSHP=XB])/, m)) {
    len = m[1] + 0; op = m[2]
    if (op ~ /[MIS=X]/) qlen += len
    if (op ~ /[MI=X]/)  alnq += len
    if (op == "S")      clip += len
    cigar = substr(cigar, RLENGTH + 1)
  }
  if (qlen <= 0) next
  aln_frac  = alnq / qlen
  clip_frac = clip / qlen
  if (aln_frac < MIN_ALN_FRAC) next
  if (clip_frac > MAX_CLIP_FRAC) next
  print
}
'''
    BASE_FILTER_FLAGS = 2308
    cmd_view = ["samtools", "view", "-h", "-F", str(BASE_FILTER_FLAGS), str(in_bam)]
    cmd_awk  = ["awk",
                f"-vMIN_MAPQ={min_mapq}",
                f"-vMIN_ALN_FRAC={min_aln_frac}",
                f"-vMAX_CLIP_FRAC={max_clip_frac}",
                f"-vEXCLUDE_FILE={str(exclude_qname_txt)}",
                awk_prog]
    cmd_sort = ["samtools", "sort", "-@", str(threads), "-T", tmp_prefix, "-o", str(out_bam)]

    print("[INFO] Running: samtools view | awk | samtools sort")
    p_view = subprocess.Popen(cmd_view, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    p_awk  = subprocess.Popen(cmd_awk,  stdin=p_view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    p_sort = subprocess.Popen(cmd_sort, stdin=p_awk.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    p_view.stdout.close(); p_awk.stdout.close()
    _, err_sort = p_sort.communicate()
    err_awk  = p_awk.stderr.read()
    err_view = p_view.stderr.read()
    rc_view, rc_awk, rc_sort = p_view.wait(), p_awk.wait(), p_sort.returncode

    SIGPIPE_RC = -signal.SIGPIPE
    if rc_view not in (0, SIGPIPE_RC):
        raise RuntimeError(f"samtools view failed (rc={rc_view}):\n{err_view}")
    if rc_awk not in (0, SIGPIPE_RC):
        raise RuntimeError(f"awk failed (rc={rc_awk}):\n{err_awk}")
    if rc_sort != 0:
        raise RuntimeError(f"samtools sort failed (rc={rc_sort}):\n{err_sort}")

    _run(["samtools", "index", "-@", str(threads), str(out_bam)])
    print("[DONE] Wrote:", out_bam)


# ============================================================================
# Driver
# ============================================================================
if __name__ == "__main__":
    # ------- Stage 00: orient -------
    total, h1_flips, bc_flips, filtered = orient_fastq_pigz(
        MERGED_RAW_FASTQ, ORIENTED_FASTQ,
        H1=H1, BCRC=BCRC, W=80,
        threads_in=16, threads_out=16, compresslevel=6,
        progress_every_seconds=5.0,
    )
    print(f"[orient] total={total:,}  h1_flips={h1_flips:,}  "
          f"bc_flips={bc_flips:,}  filtered_out={filtered:,}  "
          f"written={total - filtered:,}")

    # ------- Stage 01: primer trim (strict, 0 mismatches) -------
    count_in, count_out, n_h2, n_bc, n_drop = trim_primers_pigz(
        ORIENTED_FASTQ, HB_FASTQ,
        H2=H2, BC=BC,
        min_len=10, min_overlap=18,
        h2_mismatches=0, bc_mismatches=0,
        threads_in=16, threads_out=16, compresslevel=6,
        progress_every_seconds=5.0,
    )
    print(f"[trim]  in={count_in:,}  out={count_out:,}  "
          f"h2_trim={n_h2:,}  bc_trim={n_bc:,}  dropped={n_drop:,}")

    # ------- Stage 02: trim polyA tail -------
    # Syn1 has no native A-stretch >~10 nt; reads with >70 nt trailing A are
    # concatemers or internal-priming artifacts -> drop. Retained reads
    # have up to 35 As trimmed (99th pct of the distribution).
    _ = count_trim_A_tail_in_fastq_pigz(
        HB_FASTQ,
        output_tsv=TRAILING_A_TSV,
        trim_A=True, trim_A_limit=35, max_A_limit=70,
        output_clean_fastq_gz=FLNC_FASTQ,
        threads_in=16, threads_out=16, compresslevel=6,
        progress_every_seconds=5.0,
        store_counts_in_memory=False,
    )

    # ------- Stage 03: polyA distribution plot -------
    trailingA_df = pd.read_csv(TRAILING_A_TSV, sep="\t", header=None,
                               names=["read_idx", "nA"])
    print(trailingA_df["nA"].describe(
        percentiles=[0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))

    trailingA_counts = trailingA_df["nA"].to_numpy()
    plt.figure(figsize=(7, 4))
    x_higher = np.percentile(trailingA_counts, 99) + 2
    plt.hist(trailingA_counts, range=(0, x_higher), density=False)
    plt.xlabel("3' trailing A-run length")
    plt.ylabel("Number of reads")
    plt.title("Distribution of 3' trailing polyA length (after primer trimming)")
    plt.tight_layout()
    plt.savefig(os.path.join(INTERMEDIATE_DIR, "trailingA_distribution.pdf"))
    plt.close()

    # ------- Stage 04: mapping done externally via 01_map_sort.sh -------
    # Pause here so the user can launch mapping manually and only continue
    # to QC once the sorted BAM exists.
    IN_BAM = Path("syn1.PacBio.FLNC.sorted.bam")
    print(
        "\n" + "=" * 90 +
        f"\n[PAUSE] FLNC ready: {FLNC_FASTQ}"
        f"\n[PAUSE] Now run `bash 01_map_sort.sh` to produce {IN_BAM}."
        f"\n[PAUSE] Press ENTER once mapping is done to continue with QC, or Ctrl-C to stop.\n"
        + "=" * 90
    )
    try:
        input()
    except EOFError:
        # Non-interactive run: wait until the BAM appears on disk. Wait for the
        # .csi index too -- 01_map_sort.sh writes it only after `samtools sort`
        # has finished, so keying on the BAM alone can catch a truncated file
        # mid-write and fail in build_read_qc_table below.
        import time as _time
        _csi = Path(str(IN_BAM) + ".csi")
        while not (IN_BAM.exists() and _csi.exists()):
            print(f"[PAUSE] Waiting for {IN_BAM} (+ .csi) ...", flush=True)
            _time.sleep(30)

    # ------- Stage 05: per-read QC + HQ BAM filter -------
    OUT_BAM = Path("./syn1.PacBio.FLNC.sorted.HQ.bam")   # alongside this script after the reorg
    READ_QC_TSV       = Path("syn1.PacBio.FLNC.read_qc.tsv")
    EXCLUDE_QNAME_TXT = Path("syn1.PacBio.FLNC.exclude_concatemer_qnames.txt")
    SUMMARY_TXT       = Path("syn1.PacBio.FLNC.qc_summary.txt")

    MIN_MAPQ       = 20
    MIN_ALN_FRAC   = 0.70
    MAX_CLIP_FRAC  = 0.30
    MISMATCH_BP_MIN          = 100
    MISMATCH_FRAC_MIN        = 0.15
    HIGH_QUERY_COVERAGE_FRAC = 0.80
    MULTI_SEG_MIN            = 2
    REPEAT_LOCUS_RATIO       = 1.5
    DISTANT_GAP_BP           = 200
    MAX_QLEN_REFSPAN_DELTA   = 100   # drop reads where |qlen - primary_ref_span| > this
    THREADS = 8
    TMP_PREFIX = f"/tmp/{OUT_BAM.name}.tmp"

    assert IN_BAM.exists(), f"Input BAM not found: {IN_BAM}"
    assert shutil.which("samtools"), "samtools not on PATH"

    read_qc_df = build_read_qc_table(
        IN_BAM,
        min_mapq=MIN_MAPQ, min_aln_frac=MIN_ALN_FRAC, max_clip_frac=MAX_CLIP_FRAC,
        mismatch_bp_min=MISMATCH_BP_MIN, mismatch_frac_min=MISMATCH_FRAC_MIN,
        high_query_coverage_frac=HIGH_QUERY_COVERAGE_FRAC,
        multi_seg_min=MULTI_SEG_MIN,
        repeat_locus_ratio=REPEAT_LOCUS_RATIO,
        distant_gap_bp=DISTANT_GAP_BP,
        max_qlen_refspan_delta=MAX_QLEN_REFSPAN_DELTA,
        threads=THREADS,
    )
    read_qc_df.to_csv(READ_QC_TSV, sep="\t", index=False)

    # awk-side filter in filter_bam_hq handles MAPQ/aln_frac/clip_frac; the
    # exclude list covers the structural/concatemer/qlen-delta fails that
    # awk cannot catch.
    EXCLUDE_LABELS = {
        "FAIL_CONCATEMER_STRONG",
        "FAIL_TERMINAL_CLIP_STRUCTURE",
        "FAIL_AMBIGUOUS_STRUCTURE",
        "FAIL_QLEN_REFSPAN_DELTA",
    }
    exclude_qnames = set(
        read_qc_df.loc[read_qc_df["final_label"].isin(EXCLUDE_LABELS), "qname"].tolist()
    )
    with open(EXCLUDE_QNAME_TXT, "w") as f:
        for qn in sorted(exclude_qnames):
            f.write(qn + "\n")

    print(f"Wrote read-level QC table: {READ_QC_TSV}")
    print(f"Wrote structural exclude list: {EXCLUDE_QNAME_TXT}  (n={len(exclude_qnames):,})")
    print("\nFinal label counts:")
    print(read_qc_df["final_label"].value_counts())

    # Reuse read_qc_df for pre-filter distributions -> avoids a full BAM rescan.
    qc_before = qc_report_from_df(IN_BAM, read_qc_df, threads=THREADS)
    print_qc(qc_before, "QC BEFORE FILTERING")

    filter_bam_hq(
        IN_BAM, OUT_BAM, EXCLUDE_QNAME_TXT,
        min_mapq=MIN_MAPQ, min_aln_frac=MIN_ALN_FRAC, max_clip_frac=MAX_CLIP_FRAC,
        threads=THREADS, tmp_prefix=TMP_PREFIX,
    )

    qc_after = qc_report(OUT_BAM, sample_n=200000, seed=1, threads=THREADS)
    print_qc(qc_after, "QC AFTER FILTERING")

    before_primary = qc_before["counts"]["primary_mapped"]
    after_primary  = qc_after["counts"]["primary_mapped"]
    ret = (after_primary / before_primary) if before_primary else float("nan")
    print(f"\nRetention (primary mapped): {ret:.2%}")

    # ------- Stage 06: summary file -------
    label_counts = read_qc_df["final_label"].value_counts().sort_index()
    total_reads = len(read_qc_df)
    lines = [
        f"Input BAM: {IN_BAM}",
        f"Output BAM: {OUT_BAM}",
        "",
        "Thresholds:",
        f"  MIN_MAPQ={MIN_MAPQ}",
        f"  MIN_ALN_FRAC={MIN_ALN_FRAC}",
        f"  MAX_CLIP_FRAC={MAX_CLIP_FRAC}",
        f"  MISMATCH_BP_MIN={MISMATCH_BP_MIN}",
        f"  MISMATCH_FRAC_MIN={MISMATCH_FRAC_MIN}",
        f"  HIGH_QUERY_COVERAGE_FRAC={HIGH_QUERY_COVERAGE_FRAC}",
        f"  MULTI_SEG_MIN={MULTI_SEG_MIN}",
        f"  REPEAT_LOCUS_RATIO={REPEAT_LOCUS_RATIO}",
        f"  DISTANT_GAP_BP={DISTANT_GAP_BP}",
        f"  MAX_QLEN_REFSPAN_DELTA={MAX_QLEN_REFSPAN_DELTA}",
        "",
        f"Total reads in read-level QC table: {total_reads:,}",
        "",
    ]
    for label, count in label_counts.items():
        frac = 100.0 * count / total_reads if total_reads else 0.0
        lines.append(f"{label:30s} : {count:>10,}  ({frac:6.3f}%)")
    lines += [
        "",
        f"Structural/QC-fail qnames excluded from HQ BAM: {len(exclude_qnames):,}",
        f"HQ BAM retention among primary mapped alignments: {ret:.4%}",
    ]
    with open(SUMMARY_TXT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote summary: {SUMMARY_TXT}")

    # ------- Post-filter quick report -------
    df = pd.read_csv(READ_QC_TSV, sep="\t")
    df_pass = df[df["final_label"] == "PASS_STRICT"].copy()
    print(f"\nTotal reads: {len(df):,}")
    print(f"PASS_STRICT reads: {len(df_pass):,} ({len(df_pass)/len(df)*100:.4f}%)")

    df_pass["mismatch_bp"]   = df_pass["qlen"] - df_pass["primary_ref_span"]
    df_pass["mismatch_frac"] = df_pass["mismatch_bp"] / df_pass["qlen"]

    def _summ(arr):
        a = np.asarray(arr)
        return {
            "min": float(np.min(a)),
            "p01": float(np.quantile(a, 0.01)),
            "p05": float(np.quantile(a, 0.05)),
            "p50": float(np.quantile(a, 0.50)),
            "p95": float(np.quantile(a, 0.95)),
            "p99": float(np.quantile(a, 0.99)),
            "max": float(np.max(a)),
            "mean": float(np.mean(a)),
        }

    print("\n=== Read length (qlen) ===");              print(_summ(df_pass["qlen"]))
    print("\n=== Primary mapped span ===");             print(_summ(df_pass["primary_ref_span"]))
    print("\n=== qlen - ref_span (bp) ===");            print(_summ(df_pass["mismatch_bp"]))
    print("\n=== (qlen - ref_span) / qlen ===");        print(_summ(df_pass["mismatch_frac"]))
    print("\n=== Primary alignment fraction ===");      print(_summ(df_pass["primary_aln_frac"]))
    print("\n=== Primary clip fraction ===");           print(_summ(df_pass["primary_clip_frac"]))

    for thr in (10, 50, 100):
        pct = (df_pass["mismatch_bp"].abs() <= thr).mean() * 100
        print(f"|qlen - ref_span| <= {thr:3d} bp : {pct:.2f}%")
