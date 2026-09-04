"""Shared I/O and progress helpers for PacBio FLNC preprocessing.

Provides:
  - pigz_reader / pigz_writer: context managers wrapping `pigz` subprocesses
  - Progress: generic rate-reporting helper with arbitrary metrics
"""

from __future__ import annotations

import io
import sys
import time
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Dict


@contextmanager
def pigz_reader(path: str, threads: int = 16, bufsize: int = 1024 * 1024):
    """Decompress `path` via `pigz -dc`; yields a BufferedReader of raw bytes."""
    p = subprocess.Popen(
        ["pigz", "-dc", "-p", str(threads), path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0,
    )
    if p.stdout is None or p.stderr is None:
        if p.stdout: p.stdout.close()
        if p.stderr: p.stderr.close()
        p.kill()
        raise RuntimeError("Failed to start pigz_reader subprocess")

    try:
        yield io.BufferedReader(p.stdout, buffer_size=bufsize)
    finally:
        try:
            p.stdout.close()
        except Exception:
            pass
        err = b""
        try:
            rc = p.wait()
            err = p.stderr.read() or b""
        finally:
            try:
                p.stderr.close()
            except Exception:
                pass
        if rc != 0:
            raise RuntimeError(
                f"pigz_reader failed (code={rc}):\n{err.decode('utf-8', errors='replace')}"
            )


@contextmanager
def pigz_writer(
    path: str,
    threads: int = 16,
    bufsize: int = 1024 * 1024,
    compresslevel: int = 6,
):
    """Compress uncompressed bytes to `path` via `pigz`; yields a BufferedWriter."""
    out_f = open(path, "wb")
    p = subprocess.Popen(
        ["pigz", "-p", str(threads), f"-{int(compresslevel)}"],
        stdin=subprocess.PIPE,
        stdout=out_f,
        stderr=subprocess.PIPE,
        bufsize=0,
    )
    if p.stdin is None or p.stderr is None:
        out_f.close()
        if p.stdin: p.stdin.close()
        if p.stderr: p.stderr.close()
        p.kill()
        raise RuntimeError("Failed to start pigz_writer subprocess")

    writer = io.BufferedWriter(p.stdin, buffer_size=bufsize)
    try:
        yield writer
    finally:
        try:
            writer.flush()
        except BrokenPipeError:
            pass
        try:
            writer.close()
        except Exception:
            pass
        try:
            out_f.flush()
        except Exception:
            pass
        out_f.close()
        err = b""
        try:
            rc = p.wait()
            err = p.stderr.read() or b""
        finally:
            try:
                p.stderr.close()
            except Exception:
                pass
        if rc != 0:
            raise RuntimeError(
                f"pigz_writer failed (code={rc}):\n{err.decode('utf-8', errors='replace')}"
            )


_RC_TABLE = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

def revcomp(seq: bytes) -> bytes:
    return seq.translate(_RC_TABLE)[::-1]


@dataclass
class Progress:
    """Generic progress reporter. Pass arbitrary counters as kwargs to maybe_report/final."""
    enabled: bool = True
    every_seconds: float = 5.0
    label: str = "stage"
    stream = sys.stderr

    _t0: float = 0.0
    _t_last: float = 0.0
    _n_last: int = 0

    def start(self, n0: int = 0) -> None:
        now = time.time()
        self._t0 = now
        self._t_last = now
        self._n_last = n0

    def _fmt_metrics(self, metrics: Dict[str, int]) -> str:
        return "  ".join(f"{k}={v:,}" for k, v in metrics.items())

    def maybe_report(self, n: int, **metrics) -> None:
        if not self.enabled:
            return
        now = time.time()
        if (now - self._t_last) < self.every_seconds:
            return
        dt = now - self._t_last
        rate = ((n - self._n_last) / dt) if dt > 0 else 0.0
        total_dt = now - self._t0
        avg = (n / total_dt) if total_dt > 0 else 0.0
        print(
            f"[{self.label}] reads={n:,}  {self._fmt_metrics(metrics)}  "
            f"rate={rate:,.0f} r/s  avg={avg:,.0f} r/s",
            file=self.stream, flush=True,
        )
        self._t_last = now
        self._n_last = n

    def final(self, n: int, **metrics) -> None:
        if not self.enabled:
            return
        total_dt = time.time() - self._t0
        avg = (n / total_dt) if total_dt > 0 else 0.0
        print(
            f"[{self.label}] DONE  reads={n:,}  {self._fmt_metrics(metrics)}  "
            f"elapsed={total_dt:,.1f}s  avg={avg:,.0f} r/s",
            file=self.stream, flush=True,
        )
