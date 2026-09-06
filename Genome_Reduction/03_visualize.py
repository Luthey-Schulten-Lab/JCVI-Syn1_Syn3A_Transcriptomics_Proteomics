#!/usr/bin/env python3
"""Interactive circular visualization of the syn1 -> syn3A genome reduction.

Outer ring  : syn1 genome (CP002027.1, 1,078,809 bp)
Inner ring  : syn3A genome (CP016816.2,   543,379 bp)

Driven by aln/analysis/genome_reduction_summary.xlsx (built by
02_analyze_genome_reduction.py). Each row is one change event, classified by
the "Change Case" column into one of:

    retained_ordered      (green   - both rings)
    retained_relocated    (orange  - both rings, highlighted, with connector)
    deleted               (red     - outer ring only)
    inserted              (blue    - inner ring only)

Hover tooltips include the syn1 / syn3A gene content of each segment.

Output: aln/analysis/genome_reduction_circle.html  (open in any browser).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

# ----------------------------------------------------------------- paths
SCRIPT = Path(__file__).resolve().parent
ALN = SCRIPT / "aln"
OUT = ALN / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

EXCEL = OUT / "genome_reduction_summary.xlsx"
COORDS = ALN / "raw/syn1_vs_syn3A.coords.tsv"   # only used for genome lengths
HTML = OUT / "genome_reduction_circle.html"

for p in (EXCEL, COORDS):
    if not p.exists():
        sys.exit(f"ERROR: required input not found: {p}\n"
                 f"Run 01_align.sh and then 02_analyze.py first.")

# ----------------------------------------------------------------- load
df = pd.read_excel(EXCEL)

# Genome lengths (LENR / LENQ live in coords.tsv columns 8 / 9)
COORDS_COLS = ["S1", "E1", "S2", "E2", "LEN1", "LEN2", "PCT_IDY",
               "LENR", "LENQ", "COVR", "COVQ", "TAGR", "TAGQ"]
_coords = pd.read_csv(COORDS, sep="\t")
_coords.columns = COORDS_COLS[: len(_coords.columns)]
SYN1_LEN = int(_coords["LENR"].iloc[0])
SYN3A_LEN = int(_coords["LENQ"].iloc[0])

# ----------------------------------------------------------------- helpers
def deg(pos: float, total: int) -> float:
    return 360.0 * pos / total

def arc_center_width(start: float, end: float, total: int) -> tuple[float, float]:
    return deg((start + end) / 2.0, total), deg(end - start, total)

def fmt_genes(cell: object, max_show: int = 8) -> str:
    """Pretty-print a gene-list cell from the Excel (comma separated)."""
    if cell is None or pd.isna(cell) or str(cell) in ("", "."):
        return "(intergenic)"
    if str(cell) == "NA":
        return "NA"
    items = [s.strip() for s in str(cell).split(",") if s.strip()]
    if not items:
        return "(intergenic)"
    if len(items) <= max_show:
        return ", ".join(items)
    return ", ".join(items[:max_show]) + f" ... (+{len(items) - max_show})"

# ----------------------------------------------------------------- ring geometry
OUTER_BASE, OUTER_HEIGHT = 0.78, 0.18    # syn1 ring radial extent
INNER_BASE, INNER_HEIGHT = 0.42, 0.18    # syn3A ring radial extent

# ----------------------------------------------------------------- traces
fig = go.Figure()
retained = df[df["Change Case"].isin(["retained_ordered", "retained_relocated"])]
deletions = df[df["Change Case"] == "deleted"]
insertions = df[df["Change Case"] == "inserted"]

# --- 1. syn1 retained blocks (outer ring) -------------------------------
for case_label, case_value in (("Retained, in order", "retained_ordered"),
                               ("Retained, relocated", "retained_relocated")):
    sub = retained[retained["Change Case"] == case_value]
    if sub.empty:
        continue
    thetas, widths, hovers = [], [], []
    for _, r in sub.iterrows():
        s1_0, e1_0 = int(r["S1"]), int(r["E1"])
        s2_0, e2_0 = int(r["S2"]), int(r["E2"])
        c, w = arc_center_width(s1_0, e1_0, SYN1_LEN)
        thetas.append(c); widths.append(w)
        hovers.append(
            f"<b>syn1 block {int(r['block_index_syn1'])}</b>"
            + ("  [<b>RELOCATED</b>]" if case_value == "retained_relocated" else "")
            + f"<br>syn1: {s1_0:,}-{e1_0:,} ({int(r['LEN1']):,} bp)"
            f"<br>syn3A: {s2_0:,}-{e2_0:,} ({int(r['LEN2']):,} bp)"
            f"<br>identity: {float(r['PCT_IDY']):.2f}%"
            f"<br>syn1 genes: {fmt_genes(r['Syn1_genes'])}"
            f"<br>syn3A genes: {fmt_genes(r['Syn3A_genes'])}"
        )
    n = len(thetas)
    if case_value == "retained_relocated":
        widths = [max(w, 2.0) for w in widths]
        fig.add_trace(go.Barpolar(
            r=[OUTER_HEIGHT + 0.10] * n,
            theta=thetas, width=widths,
            base=OUTER_BASE - 0.05,
            marker_color="#ff7f0e",
            marker_line_color="black", marker_line_width=1.0,
            opacity=1.0,
            name=f"syn1 - {case_label} (lap) [n={n}]",
            hovertext=hovers, hoverinfo="text",
        ))
    else:
        fig.add_trace(go.Barpolar(
            r=[OUTER_HEIGHT] * n,
            theta=thetas, width=widths,
            base=OUTER_BASE,
            marker_color="#2ca02c",
            marker_line_color="white", marker_line_width=0.3,
            opacity=0.95,
            name=f"syn1 - {case_label} [n={n}]",
            hovertext=hovers, hoverinfo="text",
        ))

# --- 2. syn1 deletions (outer ring) -------------------------------------
thetas, widths, hovers = [], [], []
for _, r in deletions.iterrows():
    s1_0, e1_0 = int(r["S1"]), int(r["E1"])
    c, w = arc_center_width(s1_0, e1_0, SYN1_LEN)
    thetas.append(c); widths.append(w)
    hovers.append(
        f"<b>Deletion</b>"
        f"<br>syn1: {s1_0:,}-{e1_0:,} ({int(r['LEN1']):,} bp)"
        f"<br>syn1 genes: {fmt_genes(r['Syn1_genes'])}"
    )
n_del = len(thetas)
fig.add_trace(go.Barpolar(
    r=[OUTER_HEIGHT] * n_del,
    theta=thetas, width=widths,
    base=OUTER_BASE,
    marker_color="#d62728",
    marker_line_color="white", marker_line_width=0.3,
    opacity=0.85,
    name=f"syn1 - Deleted [n={n_del}]",
    hovertext=hovers, hoverinfo="text",
))

# --- 3. syn3A retained blocks (inner ring) ------------------------------
for case_label, case_value in (("Retained, in order", "retained_ordered"),
                               ("Retained, relocated", "retained_relocated")):
    sub = retained[retained["Change Case"] == case_value]
    if sub.empty:
        continue
    thetas, widths, hovers = [], [], []
    for _, r in sub.iterrows():
        s1_0, e1_0 = int(r["S1"]), int(r["E1"])
        s2_0, e2_0 = int(r["S2"]), int(r["E2"])
        c, w = arc_center_width(s2_0, e2_0, SYN3A_LEN)
        thetas.append(c); widths.append(w)
        hovers.append(
            f"<b>syn3A block {int(r['block_index_syn1'])}</b>"
            + ("  [<b>RELOCATED</b>]" if case_value == "retained_relocated" else "")
            + f"<br>syn3A: {s2_0:,}-{e2_0:,} ({int(r['LEN2']):,} bp)"
            f"<br>syn1 origin: {s1_0:,}-{e1_0:,}"
            f"<br>identity: {float(r['PCT_IDY']):.2f}%"
            f"<br>syn1 genes: {fmt_genes(r['Syn1_genes'])}"
            f"<br>syn3A genes: {fmt_genes(r['Syn3A_genes'])}"
        )
    n = len(thetas)
    if case_value == "retained_relocated":
        widths = [max(w, 2.0) for w in widths]
        fig.add_trace(go.Barpolar(
            r=[INNER_HEIGHT + 0.10] * n,
            theta=thetas, width=widths,
            base=INNER_BASE - 0.05,
            marker_color="#ff7f0e",
            marker_line_color="black", marker_line_width=1.0,
            opacity=1.0,
            name=f"syn3A - {case_label} (lap)",
            showlegend=False,
            hovertext=hovers, hoverinfo="text",
        ))
    else:
        fig.add_trace(go.Barpolar(
            r=[INNER_HEIGHT] * n,
            theta=thetas, width=widths,
            base=INNER_BASE,
            marker_color="#2ca02c",
            marker_line_color="white", marker_line_width=0.3,
            opacity=0.95,
            name=f"syn3A - {case_label}",
            showlegend=False,
            hovertext=hovers, hoverinfo="text",
        ))

# --- 4. syn3A insertions (inner ring) -----------------------------------
thetas, widths, hovers = [], [], []
for _, r in insertions.iterrows():
    s2_0, e2_0 = int(r["S2"]), int(r["E2"])
    c, w = arc_center_width(s2_0, e2_0, SYN3A_LEN)
    thetas.append(c); widths.append(w)
    note = ""
    if int(r["LEN2"]) >= 1000:
        note = "<br><i>contains JCVISYN3A_0931 (met14p) - the only novel gene</i>"
    elif int(r["LEN2"]) <= 100:
        note = "<br><i>(intergenic junction/scar)</i>"
    hovers.append(
        f"<b>Insertion</b>"
        f"<br>syn3A: {s2_0:,}-{e2_0:,} ({int(r['LEN2']):,} bp)"
        f"<br>syn3A genes: {fmt_genes(r['Syn3A_genes'])}"
        f"{note}"
    )
n_ins = len(thetas)
fig.add_trace(go.Barpolar(
    r=[INNER_HEIGHT] * n_ins,
    theta=thetas, width=widths,
    base=INNER_BASE,
    marker_color="#1f77b4",
    marker_line_color="white", marker_line_width=0.3,
    opacity=0.95,
    name=f"syn3A - Inserted [n={n_ins}]",
    hovertext=hovers, hoverinfo="text",
))

# --- 5. coordinate ticks ------------------------------------------------
tick_pos = list(range(0, SYN1_LEN, 100_000))
tick_theta = [deg(p, SYN1_LEN) for p in tick_pos]
tick_text = [f"{p // 1000}k" for p in tick_pos]

# Inner ring synthetic ticks every 50 kb
INNER_TICK_R_TEXT = INNER_BASE - 0.10
INNER_TICK_R_TIP = INNER_BASE - 0.005
INNER_TICK_R_BASE = INNER_BASE - 0.03
syn3a_tick_pos = list(range(0, SYN3A_LEN, 50_000))
syn3a_tick_theta = [deg(p, SYN3A_LEN) for p in syn3a_tick_pos]
syn3a_tick_text = [f"{p // 1000}k" for p in syn3a_tick_pos]

for th in syn3a_tick_theta:
    fig.add_trace(go.Scatterpolar(
        r=[INNER_TICK_R_BASE, INNER_TICK_R_TIP], theta=[th, th],
        mode="lines",
        line=dict(color="#444", width=1),
        hoverinfo="skip", showlegend=False,
    ))
fig.add_trace(go.Scatterpolar(
    r=[INNER_TICK_R_TEXT] * len(syn3a_tick_theta),
    theta=syn3a_tick_theta,
    mode="text",
    text=syn3a_tick_text,
    textfont=dict(size=9, color="#1f77b4"),
    hoverinfo="skip", showlegend=False,
))

# ----------------------------------------------------------------- layout
fig.update_layout(
    title=dict(
        text=f"<b>JCVI-Syn1.0 - JCVI-Syn3A genome reduction</b>"
             f"<br><sup>Outer: syn1 ({SYN1_LEN:,} bp) - Inner: syn3A ({SYN3A_LEN:,} bp) - "
             f"hover any segment for gene content</sup>"
             f"<br><sup style='color:#666'>Drag to box-zoom into a wedge - "
             f"scroll wheel zooms radially - double-click to reset</sup>",
        x=0.5, xanchor="center",
    ),
    polar=dict(
        bgcolor="white",
        radialaxis=dict(range=[0, 1.05], showticklabels=False, ticks=""),
        angularaxis=dict(
            direction="clockwise",
            rotation=90,
            tickmode="array",
            tickvals=tick_theta,
            ticktext=tick_text,
            tickfont=dict(size=10),
        ),
    ),
    legend=dict(
        orientation="v",
        x=1.05, y=0.5, yanchor="middle",
        bgcolor="rgba(255,255,255,0.8)",
    ),
    width=1100, height=950,
    margin=dict(l=40, r=240, t=110, b=40),
    dragmode="zoom",
)

plotly_config = {
    "scrollZoom": True,
    "displayModeBar": True,
    "modeBarButtonsToAdd": ["zoom2d", "pan2d", "resetScale2d"],
    "doubleClick": "reset",
    "displaylogo": False,
}

fig.write_html(HTML, include_plotlyjs="cdn", config=plotly_config)
print(f"Wrote: {HTML}")
print(f"Open in a browser:  file://{HTML}")
