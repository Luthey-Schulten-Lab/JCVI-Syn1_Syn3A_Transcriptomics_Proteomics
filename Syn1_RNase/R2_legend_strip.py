#!/usr/bin/env python
"""
R2 shared top legend = PANEL a -- born-at-size.

Molecule icons: scissors=endo (converted from the public-domain Openclipart SVG scissors_01.svg),
orange Pac-Man=exo 3'->5' (faces LEFT, chews toward 5'), blue Pac-Man=exo 5'->3' (faces RIGHT),
dark-red triangle=ribosome, black key=tmRNA/SmpB.

Emits:
  * five standalone icon PDFs (no labels, transparent) for Illustrator placement
  * R2a_legend_strip.pdf -- the 7-inch combined top strip = 4 erosion categories + 5 molecule
    icons, each icon with its enzyme names underneath (as in panel d).

Run in base env:
  /home/enguang/anaconda3/bin/python Syn1_RNase/R2_legend_strip.py
"""
import os
import re
import xml.etree.ElementTree as ET
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

mpl.rcParams.update({"font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "R2_panels")
SVG_SCISSORS = f"{OUT}/scissors_01.svg"
os.makedirs(OUT, exist_ok=True)

GREY, BLUE, RED, PURPLE = "#9e9e9e", "#3b6db3", "#c0392b", "#7a4fa3"      # erosion categories
SCIS, PAC_O, PAC_B, TRI, KEY = "#d94f2b", "#f5901f", "#1f77b4", "#a32020", "#1a1a1a"  # molecules


# ---------------------------------------------------------------- SVG path -> matplotlib
_NUM = re.compile(r"[MmLlHhVvCcSsQqTtAaZz]|-?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _parse_svg_path(d):
    toks = _NUM.findall(d); i = 0; cx = cy = sx = sy = 0.0; verts = []; codes = []; cmd = None

    def n():
        nonlocal i; v = float(toks[i]); i += 1; return v

    while i < len(toks):
        if re.match(r"[A-Za-z]", toks[i]):
            cmd = toks[i]; i += 1
        rel = cmd.islower(); C = cmd.upper(); bx, by = cx, cy
        if C == "M":
            x = n() + (bx if rel else 0); y = n() + (by if rel else 0); cx, cy = x, y; sx, sy = x, y
            verts.append((x, y)); codes.append(MplPath.MOVETO); cmd = "l" if rel else "L"
        elif C == "L":
            x = n() + (bx if rel else 0); y = n() + (by if rel else 0); cx, cy = x, y
            verts.append((x, y)); codes.append(MplPath.LINETO)
        elif C == "H":
            x = n() + (bx if rel else 0); cx = x; verts.append((x, cy)); codes.append(MplPath.LINETO)
        elif C == "V":
            y = n() + (by if rel else 0); cy = y; verts.append((cx, y)); codes.append(MplPath.LINETO)
        elif C == "C":
            pts = [(n() + (bx if rel else 0), n() + (by if rel else 0)) for _ in range(3)]
            verts += pts; codes += [MplPath.CURVE4] * 3; cx, cy = pts[2]
        elif C == "Z":
            verts.append((sx, sy)); codes.append(MplPath.CLOSEPOLY); cx, cy = sx, sy
        else:
            i += 1
    return verts, codes


def _load_svg(path):
    """Return list of (verts, codes) normalised to a centred [-1,1] box (aspect preserved)."""
    root = ET.parse(path).getroot()
    vb = [float(v) for v in root.get("viewBox").split()]; H = vb[3]
    ns = "{http://www.w3.org/2000/svg}"
    paths, allv = [], []
    for p in root.iter(ns + "path"):
        v, c = _parse_svg_path(p.get("d")); v = [(x, H - y) for x, y in v]   # SVG y is down
        paths.append((v, c)); allv += v
    xs = [p[0] for p in allv]; ys = [p[1] for p in allv]
    cx0 = (min(xs) + max(xs)) / 2; cy0 = (min(ys) + max(ys)) / 2
    half = max(max(xs) - min(xs), max(ys) - min(ys)) / 2
    return [([((x - cx0) / half, (y - cy0) / half) for x, y in v], c) for v, c in paths]


_SCISSORS = _load_svg(SVG_SCISSORS)


# ---------------------------------------------------------------- icon primitives
def icon_scissors(ax, cx, cy, s, color=SCIS):
    for v, c in _SCISSORS:
        ax.add_patch(PathPatch(MplPath([(cx + s * x, cy + s * y) for x, y in v], c),
                               facecolor=color, edgecolor="none", zorder=2))


def icon_pacman(ax, cx, cy, s, color, facing):
    t1, t2 = (32, 328) if facing == "right" else (212, 148)
    ax.add_patch(mpatches.Wedge((cx, cy), s, t1, t2, color=color, zorder=2))   # plain wedge, no eye


def icon_triangle(ax, cx, cy, s, color=TRI):
    ax.add_patch(plt.Polygon([(cx - 0.92 * s, cy - 0.78 * s), (cx + 0.92 * s, cy - 0.78 * s),
                              (cx, cy + 0.95 * s)], color=color, zorder=2))


def icon_key(ax, cx, cy, s, color=KEY):
    ax.add_patch(plt.Circle((cx - 0.60 * s, cy), 0.42 * s, color=color, zorder=2))
    ax.add_patch(plt.Circle((cx - 0.60 * s, cy), 0.17 * s, color="white", zorder=3))
    ax.add_patch(plt.Rectangle((cx - 0.28 * s, cy - 0.10 * s), 1.10 * s, 0.20 * s, color=color, zorder=2))
    for tx in (cx + 0.55 * s, cx + 0.76 * s):
        ax.add_patch(plt.Rectangle((tx, cy - 0.36 * s), 0.09 * s, 0.28 * s, color=color, zorder=2))


# exo 3'->5' faces LEFT (chews toward 5'); exo 5'->3' faces RIGHT
ICONS = {
    "endo_scissors":     lambda ax, cx, cy, s: icon_scissors(ax, cx, cy, s),
    "exo3to5_pacman":    lambda ax, cx, cy, s: icon_pacman(ax, cx, cy, s, PAC_O, "left"),
    "exo5to3_pacman":    lambda ax, cx, cy, s: icon_pacman(ax, cx, cy, s, PAC_B, "right"),
    "ribosome_triangle": lambda ax, cx, cy, s: icon_triangle(ax, cx, cy, s),
    "tmrna_key":         lambda ax, cx, cy, s: icon_key(ax, cx, cy, s),
}


def make_icons():
    for name, fn in ICONS.items():
        fig = plt.figure(figsize=(0.32, 0.32))
        ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.axis("off")
        fn(ax, 0, 0, 0.82)
        fig.savefig(f"{OUT}/R2_icon_{name}.pdf", transparent=True)
        plt.close(fig)
    print(f"[icons] {len(ICONS)} standalone icon PDFs -> {OUT}")


def make_strip():
    """7-inch combined top strip (PANEL a): 4 erosion categories + 5 molecule icons,
    each molecule icon with its enzyme names underneath."""
    H, FS, FE, SW, IC, GL = 0.70, 5.6, 4.4, 0.16, 0.18, 0.06
    fig = plt.figure(figsize=(7, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 7); ax.set_ylim(0, H); ax.axis("off")
    yc = H / 2
    fig.canvas.draw(); rend = fig.canvas.get_renderer()

    def tw(lab, fs):
        t = ax.text(0, -1, lab, fontsize=fs); w = t.get_window_extent(rend).width / fig.dpi
        t.remove(); return w

    cats = [("Unprocessed", GREY), ("5′ eroded", BLUE),
            ("3′ eroded", PAC_O), ("Both eroded", PURPLE)]   # 3' eroded matches the Exo 3'->5' icon
    mols = [("endo_scissors", "Endo", "RNase III, Y"),
            ("exo3to5_pacman", "Exo 3′→5′", "RNase R, YhaM"),
            ("exo5to3_pacman", "Exo 5′→3′", "RNase J1, J2"),
            ("ribosome_triangle", "Ribosome", ""),
            ("tmrna_key", "Rescue", "SmpB")]

    elems = [("cat", lab, col, None, SW + GL + tw(lab, FS)) for lab, col in cats]
    elems += [("div", None, None, None, 0.02)]
    elems += [("mol", role, name, enz, IC + GL + max(tw(role, FS), tw(enz, FE)))
              for name, role, enz in mols]

    content = sum(e[4] for e in elems)
    gap = max(0.08, min(0.30, (6.9 - content) / (len(elems) - 1)))
    x = (7 - (content + gap * (len(elems) - 1))) / 2
    for kind, a, b, enz, w in elems:
        if kind == "cat":
            ax.plot([x, x + SW], [yc, yc], color=b, lw=3.4, solid_capstyle="round")
            ax.text(x + SW + GL, yc, a, va="center", ha="left", fontsize=FS)
        elif kind == "div":
            ax.axvline(x + 0.01, color="#cccccc", lw=0.6)
        else:
            ICONS[b](ax, x + IC / 2, yc, 0.115)
            ylab = yc if not enz else yc + 0.095   # no enzyme sublabel (Ribosome) -> centre the role label
            ax.text(x + IC + GL, ylab, a, va="center", ha="left", fontsize=FS)
            # enzyme names NOT drawn -- space below the role label is reserved (item width still
            # accounts for the enzyme string) so they can be added in Illustrator at the same spot
        x += w + gap

    fig.savefig(f"{OUT}/R2a_legend_strip.pdf")
    plt.close(fig)
    print(f"[strip] panel-a 7in strip ({len(mols)} icons + {len(cats)} categories) "
          f"-> {OUT}/R2a_legend_strip.pdf")


if __name__ == "__main__":
    make_icons()
    make_strip()
