"""
Interactive HTML report of syn1 -> syn3A TPM/iPM changes for retained syn3A proteins.

Mimics the format of
  Syn1_Syn3A_Proteomics/proteome_common.py  (write_proteome_html)
(clickable composition bars on top -> filterable/sortable protein table below;
click Primary -> filter & hide Primary column; click Secondary or tertiary ->
filter & hide Primary+Secondary; sortable headers; CSV/Excel download), but uses
the joined RNA/protein change table as its data source. No PDF / no TSV output.

Inputs
  - Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv  (from 09; relTPM, relIPM,
    fold change, absolute change, PTR; mean-normalized per organism)
  - ../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv  (curated Primary /
    Secondary / Tertiary function + Essentiality, keyed by locus tag)

Filter
  Retained protein-coding genes only — rows whose locus_syn3a is not null
  (deleted genes are omitted).

Output
  Compare_RNA_Protein/TPM_iPM_changes_syn3A.html  (self-contained, no deps)
"""

import json
import re
from html import escape as esc

import pandas as pd

IN_TSV    = "Compare_RNA_Protein/syn1_vs_syn3a_RNA_protein.tsv"
FUNC_TSV = "../Syn1_Syn3A_Proteomics/syn3A_proteome.tsv"
OUT_HTML  = "Compare_RNA_Protein/TPM_iPM_changes_syn3A.html"

PRIM, SEC, TER = "Primary Function", "Secondary Function", "Tertiary Function"

# Compact column set shown in the table (downloads always export the full set).
TABLE_COLS = [
    "locus_syn3a", "gene_name", "gene_product",
    PRIM, SEC, TER, "essentiality",
    "relTPM_syn1", "relTPM_syn3a", "TPM_fold_change", "TPM_abs_change",
    "relIPM_syn1", "relIPM_syn3a", "iPM_fold_change", "iPM_abs_change",
    "PTR_fold_change",
]

# Primary-function palette (matches the report HTML).
PRIM_COLORS = {
    "Genetic Information Processing":       "#3b6db3",
    "Metabolism":                          "#3f9e5a",
    "Unclear":                             "#9aa0a6",
    "Cellular Processes":                  "#8e6bb1",
    "Environmental Information Processing": "#2aa6a0",
    "Exogenous":                           "#c0654e",
}


# ── Load + merge ─────────────────────────────────────────────────────────────
d = pd.read_csv(IN_TSV, sep="\t")
n_total = len(d)
d = d[d["locus_syn3a"].notna()].copy()      # omit deleted (no syn3A locus)

# Map to numeric locus for the function-annotation join.
d["locus_num"] = d["locus_syn3a"].str.extract(r"(\d+)$")[0].astype(int)
# snake_case table -> the Title-Case names used below
fa = pd.read_csv(FUNC_TSV, sep="\t").rename(columns={
    "locus_tag": "Locus Tag", "gene_name": "Gene Name",
    "gene_product": "Gene Product", "essentiality": "Essentiality",
    "primary_function": "Primary Function",
    "secondary_function": "Secondary Function",
    "tertiary_function": "Tertiary Function"})
fa["locus_num"] = fa["Locus Tag"].str.extract(r"(\d+)$")[0].astype(int)
# Pull curated annotation from the xlsx. Gene Name / Gene Product OVERRIDE the
# TSV's (09's TSV-derived names lag behind manual curation); essentiality is
# filled in where missing in the TSV.
fa = fa.rename(columns={k: v for k, v in {"Essentiality": "essentiality",
                                          "Gene Name": "gene_name",
                                          "Gene Product": "gene_product"}.items()
                        if k in fa.columns})
fa_keep = ["locus_num", PRIM, SEC, TER] + [c for c in ("gene_name", "gene_product", "essentiality")
                                            if c in fa.columns]
d = d.merge(fa[fa_keep], on="locus_num", how="left", suffixes=("", "_fa"))
# gene_name / gene_product: prefer the xlsx (curated); fall back to the TSV
for c in ("gene_name", "gene_product"):
    fc = c + "_fa"
    if fc in d.columns:
        d[c] = d[fc].fillna(d[c])
        d = d.drop(columns=[fc])
# essentiality: keep TSV value; fall back to xlsx if missing
if "essentiality_fa" in d.columns:
    d["essentiality"] = d["essentiality"].fillna(d["essentiality_fa"])
    d = d.drop(columns=["essentiality_fa"])

n_unannot = int(d[PRIM].isna().sum())
ORDER_COLS = [c for c in TABLE_COLS if c in d.columns]
print(f"Loaded {n_total} rows -> {len(d)} retained; columns in table: {len(ORDER_COLS)}; "
      f"unannotated retained: {n_unannot}")


# ── Records (NaN -> None) ────────────────────────────────────────────────────
records = json.loads(
    d[ORDER_COLS].astype(object).where(pd.notna(d[ORDER_COLS]), None)
                 .to_json(orient="records")
)
data_js = json.dumps(records).replace("</", "<\\/")
cols_js = json.dumps(ORDER_COLS).replace("</", "<\\/")


# ── Composition chart (Primary > Secondary > Tertiary; gene counts) ─────────
prim_order = d[PRIM].dropna().value_counts().index.tolist()
panels = []
for prim in prim_order:
    color = PRIM_COLORS.get(prim, "#777777")
    sub = d[d[PRIM] == prim]
    n_prim = len(sub)
    panel_body = []
    # Secondaries ordered by their own size desc; within each, tertiaries by size desc.
    secs = sub.dropna(subset=[SEC]).groupby(SEC).size().sort_values(ascending=False)
    xmax = sub.dropna(subset=[TER]).groupby(TER).size().max() if not sub.empty else 1
    for sec, n_sec in secs.items():
        panel_body.append(
            f'<div class="sec-h clk" data-p="{esc(prim, quote=True)}" '
            f'data-s="{esc(sec, quote=True)}">{esc(sec)} '
            f'<span class="sn">({n_sec})</span></div>'
        )
        tert = sub[sub[SEC] == sec].dropna(subset=[TER]).groupby(TER).size().sort_values(ascending=False)
        for ter, n_ter in tert.items():
            w = 100 * n_ter / xmax if xmax else 0
            panel_body.append(
                f'<div class="bar" data-p="{esc(prim, quote=True)}" '
                f'data-s="{esc(sec, quote=True)}" data-t="{esc(ter, quote=True)}">'
                f'<span class="lab">{esc(ter)}</span>'
                f'<span class="track"><span class="fill" style="width:{w:.1f}%"></span></span>'
                f'<span class="cnt">{n_ter}</span></div>'
            )
    panels.append(
        f'<section class="prim" style="--c:{color}">'
        f'<div class="prim-h clk" data-p="{esc(prim, quote=True)}">{esc(prim)} '
        f'<span class="pn">(n={n_prim})</span></div>'
        f'{"".join(panel_body)}</section>'
    )
chart_html = "".join(panels)


# ── HTML ─────────────────────────────────────────────────────────────────────
html_doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>syn3A retained proteins — TPM &amp; iPM changes</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
         margin: 0 auto; max-width: 1280px; padding: 20px; color:#222; }}
  h1 {{ font-size: 20px; }}
  .hint {{ color:#666; font-size:13px; margin-bottom:14px; }}
  .chart {{ display:grid; grid-template-columns: 1fr 1fr; gap: 10px 26px; }}
  @media (max-width: 800px) {{ .chart {{ grid-template-columns: 1fr; }} }}
  .prim-h {{ font-weight:700; color:var(--c); margin:10px 0 4px; font-size:15px; }}
  .pn {{ color:#888; font-weight:400; }}
  .sec-h {{ font-weight:700; color:var(--c); font-size:12px; margin:8px 0 2px 4px; opacity:.85; }}
  .sn {{ color:#999; font-weight:400; }}
  .clk {{ cursor:pointer; border-radius:4px; padding:1px 4px; }}
  .clk:hover {{ background:#f0f4fa; }}
  .bar {{ display:grid; grid-template-columns: 180px 1fr 30px; align-items:center;
          gap:6px; cursor:pointer; padding:1px 4px; border-radius:4px; }}
  .bar:hover {{ background:#f0f4fa; }}
  .bar.sel, .clk.sel {{ background:#fde9c8; }}
  .lab {{ font-size:11.5px; text-align:right; white-space:nowrap; overflow:hidden;
          text-overflow:ellipsis; }}
  .track {{ background:#eee; height:13px; border-radius:3px; }}
  .fill {{ display:block; height:100%; background:var(--c); border-radius:3px; }}
  .cnt {{ font-size:11px; color:#444; }}
  #panel {{ margin-top:26px; border-top:2px solid #ddd; padding-top:12px; }}
  #cap {{ font-weight:700; font-size:15px; }}
  #cap small {{ font-weight:400; color:#777; }}
  button {{ margin-left:10px; font-size:12px; padding:3px 8px; cursor:pointer; }}
  .tablebox {{ overflow-x:auto; margin-top:8px; max-height:75vh; }}
  table {{ border-collapse:collapse; font-size:12px; width:100%; }}
  th,td {{ border:1px solid #e2e2e2; padding:3px 7px; text-align:left;
           white-space:nowrap; vertical-align:top; }}
  th {{ position:sticky; top:0; background:#f5f5f5; cursor:pointer; user-select:none;
        white-space:nowrap; }}
  th.asc::after {{ content:" \\25B2"; font-size:9px; color:#888; }}
  th.desc::after {{ content:" \\25BC"; font-size:9px; color:#888; }}
  tr:nth-child(even) td {{ background:#fafafa; }}
  td.num {{ text-align:right; font-variant-numeric: tabular-nums; }}
  /* freeze the first three columns (Locus, Gene Name, Gene Product) when the
     table is scrolled horizontally. Widths are fixed so the cumulative left
     offsets line up. Overflow within the cell truncates with an ellipsis; the
     title attribute on each cell reveals the full text on hover. */
  th:nth-child(1), td:nth-child(1) {{ position:sticky; left:0;
       min-width:140px; max-width:140px; overflow:hidden; text-overflow:ellipsis; }}
  th:nth-child(2), td:nth-child(2) {{ position:sticky; left:140px;
       min-width:90px;  max-width:90px;  overflow:hidden; text-overflow:ellipsis; }}
  th:nth-child(3), td:nth-child(3) {{ position:sticky; left:230px;
       min-width:320px; max-width:320px; overflow:hidden; text-overflow:ellipsis;
       border-right:2px solid #ccc; }}
  td:nth-child(1), td:nth-child(2), td:nth-child(3) {{ background:white; z-index:1; }}
  tr:nth-child(even) td:nth-child(1), tr:nth-child(even) td:nth-child(2),
  tr:nth-child(even) td:nth-child(3) {{ background:#fafafa; }}
  th:nth-child(1), th:nth-child(2), th:nth-child(3) {{ background:#f5f5f5; z-index:3; }}
</style></head>
<body>
<h1>syn3A retained proteins — TPM &amp; iPM changes <small>(n = {len(d)} retained; deleted omitted)</small></h1>
<div class="hint">Click any Primary / Secondary / tertiary function to list its proteins below. Click a column header to sort. (rel = mean-normalized; FC = syn3A/syn1; Δ = syn3A − syn1; PTR_FC = iPM_FC / TPM_FC. Note: ribosomal-protein iPM is digestion-biased.)</div>
<div class="chart">{chart_html}</div>

<div id="panel">
  <span id="cap">Click a function to list its proteins</span>
  <button onclick="showAll()">Show all</button>
  <button onclick="downloadCSV()">Download CSV</button>
  <button onclick="downloadExcel()">Download Excel</button>
  <div class="tablebox"><table id="tbl"><thead><tr></tr></thead><tbody></tbody></table></div>
</div>

<script>
const DATA = {data_js};
const COLS = {cols_js};
const PRIM = "{esc(PRIM)}", SEC = "{esc(SEC)}", TER = "{esc(TER)}";
const NUMERIC = new Set(["relTPM_syn1","relTPM_syn3a","TPM_fold_change","TPM_abs_change",
                          "relIPM_syn1","relIPM_syn3a","iPM_fold_change","iPM_abs_change",
                          "PTR_fold_change"]);
const tb  = document.querySelector("#tbl tbody");
const thr = document.querySelector("#tbl thead tr");
const cap = document.getElementById("cap");
let currentRows = DATA.slice(), curCaption = "All retained proteins";
let visibleCols = COLS.slice();
let sortCol = null, sortDir = 1;            // sortCol is a column NAME

function escHTML(v) {{ return (v===null||v===undefined) ? "" :
    String(v).replace(/&/g,"&amp;").replace(/</g,"&lt;"); }}
function fmt(col, v) {{
  if (v===null || v===undefined) return "";
  if (NUMERIC.has(col) && typeof v === "number" && !Number.isInteger(v)) {{
    const a = Math.abs(v);
    return a >= 100 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toFixed(3);
  }}
  return String(v);
}}
function cellHTML(col, v) {{
  const txt = escHTML(fmt(col, v));
  const tt  = txt.replace(/"/g, '&quot;');     // hover tooltip = full value
  return NUMERIC.has(col) ? '<td class="num" title="'+tt+'">'+txt+'</td>'
                          : '<td title="'+tt+'">'+txt+'</td>';
}}
function draw() {{
  cap.innerHTML = curCaption + ' <small>(' + currentRows.length + ' proteins)</small>';
  thr.innerHTML = visibleCols.map(c => {{
    const cls = (c===sortCol) ? (sortDir>0 ? "asc" : "desc") : "";
    return '<th class="'+cls+'">'+escHTML(c)+'</th>';
  }}).join("");
  tb.innerHTML = currentRows.map(r =>
    "<tr>" + visibleCols.map(c => cellHTML(c, r[c])).join("") + "</tr>").join("");
}}
function applySort() {{
  if (sortCol !== null && visibleCols.indexOf(sortCol) !== -1) {{
    const c = sortCol;
    currentRows.sort((a, b) => {{
      let x = a[c], y = b[c];
      const xn = parseFloat(x), yn = parseFloat(y);
      const num = x!=null && y!=null && String(x).trim()!=="" && String(y).trim()!=="" &&
                  !isNaN(xn) && !isNaN(yn);
      const r = num ? (xn - yn)
                    : String(x==null?"":x).localeCompare(String(y==null?"":y));
      return r * sortDir;
    }});
  }}
  draw();
}}
function setView(rows, caption, hide) {{
  currentRows = rows.slice(); curCaption = caption;
  visibleCols = COLS.filter(c => hide.indexOf(c) === -1);
  applySort();
}}
function clearSel() {{ document.querySelectorAll(".sel").forEach(e => e.classList.remove("sel")); }}
function jump() {{ document.getElementById("panel").scrollIntoView({{behavior:"smooth", block:"start"}}); }}
function showAll() {{ clearSel(); sortCol = null; sortDir = 1; setView(DATA, "All retained proteins", []); }}

// click a Primary header -> filter, hide the Primary column
document.querySelectorAll(".prim-h.clk").forEach(el => el.addEventListener("click", () => {{
  clearSel(); el.classList.add("sel");
  const p = el.dataset.p;
  setView(DATA.filter(r => r[PRIM] === p), p, [PRIM]);
  jump();
}}));
// click a Secondary header -> filter, hide Primary + Secondary
document.querySelectorAll(".sec-h.clk").forEach(el => el.addEventListener("click", () => {{
  clearSel(); el.classList.add("sel");
  const p = el.dataset.p, s = el.dataset.s;
  setView(DATA.filter(r => r[PRIM] === p && r[SEC] === s), p + " &rsaquo; " + s, [PRIM, SEC]);
  jump();
}}));
// click a Tertiary bar -> filter, hide Primary + Secondary
document.querySelectorAll(".bar").forEach(bar => bar.addEventListener("click", () => {{
  clearSel(); bar.classList.add("sel");
  const p = bar.dataset.p, s = bar.dataset.s, t = bar.dataset.t;
  setView(DATA.filter(r => r[PRIM] === p && r[TER] === t && (!s || r[SEC] === s)),
          p + (s ? " &rsaquo; " + s : "") + " &rsaquo; " + t, [PRIM, SEC]);
  jump();
}}));
// sort on header click (column identified by name -> survives column hiding)
thr.addEventListener("click", e => {{
  const th = e.target.closest("th"); if (!th) return;
  const c = visibleCols[Array.prototype.indexOf.call(thr.children, th)];
  if (sortCol === c) {{ sortDir = -sortDir; }} else {{ sortCol = c; sortDir = 1; }}
  applySort();
}});
// downloads of the current filtered+sorted rows (always the full column set)
function dl(blob, fname) {{
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = fname; a.click();
  URL.revokeObjectURL(a.href);
}}
function csvCell(v) {{ v = (v==null) ? "" : String(v);
  return /[",\\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; }}
function downloadCSV() {{
  const lines = [COLS.map(csvCell).join(",")];
  currentRows.forEach(r => lines.push(COLS.map(c => csvCell(r[c])).join(",")));
  dl(new Blob(["\\ufeff" + lines.join("\\r\\n")], {{type:"text/csv;charset=utf-8"}}),
     "syn3A_TPM_iPM_changes.csv");
}}
function downloadExcel() {{
  let h = "<table><tr>" + COLS.map(c => "<th>" + escHTML(c) + "</th>").join("") + "</tr>";
  currentRows.forEach(r => h += "<tr>" + COLS.map(c => "<td>" + escHTML(r[c]) + "</td>").join("") + "</tr>");
  h += "</table>";
  dl(new Blob(['\\ufeff<html><head><meta charset="utf-8"></head><body>' + h + "</body></html>"],
     {{type:"application/vnd.ms-excel"}}), "syn3A_TPM_iPM_changes.xls");
}}
setView(DATA, "All retained proteins", []);
</script>
</body></html>"""

with open(OUT_HTML, "w") as fh:
    fh.write(html_doc)
print(f"Saved: {OUT_HTML}")
