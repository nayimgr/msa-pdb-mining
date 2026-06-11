"""Self-contained interactive HTML coverage viewer.

The alignment is rendered as a *wrapped* block view (à la Clustal/Jalview): the
columns are broken into blocks that fit the current viewport width, so a long
protein flows down the page as stacked blocks instead of one very wide row with
a horizontal scrollbar. The wrapping is done in the browser from an embedded
JSON payload so it re-flows live on window resize and on zoom, which a static
server-side layout cannot do (the render machine has no idea how wide the
reader's screen is).

Cells are coloured by structural-coverage depth (viridis, matching the static
heatmap); hovering a cell shows the residue position, depth, and the PDB ids
modelling it. The depth→colour classes (``.d1``…``.dN``) are emitted into the
stylesheet so the per-cell DOM stays small.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Set

import matplotlib
import matplotlib.colors as mcolors
from jinja2 import Template

from ..config import Config
from ..model import MSA, CoverageMatrix, Query
from .layout import row_label, select_rows


def _depth_colors(vmax: int) -> dict:
    cmap = matplotlib.colormaps["viridis"]
    colors = {}
    for d in range(1, vmax + 1):
        frac = (d - 1) / (vmax - 1) if vmax > 1 else 1.0
        colors[d] = mcolors.to_hex(cmap(frac))
    return colors


def _fmt_ids(ids: Set[str], limit: int = 15) -> str:
    s = sorted(ids)
    if len(s) <= limit:
        return ", ".join(s)
    return ", ".join(s[:limit]) + f" (+{len(s) - limit} more)"


def _build_payload(matrix: CoverageMatrix, msa: MSA, selected: List[int]) -> dict:
    """Compact, browser-renderable description of the alignment.

    PDB-id tooltip strings are interned into a pool and referenced by index so
    the repeated ids of a contiguous covered range cost one integer per cell
    instead of a full string.
    """
    id_pool: List[str] = []
    id_index: Dict[str, int] = {}

    def intern(ids: Set[str]) -> int:
        if not ids:
            return -1
        s = _fmt_ids(ids)
        idx = id_index.get(s)
        if idx is None:
            idx = len(id_pool)
            id_index[s] = idx
            id_pool.append(s)
        return idx

    cols = matrix.columns
    agg = matrix.column_aggregate
    aggregate = {
        "depths": [int(agg[c]) for c in range(cols)],
        "ids": [intern(matrix.column_pdb_ids[c]) for c in range(cols)],
    }

    rows = []
    for i in selected:
        seq = msa.rows[i].match_seq
        chars: List[str] = []
        depths: List[int] = []
        ids: List[int] = []
        for c in range(cols):
            ch = seq[c] if c < len(seq) else "-"
            chars.append(ch)
            if ch == "-":
                depths.append(0)
                ids.append(-1)
                continue
            d = int(matrix.depth[i, c])
            depths.append(d)
            ids.append(intern(matrix.pdb_ids[i][c]) if d else -1)
        rows.append(
            {
                "label": row_label(matrix, i),
                "isQuery": i == matrix.query_index,
                "chars": "".join(chars),
                "depths": depths,
                "ids": ids,
            }
        )

    return {
        "columns": cols,
        "querySeq": matrix.query_seq,
        "idPool": id_pool,
        "aggregate": aggregate,
        "rows": rows,
    }


_TEMPLATE = Template(
    """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ query.name }} — structural coverage</title>
<style>
 :root{--cellw:13px;--fs:11px;--labelw:190px;}
 body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;color:#1b1f23;background:#fff;}
 header{padding:14px 18px;border-bottom:1px solid #e1e4e8;position:sticky;top:0;background:#fff;z-index:5;}
 h1{font-size:16px;margin:0 0 4px;} .sub{font-size:12px;color:#586069;}
 .legend{margin-top:8px;font-size:11px;display:flex;align-items:center;gap:6px;flex-wrap:wrap;}
 .legend .sw{display:inline-block;width:14px;height:14px;border:1px solid #ccc;vertical-align:middle;}
 .toolbar{margin-top:8px;font-size:12px;display:flex;align-items:center;gap:8px;color:#586069;}
 .toolbar button{font:inherit;cursor:pointer;border:1px solid #d1d5da;background:#f6f8fa;
        border-radius:5px;padding:2px 9px;line-height:1.4;}
 .toolbar button:hover{background:#eef0f2;}
 main{padding:12px 18px 48px;}
 #viewer{font-family:'SF Mono',Menlo,Consolas,monospace;font-size:var(--fs);line-height:1.05;}
 .block{margin:0 0 20px;}
 .line{display:flex;align-items:stretch;white-space:nowrap;}
 .label{flex:0 0 var(--labelw);max-width:var(--labelw);padding-right:8px;font-size:10px;
        color:#24292e;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;
        border-right:1px solid #eee;display:flex;align-items:center;}
 .cells{white-space:nowrap;}
 .c{display:inline-block;width:var(--cellw);text-align:center;overflow:visible;}
 .tick{color:#959da5;font-size:9px;}
 .gap{color:#d1d5da;}
 .ns{background:#eef0f2;color:#9aa0a6;}   /* residue present, no structure */
 .z{background:#eef0f2;}
 .rulerline{color:#959da5;}
 .aggline .label{font-weight:600;}
 .qline .label{font-weight:700;color:#0b67d0;}
 .qline .cells{outline:1px solid rgba(11,103,208,.18);}
 @media (max-width:560px){ :root{--labelw:120px;} }
{{ depth_css }}
</style></head>
<body>
<header>
 <h1>{{ query.name }} — structural data availability across orthologs</h1>
 <div class="sub">
   {{ query.accession or '' }}{% if query.organism %} · {{ query.organism }}{% endif %}
   · {{ matrix.columns }} residues · {{ stats.covered }} covered ({{ stats.frac }})
   · {{ stats.structured_rows }} orthologs with structures · {{ stats.pdb_count }} distinct PDB entries
 </div>
 <div class="legend"><span>depth:</span>
   <span class="sw z"></span><span>0</span>
   {% for d, col in legend %}<span class="sw" style="background:{{col}}"></span><span>{{d}}</span>{% endfor %}
   <span class="sw ns" style="margin-left:10px;"></span><span>residue, no structure</span>
   <span class="sw gap"></span><span>gap</span>
 </div>
 <div class="toolbar">
   <span>zoom:</span>
   <button id="zout" type="button" title="smaller cells">&minus;</button>
   <button id="zin" type="button" title="larger cells">+</button>
   <span class="hint">alignment wraps to the window width — resize to re-flow</span>
 </div>
</header>
<main><div id="viewer"></div></main>
<script id="coverage-data" type="application/json">{{ data_json }}</script>
<script>
(function(){
 var DATA = JSON.parse(document.getElementById('coverage-data').textContent);
 var viewer = document.getElementById('viewer');
 var root = document.documentElement;
 var TICK_STEP = 10;
 var cellW = 13;

 function attr(s){
   return String(s).replace(/&/g,'&amp;').replace(/"/g,'&quot;')
                   .replace(/</g,'&lt;').replace(/>/g,'&gt;');
 }
 function cssPx(name, fallback){
   var v = parseInt(getComputedStyle(root).getPropertyValue(name), 10);
   return isNaN(v) ? fallback : v;
 }

 function ruler(start, end){
   var h = '';
   for(var c=start;c<end;c++){
     var pos = c+1;
     if(c===start || pos % TICK_STEP === 0){
       h += '<span class="c tick">'+pos+'</span>';
     } else {
       h += '<span class="c tick">&nbsp;</span>';
     }
   }
   return h;
 }
 function aggCells(start, end){
   var d = DATA.aggregate.depths, ids = DATA.aggregate.ids, h = '';
   for(var c=start;c<end;c++){
     var depth = d[c], pos = c+1;
     if(depth>0){
       var idx = ids[c], idStr = idx>=0 ? DATA.idPool[idx] : '';
       var t = 'pos '+pos+' ('+DATA.querySeq[c]+') · '+depth+' structures'+(idStr?' · '+idStr:'');
       h += '<span class="c d'+depth+'" title="'+attr(t)+'">&nbsp;</span>';
     } else {
       h += '<span class="c z">&nbsp;</span>';
     }
   }
   return h;
 }
 function rowCells(row, start, end){
   var chars = row.chars, depths = row.depths, ids = row.ids, h = '';
   for(var c=start;c<end;c++){
     var ch = chars.charAt(c), pos = c+1;
     if(ch==='-'){ h += '<span class="c gap">-</span>'; continue; }
     var depth = depths[c];
     if(depth>0){
       var idx = ids[c], idStr = idx>=0 ? DATA.idPool[idx] : '';
       var t = ch+pos+' · depth '+depth+(idStr?' · '+idStr:'');
       h += '<span class="c d'+depth+'" title="'+attr(t)+'">'+ch+'</span>';
     } else {
       h += '<span class="c ns">'+ch+'</span>';
     }
   }
   return h;
 }

 function render(){
   var labelW = cssPx('--labelw', 190);
   var avail = viewer.clientWidth;
   var perBlock = Math.max(1, Math.floor((avail - labelW) / cellW));
   var html = '';
   for(var start=0; start<DATA.columns; start+=perBlock){
     var end = Math.min(start+perBlock, DATA.columns);
     html += '<div class="block">';
     html += '<div class="line rulerline"><div class="label"></div><div class="cells">'+ruler(start,end)+'</div></div>';
     html += '<div class="line aggline"><div class="label">all orthologs (union)</div><div class="cells">'+aggCells(start,end)+'</div></div>';
     for(var r=0;r<DATA.rows.length;r++){
       var row = DATA.rows[r];
       html += '<div class="line'+(row.isQuery?' qline':'')+'">'
            +  '<div class="label" title="'+attr(row.label)+'">'+attr(row.label)+'</div>'
            +  '<div class="cells">'+rowCells(row,start,end)+'</div></div>';
     }
     html += '</div>';
   }
   viewer.innerHTML = html;
 }

 function setCell(w){
   cellW = Math.max(7, Math.min(24, w));
   root.style.setProperty('--cellw', cellW+'px');
   root.style.setProperty('--fs', Math.max(8, Math.round(cellW*0.82))+'px');
   render();
 }

 document.getElementById('zin').addEventListener('click', function(){ setCell(cellW+2); });
 document.getElementById('zout').addEventListener('click', function(){ setCell(cellW-2); });
 var timer;
 window.addEventListener('resize', function(){
   clearTimeout(timer); timer = setTimeout(render, 120);
 });
 setCell(cellW);
})();
</script>
</body></html>
"""
)


def write_html(
    out_dir: Path, matrix: CoverageMatrix, msa: MSA, query: Query, config: Config
) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = select_rows(matrix, config.max_rows_render)
    vmax = max(int(matrix.depth.max()), int(matrix.column_aggregate.max()), 1)
    colors = _depth_colors(vmax)
    depth_css = "\n".join(f" .d{d}{{background:{c};color:#fff;}}" for d, c in colors.items())

    payload = _build_payload(matrix, msa, selected)
    # Embedded as a JSON <script>; escape '<' so a stray '</script>' in the data
    # (or the '<' in any tooltip) can never terminate the block early.
    data_json = json.dumps(payload, separators=(",", ":")).replace("<", "\\u003c")

    agg = matrix.column_aggregate
    covered = int((agg > 0).sum())
    all_ids: Set[str] = set()
    for ids in matrix.column_pdb_ids:
        all_ids |= ids
    stats = {
        "covered": covered,
        "frac": f"{(covered / matrix.columns if matrix.columns else 0):.0%}",
        "structured_rows": sum(1 for s in matrix.row_status if s == "structured"),
        "pdb_count": len(all_ids),
    }

    html = _TEMPLATE.render(
        query=query,
        matrix=matrix,
        depth_css=depth_css,
        legend=sorted(colors.items()),
        data_json=data_json,
        stats=stats,
    )
    path = out_dir / "alignment.html"
    path.write_text(html)
    return [path]
