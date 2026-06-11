"""Self-contained interactive HTML coverage viewer."""

from __future__ import annotations

from pathlib import Path
from typing import List, Set

import matplotlib
import matplotlib.colors as mcolors
from jinja2 import Template

from ..config import Config
from ..model import MSA, CoverageMatrix, Query
from .layout import row_label, select_rows

_CELL_W = 13  # px; keeps ruler ticks aligned to columns


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


def _cells_for_aggregate(matrix: CoverageMatrix, step: int) -> str:
    parts: List[str] = []
    agg = matrix.column_aggregate
    for col in range(matrix.columns):
        d = int(agg[col])
        pos = col + 1
        if d:
            title = f"pos {pos} ({matrix.query_seq[col]}) · {d} structures · {_fmt_ids(matrix.column_pdb_ids[col])}"
            parts.append(f'<span class="c d{d}" title="{title}">&nbsp;</span>')
        else:
            parts.append('<span class="c z">&nbsp;</span>')
    return "".join(parts)


def _cells_for_row(matrix: CoverageMatrix, msa: MSA, i: int) -> str:
    parts: List[str] = []
    seq = msa.rows[i].match_seq
    for col in range(matrix.columns):
        ch = seq[col] if col < len(seq) else "-"
        if ch == "-":
            parts.append('<span class="c gap">-</span>')
            continue
        d = int(matrix.depth[i, col])
        if d:
            ids = matrix.pdb_ids[i][col]
            title = f"{ch}{col + 1} · depth {d} · {_fmt_ids(ids)}"
            parts.append(f'<span class="c d{d}" title="{title}">{ch}</span>')
        else:
            parts.append(f'<span class="c ns">{ch}</span>')
    return "".join(parts)


def _ruler(n_cols: int, step: int) -> str:
    parts: List[str] = []
    for col in range(n_cols):
        if (col + 1) % step == 0 or col == 0:
            parts.append(f'<span class="c tick">{col + 1}</span>')
        else:
            parts.append('<span class="c tick">&nbsp;</span>')
    return "".join(parts)


_TEMPLATE = Template(
    """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>{{ query.name }} — structural coverage</title>
<style>
 body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;color:#1b1f23;background:#fff;}
 header{padding:14px 18px;border-bottom:1px solid #e1e4e8;position:sticky;top:0;background:#fff;z-index:5;}
 h1{font-size:16px;margin:0 0 4px;} .sub{font-size:12px;color:#586069;}
 .legend{margin-top:8px;font-size:11px;display:flex;align-items:center;gap:6px;flex-wrap:wrap;}
 .legend .sw{display:inline-block;width:14px;height:14px;border:1px solid #ccc;vertical-align:middle;}
 .scroll{overflow:auto;padding:10px 18px 40px;}
 .grid{font-family:'SF Mono',Menlo,Consolas,monospace;font-size:11px;white-space:nowrap;line-height:1.0;}
 .line{display:flex;align-items:stretch;}
 .label{position:sticky;left:0;background:#fff;z-index:3;min-width:230px;max-width:230px;
        padding-right:8px;font-size:10px;color:#24292e;overflow:hidden;text-overflow:ellipsis;
        border-right:1px solid #eee;display:flex;align-items:center;}
 .cells{display:inline-block;}
 .c{display:inline-block;width:{{cw}}px;text-align:center;overflow:visible;}
 .tick{color:#959da5;font-size:9px;}
 .gap{color:#d1d5da;}
 .ns{background:#eef0f2;color:#9aa0a6;}   /* residue present, no structure */
 .z{background:#eef0f2;}
 .queryline{position:sticky;top:0;}
 .aggline .label{font-weight:600;}
 .qline .label{font-weight:600;color:#0b67d0;}
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
</header>
<div class="scroll"><div class="grid">
 <div class="line"><div class="label">position</div><div class="cells">{{ ruler }}</div></div>
 <div class="line aggline"><div class="label">all orthologs (union)</div><div class="cells">{{ aggregate }}</div></div>
 {% for r in rows %}
 <div class="line {{ 'qline' if r.is_query else '' }}"><div class="label" title="{{ r.label }}">{{ r.label }}</div><div class="cells">{{ r.cells }}</div></div>
 {% endfor %}
</div></div>
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
    step = max(10, (matrix.columns // 25 // 10 + 1) * 10)

    rows = []
    for i in selected:
        rows.append(
            {
                "label": row_label(matrix, i),
                "cells": _cells_for_row(matrix, msa, i),
                "is_query": i == matrix.query_index,
            }
        )

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
        cw=_CELL_W,
        depth_css=depth_css,
        legend=sorted(colors.items()),
        ruler=_ruler(matrix.columns, step),
        aggregate=_cells_for_aggregate(matrix, step),
        rows=rows,
        stats=stats,
    )
    path = out_dir / "alignment.html"
    path.write_text(html)
    return [path]
