"""Static publication-style coverage heatmap (PNG + SVG) via matplotlib."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch, Polygon  # noqa: E402

from ..config import Config  # noqa: E402
from ..model import CoverageMatrix, Query  # noqa: E402
from ..secondary_structure import HELIX, STRAND  # noqa: E402
from .layout import row_label, select_rows, ss_runs  # noqa: E402
from .palette import SS_HELIX_COLOR, SS_STRAND_COLOR, depth_cmap  # noqa: E402

_BG = "#eef0f2"  # depth-zero cells / disordered


def write_image(
    out_dir: Path,
    matrix: CoverageMatrix,
    query: Query,
    config: Config,
    formats: Sequence[str] = ("png", "svg"),
) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = select_rows(matrix, config.max_rows_render)
    sub = matrix.depth[rows, :]
    aggregate = matrix.column_aggregate.reshape(1, -1)
    vmax = max(int(sub.max()), int(aggregate.max()), 1)

    n_cols = matrix.columns
    n_rows = len(rows)
    has_ss = matrix.has_secondary_structure
    width = min(40.0, max(8.0, n_cols * 0.05))
    height = min(30.0, max(3.5, (n_rows + 3 + (1 if has_ss else 0)) * 0.26))

    fig = plt.figure(figsize=(width, height))
    # Optional secondary-structure cartoon row on top, then the aggregate track,
    # then the per-row heatmap; a shared colorbar spans the right column.
    height_ratios = ([0.8] if has_ss else []) + [1, max(n_rows, 1)]
    gs = GridSpec(
        len(height_ratios), 2, height_ratios=height_ratios, width_ratios=[40, 1],
        hspace=0.06, wspace=0.02,
    )
    cmap = depth_cmap().copy()
    cmap.set_bad(_BG)

    heat_row = len(height_ratios) - 1
    agg_row = heat_row - 1

    # --- aggregate track ---
    ax_top = fig.add_subplot(gs[agg_row, 0])
    agg_masked = np.ma.masked_equal(aggregate, 0)
    ax_top.imshow(agg_masked, aspect="auto", cmap=cmap, vmin=1, vmax=vmax,
                  interpolation="nearest")
    ax_top.set_yticks([0])
    ax_top.set_yticklabels(["all orthologs"], fontsize=8)
    ax_top.set_xticks([])

    title = _title(matrix, query)
    top_ax = ax_top
    # --- secondary-structure cartoon ---
    if has_ss:
        ax_ss = fig.add_subplot(gs[0, 0], sharex=ax_top)
        _draw_ss_track(ax_ss, matrix.query_ss, n_cols)
        top_ax = ax_ss
    top_ax.set_title(title, fontsize=11, loc="left")

    # --- per-row heatmap ---
    ax = fig.add_subplot(gs[heat_row, 0], sharex=ax_top)
    masked = np.ma.masked_equal(sub, 0)
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=1, vmax=vmax,
                   interpolation="nearest")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([_short(row_label(matrix, i)) for i in rows], fontsize=7)
    ax.set_xlabel(f"{query.name} residue position")
    step = max(1, n_cols // 25)
    xticks = list(range(0, n_cols, step))
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x + 1) for x in xticks], fontsize=7)

    cax = fig.add_subplot(gs[:, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("# structures modelling residue", fontsize=8)

    paths: List[Path] = []
    for fmt in formats:
        path = out_dir / f"alignment.{fmt}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def _draw_ss_track(ax, track: Sequence[Optional[str]], n_cols: int) -> None:
    """Draw the query topology cartoon: helices as cylinders, strands as arrows.

    Columns are in imshow data coordinates (column *c* spans ``c-0.5 .. c+0.5``),
    so the cartoon lines up exactly with the heatmap beneath it.
    """
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.5])
    ax.set_yticklabels(["query SS"], fontsize=8)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    # backbone / coil baseline
    ax.plot([-0.5, n_cols - 0.5], [0.5, 0.5], color="#b0b4b8", lw=1.0, zorder=1)

    for element, c0, c1 in ss_runs(track):
        x0, x1 = c0 - 0.5, c1 + 0.5
        if element == HELIX:
            ax.add_patch(FancyBboxPatch(
                (x0, 0.28), x1 - x0, 0.44,
                boxstyle="round,pad=0,rounding_size=0.18",
                mutation_aspect=0.5, facecolor=SS_HELIX_COLOR,
                edgecolor="#7a1f1f", lw=0.6, zorder=2,
            ))
        elif element == STRAND:
            head = min(2.0, x1 - x0)
            body_end = x1 - head
            pts = [
                (x0, 0.40), (body_end, 0.40), (body_end, 0.28),
                (x1, 0.50), (body_end, 0.72), (body_end, 0.60), (x0, 0.60),
            ]
            ax.add_patch(Polygon(
                pts, closed=True, facecolor=SS_STRAND_COLOR,
                edgecolor="#9c7600", lw=0.6, zorder=2,
            ))
    # Legend in the lower band of the track: the cartoon lives in the middle
    # band (y 0.28-0.72), so this never overlaps a helix/strand and stays clear
    # of the title drawn above the axes.
    ax.legend(
        handles=[
            Patch(facecolor=SS_HELIX_COLOR, edgecolor="#7a1f1f", label="α-helix"),
            Patch(facecolor=SS_STRAND_COLOR, edgecolor="#9c7600", label="β-strand"),
        ],
        loc="lower right", ncol=2, fontsize=7, frameon=False,
        handlelength=1.2, handleheight=0.7, columnspacing=1.0, borderaxespad=0.2,
    )


def _short(text: str, limit: int = 34) -> str:
    return text if len(text) <= limit else text[: limit - 1] + "…"


def _title(matrix: CoverageMatrix, query: Query) -> str:
    agg = matrix.column_aggregate
    covered = int((agg > 0).sum())
    frac = covered / matrix.columns if matrix.columns else 0.0
    return (
        f"Structural coverage of {query.name} across orthologs  "
        f"({covered}/{matrix.columns} residues covered, {frac:.0%})"
    )
