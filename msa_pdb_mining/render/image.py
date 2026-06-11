"""Static publication-style coverage heatmap (PNG + SVG) via matplotlib."""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from ..config import Config  # noqa: E402
from ..model import CoverageMatrix, Query  # noqa: E402
from .layout import row_label, select_rows  # noqa: E402

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
    width = min(40.0, max(8.0, n_cols * 0.05))
    height = min(30.0, max(3.5, (n_rows + 3) * 0.26))

    fig = plt.figure(figsize=(width, height))
    gs = GridSpec(
        2, 2, height_ratios=[1, max(n_rows, 1)], width_ratios=[40, 1],
        hspace=0.06, wspace=0.02,
    )
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(_BG)

    # --- aggregate track ---
    ax_top = fig.add_subplot(gs[0, 0])
    agg_masked = np.ma.masked_equal(aggregate, 0)
    ax_top.imshow(agg_masked, aspect="auto", cmap=cmap, vmin=1, vmax=vmax,
                  interpolation="nearest")
    ax_top.set_yticks([0])
    ax_top.set_yticklabels(["all orthologs"], fontsize=8)
    ax_top.set_xticks([])
    ax_top.set_title(_title(matrix, query), fontsize=11, loc="left")

    # --- per-row heatmap ---
    ax = fig.add_subplot(gs[1, 0], sharex=ax_top)
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
