"""Shared helpers for selecting/ordering rows for visual renderers."""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..model import CoverageMatrix


def select_rows(matrix: CoverageMatrix, max_rows: int) -> List[int]:
    """Query first, then rows with structural coverage by descending total depth.

    Rows with no covered residues are omitted from the visuals (they carry no
    structural signal); the full data is preserved in the CSV/JSON outputs.
    """
    totals = matrix.depth.sum(axis=1)
    others = [i for i in range(matrix.n_rows) if i != matrix.query_index and totals[i] > 0]
    others.sort(key=lambda i: totals[i], reverse=True)
    ordered = [matrix.query_index] + others
    return ordered[: max(1, max_rows)]


def ss_runs(track: Sequence[Optional[str]]) -> List[Tuple[str, int, int]]:
    """Collapse a per-column SS track into ``(element, start_col, end_col)`` runs.

    Consecutive columns with the same non-``None`` element become one segment —
    what the renderers draw as a single helix cylinder or strand arrow.
    """
    runs: List[Tuple[str, int, int]] = []
    start: Optional[int] = None
    cur: Optional[str] = None
    for col, el in enumerate(track):
        if el != cur:
            if start is not None and cur is not None:
                runs.append((cur, start, col - 1))
            start = col if el is not None else None
            cur = el
    if start is not None and cur is not None:
        runs.append((cur, start, len(track) - 1))
    return runs


def row_label(matrix: CoverageMatrix, i: int) -> str:
    label = matrix.row_accessions[i] or matrix.row_labels[i]
    gene = matrix.row_genes[i]
    if gene:
        label = f"{label} {gene}"
    org = matrix.row_organisms[i]
    if org:
        label = f"{label} · {org}"
    return label
