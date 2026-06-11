"""Shared helpers for selecting/ordering rows for visual renderers."""

from __future__ import annotations

from typing import List

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


def row_label(matrix: CoverageMatrix, i: int) -> str:
    label = matrix.row_accessions[i] or matrix.row_labels[i]
    gene = matrix.row_genes[i]
    if gene:
        label = f"{label} {gene}"
    org = matrix.row_organisms[i]
    if org:
        label = f"{label} · {org}"
    return label
