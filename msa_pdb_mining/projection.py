"""Project per-accession structural depth onto query-anchored MSA columns."""

from __future__ import annotations

from typing import Dict, List, Optional, Set

import numpy as np

from .model import (
    STATUS_NO_STRUCTURES,
    STATUS_QUERY,
    STATUS_STRUCTURED,
    STATUS_UNKNOWN,
    STATUS_UNMAPPED,
    CoverageMatrix,
    MSA,
    StructureCoverage,
)
from .msa.a3m import find_residue_offset, iter_match_residues
from .secondary_structure import SecondaryStructure


def _project_query_ss(
    msa: MSA, n_cols: int, query_ss: Optional[SecondaryStructure]
) -> Optional[List[Optional[str]]]:
    """Map the query's per-residue SS elements onto query-anchored columns.

    Uses the same residue/column walk as the coverage projection: locate the
    query within the SS endpoint's full sequence, then place each consensus
    element at its column. Returns ``None`` if there is nothing to draw.
    """
    if query_ss is None or not query_ss.has_elements or not query_ss.full_sequence:
        return None
    query_row = msa.rows[0]
    offset = find_residue_offset(query_row.a3m_seq, query_ss.full_sequence)
    if offset is None:
        return None
    track: List[Optional[str]] = [None] * n_cols
    placed = False
    for col, resnum in iter_match_residues(query_row.a3m_seq, offset):
        if col >= n_cols:
            break
        element = query_ss.element(resnum)
        if element is not None:
            track[col] = element
            placed = True
    return track if placed else None


def build_coverage_matrix(
    msa: MSA,
    coverage_by_acc: Dict[str, StructureCoverage],
    query_ss: Optional[SecondaryStructure] = None,
) -> CoverageMatrix:
    n_cols = msa.length
    n_rows = len(msa.rows)

    depth = np.zeros((n_rows, n_cols), dtype=int)
    pdb_ids: List[List[Set[str]]] = [
        [set() for _ in range(n_cols)] for _ in range(n_rows)
    ]
    column_pdb_ids: List[Set[str]] = [set() for _ in range(n_cols)]

    row_labels: List[str] = []
    row_accessions: List[Optional[str]] = []
    row_organisms: List[Optional[str]] = []
    row_genes: List[Optional[str]] = []
    row_status: List[str] = []

    for i, row in enumerate(msa.rows):
        acc = row.accession
        cov = coverage_by_acc.get(acc) if acc else None

        if row.is_query:
            status = STATUS_QUERY
        elif acc is None:
            status = STATUS_UNKNOWN
        elif cov is None or not cov.has_structures:
            status = STATUS_NO_STRUCTURES
        else:
            status = STATUS_STRUCTURED

        if cov is not None and cov.has_structures and cov.full_sequence:
            offset = find_residue_offset(row.a3m_seq, cov.full_sequence)
            if offset is None:
                if not row.is_query:
                    status = STATUS_UNMAPPED
            else:
                for col, resnum in iter_match_residues(row.a3m_seq, offset):
                    if col >= n_cols:  # guard against a malformed over-long row
                        break
                    d = cov.depth(resnum)
                    if d:
                        depth[i, col] = d
                        ids = cov.pdb_ids(resnum)
                        pdb_ids[i][col] = ids
                        column_pdb_ids[col] |= ids

        row_labels.append(row.row_id)
        row_accessions.append(acc)
        row_organisms.append(row.organism)
        row_genes.append(row.gene)
        row_status.append(status)

    return CoverageMatrix(
        columns=n_cols,
        query_seq=msa.query.sequence,
        row_labels=row_labels,
        row_accessions=row_accessions,
        row_organisms=row_organisms,
        row_genes=row_genes,
        row_status=row_status,
        depth=depth,
        pdb_ids=pdb_ids,
        column_pdb_ids=column_pdb_ids,
        query_index=0,
        query_ss=_project_query_ss(msa, n_cols, query_ss),
    )
