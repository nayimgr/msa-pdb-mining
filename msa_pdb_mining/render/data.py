"""Machine-readable outputs: raw MSA, per-column summary, long coverage table, JSON."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import List

import numpy as np

from ..model import MSA, CoverageMatrix, Query


def _join_ids(ids) -> str:
    return ";".join(sorted(ids))


def write_raw_msa(out_dir: Path, msa: MSA, a3m_text: str) -> List[Path]:
    a3m_path = out_dir / "msa.a3m"
    a3m_path.write_text(a3m_text)

    fasta_path = out_dir / "msa.query_anchored.fasta"
    lines = []
    for row in msa.rows:
        lines.append(f">{row.row_id}")
        lines.append(row.match_seq)
    fasta_path.write_text("\n".join(lines) + "\n")
    return [a3m_path, fasta_path]


def write_column_summary(out_dir: Path, matrix: CoverageMatrix) -> Path:
    path = out_dir / "column_summary.csv"
    aggregate = matrix.column_aggregate
    n_orthologs = (matrix.depth > 0).sum(axis=0)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["query_resnum", "query_aa", "n_structures", "n_orthologs_covered", "pdb_ids"]
        )
        for col in range(matrix.columns):
            w.writerow(
                [
                    col + 1,
                    matrix.query_seq[col],
                    int(aggregate[col]),
                    int(n_orthologs[col]),
                    _join_ids(matrix.column_pdb_ids[col]),
                ]
            )
    return path


def write_coverage_long(out_dir: Path, matrix: CoverageMatrix) -> Path:
    """One row per (alignment row, column) cell that has >=1 observing structure."""
    path = out_dir / "coverage_long.csv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["row_label", "accession", "gene", "organism", "status",
             "query_resnum", "depth", "pdb_ids"]
        )
        nz_rows, nz_cols = np.nonzero(matrix.depth)
        for i, col in zip(nz_rows.tolist(), nz_cols.tolist()):
            w.writerow(
                [
                    matrix.row_labels[i],
                    matrix.row_accessions[i] or "",
                    matrix.row_genes[i] or "",
                    matrix.row_organisms[i] or "",
                    matrix.row_status[i],
                    col + 1,
                    int(matrix.depth[i, col]),
                    _join_ids(matrix.pdb_ids[i][col]),
                ]
            )
    return path


def write_summary_json(
    out_dir: Path, matrix: CoverageMatrix, query: Query, params: dict
) -> Path:
    path = out_dir / "summary.json"
    aggregate = matrix.column_aggregate
    covered = int((aggregate > 0).sum())
    structured_rows = int(sum(1 for s in matrix.row_status if s == "structured"))
    all_ids = set()
    for ids in matrix.column_pdb_ids:
        all_ids |= ids

    summary = {
        "query": {
            "name": query.name,
            "accession": query.accession,
            "gene": query.gene,
            "organism": query.organism,
            "tax_id": query.tax_id,
            "length": query.length,
        },
        "params": params,
        "n_rows": matrix.n_rows,
        "n_structured_rows": structured_rows,
        "n_distinct_pdb_entries": len(all_ids),
        "columns_with_any_structure": covered,
        "columns_total": matrix.columns,
        "fraction_query_covered": round(covered / matrix.columns, 4) if matrix.columns else 0.0,
        "max_column_depth": int(aggregate.max()) if matrix.columns else 0,
        "column_aggregate": [int(x) for x in aggregate],
    }
    path.write_text(json.dumps(summary, indent=2))
    return path


def write_data_outputs(
    out_dir: Path, msa: MSA, matrix: CoverageMatrix, query: Query, a3m_text: str, params: dict
) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    paths += write_raw_msa(out_dir, msa, a3m_text)
    paths.append(write_column_summary(out_dir, matrix))
    paths.append(write_coverage_long(out_dir, matrix))
    paths.append(write_summary_json(out_dir, matrix, query, params))
    return paths
