"""Core data structures shared across the pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

import numpy as np


@dataclass
class Query:
    """The protein the run is anchored on."""

    sequence: str
    accession: Optional[str] = None
    gene: Optional[str] = None
    organism: Optional[str] = None
    tax_id: Optional[int] = None
    name: str = "query"

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class MSARow:
    """One sequence of the alignment, anchored to the query columns.

    ``a3m_seq`` is the raw A3M string (uppercase = aligned/match residue,
    lowercase = insertion relative to query, ``-`` = deletion). ``match_seq`` is
    the query-anchored projection: uppercase + ``-`` only, so ``len(match_seq)``
    equals the query length for every row.
    """

    row_id: str
    raw_header: str
    a3m_seq: str
    match_seq: str
    accession: Optional[str] = None
    organism: Optional[str] = None
    gene: Optional[str] = None
    reviewed: bool = False
    is_query: bool = False


@dataclass
class MSA:
    query: Query
    rows: List[MSARow]  # rows[0] is the query

    @property
    def length(self) -> int:
        return self.query.length


@dataclass
class StructureCoverage:
    """Per-residue experimental-structure coverage for one UniProt accession.

    ``pdb_ids_per_res`` maps a 1-based UniProt residue number to the set of PDB
    entries that *observe* (model) it. Depth at a residue is the size of that set.
    """

    accession: str
    full_sequence: str
    pdb_ids_per_res: Dict[int, Set[str]] = field(default_factory=dict)
    method: Dict[str, str] = field(default_factory=dict)  # pdb_id -> experiment
    resolution: Dict[str, Optional[float]] = field(default_factory=dict)

    @property
    def has_structures(self) -> bool:
        return bool(self.pdb_ids_per_res)

    @property
    def n_structures(self) -> int:
        seen: Set[str] = set()
        for ids in self.pdb_ids_per_res.values():
            seen |= ids
        return len(seen)

    def depth(self, res: int) -> int:
        return len(self.pdb_ids_per_res.get(res, ()))

    def pdb_ids(self, res: int) -> Set[str]:
        return self.pdb_ids_per_res.get(res, set())


# Row-level structural status, used for colouring distinct cases.
STATUS_QUERY = "query"
STATUS_STRUCTURED = "structured"  # accession resolved, >=1 PDB
STATUS_NO_STRUCTURES = "no_structures"  # accession resolved, 0 PDB
STATUS_UNKNOWN = "unknown"  # no UniProt accession resolved
STATUS_UNMAPPED = "unmapped"  # accession known but residues couldn't be located


@dataclass
class CoverageMatrix:
    """Structural-coverage depth projected onto query-anchored alignment columns."""

    columns: int  # == query length
    query_seq: str
    row_labels: List[str]
    row_accessions: List[Optional[str]]
    row_organisms: List[Optional[str]]
    row_genes: List[Optional[str]]
    row_status: List[str]
    depth: np.ndarray  # shape [n_rows, columns], int; structures observing each cell
    pdb_ids: List[List[Set[str]]]  # [n_rows][columns] sets of pdb ids
    column_pdb_ids: List[Set[str]]  # union of pdb ids across all rows, per column
    query_index: int = 0

    @property
    def n_rows(self) -> int:
        return self.depth.shape[0]

    @property
    def column_aggregate(self) -> np.ndarray:
        """Distinct structures (across all orthologs) covering each query column."""
        return np.array([len(ids) for ids in self.column_pdb_ids], dtype=int)
