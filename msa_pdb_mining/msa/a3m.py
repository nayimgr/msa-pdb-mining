"""A3M parsing and query-anchored projection.

A3M convention (as produced by the ColabFold MMseqs2 server):
  * the first record is the query; it is all-uppercase, no gaps, and defines the
    match columns;
  * for every other record, uppercase letters are residues aligned to a query
    column, lowercase letters are insertions relative to the query, and ``-`` is
    a deletion. ``.`` is sometimes used as insertion padding.

The "query-anchored" view keeps only match columns (uppercase + ``-``), dropping
insertions, so every row has exactly ``len(query)`` columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Optional, Tuple


@dataclass
class A3MRecord:
    header: str  # text after '>' on the header line
    seq: str  # raw A3M aligned sequence (may contain upper/lower/-/.)

    @property
    def first_token(self) -> str:
        return self.header.split()[0] if self.header.split() else self.header


def parse_a3m(text: str) -> List[A3MRecord]:
    records: List[A3MRecord] = []
    header: Optional[str] = None
    chunks: List[str] = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(A3MRecord(header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        records.append(A3MRecord(header, "".join(chunks)))
    return records


def match_columns(a3m_seq: str) -> str:
    """Return the query-anchored projection: uppercase residues + ``-`` only."""
    return "".join(c for c in a3m_seq if c == "-" or (c.isalpha() and c.isupper()))


def ungapped(a3m_seq: str) -> str:
    """All residues (upper- and lowercase) as the contiguous hit subsequence."""
    return "".join(c.upper() for c in a3m_seq if c.isalpha())


def iter_match_residues(
    a3m_seq: str, first_residue: int = 1
) -> Iterator[Tuple[int, int]]:
    """Yield ``(column, residue_number)`` for each match (uppercase) position.

    ``column`` is the 0-based index into the query-anchored alignment.
    ``residue_number`` is the 1-based position in the hit's *full* sequence,
    where the first residue consumed is ``first_residue``. Insertions (lowercase)
    consume a hit residue but produce no column; deletions (``-``) advance the
    column but consume no residue.
    """
    col = -1
    resnum = first_residue - 1
    for c in a3m_seq:
        if c == "-":
            col += 1
        elif c == ".":
            continue  # insertion padding: neither column nor residue
        elif c.isupper() and c.isalpha():
            col += 1
            resnum += 1
            yield col, resnum
        elif c.islower():
            resnum += 1  # insertion residue: consumes a residue, no column


def find_residue_offset(a3m_seq: str, full_sequence: str) -> Optional[int]:
    """1-based residue number of the hit's first residue within ``full_sequence``.

    The hit residues (``ungapped``) are a contiguous slice of the full UniProt
    sequence, so a substring search recovers the offset. Returns ``None`` if the
    slice can't be located exactly (e.g. sequence-version mismatch).
    """
    hit = ungapped(a3m_seq)
    if not hit:
        return None
    idx = full_sequence.find(hit)
    if idx >= 0:
        return idx + 1
    return None
