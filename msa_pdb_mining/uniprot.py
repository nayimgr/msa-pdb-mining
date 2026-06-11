"""Batched UniProt metadata lookup for MSA hit accessions.

One TSV request to ``/uniprotkb/accessions`` returns, for up to ~100 accessions,
their review status (Swiss-Prot vs TrEMBL), organism, primary gene name, and the
list of cross-referenced PDB ids. We use this to (a) label rows with organism/gene,
(b) restrict to reviewed entries when asked, and (c) pre-filter which accessions
actually have experimental structures — so the expensive per-residue PDBe lookup
runs only where there is something to find.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set

import httpx

from .cache import DiskCache
from .config import Config
from .net import request

log = logging.getLogger(__name__)

_FIELDS = "accession,reviewed,organism_name,gene_primary,xref_pdb"
_CHUNK = 100


@dataclass
class AccMeta:
    reviewed: bool = False
    organism: Optional[str] = None
    gene: Optional[str] = None
    pdb_ids: Set[str] = field(default_factory=set)
    found: bool = False

    @property
    def has_pdb(self) -> bool:
        return bool(self.pdb_ids)

    def to_dict(self) -> dict:
        return {
            "reviewed": self.reviewed,
            "organism": self.organism,
            "gene": self.gene,
            "pdb_ids": sorted(self.pdb_ids),
            "found": self.found,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "AccMeta":
        return cls(
            reviewed=bool(d.get("reviewed")),
            organism=d.get("organism"),
            gene=d.get("gene"),
            pdb_ids=set(d.get("pdb_ids", [])),
            found=bool(d.get("found", True)),
        )


def _clean_organism(text: str) -> Optional[str]:
    """Scientific name only: 'Homo sapiens (Human)' -> 'Homo sapiens'."""
    s = text.strip()
    if not s:
        return None
    i = s.find(" (")
    return s[:i] if i > 0 else s


def parse_tsv(text: str) -> Dict[str, AccMeta]:
    """Parse the accessions TSV (header: Entry, Reviewed, Organism, Gene, PDB)."""
    out: Dict[str, AccMeta] = {}
    lines = text.splitlines()
    for line in lines[1:]:  # skip header row
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 5:
            cols += [""] * (5 - len(cols))
        acc = cols[0].strip()
        if not acc:
            continue
        out[acc] = AccMeta(
            reviewed=cols[1].strip().lower() == "reviewed",
            organism=_clean_organism(cols[2]),
            gene=cols[3].strip() or None,
            pdb_ids={p.strip().lower() for p in cols[4].split(";") if p.strip()},
            found=True,
        )
    return out


def _chunks(items: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _fetch_chunk(client: httpx.Client, config: Config, chunk: List[str]) -> Dict[str, AccMeta]:
    resp = request(
        client,
        "GET",
        f"{config.uniprot_rest}/uniprotkb/accessions",
        params={"accessions": ",".join(chunk), "fields": _FIELDS, "format": "tsv"},
        headers={"Accept": "text/plain"},
    )
    if resp is None or resp.status_code >= 400:
        return {}
    return parse_tsv(resp.text)


def batch_lookup(
    client: httpx.Client, cache: DiskCache, config: Config, accessions: List[str]
) -> Dict[str, AccMeta]:
    """Return metadata for each accession, using a per-accession cache + batched fetch."""
    result: Dict[str, AccMeta] = {}
    misses: List[str] = []
    for acc in accessions:
        hit = cache.get_json("uniprot_meta", acc)
        if hit is not None:
            result[acc] = AccMeta.from_dict(hit)
        else:
            misses.append(acc)

    for chunk in _chunks(misses, _CHUNK):
        parsed = _fetch_chunk(client, config, chunk)
        for acc in chunk:
            meta = parsed.get(acc, AccMeta(found=False))
            cache.set_json("uniprot_meta", acc, meta.to_dict())
            result[acc] = meta

    n_pdb = sum(1 for m in result.values() if m.has_pdb)
    n_rev = sum(1 for m in result.values() if m.reviewed)
    log.info(
        "UniProt metadata: %d accessions (%d reviewed, %d with PDB)",
        len(result), n_rev, n_pdb,
    )
    return result
