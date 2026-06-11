"""Per-accession experimental-structure coverage from PDBe SIFTS (unipdb).

``graph-api/uniprot/unipdb/<acc>`` returns, per PDB structure, ranges of UniProt
residues tagged ``observed: "Y"`` (modelled) or ``"N"`` (disordered/missing). We
aggregate the ``observed:"Y"`` ranges into a per-residue set of PDB ids; depth at
a residue is the number of distinct structures observing it. Because the ranges
are in UniProt coordinates, affinity tags / fusions (which have no UniProt
mapping) never appear here, and PDB author numbering is irrelevant.
"""

from __future__ import annotations

import logging
from typing import Any, Optional

import httpx

from .cache import DiskCache
from .config import Config
from .model import StructureCoverage
from .net import request

log = logging.getLogger(__name__)


def build_coverage(accession: str, raw: Optional[Any]) -> StructureCoverage:
    """Construct a StructureCoverage from a raw unipdb JSON response (pure)."""
    cov = StructureCoverage(accession=accession, full_sequence="")
    if not raw:
        return cov

    entry = raw.get(accession) if isinstance(raw, dict) else None
    if entry is None and isinstance(raw, dict) and raw:
        entry = next(iter(raw.values()))
    if not isinstance(entry, dict):
        return cov

    cov.full_sequence = entry.get("sequence", "") or ""
    for struct in entry.get("data", []) or []:
        pdb_id = (struct.get("accession") or struct.get("name") or "").lower()
        if not pdb_id:
            continue
        add = struct.get("additionalData") or {}
        cov.method[pdb_id] = add.get("experiment")
        cov.resolution[pdb_id] = add.get("resolution")
        for r in struct.get("residues", []) or []:
            if r.get("indexType") != "UNIPROT" or r.get("observed") != "Y":
                continue
            start, end = r.get("startIndex"), r.get("endIndex")
            if start is None or end is None:
                continue
            for res in range(int(start), int(end) + 1):
                cov.pdb_ids_per_res.setdefault(res, set()).add(pdb_id)
    return cov


def fetch_coverage(
    client: httpx.Client, cache: DiskCache, config: Config, accession: str
) -> StructureCoverage:
    """Fetch (cached) the unipdb response for an accession and build coverage.

    A definitive 404 (no SIFTS mapping = no structures) is cached as ``{}`` so we
    never re-query it; transient failures return ``None`` and are *not* cached.
    """
    url = f"{config.pdbe_graph_api}/uniprot/unipdb/{accession}"

    def _fetch():
        resp = request(client, "GET", url)
        if resp is None:
            return None  # transient (connection / exhausted retries): don't cache
        if resp.status_code == 404:
            return {}  # definitive: accession has no experimental structures
        if resp.status_code >= 400:
            return None
        try:
            return resp.json()
        except ValueError:
            return None

    raw = cache.cached_json("unipdb", accession, _fetch)
    cov = build_coverage(accession, raw)
    if cov.has_structures:
        log.info("%s: %d structures", accession, cov.n_structures)
    return cov
