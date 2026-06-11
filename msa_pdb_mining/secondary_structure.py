"""Query secondary structure (helices/strands) from PDBe SIFTS, in UniProt coords.

``graph-api/uniprot/secondary_structures/<acc>`` returns, across every PDB
structure of an accession, the residue ranges assigned to a ``Helix`` or a
``Strand`` — each range already mapped to **UniProt** residue numbers
(``indexType == "UNIPROT"``). Because the tool's whole coordinate system is
UniProt numbering (see CLAUDE.md), these ranges drop straight onto the query
columns of the alignment, letting the renderers draw a topology cartoon (helices
and strand-arrows) on top of the coverage heatmap.

The same residue can be a helix in one deposition and a strand in another, so we
keep a per-residue tally of each assignment and take the **majority** (ties ->
helix). The result is a single consensus element per residue: ``"H"`` or ``"E"``.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import httpx

from .cache import DiskCache
from .config import Config
from .net import request

log = logging.getLogger(__name__)

# Single-letter element codes (DSSP-style: H = helix, E = extended/strand).
HELIX = "H"
STRAND = "E"


@dataclass
class SecondaryStructure:
    """Consensus secondary-structure element per UniProt residue for one accession.

    ``element_per_res`` maps a 1-based UniProt residue number to ``HELIX`` or
    ``STRAND``; residues with no assignment (coil / not in any structure) are
    simply absent.
    """

    accession: str
    full_sequence: str = ""
    element_per_res: Dict[int, str] = field(default_factory=dict)

    @property
    def has_elements(self) -> bool:
        return bool(self.element_per_res)

    def element(self, res: int) -> Optional[str]:
        return self.element_per_res.get(res)


def build_secondary_structure(
    accession: str, raw: Optional[Any]
) -> SecondaryStructure:
    """Construct a SecondaryStructure from a raw secondary_structures response (pure)."""
    ss = SecondaryStructure(accession=accession)
    if not raw:
        return ss

    entry = raw.get(accession) if isinstance(raw, dict) else None
    if entry is None and isinstance(raw, dict) and raw:
        entry = next(iter(raw.values()))
    if not isinstance(entry, dict):
        return ss

    ss.full_sequence = entry.get("sequence", "") or ""

    # Tally each residue's helix vs strand assignments across all depositions,
    # then resolve to a single consensus element (ties favour helix).
    helix: Dict[int, int] = defaultdict(int)
    strand: Dict[int, int] = defaultdict(int)
    for block in entry.get("data", []) or []:
        name = (block.get("name") or block.get("dataType") or "").lower()
        if "helix" in name:
            tally = helix
        elif "strand" in name:
            tally = strand
        else:
            continue  # ignore turns / other element kinds
        for r in block.get("residues", []) or []:
            if r.get("indexType") != "UNIPROT":
                continue
            start, end = r.get("startIndex"), r.get("endIndex")
            if start is None or end is None:
                continue
            for res in range(int(start), int(end) + 1):
                tally[res] += 1

    for res in set(helix) | set(strand):
        ss.element_per_res[res] = HELIX if helix[res] >= strand[res] else STRAND
    return ss


def fetch_secondary_structure(
    client: httpx.Client, cache: DiskCache, config: Config, accession: str
) -> SecondaryStructure:
    """Fetch (cached) the secondary_structures response for an accession.

    Caching follows the same contract as ``structures.fetch_coverage``: a
    definitive 404 (no SS mapping) is cached as ``{}``; transient failures return
    ``None`` and are not cached, so they retry on the next run.
    """
    url = f"{config.pdbe_graph_api}/uniprot/secondary_structures/{accession}"

    def _fetch():
        resp = request(client, "GET", url)
        if resp is None:
            return None  # transient: don't cache
        if resp.status_code == 404:
            return {}  # definitive: no secondary-structure mapping
        if resp.status_code >= 400:
            return None
        try:
            return resp.json()
        except ValueError:
            return None

    raw = cache.cached_json("unipdb_ss", accession, _fetch)
    ss = build_secondary_structure(accession, raw)
    if ss.has_elements:
        n_h = sum(1 for v in ss.element_per_res.values() if v == HELIX)
        n_e = len(ss.element_per_res) - n_h
        log.info("%s: secondary structure %d helix / %d strand residues",
                 accession, n_h, n_e)
    return ss
