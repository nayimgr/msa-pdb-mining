"""Curated-ortholog MSA source.

The ColabFold folding-MSA filters out redundant near-identical orthologs and
returns UniRef cluster representatives, so curated Swiss-Prot orthologs (mouse,
rat, ...) are largely absent. This module instead pulls the reviewed members of
the query's protein family directly from UniProt and aligns them, guaranteeing
the curated orthologs appear. (Close paralogs sharing the family — e.g. p63/p73
for p53 — are included too, distinguishable by gene name.)

Two alignment engines, selected by ``Config.ortholog_aligner``:
  * ``"famsa"`` (default) — one **true multiple alignment** of the query plus
    every ortholog in a single FAMSA pass (in-process via pyfamsa), so columns
    are family-aware. The MSA is then projected back to the query-anchored A3M
    the rest of the pipeline expects (``a3m_from_aligned``).
  * ``"pairwise"`` — the original **star** alignment: each ortholog aligned to
    the query independently (``a3m_from_pairwise``).

Either way the result is a query-anchored A3M whose every row's ungapped
sequence equals the full UniProt sequence, so the downstream offset/depth
projection is identical for both engines.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import httpx
from Bio.Align import PairwiseAligner, substitution_matrices
from tqdm import tqdm

from .cache import DiskCache
from .config import Config
from .model import MSA, MSARow, Query
from .msa.a3m import match_columns
from .net import request
from .uniprot import AccMeta, _clean_organism

log = logging.getLogger(__name__)

_FIELDS = "accession,reviewed,organism_name,gene_primary,sequence,xref_pdb"


@dataclass
class Ortholog:
    accession: str
    organism: Optional[str]
    gene: Optional[str]
    sequence: str
    pdb_ids: Set[str]


# --- family resolution ---
def fetch_query_family(
    client: httpx.Client, cache: DiskCache, config: Config, accession: str
) -> Optional[str]:
    """Most-specific family from the query's SIMILARITY comment.

    'Belongs to the protein kinase superfamily. Tyr ... family. INSR subfamily.'
    -> 'INSR subfamily' (tightest grouping = closest orthologs).
    """
    cached = cache.get_json("uniprot_family", accession)
    if cached is not None:
        return cached.get("family")

    family = None
    resp = request(
        client, "GET", f"{config.uniprot_rest}/uniprotkb/{accession}.json",
        params={"fields": "cc_similarity"}, headers={"Accept": "application/json"},
    )
    if resp is not None and resp.status_code < 400:
        try:
            for c in resp.json().get("comments", []):
                if c.get("commentType") != "SIMILARITY":
                    continue
                for t in c.get("texts", []):
                    value = t.get("value", "")
                    if value.startswith("Belongs to the "):
                        body = value[len("Belongs to the "):].rstrip(".")
                        parts = [p.strip() for p in body.split(".") if p.strip()]
                        if parts:
                            family = parts[-1]
        except ValueError:
            pass
    cache.set_json("uniprot_family", accession, {"family": family})
    return family


def parse_orthologs(text: str, exclude: Optional[str]) -> List[Ortholog]:
    out: List[Ortholog] = []
    for line in text.splitlines()[1:]:  # skip header
        cols = line.split("\t")
        if len(cols) < 6:
            cols += [""] * (6 - len(cols))
        acc = cols[0].strip()
        seq = cols[4].strip().upper()
        if not acc or acc == exclude or not seq:
            continue
        out.append(
            Ortholog(
                accession=acc,
                organism=_clean_organism(cols[2]),
                gene=cols[3].strip() or None,
                sequence=seq,
                pdb_ids={p.strip().lower() for p in cols[5].split(";") if p.strip()},
            )
        )
    return out


def fetch_reviewed_orthologs(
    client: httpx.Client, cache: DiskCache, config: Config, query: Query
) -> List[Ortholog]:
    family = fetch_query_family(client, cache, config, query.accession) if query.accession else None
    if family:
        uq = f'family:"{family}" AND reviewed:true'
    elif query.gene:
        uq = f"gene:{query.gene} AND reviewed:true"
    else:
        raise ValueError("Ortholog mode needs a UniProt accession or gene to define the family")
    log.info("Ortholog UniProt query: %s", uq)

    text = cache.get_text("orthologs", uq)
    if text is None:
        resp = request(
            client, "GET", f"{config.uniprot_rest}/uniprotkb/search",
            params={"query": uq, "fields": _FIELDS, "format": "tsv", "size": str(config.max_orthologs)},
            headers={"Accept": "text/plain"},
        )
        if resp is None or resp.status_code >= 400 or not resp.text.strip():
            raise ValueError(f"UniProt ortholog query returned nothing: {uq}")
        text = resp.text
        cache.set_text("orthologs", uq, text)
    return parse_orthologs(text, exclude=query.accession)


# --- pairwise alignment to the query ---
def _make_aligner() -> Tuple[PairwiseAligner, Set[str]]:
    m = substitution_matrices.load("BLOSUM62")
    a = PairwiseAligner()
    a.substitution_matrix = m
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    return a, set(m.alphabet)


def _sanitize(seq: str, alphabet: Set[str]) -> str:
    return "".join(c if c in alphabet else "X" for c in seq)


# FAMSA validates its input against the standard amino-acid alphabet plus the
# ambiguity codes B/Z/X; anything else (selenocysteine U, pyrrolysine O, ...) is
# rejected and must be substituted before alignment. We restore the real residues
# from the original sequence afterwards, so this only affects scoring/gap
# placement, never the residues emitted into the A3M.
_FAMSA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZX")


def a3m_from_pairwise(query_aln: str, ortho_aln: str, ortho_seq: str) -> str:
    """Build a query-anchored A3M row from a pairwise alignment.

    Uses the *original* ortholog residues (``ortho_seq``) so the row's ungapped
    sequence equals the full ortholog sequence (offset detection stays exact even
    when scoring used sanitised residues). Uppercase = match column, ``-`` =
    deletion, lowercase = insertion relative to the query.
    """
    parts: List[str] = []
    oi = 0
    for qc, oc in zip(query_aln, ortho_aln):
        real = None
        if oc != "-":
            real = ortho_seq[oi]
            oi += 1
        if qc != "-" and real is not None:
            parts.append(real.upper())
        elif qc != "-" and real is None:
            parts.append("-")
        elif qc == "-" and real is not None:
            parts.append(real.lower())
    return "".join(parts)


# --- true multiple alignment (FAMSA) ---
def a3m_from_aligned(query_aln: str, row_aln: str, row_seq: str) -> str:
    """Project one row of a true MSA into a query-anchored A3M string.

    The inverse of :func:`match_columns`. ``query_aln`` / ``row_aln`` are aligned
    (equal-length, ``-``-gapped) rows from a multiple alignment; ``row_seq`` is the
    row's *original* ungapped sequence. A column where the query has a residue is a
    **match column** — emit the row's residue uppercase, or ``-`` if the row is
    gapped there; a column where the query is gapped is an **insertion relative to
    the query** — emit the row's residue lowercase (nothing if both are gapped).

    Residue identities come from ``row_seq``, not from ``row_aln``, so the row's
    ungapped sequence stays byte-identical to the UniProt sequence even if the
    aligner normalised any residue — keeping offset detection exact.
    """
    parts: List[str] = []
    oi = 0
    for qc, rc in zip(query_aln, row_aln):
        real = None
        if rc != "-":
            real = row_seq[oi]
            oi += 1
        if qc != "-":  # match column (query has a residue here)
            parts.append(real.upper() if real is not None else "-")
        elif real is not None:  # insertion relative to the query
            parts.append(real.lower())
        # else: both gapped -> contributes nothing
    return "".join(parts)


def _run_famsa(
    query_seq: str, ortholog_seqs: List[str]
) -> Tuple[str, List[Optional[str]]]:
    """Multiple-align the query + orthologs in one FAMSA pass (in-process).

    Returns the gapped query string and a list of gapped ortholog strings in the
    same order as ``ortholog_seqs`` (``None`` for any sequence FAMSA failed to
    return). Only the gap patterns matter here — residue identities are pulled
    from the originals by :func:`a3m_from_aligned`.
    """
    # Imported lazily: only the live FAMSA path needs the native extension, so the
    # pure helpers above (and their offline tests) don't require pyfamsa installed.
    from pyfamsa import Aligner, Sequence

    # Sanitise only for FAMSA's alphabet check; a3m_from_aligned restores the real
    # residues from the originals, so the emitted rows are unaffected.
    inputs = [Sequence(b"q", _sanitize(query_seq, _FAMSA_ALPHABET).encode())]
    inputs += [
        Sequence(f"o{i}".encode(), _sanitize(seq, _FAMSA_ALPHABET).encode())
        for i, seq in enumerate(ortholog_seqs)
    ]
    # keep_duplicates=True guarantees every input reappears in the output, so an
    # ortholog identical to another never silently vanishes.
    alignment = Aligner(keep_duplicates=True).align(inputs)
    by_id = {seq.id.decode(): seq.sequence.decode() for seq in alignment}
    query_aln = by_id["q"]
    ortho_alns = [by_id.get(f"o{i}") for i in range(len(ortholog_seqs))]
    return query_aln, ortho_alns


def _align_famsa(query: Query, orthologs: List[Ortholog]) -> List[Optional[str]]:
    """Query-anchored A3M per ortholog from a single FAMSA multiple alignment."""
    if not orthologs:
        return []
    log.info(
        "FAMSA: aligning query + %d orthologs in one multiple alignment", len(orthologs)
    )
    query_aln, ortho_alns = _run_famsa(query.sequence, [o.sequence for o in orthologs])
    out: List[Optional[str]] = []
    for o, oa in zip(orthologs, ortho_alns):
        if oa is None:  # pragma: no cover - FAMSA should return every input
            log.warning("FAMSA returned no row for %s; skipping", o.accession)
            out.append(None)
        else:
            out.append(a3m_from_aligned(query_aln, oa, o.sequence))
    return out


def _align_pairwise(query: Query, orthologs: List[Ortholog]) -> List[Optional[str]]:
    """Query-anchored A3M per ortholog from independent pairwise alignments (star)."""
    aligner, alphabet = _make_aligner()
    q_san = _sanitize(query.sequence, alphabet)
    out: List[Optional[str]] = []
    for o in tqdm(orthologs, desc="Aligning orthologs (pairwise)", unit="seq"):
        try:
            aln = aligner.align(q_san, _sanitize(o.sequence, alphabet))[0]
        except Exception as exc:  # pragma: no cover - defensive
            log.warning("alignment failed for %s: %s", o.accession, exc)
            out.append(None)
            continue
        out.append(a3m_from_pairwise(aln[0], aln[1], o.sequence))
    return out


def build_ortholog_msa(
    client: httpx.Client, cache: DiskCache, config: Config, query: Query
) -> Tuple[MSA, Dict[str, AccMeta], str]:
    orthologs = fetch_reviewed_orthologs(client, cache, config, query)
    log.info("Fetched %d reviewed family members (orthologs + paralogs)", len(orthologs))

    if config.ortholog_aligner == "pairwise":
        a3ms = _align_pairwise(query, orthologs)
    else:
        a3ms = _align_famsa(query, orthologs)

    rows: List[MSARow] = [
        MSARow(
            row_id=query.name, raw_header="query",
            a3m_seq=query.sequence, match_seq=query.sequence,
            accession=query.accession, organism=query.organism, gene=query.gene,
            reviewed=True, is_query=True,
        )
    ]
    meta: Dict[str, AccMeta] = {}
    fasta_chunks = [f">{query.name}\n{query.sequence}"]

    for o, a3m in zip(orthologs, a3ms):
        if a3m is None:  # alignment failed / dropped for this ortholog
            continue
        rows.append(
            MSARow(
                row_id=o.accession, raw_header=o.accession,
                a3m_seq=a3m, match_seq=match_columns(a3m),
                accession=o.accession, organism=o.organism, gene=o.gene, reviewed=True,
            )
        )
        meta[o.accession] = AccMeta(
            reviewed=True, organism=o.organism, gene=o.gene, pdb_ids=o.pdb_ids, found=True
        )
        fasta_chunks.append(f">{o.accession}\n{a3m}")

    a3m_text = "\n".join(fasta_chunks) + "\n"
    return MSA(query=query, rows=rows), meta, a3m_text
