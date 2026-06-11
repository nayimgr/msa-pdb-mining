"""End-to-end orchestration: input -> MSA -> structures -> projection -> render."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

from tqdm import tqdm

from .accession_map import extract_accession
from .cache import DiskCache
from .config import Config
from .input_resolver import resolve_query
from .model import MSA, MSARow, Query, StructureCoverage
from .msa.a3m import match_columns, parse_a3m
from .msa.colabfold import ColabFoldMMseqs2Backend
from .net import make_client
from .projection import build_coverage_matrix
from .render import write_data_outputs, write_html, write_image
from .orthologs import build_ortholog_msa
from .secondary_structure import fetch_secondary_structure
from .structures import fetch_coverage
from .uniprot import AccMeta, batch_lookup

log = logging.getLogger(__name__)

_ALL_FORMATS = {"html", "image", "data"}


def parse_formats(spec: str) -> Set[str]:
    if not spec or spec.strip() == "all":
        return set(_ALL_FORMATS)
    out: Set[str] = set()
    for tok in spec.split(","):
        t = tok.strip().lower()
        if t in ("png", "svg", "image"):
            out.add("image")
        elif t in ("csv", "json", "data"):
            out.add("data")
        elif t == "html":
            out.add("html")
        elif t == "all":
            out |= _ALL_FORMATS
    return out or set(_ALL_FORMATS)


def build_msa(query: Query, a3m_text: str) -> MSA:
    records = parse_a3m(a3m_text)
    if not records:
        raise ValueError("Empty MSA returned by backend")

    # Record 0 is the query (header '101'); it defines the match columns.
    q_seq = records[0].seq
    rows: List[MSARow] = [
        MSARow(
            row_id=query.name,
            raw_header=records[0].header,
            a3m_seq=q_seq,
            match_seq=match_columns(q_seq),
            accession=query.accession,
            organism=query.organism,
            gene=query.gene,
            is_query=True,
        )
    ]
    seen: Set[str] = {rows[0].match_seq}
    for rec in records[1:]:
        match = match_columns(rec.seq)
        if match in seen:  # drop exact-duplicate alignments
            continue
        seen.add(match)
        rows.append(
            MSARow(
                row_id=rec.first_token,
                raw_header=rec.header,
                a3m_seq=rec.seq,
                match_seq=match,
                accession=extract_accession(rec.header),
            )
        )
    log.info("MSA: %d unique rows (%d columns)", len(rows), query.length)
    return MSA(query=query, rows=rows)


def unique_accessions(msa: MSA, limit: int) -> List[str]:
    ordered: List[str] = []
    seen: Set[str] = set()
    for row in msa.rows:
        acc = row.accession
        if acc and acc not in seen:
            seen.add(acc)
            ordered.append(acc)
    return ordered[:limit]


def annotate_rows(msa: MSA, meta: Dict[str, AccMeta]) -> None:
    """Attach organism / gene / review status from UniProt metadata to each row."""
    for row in msa.rows:
        if row.is_query:
            continue
        m = meta.get(row.accession or "")
        if m and m.found:
            row.organism = m.organism or row.organism
            row.gene = m.gene or row.gene
            row.reviewed = m.reviewed


def filter_reviewed(msa: MSA, meta: Dict[str, AccMeta]) -> MSA:
    """Keep the query plus reviewed (Swiss-Prot) hits; drop TrEMBL/unknown noise."""
    kept = []
    for row in msa.rows:
        if row.is_query:
            kept.append(row)
            continue
        m = meta.get(row.accession or "")
        if m and m.reviewed:
            kept.append(row)
    log.info(
        "reviewed-only: kept %d of %d rows (dropped %d unreviewed/unknown)",
        len(kept), len(msa.rows), len(msa.rows) - len(kept),
    )
    return MSA(query=msa.query, rows=kept)


def collect_structures(
    client, cache, config: Config, msa: MSA, meta: Dict[str, AccMeta]
) -> Dict[str, StructureCoverage]:
    """Fetch observed-residue depth, only for the query and accessions known (from
    UniProt) to carry PDB cross-references — avoiding blind PDBe calls."""
    candidates: List[str] = []
    seen: Set[str] = set()
    if msa.query.accession:
        candidates.append(msa.query.accession)
        seen.add(msa.query.accession)
    for row in msa.rows:
        acc = row.accession
        if not acc or acc in seen:
            continue
        m = meta.get(acc)
        if m and m.has_pdb:
            seen.add(acc)
            candidates.append(acc)

    capped = candidates[: config.max_structured_rows + (1 if msa.query.accession else 0)]
    if len(candidates) > len(capped):
        log.info("Capping observed-residue lookups at %d (of %d with PDB)",
                 len(capped), len(candidates))

    coverage: Dict[str, StructureCoverage] = {}
    for acc in tqdm(capped, desc="PDB observed residues", unit="acc"):
        coverage[acc] = fetch_coverage(client, cache, config, acc)
    n_struct = sum(1 for c in coverage.values() if c.has_structures)
    log.info("%d/%d accessions returned observed structures", n_struct, len(capped))
    return coverage


def run_pipeline(
    config: Config,
    out_dir: Path,
    formats: Set[str],
    *,
    uniprot: Optional[str] = None,
    gene: Optional[str] = None,
    organism: Optional[str] = None,
    sequence: Optional[str] = None,
    name: Optional[str] = None,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = DiskCache(config.cache_dir, enabled=config.use_cache)

    with make_client(config) as client:
        log.info("Resolving query...")
        query = resolve_query(
            client, cache, config,
            uniprot=uniprot, gene=gene, organism=organism, sequence=sequence, name=name,
        )
        log.info("Query: %s (%d aa)%s", query.name, query.length,
                 f" — {query.organism}" if query.organism else "")

        if config.source == "orthologs":
            log.info("Building curated-ortholog alignment from UniProt...")
            msa, meta, a3m_text = build_ortholog_msa(client, cache, config, query)
            source_label = "uniprot-orthologs"
        else:
            log.info("Building homology MSA via ColabFold MMseqs2...")
            backend = ColabFoldMMseqs2Backend(client, cache, config)
            a3m_text = backend.run(query.sequence)
            msa = build_msa(query, a3m_text)
            accs = unique_accessions(msa, config.max_lookup_accessions)
            log.info("Looking up UniProt metadata for %d accessions...", len(accs))
            meta = batch_lookup(client, cache, config, accs)
            annotate_rows(msa, meta)
            if config.reviewed_only:
                msa = filter_reviewed(msa, meta)
            source_label = backend.name

        coverage = collect_structures(client, cache, config, msa, meta)

        # Query secondary structure (helices/strands, in UniProt coords) for the
        # topology cartoon drawn on top of the alignment.
        query_ss = None
        if config.show_secondary_structure and query.accession:
            log.info("Fetching query secondary structure for %s ...", query.accession)
            query_ss = fetch_secondary_structure(client, cache, config, query.accession)

        matrix = build_coverage_matrix(msa, coverage, query_ss=query_ss)

    params = {
        "source": source_label,
        "ortholog_aligner": config.ortholog_aligner if config.source == "orthologs" else None,
        "colabfold_mode": config.colabfold_mode if config.source == "msa" else None,
        "reviewed_only": config.reviewed_only,
        "max_structured_rows": config.max_structured_rows,
        "n_msa_rows": len(msa.rows),
        "secondary_structure": matrix.has_secondary_structure,
    }

    outputs: List[Path] = []
    if "data" in formats:
        outputs += write_data_outputs(out_dir, msa, matrix, query, a3m_text, params)
    if "image" in formats:
        outputs += write_image(out_dir, matrix, query, config)
    if "html" in formats:
        outputs += write_html(out_dir, matrix, msa, query, config)

    return {"query": query, "msa": msa, "matrix": matrix, "outputs": outputs}
