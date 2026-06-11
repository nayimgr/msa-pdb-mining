"""Resolve user input (UniProt accession / gene+organism / raw sequence) to a Query."""

from __future__ import annotations

import logging
import re
from typing import Optional, Tuple
from urllib.parse import quote

import httpx

from .cache import DiskCache
from .config import Config
from .model import Query
from .net import request

log = logging.getLogger(__name__)

_HDR_FIELD = {
    "OS": "organism",
    "OX": "tax_id",
    "GN": "gene",
}
# Tokens that introduce header key=value fields in a UniProt FASTA description.
_FIELD_RE = re.compile(r"\b(OS|OX|GN|PE|SV)=")


def _get_text(
    client: httpx.Client, cache: DiskCache, namespace: str, key: str, url: str
) -> Optional[str]:
    hit = cache.get_text(namespace, key)
    if hit is not None:
        return hit
    resp = request(client, "GET", url, headers={"Accept": "text/plain"})
    if resp is None or resp.status_code >= 400 or not resp.text.strip():
        return None
    cache.set_text(namespace, key, resp.text)
    return resp.text


def parse_uniprot_fasta(text: str) -> Tuple[Optional[str], str, dict]:
    """Return (accession, sequence, fields) from a single-entry UniProt FASTA."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith(">"):
        return None, "".join(lines), {}
    header = lines[0][1:]
    seq = "".join(lines[1:])

    accession: Optional[str] = None
    if header[:3] in ("sp|", "tr|") and header.count("|") >= 2:
        accession = header.split("|")[1]

    fields: dict = {}
    matches = list(_FIELD_RE.finditer(header))
    for i, m in enumerate(matches):
        key = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(header)
        value = header[start:end].strip()
        if key in _HDR_FIELD:
            fields[_HDR_FIELD[key]] = value
    return accession, seq, fields


def _query_from_fasta(text: str, fallback_name: str) -> Query:
    accession, seq, fields = parse_uniprot_fasta(text)
    tax_id = None
    if fields.get("tax_id", "").isdigit():
        tax_id = int(fields["tax_id"])
    return Query(
        sequence=seq.upper(),
        accession=accession,
        gene=fields.get("gene"),
        organism=fields.get("organism"),
        tax_id=tax_id,
        name=accession or fallback_name,
    )


def resolve_by_accession(
    client: httpx.Client, cache: DiskCache, config: Config, accession: str
) -> Query:
    url = f"{config.uniprot_rest}/uniprotkb/{accession}.fasta"
    text = _get_text(client, cache, "uniprot_fasta", accession, url)
    if not text:
        raise ValueError(f"Could not fetch UniProt entry for accession '{accession}'")
    return _query_from_fasta(text, accession)


def resolve_by_gene(
    client: httpx.Client,
    cache: DiskCache,
    config: Config,
    gene: str,
    organism: Optional[str],
) -> Query:
    if organism is None:
        log.warning("No --organism given for gene '%s'; defaulting to human (9606)", gene)
        organism = "9606"

    org_clause = (
        f"organism_id:{organism}"
        if str(organism).isdigit()
        else f'organism_name:"{organism}"'
    )
    base = f"{config.uniprot_rest}/uniprotkb/search"
    for reviewed in (" AND reviewed:true", ""):
        query = f"gene:{gene} AND {org_clause}{reviewed}"
        url = f"{base}?query={quote(query)}&format=fasta&size=1"
        key = f"{gene}__{organism}__{'rev' if reviewed else 'all'}"
        text = _get_text(client, cache, "uniprot_gene", key, url)
        if text and text.lstrip().startswith(">"):
            q = _query_from_fasta(text, f"{gene}_{organism}")
            if not q.gene:
                q.gene = gene
            return q
    raise ValueError(f"No UniProt entry found for gene '{gene}' in organism '{organism}'")


def resolve_query(
    client: httpx.Client,
    cache: DiskCache,
    config: Config,
    *,
    uniprot: Optional[str] = None,
    gene: Optional[str] = None,
    organism: Optional[str] = None,
    sequence: Optional[str] = None,
    name: Optional[str] = None,
) -> Query:
    if uniprot:
        return resolve_by_accession(client, cache, config, uniprot.strip())
    if gene:
        return resolve_by_gene(client, cache, config, gene.strip(), organism)
    if sequence:
        clean = "".join(sequence.split()).upper()
        return Query(sequence=clean, name=name or "query")
    raise ValueError("Provide one of: --uniprot, --gene, or --sequence")
