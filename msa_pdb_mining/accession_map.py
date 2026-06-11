"""Resolve MSA hit headers to UniProt accessions.

ColabFold's ``uniref.a3m`` headers carry the target identifier as the first
token, typically ``UniRef100_<accession>`` (occasionally a bare accession or an
``sp|``/``tr|`` style id), sometimes with a ``/start-end`` range suffix. We strip
those decorations and validate against the UniProt accession grammar. Identifiers
that don't resolve (e.g. UniParc ``UPI...`` representatives, environmental hits)
are returned as ``None`` and carried through as "unknown" rows.
"""

from __future__ import annotations

import re
from typing import Optional

# Official UniProtKB accession pattern.
_UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
_UNIREF_PREFIX = re.compile(r"^UniRef(?:100|90|50)_")


def extract_accession(header: str) -> Optional[str]:
    """Best-effort UniProt accession from an A3M/FASTA header. ``None`` if none."""
    token = header.split()[0] if header.split() else header

    # sp|ACC|NAME or tr|ACC|NAME
    if token[:3] in ("sp|", "tr|") and token.count("|") >= 2:
        token = token.split("|")[1]
    else:
        token = _UNIREF_PREFIX.sub("", token)

    # drop a trailing /start-end range
    token = token.split("/")[0]
    # drop an isoform suffix (P04637-2 -> P04637); unipdb keys on the canonical
    candidate = token.split("-")[0]

    return candidate if _UNIPROT_RE.match(candidate) else None
