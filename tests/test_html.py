"""Offline checks for the wrapped HTML coverage viewer.

The visual wrapping itself happens in the browser, so these tests exercise what
can be verified without one: the embedded JSON payload is well-formed and
correctly escaped, has the right shape, and the page carries the responsive
render machinery.
"""

import json
import re

from msa_pdb_mining.config import Config
from msa_pdb_mining.model import MSA, MSARow, Query, StructureCoverage
from msa_pdb_mining.msa.a3m import match_columns
from msa_pdb_mining.projection import build_coverage_matrix
from msa_pdb_mining.render.html import write_html

QUERY_SEQ = "ACDEFGHIK"


def _matrix_and_msa():
    query = Query(sequence=QUERY_SEQ, accession="QUERY", name="q")
    rows = [
        MSARow("q", "101", QUERY_SEQ, QUERY_SEQ, accession="QUERY", is_query=True),
        MSARow("hit", "HIT", QUERY_SEQ, match_columns(QUERY_SEQ), accession="HIT"),
    ]
    msa = MSA(query=query, rows=rows)
    cov = {
        "QUERY": StructureCoverage("QUERY", QUERY_SEQ, {r: {"1abc", "2xyz"} for r in range(1, 6)}),
        "HIT": StructureCoverage("HIT", QUERY_SEQ, {3: {"3hit"}, 4: {"3hit"}}),
    }
    return build_coverage_matrix(msa, cov), msa, query


def _payload(html: str) -> dict:
    m = re.search(
        r'<script id="coverage-data" type="application/json">(.*?)</script>', html, re.S
    )
    assert m, "embedded JSON payload not found"
    raw = m.group(1)
    # The '<' escaping is what prevents a stray '</script>' or '<' in tooltips
    # from terminating the data block early.
    assert "<" not in raw, "payload contains an unescaped '<'"
    return json.loads(raw)


def test_payload_shape(tmp_path):
    matrix, msa, query = _matrix_and_msa()
    html = write_html(tmp_path, matrix, msa, query, Config())[0].read_text()
    data = _payload(html)

    assert data["columns"] == len(QUERY_SEQ)
    assert data["querySeq"] == QUERY_SEQ
    assert len(data["rows"]) == 2
    assert data["rows"][0]["isQuery"] is True
    assert data["rows"][1]["isQuery"] is False
    for row in data["rows"]:
        assert len(row["chars"]) == data["columns"]
        assert len(row["depths"]) == data["columns"]
        assert len(row["ids"]) == data["columns"]
    # Aggregate = distinct PDB entries across all rows per column. Cols 0-4 are
    # covered by QUERY (ids 1abc,2xyz); cols 2,3 additionally by HIT (3hit), so
    # col 2 sees 3 distinct entries. Col 5 (residue 6) is covered by neither.
    assert data["aggregate"]["depths"][2] == 3
    assert data["aggregate"]["depths"][5] == 0


def test_ids_interned_and_referenced(tmp_path):
    matrix, msa, query = _matrix_and_msa()
    html = write_html(tmp_path, matrix, msa, query, Config())[0].read_text()
    data = _payload(html)

    pool = data["idPool"]
    assert "1abc, 2xyz" in pool  # sorted, joined
    # Every non-negative id index in a covered cell must resolve into the pool.
    for row in data["rows"]:
        for idx in row["ids"]:
            assert idx == -1 or 0 <= idx < len(pool)
    # Gap/uncovered cells carry no id reference.
    for row in data["rows"]:
        for ch, d, idx in zip(row["chars"], row["depths"], row["ids"]):
            if ch == "-" or d == 0:
                assert idx == -1


def test_page_has_responsive_machinery(tmp_path):
    matrix, msa, query = _matrix_and_msa()
    html = write_html(tmp_path, matrix, msa, query, Config())[0].read_text()
    for token in (
        'id="viewer"',
        "function render",
        "--labelw",
        "--cellw",
        "addEventListener('resize'",
        'all orthologs (union)',
    ):
        assert token in html, token
