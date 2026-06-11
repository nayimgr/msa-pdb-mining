from msa_pdb_mining.model import (
    STATUS_QUERY,
    STATUS_STRUCTURED,
    MSA,
    MSARow,
    Query,
    StructureCoverage,
)
from msa_pdb_mining.msa.a3m import match_columns
from msa_pdb_mining.projection import build_coverage_matrix

QUERY_SEQ = "ACDEFGHIK"
HIT_A3M = "ACDefEFG-IK"  # ungapped: ACDEFEFGIK


def _make_msa():
    query = Query(sequence=QUERY_SEQ, accession="QUERY", name="q")
    rows = [
        MSARow("q", "101", QUERY_SEQ, QUERY_SEQ, accession="QUERY", is_query=True),
        MSARow("hit", "HIT", HIT_A3M, match_columns(HIT_A3M), accession="HIT"),
    ]
    return MSA(query=query, rows=rows)


def _coverage():
    query_cov = StructureCoverage(
        accession="QUERY",
        full_sequence=QUERY_SEQ,
        pdb_ids_per_res={r: {"q1"} for r in range(1, 6)},  # residues 1-5
    )
    hit_cov = StructureCoverage(
        accession="HIT",
        full_sequence="ACDEFEFGIK",
        pdb_ids_per_res={6: {"h1"}, 7: {"h1"}},  # hit residues 6,7 == query cols 3,4
    )
    return {"QUERY": query_cov, "HIT": hit_cov}


def test_projection_places_depth_on_correct_columns():
    matrix = build_coverage_matrix(_make_msa(), _coverage())

    # query row: residues 1-5 covered -> columns 0-4 depth 1, rest 0
    assert list(matrix.depth[0, :5]) == [1, 1, 1, 1, 1]
    assert list(matrix.depth[0, 5:]) == [0, 0, 0, 0]

    # hit row: hit residues 6,7 land on query columns 3,4
    assert matrix.depth[1, 3] == 1
    assert matrix.depth[1, 4] == 1
    assert matrix.depth[1, 0] == 0  # hit residue 1 not in a structure


def test_column_aggregate_unions_pdb_ids():
    matrix = build_coverage_matrix(_make_msa(), _coverage())
    agg = matrix.column_aggregate
    assert agg[0] == 1  # only query's q1
    assert agg[3] == 2  # query q1 + hit h1
    assert agg[4] == 2
    assert agg[5] == 0  # neither covers query position 6


def test_row_status():
    matrix = build_coverage_matrix(_make_msa(), _coverage())
    assert matrix.row_status[0] == STATUS_QUERY
    assert matrix.row_status[1] == STATUS_STRUCTURED
