from msa_pdb_mining.model import MSA, MSARow, Query, StructureCoverage
from msa_pdb_mining.projection import build_coverage_matrix
from msa_pdb_mining.secondary_structure import build_secondary_structure

# Two helix blocks and one strand block, plus a conflicting assignment on
# residue 12 (helix in two depositions, strand in one -> consensus helix) and a
# non-UNIPROT range that must be ignored.
RAW = {
    "P0TEST": {
        "sequence": "M" * 60,
        "length": 60,
        "dataType": "SECONDARY STRUCTURES",
        "data": [
            {"name": "Helix", "residues": [
                {"startIndex": 10, "endIndex": 12, "indexType": "UNIPROT"},
            ]},
            {"name": "Helix", "residues": [
                {"startIndex": 11, "endIndex": 12, "indexType": "UNIPROT"},
            ]},
            {"name": "Strand", "residues": [
                {"startIndex": 12, "endIndex": 12, "indexType": "UNIPROT"},
                {"startIndex": 20, "endIndex": 23, "indexType": "UNIPROT"},
                {"startIndex": 30, "endIndex": 31, "indexType": "PDB"},  # ignored
            ]},
        ],
    }
}


def test_build_secondary_structure_consensus_and_ranges():
    ss = build_secondary_structure("P0TEST", RAW)
    assert ss.has_elements
    assert ss.full_sequence == "M" * 60
    # residues 10,11 helix only; 12 helix (2 votes) beats strand (1 vote)
    assert ss.element(10) == "H"
    assert ss.element(11) == "H"
    assert ss.element(12) == "H"
    # strand range 20-23
    assert [ss.element(r) for r in range(20, 24)] == ["E", "E", "E", "E"]
    # non-UNIPROT range ignored; coil residues absent
    assert ss.element(30) is None
    assert ss.element(5) is None


def test_empty_secondary_structure():
    ss = build_secondary_structure("P0NONE", None)
    assert not ss.has_elements
    assert build_secondary_structure("P0NONE", {}).has_elements is False


def _make_msa(query_seq: str):
    query = Query(sequence=query_seq, accession="P0TEST", name="q")
    rows = [MSARow("q", "101", query_seq, query_seq, accession="P0TEST", is_query=True)]
    return MSA(query=query, rows=rows)


def test_projection_places_ss_on_query_columns():
    seq = "M" * 60
    ss = build_secondary_structure("P0TEST", RAW)
    cov = {"P0TEST": StructureCoverage("P0TEST", seq)}
    matrix = build_coverage_matrix(_make_msa(seq), cov, query_ss=ss)

    assert matrix.has_secondary_structure
    track = matrix.query_ss
    assert len(track) == 60
    # residue r (1-based) -> column r-1
    assert track[9] == "H" and track[11] == "H"   # residues 10, 12
    assert track[19] == "E" and track[22] == "E"  # residues 20, 23
    assert track[0] is None                       # residue 1: coil


def test_projection_without_ss_leaves_track_none():
    seq = "M" * 60
    cov = {"P0TEST": StructureCoverage("P0TEST", seq)}
    matrix = build_coverage_matrix(_make_msa(seq), cov, query_ss=None)
    assert matrix.query_ss is None
    assert matrix.has_secondary_structure is False
