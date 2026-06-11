from msa_pdb_mining.msa.a3m import (
    find_residue_offset,
    iter_match_residues,
    match_columns,
    parse_a3m,
    ungapped,
)

# query columns:  A C D E F G H I K  (9)
# hit:            ACD ef EF G - IK  -> insertion "ef", deletion at col 6 (G/-)
QUERY = "ACDEFGHIK"
HIT_A3M = "ACDefEFG-IK"


def test_parse_a3m_roundtrip():
    text = f">101\n{QUERY}\n>UniRef100_P12345 desc here\n{HIT_A3M}\n"
    recs = parse_a3m(text)
    assert len(recs) == 2
    assert recs[0].seq == QUERY
    assert recs[1].first_token == "UniRef100_P12345"


def test_match_columns_is_query_length():
    assert match_columns(QUERY) == QUERY
    mc = match_columns(HIT_A3M)
    assert mc == "ACDEFG-IK"
    assert len(mc) == len(QUERY)  # insertions dropped, gaps kept


def test_ungapped_includes_insertions():
    assert ungapped(HIT_A3M) == "ACDEFEFGIK"


def test_iter_match_residues_maps_columns_to_residues():
    pairs = list(iter_match_residues(HIT_A3M, first_residue=1))
    assert pairs == [(0, 1), (1, 2), (2, 3), (3, 6), (4, 7), (5, 8), (7, 9), (8, 10)]
    # column 6 (query position 7) is a deletion in the hit -> no residue mapped
    cols = {c for c, _ in pairs}
    assert 6 not in cols


def test_find_residue_offset():
    full = "XX" + ungapped(HIT_A3M) + "YY"  # hit slice starts at residue 3
    assert find_residue_offset(HIT_A3M, full) == 3
    # apply the offset: first match column now maps to residue 3
    first_col, first_res = next(iter(iter_match_residues(HIT_A3M, 3)))
    assert (first_col, first_res) == (0, 3)


def test_find_residue_offset_missing():
    assert find_residue_offset(HIT_A3M, "TOTALLYDIFFERENT") is None
