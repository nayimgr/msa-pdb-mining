from msa_pdb_mining.msa.a3m import match_columns, ungapped
from msa_pdb_mining.orthologs import (
    a3m_from_aligned,
    a3m_from_pairwise,
    build_group_query,
    extract_ortholog_groups,
    parse_orthologs,
)

QUERY = "ACDEFGHIK"


def test_a3m_from_pairwise_handles_insertion_and_deletion():
    # query ACDEFGHIK aligned to ortholog ACDQQEFGHK (QQ inserted, I deleted)
    q_aln = "ACD--EFGHIK"
    o_aln = "ACDQQEFGH-K"
    o_seq = "ACDQQEFGHK"
    a3m = a3m_from_pairwise(q_aln, o_aln, o_seq)
    assert a3m == "ACDqqEFGH-K"
    # query-anchored projection has exactly len(query) columns
    assert len(match_columns(a3m)) == len(QUERY)
    # ungapped row equals the full ortholog sequence -> offset detection stays exact
    assert ungapped(a3m) == o_seq


def test_a3m_from_aligned_handles_insertion_and_deletion():
    # one row of a true MSA: query gapped where the ortholog inserts QQ; ortholog
    # gapped (deletion) where the query has I.
    q_aln = "ACD--EFGHIK"
    o_aln = "ACDQQEFGH-K"
    o_seq = "ACDQQEFGHK"
    a3m = a3m_from_aligned(q_aln, o_aln, o_seq)
    # same query-anchored A3M the pairwise engine would yield for this block
    assert a3m == "ACDqqEFGH-K"
    assert len(match_columns(a3m)) == len(QUERY)
    assert ungapped(a3m) == o_seq


def test_a3m_from_aligned_takes_residues_from_original_not_aligner():
    # the aligner may normalise a non-standard residue (U -> X); the gap pattern is
    # what we trust, but each residue must come from the original sequence so the
    # row's ungapped form still matches the UniProt sequence.
    q_aln = "ACDEF"
    o_aln = "ACXEF"  # aligner emitted X at the non-standard position
    o_seq = "ACUEF"  # real residue is selenocysteine U
    a3m = a3m_from_aligned(q_aln, o_aln, o_seq)
    assert a3m == "ACUEF"
    assert ungapped(a3m) == o_seq


def test_a3m_from_aligned_drops_columns_where_both_are_gapped():
    # a multiple alignment can leave columns where neither this row nor the query
    # has a residue (another family member's insertion); they contribute nothing.
    q_aln = "AC--DEF"
    o_aln = "AG-WDEF"
    o_seq = "AGWDEF"
    a3m = a3m_from_aligned(q_aln, o_aln, o_seq)
    # col0 A/A match -> 'A'; col1 C/G match -> 'G'; col2 -/- both gapped -> dropped;
    # col3 -/W query-gap insertion -> 'w'; cols D/E/F match -> 'DEF'
    assert a3m == "AGwDEF"
    assert len(match_columns(a3m)) == len("ACDEF")
    assert ungapped(a3m) == o_seq


def test_a3m_uses_original_residues_not_sanitised():
    # ortholog has a non-standard residue 'U'; scoring may sanitise to X, but the
    # row must keep the real residue so it matches the unipdb full sequence.
    q_aln = "ACDEF"
    o_aln = "ACUEF"
    o_seq = "ACUEF"
    a3m = a3m_from_pairwise(q_aln, o_aln, o_seq)
    assert a3m == "ACUEF"
    assert ungapped(a3m) == "ACUEF"


TSV = (
    "Entry\tReviewed\tOrganism\tGene Names (primary)\tSequence\tPDB\n"
    "P10361\treviewed\tRattus norvegicus (Rat)\tTp53\tMEDSQSD\t\n"
    "G0RYI5\tunreviewed\tThermochaetoides thermophila\t\tMEDSQTD\t5FM6;\n"
    "P04637\treviewed\tHomo sapiens (Human)\tTP53\tMEEPQSD\t1TUP;\n"
)


def test_parse_orthologs_excludes_query_and_keeps_metadata():
    orths = parse_orthologs(TSV, exclude="P04637")
    assert len(orths) == 2
    o = next(o for o in orths if o.accession == "P10361")
    assert o.organism == "Rattus norvegicus"
    assert o.gene == "Tp53"
    assert o.sequence == "MEDSQSD"
    assert o.pdb_ids == set()
    assert o.reviewed is True


def test_parse_orthologs_captures_unreviewed_with_pdb():
    # an unreviewed (TrEMBL) ortholog that has a PDB structure is kept and flagged
    # unreviewed (e.g. the Chaetomium RUVBL1 of 5FM6).
    orths = parse_orthologs(TSV, exclude="P04637")
    o = next(o for o in orths if o.accession == "G0RYI5")
    assert o.reviewed is False
    assert o.pdb_ids == {"5fm6"}
    assert o.gene is None


# A trimmed UniProt entry JSON (CTC1, Q2NKJ3) carrying the orthology cross-refs
# plus an unrelated one we must ignore.
ENTRY = {
    "uniProtKBCrossReferences": [
        {"database": "OrthoDB", "id": "2314520at2759"},
        {"database": "eggNOG", "id": "ENOG502RBD3"},
        {"database": "PANTHER", "id": "PTHR14865"},
        {"database": "PANTHER", "id": "PTHR14865:SF2"},  # subfamily -> family id
        {"database": "GeneTree", "id": "ENSGT00390000011553"},
        {"database": "OMA", "id": "HTDYTPT"},
        {"database": "PDB", "id": "6XYZ"},  # not an orthology db -> ignored
    ]
}
DBS = ("OrthoDB", "eggNOG", "PANTHER", "GeneTree", "OMA")


def test_extract_ortholog_groups_maps_tokens_and_strips_panther_subfamily():
    groups = extract_ortholog_groups(ENTRY, DBS)
    # PDB ignored; PANTHER family + subfamily collapse to one PTHR14865; tokens
    # lower-cased; order preserved.
    assert groups == [
        ("orthodb", "2314520at2759"),
        ("eggnog", "ENOG502RBD3"),
        ("panther", "PTHR14865"),
        ("genetree", "ENSGT00390000011553"),
        ("oma", "HTDYTPT"),
    ]


def test_extract_ortholog_groups_respects_db_selection():
    groups = extract_ortholog_groups(ENTRY, ("OrthoDB",))
    assert groups == [("orthodb", "2314520at2759")]


def test_extract_ortholog_groups_empty_when_no_orthology_xrefs():
    assert extract_ortholog_groups({"uniProtKBCrossReferences": []}, DBS) == []
    assert extract_ortholog_groups({}, DBS) == []


def test_build_group_query_ors_clauses_and_keeps_reviewed_or_structured():
    q = build_group_query([("orthodb", "2314520at2759"), ("panther", "PTHR14865")])
    assert q == (
        "(xref:orthodb-2314520at2759 OR xref:panther-PTHR14865) "
        "AND (reviewed:true OR database:pdb)"
    )
