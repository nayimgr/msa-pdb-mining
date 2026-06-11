from msa_pdb_mining.msa.a3m import match_columns, ungapped
from msa_pdb_mining.orthologs import a3m_from_pairwise, parse_orthologs

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
    "P04637\treviewed\tHomo sapiens (Human)\tTP53\tMEEPQSD\t1TUP;\n"
)


def test_parse_orthologs_excludes_query_and_keeps_metadata():
    orths = parse_orthologs(TSV, exclude="P04637")
    assert len(orths) == 1
    o = orths[0]
    assert o.accession == "P10361"
    assert o.organism == "Rattus norvegicus"
    assert o.gene == "Tp53"
    assert o.sequence == "MEDSQSD"
    assert o.pdb_ids == set()
