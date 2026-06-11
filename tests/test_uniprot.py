from msa_pdb_mining.uniprot import AccMeta, parse_tsv

TSV = (
    "Entry\tReviewed\tOrganism\tGene Names (primary)\tPDB\n"
    "P04637\treviewed\tHomo sapiens (Human)\tTP53\t1TUP;1TSR;\n"
    "A0A2K5XT84\tunreviewed\tMandrillus leucophaeus (Drill) (Papio leucophaeus)\t\t\n"
    "Q00987\treviewed\tHomo sapiens (Human)\tMDM2\t1RV1;1T4E\n"
)


def test_parse_tsv_reviewed_with_structures():
    d = parse_tsv(TSV)
    p53 = d["P04637"]
    assert p53.reviewed is True
    assert p53.organism == "Homo sapiens"  # common name stripped
    assert p53.gene == "TP53"
    assert p53.pdb_ids == {"1tup", "1tsr"}  # lowercased, trailing ';' ignored
    assert p53.has_pdb


def test_parse_tsv_unreviewed_no_structures():
    m = parse_tsv(TSV)["A0A2K5XT84"]
    assert m.reviewed is False
    assert m.organism == "Mandrillus leucophaeus"
    assert m.gene is None
    assert not m.has_pdb


def test_accmeta_roundtrip():
    m = AccMeta(reviewed=True, organism="Homo sapiens", gene="TP53", pdb_ids={"1tup"}, found=True)
    assert AccMeta.from_dict(m.to_dict()) == m
