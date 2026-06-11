from msa_pdb_mining.structures import build_coverage

RAW = {
    "P0TEST": {
        "sequence": "M" * 300,
        "length": 300,
        "dataType": "UNIPDB",
        "data": [
            {
                "name": "1abc",
                "accession": "1abc",
                "additionalData": {"experiment": "X-ray diffraction", "resolution": 2.0},
                "residues": [
                    {"startIndex": 93, "endIndex": 94, "indexType": "UNIPROT", "observed": "N"},
                    {"startIndex": 95, "endIndex": 291, "indexType": "UNIPROT", "observed": "Y"},
                ],
            },
            {
                "name": "2xyz",
                "accession": "2xyz",
                "additionalData": {"experiment": "Electron Microscopy", "resolution": 3.2},
                "residues": [
                    {"startIndex": 100, "endIndex": 150, "indexType": "UNIPROT", "observed": "Y"},
                ],
            },
        ],
    }
}


def test_build_coverage_depth_counts_distinct_structures():
    cov = build_coverage("P0TEST", RAW)
    assert cov.has_structures
    assert cov.n_structures == 2
    assert cov.full_sequence == "M" * 300
    # residue 100 is observed by both structures
    assert cov.depth(100) == 2
    assert cov.pdb_ids(100) == {"1abc", "2xyz"}
    # residue 95 only by 1abc; 291 edge of range
    assert cov.depth(95) == 1
    assert cov.depth(291) == 1


def test_unobserved_and_outside_ranges_are_zero():
    cov = build_coverage("P0TEST", RAW)
    assert cov.depth(93) == 0  # observed == "N" (disordered)
    assert cov.depth(94) == 0
    assert cov.depth(292) == 0  # outside any range
    assert cov.depth(1) == 0


def test_methods_recorded():
    cov = build_coverage("P0TEST", RAW)
    assert cov.method["1abc"] == "X-ray diffraction"
    assert cov.resolution["2xyz"] == 3.2


def test_empty_response():
    cov = build_coverage("P0NONE", None)
    assert not cov.has_structures
    assert cov.n_structures == 0
