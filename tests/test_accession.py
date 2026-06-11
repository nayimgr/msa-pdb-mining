import pytest

from msa_pdb_mining.accession_map import extract_accession


@pytest.mark.parametrize(
    "header,expected",
    [
        ("UniRef100_P04637", "P04637"),
        ("UniRef100_A0A2K5XT84/12-340", "A0A2K5XT84"),
        ("UniRef90_P04637 some description", "P04637"),
        ("sp|P04637|P53_HUMAN Cellular tumor antigen", "P04637"),
        ("tr|notanacc|FOO_XENLA", None),  # malformed accession
        ("P04637-2", "P04637"),  # isoform -> canonical
        ("UniRef100_UPI000123ABCD", None),  # UniParc representative
        ("randomjunk", None),
    ],
)
def test_extract_accession(header, expected):
    assert extract_accession(header) == expected
