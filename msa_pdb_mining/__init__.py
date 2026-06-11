"""msa-pdb-mining: overview structural data availability across orthologs.

Given a UniProt accession / gene name (or a sequence), build a multiple sequence
alignment of homologs, find every experimental PDB structure for those sequences,
determine which residues are actually *modelled* (observed, not just present in
SEQRES), and render the MSA coloured by structural-coverage depth.
"""

__version__ = "0.1.0"
