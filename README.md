# msa-pdb-mining

Get a fast overview of **what experimental structural data already exists for specific
regions of a protein, across all its orthologs**.

Given a UniProt accession (or gene name, or raw sequence), the tool:

1. builds a multiple sequence alignment of homologs with the **ColabFold MMseqs2** server;
2. finds every **experimental PDB** structure for each aligned sequence (PDBe SIFTS);
3. determines which residues are actually **modelled** — observed in the density, not just
   present in SEQRES (disordered/missing residues don't count);
4. renders the alignment **coloured by structural-coverage depth** — how many structures
   model each residue — as an interactive HTML viewer, a static figure, and data tables.

Glance at the result and immediately see which parts of your protein are structurally
characterised (in itself or any ortholog) and which are blind spots.

## Why the coordinates are trustworthy

All structural information is consumed in **UniProt coordinates** (from SIFTS), never PDB
author numbering — so inconsistent residue numbering across depositions is irrelevant, and
the **MSA column is the absolute, cross-ortholog index**. Affinity tags, fusion partners and
cloning artefacts have no UniProt mapping and therefore can never contribute to coverage.

## Install

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

## Usage

```bash
# by UniProt accession
msa-pdb-mining --uniprot P04637 --out results/p53

# by gene + organism (NCBI taxid or name)
msa-pdb-mining --gene TP53 --organism 9606 --out results/p53

# by raw sequence
msa-pdb-mining --sequence MEEPQSDPSV... --name myprotein --out results/x
```

Useful flags:

| flag | meaning |
|---|---|
| `--formats html,image,data` | which outputs to write (default: all) |
| `--email you@example.org` | contact sent to EBI/ColabFold (etiquette; recommended) |
| `--reviewed-only` | keep only reviewed (Swiss-Prot) hits — drops uncharacterised TrEMBL noise |
| `--max-structured-rows N` | cap accessions resolved to observed PDB residues (default 200) |
| `--max-lookup-accessions N` | cap accessions fetched for UniProt metadata (batched; default 2000) |
| `--max-rows-render N` | rows drawn in the image/HTML (default 60) |
| `--include-env` | also include environmental/metagenomic hits (no PDB) |
| `--cache-dir DIR` / `--no-cache` | control the on-disk response cache |
| `-v` | verbose progress logging |

## Outputs (in `--out`)

| file | contents |
|---|---|
| `alignment.html` | interactive viewer; cells coloured by depth, hover for PDB ids + method/resolution |
| `alignment.png` / `alignment.svg` | static coverage heatmap (query + aggregate track) |
| `column_summary.csv` | per query residue: #structures, #orthologs covered, PDB ids |
| `coverage_long.csv` | one row per covered (ortholog, residue) cell |
| `summary.json` | run parameters + per-column aggregate coverage vector |
| `msa.a3m`, `msa.query_anchored.fasta` | the raw alignment |

## How it works

```
input ─▶ resolve query (UniProt) ─▶ ColabFold MMseqs2 MSA ─▶ query-anchored alignment
      ─▶ accession per row ─▶ PDBe SIFTS observed residues ─▶ depth per residue
      ─▶ project onto query columns ─▶ HTML / image / data
```

The MSA backend is abstracted (`msa_pdb_mining/msa/base.py`), so a local jackhmmer / MMseqs2
engine can be added later without touching the rest of the pipeline.

## Develop / test

```bash
.venv/bin/pytest -q          # offline unit tests (a3m, structures, projection, accession)
```

## Status

Implemented: `--uniprot` / `--gene` / `--sequence` → ColabFold MSA → experimental-PDB depth →
HTML + image + data. Caching of all API responses.

Planned: species-list ortholog mode (OrthoDB/OMA + alignment); pluggable local MSA backends;
optional AlphaFold/computed-model coverage; sequence-search fallback for non-UniProt hits.
