# msa-pdb-mining

Get a fast overview of **what experimental structural data already exists for specific
regions of a protein, across all its orthologs**.

Given a UniProt accession (or gene name, or raw sequence), the tool:

1. assembles the relevant sequences — **by default the curated reviewed orthologs** of the
   query, taken from its UniProt **orthology groups** (OrthoDB, eggNOG, PANTHER, GeneTree,
   OMA), combined into a true **multiple sequence alignment** (FAMSA, in-process); or,
   optionally, a broad **ColabFold MMseqs2** homology MSA (`--source msa`);
2. finds every **experimental PDB** structure for each sequence (PDBe SIFTS);
3. determines which residues are actually **modelled** — observed in the density, not just
   present in SEQRES (disordered/missing residues don't count);
4. renders the alignment **coloured by structural-coverage depth** — how many structures
   model each residue (seaborn's *mako* palette) — as an interactive HTML viewer, a static
   figure, and data tables.

When the query itself has experimental structures, its **secondary structure** (helices and
strands, in UniProt coordinates from PDBe) is drawn as a topology cartoon — α-helix cylinders
and β-strand arrows — on top of the alignment, so you can read coverage against fold.

> **Why curated orthologs by default?** A folding-style MMseqs2 MSA searches a *clustered*
> database and deliberately filters out redundant near-identical sequences, so curated
> Swiss-Prot orthologs (mouse, rat, …) are largely absent. For "what structural data exists
> across the orthologs", you want completeness per species — so the default pulls the reviewed
> members of the query's **orthology groups** straight from UniProt. Crucially this is true
> orthology (OrthoDB/eggNOG/PANTHER/…), **not** the free-text "Belongs to the … family"
> annotation: the family text conflates evolutionary orthologs with same-family *interactors*
> and unrelated members — e.g. it files several species' **STN1** under the *CTC1* family, so a
> family search for CTC1 wrongly drags in its CST-complex partner — whereas no orthology group
> does. The members kept are reviewed Swiss-Prot entries **plus** any unreviewed (TrEMBL)
> ortholog that has a PDB structure — so structural-biology workhorses that happen to be
> non-model organisms still count (e.g. the *Thermochaetoides/Chaetomium* RUVBL1 of PDB 5FM6),
> while the thousands of sequence-only TrEMBL entries per group are dropped. Use `--source msa`
> for broad/distant homology instead.

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

# broad homology MSA instead of curated orthologs
msa-pdb-mining --uniprot P04637 --source msa --out results/p53_msa
```

Useful flags:

| flag | meaning |
|---|---|
| `--source orthologs\|msa` | row source: curated UniProt orthologs (default) or ColabFold MMseqs2 |
| `--max-orthologs N` | max orthology-group members to fetch in ortholog mode (default 500) |
| `--ortholog-aligner famsa\|pairwise` | ortholog-mode alignment: one true FAMSA MSA (default) or star (each aligned to the query) |
| `--formats html,image,data` | which outputs to write (default: all) |
| `--email you@example.org` | contact sent to EBI/ColabFold (etiquette; recommended) |
| `--reviewed-only` | keep only reviewed (Swiss-Prot) hits — drops uncharacterised TrEMBL noise |
| `--max-structured-rows N` | cap accessions resolved to observed PDB residues (default 200) |
| `--max-lookup-accessions N` | cap accessions fetched for UniProt metadata (batched; default 2000) |
| `--max-rows-render N` | rows drawn in the image/HTML (default 60) |
| `--no-secondary-structure` | don't draw the query's secondary-structure cartoon |
| `--include-env` | also include environmental/metagenomic hits (no PDB) |
| `--cache-dir DIR` / `--no-cache` | control the on-disk response cache |
| `-v` | verbose progress logging |

## Outputs (in `--out`)

| file | contents |
|---|---|
| `alignment.html` | interactive viewer; cells coloured by depth (mako), query SS cartoon on top, hover for PDB ids + method/resolution |
| `alignment.png` / `alignment.svg` | static coverage heatmap (query SS cartoon + aggregate track) |
| `column_summary.csv` | per query residue: #structures, #orthologs covered, PDB ids |
| `coverage_long.csv` | one row per covered (ortholog, residue) cell |
| `summary.json` | run parameters + per-column aggregate coverage vector |
| `msa.a3m`, `msa.query_anchored.fasta` | the raw alignment |

## How it works

```
                ┌─ default: UniProt orthology-group orthologs ─▶ FAMSA multiple alignment
input ─▶ query ─┤
                └─ --source msa: ColabFold MMseqs2 homology search ─▶ query-anchored a3m
      ─▶ UniProt metadata (organism/gene/PDB existence, batched)
      ─▶ PDBe SIFTS observed residues ─▶ depth per residue
      ─▶ project onto query columns ─▶ HTML / image / data
```

Both sources produce the same query-anchored `MSA`, so the structure-mapping and rendering
stages are shared. The MMseqs2 path is abstracted behind `msa_pdb_mining/msa/base.py` for
future local engines; the ortholog path lives in `msa_pdb_mining/orthologs.py`.

## Develop / test

```bash
.venv/bin/pytest -q          # offline unit tests (a3m, structures, projection, accession)
```

## Status

Implemented: `--uniprot` / `--gene` / `--sequence`; **curated UniProt ortholog source
(default)** — orthologs taken from the query's OrthoDB/eggNOG/PANTHER/GeneTree/OMA orthology
groups (true orthology, not the SIMILARITY family text) — and ColabFold MMseqs2 source
(`--source msa`); experimental-PDB observed-residue depth; HTML + image + data outputs;
caching of all API responses.

Planned: explicit species lists; pluggable local MSA backends; optional AlphaFold/computed-model
coverage; sequence-search fallback for hits without a UniProt accession.
