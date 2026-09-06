# CLAUDE.md — msa-pdb-mining

Developer/agent onboarding for this repo. The user-facing pitch and CLI usage live in
`README.md`; this file is the internal map — architecture, invariants, conventions, and the
non-obvious decisions that aren't visible in the code. Read this first, skim `README.md` for
the product framing, then dive into modules as needed.

## What this tool does (one paragraph)

Given a protein (UniProt accession, gene+organism, or raw sequence), it builds a multiple
sequence alignment of homologs via the **ColabFold MMseqs2** server, finds every
**experimental PDB** structure for each aligned sequence through **PDBe SIFTS**, determines
which residues are actually **modelled** (observed in the density, not just present in
SEQRES), and renders the alignment **coloured by structural-coverage depth** (how many
structures model each residue) as interactive HTML, a static heatmap, and data tables. The
output answers: *which regions of my protein — in itself or any ortholog — already have
experimental structural data, and which are blind spots?*

## The central invariant: everything is in UniProt coordinates

This is the design decision the whole tool hangs on. All structural data is consumed in
**UniProt residue numbering** (via SIFTS), never PDB author numbering. Consequences:
- Inconsistent residue numbering across PDB depositions is irrelevant.
- The **MSA column is the absolute, cross-ortholog index** — column *N* means the same thing
  for every row.
- Affinity tags, fusion partners, and cloning artefacts have **no** UniProt mapping, so they
  can *never* contribute coverage. This is a correctness feature, not a limitation.

"Depth" at a residue = number of **distinct PDB entries** that observe (model) it. Disordered/
missing residues (`observed: "N"` in SIFTS) do not count.

## Pipeline / data flow

```
input ─▶ resolve query (UniProt) ─▶ ColabFold MMseqs2 MSA ─▶ query-anchored alignment
      ─▶ accession per row ─▶ PDBe SIFTS observed residues ─▶ depth per residue
      ─▶ project onto query columns ─▶ HTML / image / data
```

Orchestrated end-to-end by `pipeline.py:run_pipeline`, which is the best entry point for
understanding execution order:
1. `input_resolver.resolve_query` → a `Query` (sequence + metadata).
2. `msa/colabfold.ColabFoldMMseqs2Backend.run` → an A3M string.
3. `pipeline.build_msa` → an `MSA` (parses A3M, dedups rows, extracts accessions).
4. `pipeline.collect_structures` → `{accession: StructureCoverage}` from PDBe.
5. `secondary_structure.fetch_secondary_structure` → the query's `SecondaryStructure`
   (helices/strands in UniProt coords), only when the query has an accession.
6. `projection.build_coverage_matrix` → a `CoverageMatrix` (depth projected onto query
   columns; the query SS is projected onto the same columns into `matrix.query_ss`).
7. `render/*` → data files, image, HTML.

## Module map (`msa_pdb_mining/`)

| Module | Responsibility |
|---|---|
| `cli.py` | argparse entry point (`msa-pdb-mining` console script → `main`). Builds `Config`, calls `run_pipeline`, prints a summary. |
| `config.py` | `Config` dataclass: all tunable knobs + API endpoints (ColabFold, PDBe REST + graph-api, UniProt REST). Default cache dir via `platformdirs`. |
| `model.py` | Core dataclasses: `Query`, `MSARow`, `MSA`, `StructureCoverage`, `CoverageMatrix`, and the `STATUS_*` row-status constants. **No logic beyond derived properties** — read this to learn the vocabulary. |
| `input_resolver.py` | Resolve `--uniprot` / `--gene`+`--organism` / `--sequence` to a `Query`. Parses UniProt FASTA headers (OS/OX/GN fields). Gene lookup tries `reviewed:true` first, then any. |
| `net.py` | Shared `httpx.Client` factory + `request()` with exponential-backoff retry on `{429,500,502,503,504}`. Returns `None` on give-up so callers distinguish transient failure. |
| `cache.py` | `DiskCache`: tiny namespace+key on-disk cache (JSON and text). `cached_json` only writes non-`None` values — see caching contract below. |
| `accession_map.py` | `extract_accession(header)`: best-effort UniProt accession from an A3M header (strips `UniRef100_`, `sp|`/`tr|`, `/start-end`, isoform `-2`). Returns `None` for UniParc/env hits → carried as "unknown" rows. |
| `structures.py` | `fetch_coverage` / `build_coverage`: PDBe graph-api `uniprot/unipdb/<acc>` → `StructureCoverage`. Only `indexType==UNIPROT` + `observed=="Y"` ranges count. |
| `secondary_structure.py` | `fetch_secondary_structure` / `build_secondary_structure`: PDBe graph-api `uniprot/secondary_structures/<acc>` → `SecondaryStructure`. Per-residue **consensus** element (`H`/`E`) across all depositions (ties → helix); only the query's is fetched. |
| `projection.py` | `build_coverage_matrix`: for each row, locate its hit-subsequence offset in the full UniProt sequence, then map each match column to a residue number and pull depth. Assigns `STATUS_*`. |
| `msa/base.py` | `MSABackend` ABC. The homology search is abstracted so the remote ColabFold engine can later be swapped for a local jackhmmer/MMseqs2 without touching the rest. A backend turns a query sequence into a query-anchored A3M. |
| `msa/a3m.py` | A3M parsing + the column/residue math: `parse_a3m`, `match_columns`, `ungapped`, `iter_match_residues`, `find_residue_offset`. **The trickiest, most important file — see A3M conventions below.** |
| `msa/colabfold.py` | The remote backend: POST to `ticket/msa`, poll `ticket/<id>` to `COMPLETE`, download+untar, extract `uniref.a3m` (optionally merge env member). Result A3M is cached by sequence hash. |
| `render/data.py` | `write_data_outputs`: `msa.a3m`, `msa.query_anchored.fasta`, `column_summary.csv`, `coverage_long.csv`, `summary.json`. |
| `render/image.py` | `write_image`: matplotlib **mako** heatmap (PNG+SVG), aggregate "all orthologs" track, and (when present) the query secondary-structure cartoon on top. Headless (`Agg`). |
| `render/html.py` | `write_html`: self-contained interactive viewer. Embeds a compact JSON payload (`_build_payload`) and re-flows the alignment **in the browser** into width-fitted blocks (wrapped MSA view) — no horizontal scroll; re-wraps on resize/zoom. Cells coloured by depth (mako, luminance-aware text, `.dN` classes in the inline stylesheet); the query secondary structure is drawn as an inline SVG cartoon per block; hover shows PDB ids + position. Jinja2 template inline. |
| `render/layout.py` | Shared row selection/ordering for the visual renderers (query first, then by descending total depth), plus `ss_runs` (collapse the per-column SS track into helix/strand segments). |
| `render/palette.py` | Shared colours: `depth_cmap()` (seaborn mako, lazily imported + cached) and the SS cartoon colours, so image and HTML stay in lockstep. |

## A3M conventions (read before touching `msa/a3m.py` or `projection.py`)

The ColabFold MMseqs2 server returns A3M where:
- **Record 0 is the query** (header `101`), all-uppercase, no gaps — it defines the match columns.
- For other records: **UPPERCASE** = residue aligned to a query column; **lowercase** =
  insertion relative to the query; **`-`** = deletion; **`.`** = insertion padding.

Two projections of an A3M row:
- `match_columns(seq)` — the **query-anchored** view: uppercase + `-` only, so every row has
  exactly `len(query)` columns. Used for display and the FASTA output.
- `ungapped(seq)` — the **contiguous hit subsequence** (upper+lowercase, no gaps). Used to
  locate the row within its full UniProt sequence.

How residue numbers are recovered (this is the clever bit): `find_residue_offset` finds the
hit's `ungapped` subsequence inside the full UniProt sequence (substring search) → 1-based
offset of the first hit residue. Then `iter_match_residues(seq, offset)` walks the A3M
yielding `(column, residue_number)`: deletions advance the column but consume no residue,
insertions consume a residue but produce no column, matches advance both. If the offset
can't be found exactly (e.g. UniProt sequence-version mismatch), the row is marked
`STATUS_UNMAPPED` and contributes no coverage.

## Row status vocabulary (`model.py`)

Each alignment row gets exactly one `STATUS_*` after projection:
- `query` — the anchor row.
- `structured` — accession resolved and has ≥1 PDB entry.
- `no_structures` — accession resolved but 0 PDB entries.
- `unknown` — no UniProt accession could be extracted from the header (UniParc/env hit).
- `unmapped` — accession known and has structures, but the hit subsequence couldn't be
  located in the full UniProt sequence, so residues couldn't be placed.

## External services (all public, no auth)

| Service | Base | Used for |
|---|---|---|
| ColabFold MMseqs2 | `https://api.colabfold.com` | the MSA (ticket submit/poll/download) |
| PDBe graph-api | `https://www.ebi.ac.uk/pdbe/graph-api` | `uniprot/unipdb/<acc>` — observed residue ranges; `uniprot/secondary_structures/<acc>` — query helices/strands |
| UniProt REST | `https://rest.uniprot.org` | FASTA by accession; gene+organism search |

Network etiquette: pass `--email` to set a contact in the User-Agent (recommended, not
required). `request_timeout` default 30s; ColabFold poll interval 5s, max 1800s.

## Caching contract (important for offline work and tests)

`DiskCache` lives at `platformdirs.user_cache_dir("msa-pdb-mining")` by default; override with
`--cache-dir` or disable with `--no-cache`. Namespaces in use: `a3m`, `unipdb`, `unipdb_ss`,
`uniprot_fasta`, `uniprot_gene`, `uniprot_meta`.

The subtle rule (`cache.cached_json` + `structures.fetch_coverage`): a **definitive 404**
(accession has no SIFTS/PDB mapping) is cached as `{}` so it's never re-queried; a **transient
failure** returns `None` and is **not** cached, so it will be retried on the next run. Never
cache `None`. Preserve this distinction if you touch the caching path.

## Build / install / test

```bash
python3 -m venv .venv
.venv/bin/pip install -e .          # editable install; console script: msa-pdb-mining
.venv/bin/pip install -e '.[dev]'   # + pytest
.venv/bin/pytest -q                 # offline unit tests, no network
```

The test suite (`tests/`) is **fully offline** — it exercises pure functions with synthetic
fixtures (`test_a3m`, `test_structures`, `test_secondary_structure`, `test_projection`,
`test_accession`). There is no
integration test that hits the live services. When adding logic, prefer factoring a pure
function (like `build_coverage`, `build_coverage_matrix`, `extract_accession`) that can be
unit-tested without network, mirroring the existing split between `fetch_*` (I/O) and
`build_*` (pure).

## Conventions / house style

- Python ≥3.9, `from __future__ import annotations` at the top of every module.
- Type hints throughout; dataclasses for all data structures.
- Module docstrings explain the *why* and any external-format quirks — match that density.
- Logging via `logging.getLogger(__name__)`; `-v` flips level INFO↔WARNING. No prints except
  the final CLI summary.
- Pure/I-O split: keep network in `fetch_*`/backend methods, keep transformation pure and
  testable.
- numpy for the depth matrix; pandas/matplotlib/seaborn/jinja2 are dependencies but used
  narrowly (seaborn only supplies the mako colormap, via `render/palette.depth_cmap`).

## Outputs written to `--out` (default `results/`)

`alignment.html` (interactive), `alignment.png` + `alignment.svg` (static heatmap),
`column_summary.csv` (per query residue: #structures, #orthologs, PDB ids), `coverage_long.csv`
(one row per covered ortholog×residue cell), `summary.json` (run params + per-column aggregate
vector), `msa.a3m` + `msa.query_anchored.fasta` (raw alignment). The depth heatmap uses the
mako palette; when the query has experimental structures, both visual renderers draw its
secondary structure (helix cylinders / strand arrows) as a track on top — suppress with
`--no-secondary-structure`. Visual renderers cap rows at `--max-rows-render` (default 60);
data files keep **all** rows. `results/` is gitignored; `results/p53`, `results/ruvbl1`,
`results/ruvbl2` are committed sample runs (predate these visual changes — regenerate to refresh).

## Repo facts

- Git repository, hosted at `nayimgr/msa-pdb-mining` on GitHub. No CI configured.
- `.gitignore` excludes `.venv/`, `__pycache__/`, `*.egg-info/`, `results/`, `.pytest_cache/`.
- `.claude/settings.local.json` pre-allows the `pip install`, `pytest -q`, and
  `msa-pdb-mining` Bash invocations.
- Version is single-sourced in `msa_pdb_mining/__init__.py` (`__version__`) and `pyproject.toml`.

## Status & roadmap (from README)

**Implemented:** `--uniprot`/`--gene`/`--sequence` → ColabFold MSA → experimental-PDB depth →
HTML + image + data (mako palette, query secondary-structure overlay), with full API-response
caching.

**Planned:** species-list ortholog mode (OrthoDB/OMA + alignment); pluggable local MSA
backends (the `MSABackend` ABC exists for exactly this); optional AlphaFold/computed-model
coverage; sequence-search fallback for hits without a UniProt accession.
