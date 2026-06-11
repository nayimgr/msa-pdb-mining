"""Command-line entry point."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import Config
from .pipeline import parse_formats, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="msa-pdb-mining",
        description="Build an MSA for a protein and colour it by experimental-structure "
        "(PDB) coverage depth across orthologs.",
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--uniprot", help="UniProt accession, e.g. P04637")
    src.add_argument("--gene", help="Gene name (use with --organism)")
    src.add_argument("--sequence", help="Raw amino-acid sequence (no DB lookup for the query)")

    p.add_argument("--organism", help="Organism for --gene: NCBI taxid or name (default: 9606)")
    p.add_argument("--name", help="Display name when using --sequence")
    p.add_argument("--out", "-o", default="results", help="Output directory (default: results)")
    p.add_argument(
        "--formats", default="all",
        help="Comma list of html,image,data (or 'all'; default: all)",
    )
    p.add_argument("--email", help="Contact email sent to EBI/ColabFold (polite, recommended)")
    p.add_argument("--cache-dir", help="Override cache directory")
    p.add_argument("--no-cache", action="store_true", help="Disable the on-disk cache")
    p.add_argument(
        "--mode", default="env",
        help="ColabFold MMseqs2 mode (default: env)",
    )
    p.add_argument(
        "--include-env", action="store_true",
        help="Also include environmental/metagenomic hits (no PDB; off by default)",
    )
    p.add_argument(
        "--reviewed-only", action="store_true",
        help="Keep only reviewed (Swiss-Prot) hits — drops uncharacterised TrEMBL noise",
    )
    p.add_argument(
        "--max-structured-rows", type=int, default=200,
        help="Max distinct accessions resolved to observed PDB residues (default: 200)",
    )
    p.add_argument(
        "--max-lookup-accessions", type=int, default=2000,
        help="Max accessions to fetch UniProt metadata for (batched; default: 2000)",
    )
    p.add_argument(
        "--max-rows-render", type=int, default=60,
        help="Max alignment rows drawn in the image/HTML (default: 60)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )

    config = Config(
        email=args.email,
        use_cache=not args.no_cache,
        colabfold_mode=args.mode,
        include_env_hits=args.include_env,
        reviewed_only=args.reviewed_only,
        max_structured_rows=args.max_structured_rows,
        max_lookup_accessions=args.max_lookup_accessions,
        max_rows_render=args.max_rows_render,
    )
    if args.cache_dir:
        config.cache_dir = Path(args.cache_dir)

    formats = parse_formats(args.formats)
    out_dir = Path(args.out)

    try:
        result = run_pipeline(
            config, out_dir, formats,
            uniprot=args.uniprot, gene=args.gene, organism=args.organism,
            sequence=args.sequence, name=args.name,
        )
    except (ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    matrix = result["matrix"]
    agg = matrix.column_aggregate
    covered = int((agg > 0).sum())
    print(f"\nQuery: {result['query'].name} ({matrix.columns} aa)")
    print(f"MSA rows: {matrix.n_rows}  |  rows with structures: "
          f"{sum(1 for s in matrix.row_status if s == 'structured')}")
    print(f"Query residues covered by >=1 structure: {covered}/{matrix.columns} "
          f"({covered / matrix.columns:.0%})" if matrix.columns else "")
    print(f"\nWrote {len(result['outputs'])} files to {out_dir}/:")
    for path in result["outputs"]:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
