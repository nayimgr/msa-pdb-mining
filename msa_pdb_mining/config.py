"""Runtime configuration: endpoints, cache location, network etiquette."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import platformdirs

COLABFOLD_HOST = "https://api.colabfold.com"
PDBE_API = "https://www.ebi.ac.uk/pdbe/api"
PDBE_GRAPH_API = "https://www.ebi.ac.uk/pdbe/graph-api"
UNIPROT_REST = "https://rest.uniprot.org"


def _default_cache_dir() -> Path:
    return Path(platformdirs.user_cache_dir("msa-pdb-mining"))


@dataclass
class Config:
    """All tunable knobs for a run. Sensible defaults; overridden by the CLI."""

    # --- network etiquette ---
    email: str | None = None
    request_timeout: float = 30.0

    # --- caching ---
    cache_dir: Path = field(default_factory=_default_cache_dir)
    use_cache: bool = True

    # --- MSA backend (ColabFold MMseqs2) ---
    colabfold_host: str = COLABFOLD_HOST
    # "env" searches UniRef30 + environmental DBs. We consume only the UniRef
    # hits by default (environmental/metagenomic hits never have PDB entries).
    colabfold_mode: str = "env"
    include_env_hits: bool = False
    poll_interval: float = 5.0
    max_poll_seconds: float = 1800.0

    # --- structural mapping ---
    pdbe_api: str = PDBE_API
    pdbe_graph_api: str = PDBE_GRAPH_API
    uniprot_rest: str = UNIPROT_REST
    # Cap how many distinct accessions we resolve to PDB (highest-identity first).
    max_structured_rows: int = 200
    # Cap how many accessions we look up metadata for (batched, cheap).
    max_lookup_accessions: int = 2000
    # Keep only reviewed (Swiss-Prot) hits, dropping uncharacterised TrEMBL noise.
    reviewed_only: bool = False

    # --- rendering ---
    max_rows_render: int = 60  # rows shown in the image/HTML (data files keep all)

    def user_agent(self) -> str:
        contact = f"; mailto:{self.email}" if self.email else ""
        return f"msa-pdb-mining/0.1 (+https://github.com/msa-pdb-mining{contact})"
