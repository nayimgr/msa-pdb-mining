"""Pluggable MSA backend interface.

The homology search is abstracted behind ``MSABackend`` so that the remote
ColabFold MMseqs2 engine used today can later be swapped for a local jackhmmer /
MMseqs2 engine without touching the rest of the pipeline. A backend's only job is
to turn a query sequence into an A3M string anchored on that query.
"""

from __future__ import annotations

import abc


class MSABackend(abc.ABC):
    name: str = "base"

    @abc.abstractmethod
    def run(self, query_sequence: str) -> str:
        """Return an A3M string whose first record is the query."""
        raise NotImplementedError
