"""ColabFold MMseqs2 remote MSA backend.

Protocol (https://api.colabfold.com): POST the query FASTA to ``ticket/msa``,
poll ``ticket/<id>`` until ``COMPLETE``, then download ``result/download/<id>``
(a gzipped tar). We consume ``uniref.a3m`` by default (environmental hits never
have PDB entries); ``bfd.mgnify30.metaeuk30.smag30.a3m`` is merged in optionally.
"""

from __future__ import annotations

import hashlib
import io
import logging
import tarfile
import time
from typing import List, Optional

import httpx

from ..cache import DiskCache
from ..config import Config
from .a3m import parse_a3m
from .base import MSABackend

log = logging.getLogger(__name__)

_UNIREF_MEMBER = "uniref.a3m"
_ENV_MEMBER = "bfd.mgnify30.metaeuk30.smag30.a3m"
_TERMINAL_OK = {"COMPLETE"}
_RUNNING = {"PENDING", "RUNNING", "UNKNOWN"}
_RETRY_SUBMIT = {"RATELIMIT", "MAINTENANCE"}


class ColabFoldError(RuntimeError):
    pass


class ColabFoldMMseqs2Backend(MSABackend):
    name = "colabfold"

    def __init__(self, client: httpx.Client, cache: DiskCache, config: Config) -> None:
        self.client = client
        self.cache = cache
        self.config = config

    def _cache_key(self, query_sequence: str) -> str:
        payload = f"{query_sequence}|{self.config.colabfold_mode}|env={self.config.include_env_hits}"
        return hashlib.sha1(payload.encode()).hexdigest()

    def run(self, query_sequence: str) -> str:
        key = self._cache_key(query_sequence)
        cached = self.cache.get_text("a3m", key)
        if cached is not None:
            log.info("Using cached MSA (%s)", key[:8])
            return cached

        archive = self._submit_and_download(query_sequence)
        a3m = self._extract_a3m(archive)
        self.cache.set_text("a3m", key, a3m)
        return a3m

    # --- HTTP workflow ---
    def _submit_and_download(self, query_sequence: str) -> bytes:
        host = self.config.colabfold_host
        fasta = f">101\n{query_sequence}\n"
        ticket = self._submit(host, fasta)
        ticket_id = ticket["id"]
        log.info("ColabFold ticket %s submitted", ticket_id)
        self._poll(host, ticket_id)
        resp = self.client.get(f"{host}/result/download/{ticket_id}")
        resp.raise_for_status()
        return resp.content

    def _submit(self, host: str, fasta: str) -> dict:
        deadline = time.time() + self.config.max_poll_seconds
        while True:
            resp = self.client.post(
                f"{host}/ticket/msa",
                data={"q": fasta, "mode": self.config.colabfold_mode},
            )
            resp.raise_for_status()
            ticket = resp.json()
            status = ticket.get("status", "UNKNOWN")
            if status in _RETRY_SUBMIT:
                if time.time() > deadline:
                    raise ColabFoldError(f"submission kept returning {status}")
                log.warning("ColabFold %s; waiting before resubmit", status)
                time.sleep(max(self.config.poll_interval, 15.0))
                continue
            if "id" not in ticket:
                raise ColabFoldError(f"unexpected submit response: {ticket}")
            return ticket

    def _poll(self, host: str, ticket_id: str) -> None:
        deadline = time.time() + self.config.max_poll_seconds
        while True:
            resp = self.client.get(f"{host}/ticket/{ticket_id}")
            resp.raise_for_status()
            status = resp.json().get("status", "UNKNOWN")
            if status in _TERMINAL_OK:
                return
            if status not in _RUNNING:
                raise ColabFoldError(f"ColabFold job {ticket_id} failed: {status}")
            if time.time() > deadline:
                raise ColabFoldError(f"ColabFold job {ticket_id} timed out")
            time.sleep(self.config.poll_interval)

    # --- archive handling ---
    def _extract_a3m(self, archive: bytes) -> str:
        members = self._read_members(archive)
        if _UNIREF_MEMBER not in members:
            raise ColabFoldError(f"{_UNIREF_MEMBER} missing from result archive")
        a3m = members[_UNIREF_MEMBER]
        if self.config.include_env_hits and _ENV_MEMBER in members:
            a3m = _merge_a3m(a3m, members[_ENV_MEMBER])
        return a3m

    @staticmethod
    def _read_members(archive: bytes) -> dict:
        out: dict = {}
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:*") as tar:
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                out[member.name.split("/")[-1]] = fh.read().decode("utf-8", "replace")
        return out


def _merge_a3m(primary: str, secondary: str) -> str:
    """Append the hit records of ``secondary`` (dropping its query) onto ``primary``."""
    extra: List[str] = []
    records = parse_a3m(secondary)
    for rec in records[1:]:  # skip the duplicated query record
        extra.append(f">{rec.header}\n{rec.seq}")
    if not extra:
        return primary
    return primary.rstrip("\n") + "\n" + "\n".join(extra) + "\n"
