"""Tiny on-disk cache for API responses, keyed by namespace + key.

Keeps re-runs cheap and offline-friendly, and avoids hammering the public
EBI / ColabFold services during iterative development.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Callable, Optional

_SAFE = re.compile(r"[^A-Za-z0-9._-]")


def _safe_name(key: str) -> str:
    cleaned = _SAFE.sub("_", key)
    if len(cleaned) > 100:  # keep filenames sane; disambiguate with a hash
        digest = hashlib.sha1(key.encode()).hexdigest()[:12]
        cleaned = cleaned[:80] + "_" + digest
    return cleaned


class DiskCache:
    def __init__(self, root: Path, enabled: bool = True) -> None:
        self.root = Path(root)
        self.enabled = enabled

    def _path(self, namespace: str, key: str, ext: str) -> Path:
        return self.root / namespace / f"{_safe_name(key)}.{ext}"

    # --- JSON ---
    def get_json(self, namespace: str, key: str) -> Optional[Any]:
        if not self.enabled:
            return None
        path = self._path(namespace, key, "json")
        if path.exists():
            try:
                return json.loads(path.read_text())
            except (ValueError, OSError):
                return None
        return None

    def set_json(self, namespace: str, key: str, value: Any) -> None:
        if not self.enabled:
            return
        path = self._path(namespace, key, "json")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))

    def cached_json(
        self, namespace: str, key: str, fetch: Callable[[], Optional[Any]]
    ) -> Optional[Any]:
        hit = self.get_json(namespace, key)
        if hit is not None:
            return hit
        value = fetch()
        if value is not None:
            self.set_json(namespace, key, value)
        return value

    # --- text (e.g. a3m) ---
    def get_text(self, namespace: str, key: str) -> Optional[str]:
        if not self.enabled:
            return None
        path = self._path(namespace, key, "txt")
        if path.exists():
            try:
                return path.read_text()
            except OSError:
                return None
        return None

    def set_text(self, namespace: str, key: str, value: str) -> None:
        if not self.enabled:
            return
        path = self._path(namespace, key, "txt")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value)
