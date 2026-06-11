"""Shared HTTP client with retry/backoff and polite headers."""

from __future__ import annotations

import time
from typing import Any, Optional

import httpx

from .config import Config

# Status codes worth retrying (transient server / rate-limit conditions).
_RETRY_STATUS = {429, 500, 502, 503, 504}


def make_client(config: Config) -> httpx.Client:
    transport = httpx.HTTPTransport(retries=2)  # connection-level retries
    return httpx.Client(
        headers={"User-Agent": config.user_agent(), "Accept": "application/json"},
        timeout=config.request_timeout,
        transport=transport,
        follow_redirects=True,
    )


def request(
    client: httpx.Client,
    method: str,
    url: str,
    *,
    retries: int = 4,
    backoff: float = 2.0,
    **kwargs: Any,
) -> Optional[httpx.Response]:
    """Issue a request, retrying transient failures with exponential backoff.

    Returns the response (even for non-retryable 4xx so callers can inspect the
    status), or ``None`` if every attempt failed to connect / kept being throttled.
    """
    for attempt in range(retries + 1):
        try:
            resp = client.request(method, url, **kwargs)
        except httpx.HTTPError:  # connection/read errors
            resp = None
        else:
            if resp.status_code not in _RETRY_STATUS:
                return resp
        if attempt < retries:
            time.sleep(backoff * (2 ** attempt))
    return None


def get_json(client: httpx.Client, url: str, **kwargs: Any) -> Optional[Any]:
    """GET a URL and parse JSON. Returns ``None`` on 404 / failure / empty body."""
    resp = request(client, "GET", url, **kwargs)
    if resp is None or resp.status_code == 404:
        return None
    if resp.status_code >= 400:
        return None
    try:
        return resp.json()
    except ValueError:
        return None
