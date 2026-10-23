"""Small HTTP reads (ENA reports, NCBI listings) retried through transient failures.

ENA's portal and NCBI's FTP-over-HTTPS both drop the odd request. A timeout,
reset or 5xx is retried with back-off; a 4xx other than 429 is an answer
(wrong accession, missing file) and is raised at once for the caller to handle.
"""
from __future__ import annotations

import time
import urllib.error
import urllib.request
from http.client import HTTPException

TIMEOUT = 120
BACKOFF = [5, 15, 45]


def _transient(exc: Exception) -> bool:
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code == 429 or exc.code >= 500
    return isinstance(exc, (urllib.error.URLError, HTTPException, OSError))


def get_text(url: str, *, timeout: float = TIMEOUT, backoff: list[float] = BACKOFF,
             sleep=time.sleep, opener=urllib.request.urlopen) -> str:
    """The body of `url` as UTF-8 text, tried len(backoff) + 1 times."""
    for attempt in range(len(backoff) + 1):
        try:
            with opener(url, timeout=timeout) as fh:  # noqa: S310 - fixed EBI / NCBI hosts
                return fh.read().decode("utf-8")
        except (urllib.error.URLError, HTTPException, OSError) as exc:
            if not _transient(exc) or attempt == len(backoff):
                raise
            sleep(backoff[attempt])
    raise AssertionError("unreachable")
