"""Retried HTTP reads (src/whale/net.py). No network access."""
from __future__ import annotations

import io
import urllib.error

import pytest

from whale.net import get_text


def _opener(outcomes):
    """Each call takes the next outcome: an exception to raise, or a body to return."""
    calls = []

    def opener(url, timeout):
        calls.append(url)
        out = outcomes[len(calls) - 1]
        if isinstance(out, Exception):
            raise out
        return io.BytesIO(out.encode())
    return opener, calls


def _http(code):
    return urllib.error.HTTPError("https://x", code, "err", {}, None)


def test_transient_failures_are_retried():
    opener, calls = _opener([urllib.error.URLError("reset"), TimeoutError("slow"), _http(503), "ok"])
    slept = []
    assert get_text("https://x", opener=opener, sleep=slept.append) == "ok"
    assert len(calls) == 4 and slept == [5, 15, 45]


def test_gives_up_after_the_last_retry():
    opener, calls = _opener([TimeoutError("slow")] * 4)
    with pytest.raises(TimeoutError):
        get_text("https://x", opener=opener, sleep=lambda s: None)
    assert len(calls) == 4


@pytest.mark.parametrize("code", [400, 404])
def test_client_errors_are_raised_at_once(code):
    opener, calls = _opener([_http(code), "never"])
    with pytest.raises(urllib.error.HTTPError):
        get_text("https://x", opener=opener, sleep=lambda s: None)
    assert len(calls) == 1


def test_rate_limit_is_retried():
    opener, calls = _opener([_http(429), "ok"])
    assert get_text("https://x", opener=opener, sleep=lambda s: None) == "ok"
