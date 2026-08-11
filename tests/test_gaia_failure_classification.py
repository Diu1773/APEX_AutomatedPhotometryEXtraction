"""The failure label has to name the cause, not the first word that matched.

Ordering is the whole content of these classifiers, and it was wrong twice:

* a bare "timed out" was tested before any connection test, so "Connection
  timed out" — a host that cannot be reached — was filed as TIMEOUT;
* SSL was tested last, after the generic "connection" test, so an untrusted
  certificate (which surfaces as a connection error) was filed as
  NETWORK_ERROR — the one label that suggests waiting rather than fixing.

Neither is cosmetic. A Phase 3 run reported TIMEOUT eight times out of eight
and the label was the only thing kept, so "the server is slow" and "this
machine could not reach it" were indistinguishable afterwards; ESA answered in
4.8 s when retried by hand. TIMEOUT and NETWORK_ERROR route to the same VizieR
fallback, so what changed is the record, not the behaviour — and the record is
what the next investigation reads.
"""

from __future__ import annotations

import ssl
from types import SimpleNamespace

import pytest

from apex.utils.gaia_catalog_service import (
    GaiaCatalogService,
    gaia_failure_reason,
)


@pytest.fixture()
def classify(tmp_path):
    svc = GaiaCatalogService(SimpleNamespace(P=SimpleNamespace()), tmp_path)
    return svc._classify_query_failure


@pytest.mark.parametrize("message, label", [
    # The regression: transport faults that contain the word "timeout".
    ("Connection timed out", "NETWORK_ERROR"),
    ("[Errno 110] Connection timed out", "NETWORK_ERROR"),
    ("HTTPSConnectionPool: Failed to establish a new connection", "NETWORK_ERROR"),
    ("[Errno 111] Connection refused", "NETWORK_ERROR"),
    ("Network is unreachable", "NETWORK_ERROR"),
    ("getaddrinfo failed", "NETWORK_ERROR"),
    ("Connection reset by peer", "NETWORK_ERROR"),
    # Genuine timeouts stay timeouts.
    ("TAP query timed out (hard deadline 60s)", "TIMEOUT"),
    ("Read timed out", "TIMEOUT"),
    ("The read operation timed out", "TIMEOUT"),
    # The specific causes keep priority over both.
    ("503 Service Unavailable", "SERVER_DOWN"),
    ("404 job not found", "SERVER_JOB_LOST"),
    ("Your IP has been blocked due to heavy usage", "IP_BANNED"),
])
def test_classification_names_the_cause(classify, message, label):
    assert classify(RuntimeError(message)) == label


def test_our_own_deadline_wins_over_a_connection_word(classify):
    """APEX raised it, so it is a timeout however the text reads."""
    exc = TimeoutError(
        "TAP query timed out (hard deadline 60s) while holding a connection")
    assert classify(exc) == "TIMEOUT"


def test_ssl_is_recognised_even_though_it_reads_as_a_connection_error(classify):
    """The second regression: a fixable local trust problem, not a network one."""
    exc = ssl.SSLCertVerificationError(
        "[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: "
        "unable to get local issuer certificate (_ssl.c:1006)")
    assert classify(exc) == "SSL_CERTIFICATE_VERIFY_FAILED"


def test_unknown_stays_unknown(classify):
    assert classify(RuntimeError("something entirely new")) == "UNKNOWN"


@pytest.mark.parametrize("message, reason", [
    ("Connection timed out", "network_error"),
    ("TAP query timed out (hard deadline 60s)", "timeout"),
    ("Read timed out", "timeout"),
    ("[Errno 111] Connection refused", "network_error"),
    ("503 Service Unavailable", "server_down"),
])
def test_module_level_reason_uses_the_same_ordering(message, reason):
    """`gaia_failure_reason` feeds the same logs and must not disagree."""
    assert gaia_failure_reason(RuntimeError(message)) == reason


def test_reason_still_reports_timeout_for_the_deadline_error():
    """Pins the contract test_gaia_deadline.py depends on."""
    assert gaia_failure_reason(
        TimeoutError("TAP query timed out (hard deadline 60s)")) == "timeout"
