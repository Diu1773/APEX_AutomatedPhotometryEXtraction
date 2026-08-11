"""A Gaia query failure must leave the exception itself in the log.

The classification alone is not enough to act on, and this was found the hard
way. A Phase 3 run reported "ESA TAP failed [TIMEOUT]" eight times out of eight
and the raw exception was never written anywhere: on the common path — where
the VizieR fallback then succeeds — the original exception was caught, labelled
and discarded. So there was no way afterwards to tell a slow server from a
refused connection. (ESA answered in 4.8 s when retried by hand the next
morning, so the label was at best incomplete.)

Two things make the label checkable, and these tests pin both:

* the exception text, verbatim, so the actual failure can be read; and
* the elapsed time, because the label alone cannot say how long APEX waited.
  A hard deadline fires at ~max(2*timeout_s, 60) s; a refused connection
  returns in milliseconds. The duration separates a server that was slow from
  one that was never reached, and it is also how the settings in force
  (timeout_s, hard_deadline_s) get checked against reality.
"""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest
from astropy.coordinates import SkyCoord

import apex.utils.gaia_catalog_service as gcs
from apex.utils.gaia_catalog_service import GaiaCatalogService


@pytest.fixture()
def service_and_log(tmp_path, monkeypatch):
    """A service whose ESA query always fails, with its log captured."""
    lines: list[str] = []
    params = SimpleNamespace(P=SimpleNamespace(
        gaia_mag_max=20.0, gaia_wcs_mag_max=18.0, gaia_retry=1,
        gaia_backoff_s=0.0, gaia_allow_no_cache=True, gaia_timeout_s=30.0,
    ))
    svc = GaiaCatalogService(params, tmp_path, log_fn=lines.append)
    monkeypatch.setattr(gcs, "_HAS_GAIA", True)
    monkeypatch.setattr(gcs, "Gaia", SimpleNamespace(ROW_LIMIT=-1))
    monkeypatch.setattr(gcs, "_get_gaia", lambda: gcs.Gaia)
    return svc, lines


def _fail_esa_with(monkeypatch, exc: Exception):
    def _boom(*_a, **_k):
        raise exc
    monkeypatch.setattr(gcs, "tap_query_with_deadline", _boom)


def test_exception_text_survives_a_successful_fallback(
        service_and_log, monkeypatch):
    """The case that lost the evidence: fallback works, so nothing raises."""
    svc, lines = service_and_log
    _fail_esa_with(monkeypatch, TimeoutError("TAP query timed out (hard deadline 60s)"))
    monkeypatch.setattr(
        GaiaCatalogService, "_query_gaia_vizier",
        lambda self, *a, **k: pd.DataFrame(
            {"source_id": [1], "ra": [10.0], "dec": [20.0],
             "phot_g_mean_mag": [15.0], "ruwe": [1.0]}),
    )

    df = svc._query_gaia(SkyCoord(10.0, 20.0, unit="deg"), 0.1, 18.0)
    assert df.attrs["gaia_source"] == "vizier_fallback", "fallback must succeed"

    blob = "\n".join(lines)
    assert "hard deadline 60s" in blob, "the exception text must be logged verbatim"
    assert "TimeoutError" in blob, "the exception type must be logged"
    assert "[TIMEOUT]" in blob, "the classification is still useful, keep it"


def test_elapsed_time_and_settings_are_logged(service_and_log, monkeypatch):
    """Duration is what separates a real timeout from an instant refusal."""
    svc, lines = service_and_log
    _fail_esa_with(monkeypatch, ConnectionRefusedError("Connection refused"))
    monkeypatch.setattr(
        GaiaCatalogService, "_query_gaia_vizier",
        lambda self, *a, **k: pd.DataFrame(
            {"source_id": [1], "ra": [10.0], "dec": [20.0],
             "phot_g_mean_mag": [15.0], "ruwe": [1.0]}),
    )

    svc._query_gaia(SkyCoord(10.0, 20.0, unit="deg"), 0.1, 18.0)

    blob = "\n".join(lines)
    assert "Connection refused" in blob
    assert "after 0." in blob, "an instant failure must be visibly instant"
    assert "timeout_s=30" in blob and "hard_deadline_s=60" in blob, (
        "the settings in force must be logged next to the duration, or the "
        "duration cannot be judged"
    )


def test_connection_timeout_is_labelled_network_error(
        service_and_log, monkeypatch):
    """This assertion used to read `[TIMEOUT]`, and that was the bug.

    When the raw text was first added to the log, the classifier still tested
    a bare "timed out" before any connection test, so an unreachable host was
    labelled TIMEOUT and only the text gave it away. The ordering has since
    been fixed (see test_gaia_failure_classification.py), so the label now
    agrees with the text instead of contradicting it — and the text is still
    there, which is what makes the label checkable at all.
    """
    svc, lines = service_and_log
    _fail_esa_with(monkeypatch, OSError("Connection timed out"))
    monkeypatch.setattr(
        GaiaCatalogService, "_query_gaia_vizier",
        lambda self, *a, **k: pd.DataFrame(
            {"source_id": [1], "ra": [10.0], "dec": [20.0],
             "phot_g_mean_mag": [15.0], "ruwe": [1.0]}),
    )

    svc._query_gaia(SkyCoord(10.0, 20.0, unit="deg"), 0.1, 18.0)

    blob = "\n".join(lines)
    assert "[NETWORK_ERROR]" in blob, "an unreachable host is not a slow server"
    assert "Connection timed out" in blob, "the text must still be verbatim"


def test_fallback_failure_also_logs_its_own_exception(
        service_and_log, monkeypatch):
    svc, lines = service_and_log
    _fail_esa_with(monkeypatch, TimeoutError("esa timed out"))

    def _viz_boom(self, *_a, **_k):
        raise RuntimeError("VizieR 500 Internal Server Error")
    monkeypatch.setattr(GaiaCatalogService, "_query_gaia_vizier", _viz_boom)

    with pytest.raises(RuntimeError):
        svc._query_gaia(SkyCoord(10.0, 20.0, unit="deg"), 0.1, 18.0)

    blob = "\n".join(lines)
    assert "esa timed out" in blob
    assert "VizieR 500 Internal Server Error" in blob
