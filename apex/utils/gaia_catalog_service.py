"""Shared Gaia catalog loading/query service for Step 5 WCS workflows."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table

from apex.utils.gaia_columns import (
    canonical_rename_map,
    esa_select_clause,
    missing_columns,
    vizier_select_clause,
)
from apex.utils.io_utils import coerce_int64_source_id
from apex.utils.ssl_certificates import configure_ssl_certificates, ssl_error_reason
from apex.utils.step_paths import step5_wcs_dir

_SSL_CERT_OK, _SSL_CERT_DETAIL = configure_ssl_certificates()

# astroquery.gaia pulls in TAP/VOTable/SAMP/astropy.io.votable on import — in a
# PyInstaller frozen build that resolution can take 5–30 s on the first call.
# Importing it at module scope used to freeze the main window for the whole
# duration the first time Step 5 was opened, because step5 imports this module
# at its own top level. Defer the import to first use via _get_gaia(), and let
# the rest of the module reference Gaia through that accessor.
Gaia = None
_HAS_GAIA: bool | None = None  # tri-state: None = not yet probed


def _get_gaia():
    """Lazy import of astroquery.gaia.Gaia. Returns the module or None."""
    global Gaia, _HAS_GAIA
    if _HAS_GAIA is not None:
        return Gaia
    try:
        from astroquery.gaia import Gaia as _Gaia
        Gaia = _Gaia
        _HAS_GAIA = True
    except Exception:
        Gaia = None
        _HAS_GAIA = False
    return Gaia


LogFn = Callable[[str], None]
StopFn = Callable[[], bool]


def _tail_text(value: object, limit: int = 800, max_lines: int = 8) -> str:
    if value is None:
        return ""
    s = str(value).replace("\r\n", "\n").replace("\r", "\n").strip()
    if not s:
        return ""
    lines = [ln.strip() for ln in s.split("\n") if ln.strip()]
    if max_lines > 0 and len(lines) > max_lines:
        lines = lines[-max_lines:]
    one_line = " | ".join(lines)
    if len(one_line) > limit:
        one_line = "..." + one_line[-limit:]
    return one_line


def _exc_brief(exc: Exception, limit: int = 260) -> str:
    return _tail_text(f"{type(exc).__name__}: {exc}", limit=limit, max_lines=4)


# A transport-level fault that happens to contain the word "timeout". These
# have to be recognised BEFORE a bare "timed out" test, or a host that cannot
# be reached at all is filed as a slow server. That mislabel is not academic:
# a Phase 3 run reported TIMEOUT eight times out of eight, and the label was
# the only thing recorded, so "the ESA service is slow" and "this machine
# could not reach it" stayed indistinguishable until ESA was re-queried by
# hand the next morning and answered in 4.8 s.
_CONNECTION_FAULT_PHRASES = (
    "connection refused",
    "connection reset",
    "connection aborted",
    "connection timed out",
    "connection timeout",
    "failed to establish a new connection",
    "cannot connect",
    "unreachable",
    "getaddrinfo",
    "name or service not known",
    "nodename nor servname",
    "temporary failure in name resolution",
)

# Text APEX itself puts in the exception it raises when its own wall-clock
# deadline fires. That one is a timeout by construction, whatever else the
# message says.
_OWN_DEADLINE_MARKER = "hard deadline"


def gaia_failure_reason(exc: Exception | None) -> str:
    if exc is None:
        return "unknown"
    ssl_reason = ssl_error_reason(exc)
    if ssl_reason:
        return ssl_reason
    text = _exc_brief(exc, limit=180)
    low = text.lower()
    if "server_down" in low or "503" in low or "502" in low or "500" in low:
        return "server_down"
    if _OWN_DEADLINE_MARKER in low:
        return "timeout"
    if any(p in low for p in _CONNECTION_FAULT_PHRASES):
        return "network_error"
    if "timeout" in low or "timed out" in low:
        return "timeout"
    if "connection" in low or "refused" in low:
        return "network_error"
    if "vizier fallback" in low:
        return "vizier_fallback_failed"
    return text or "unknown"


def tap_launch_job_async(service, adql: str, *, timeout_s: float | None = None, **kwargs):
    if timeout_s is not None and np.isfinite(timeout_s) and timeout_s > 0:
        try:
            setattr(service, "TIMEOUT", float(timeout_s))
        except Exception:
            pass
        try:
            return service.launch_job_async(adql, timeout=float(timeout_s), **kwargs)
        except TypeError as exc:
            if "timeout" not in str(exc).lower():
                raise
    return service.launch_job_async(adql, **kwargs)


def tap_query_with_deadline(service, adql: str, *, timeout_s: float | None = None,
                            deadline_s: float | None = None, **kwargs):
    """Run an async TAP job and fetch its results under a HARD wall-clock deadline.

    astroquery's ``job.get_results()`` polls the server and does NOT honour the
    requested timeout, so a stuck ESA/VizieR queue blocks the caller forever
    (the "10-minute hang" seen in practice). Run launch + fetch on a worker
    thread and abandon it if it exceeds ``deadline_s``, raising a ``TimeoutError``
    (message contains "timed out" so ``_classify_query_failure`` maps it to
    TIMEOUT and callers fall back to another mirror).

    The abandoned query is left to finish on a **daemon** thread. That detail
    is the whole point: an earlier version used a ``ThreadPoolExecutor`` and
    dropped it with ``shutdown(wait=False)``, which frees the *caller* but not
    the *process* — futures' worker threads are non-daemon and
    ``concurrent.futures`` installs an ``atexit`` hook that joins them, so the
    interpreter sat waiting for a query nobody was listening to any more.
    Measured on a benchmark run whose pipeline finished in 226 s and whose
    process took 2,424 s to exit; three of six sweep points were inflated the
    same way, and the giveaway was "Query finished" printing *after* the run
    summary. A daemon thread is killed at exit instead.
    """
    import queue as _queue
    import threading as _threading

    if not (deadline_s is not None and np.isfinite(deadline_s) and deadline_s > 0):
        base = float(timeout_s) if (timeout_s and np.isfinite(timeout_s) and timeout_s > 0) else 30.0
        deadline_s = max(base * 2.0, 60.0)

    outcome: _queue.Queue = _queue.Queue(maxsize=1)

    def _run():
        try:
            job = tap_launch_job_async(service, adql, timeout_s=timeout_s, **kwargs)
            outcome.put(("ok", job.get_results()))
        except BaseException as exc:  # noqa: BLE001 - relayed to the caller
            outcome.put(("error", exc))

    _threading.Thread(target=_run, name="apex-tap-query", daemon=True).start()
    try:
        status, payload = outcome.get(timeout=float(deadline_s))
    except _queue.Empty as exc:
        raise TimeoutError(
            f"TAP query timed out (hard deadline {deadline_s:.0f}s)"
        ) from exc
    if status == "error":
        raise payload
    return payload


def gaia_runtime_available() -> tuple[bool, str]:
    ssl_detail = "" if _SSL_CERT_OK else f" SSL setup issue: {_SSL_CERT_DETAIL}"
    _get_gaia()  # probe lazily; safe to call repeatedly
    if not _HAS_GAIA:
        return (
            False,
            "astroquery.gaia is not importable. Gaia catalog attach, WCS refine, "
            f"resid_med, and Step 6 Gaia stats require astroquery in the packaged build.{ssl_detail}",
        )
    try:
        from astroquery.utils.tap.core import TapPlus  # noqa: F401
    except Exception as exc:
        return (
            True,
            f"astroquery.gaia is importable, but VizieR TAP fallback is unavailable: {_exc_brief(exc)}.{ssl_detail}",
        )
    return True, f"astroquery.gaia and TAP fallback imports are available.{ssl_detail}"


def coerce_source_id_int64(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize Gaia source_id as signed int64 without float precision loss."""
    if df is None or df.empty or "source_id" not in df.columns:
        return df
    out = df.copy()
    sid = coerce_int64_source_id(out["source_id"])
    valid = sid.notna()
    out = out.loc[valid].copy()
    out["source_id"] = sid.loc[valid].astype("int64")
    return out


def _normalise_gaia_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename either service's column names to the canonical set.

    The mapping is generated from apex.utils.gaia_columns, so a column added to
    the contract is understood here without a second edit.
    """
    out = df.copy()
    canonical = canonical_rename_map()
    rename = {c: canonical[str(c).strip().lower()]
              for c in out.columns if str(c).strip().lower() in canonical}
    out = out.rename(columns=rename)
    out.columns = [str(c).strip().lower() for c in out.columns]
    if "phot_variable_flag" not in out.columns:
        out["phot_variable_flag"] = ""
    if "bp_rp" not in out.columns and {"phot_bp_mean_mag", "phot_rp_mean_mag"} <= set(out.columns):
        bp = pd.to_numeric(out["phot_bp_mean_mag"], errors="coerce")
        rp = pd.to_numeric(out["phot_rp_mean_mag"], errors="coerce")
        out["bp_rp"] = bp - rp
    return coerce_source_id_int64(out)


class GaiaCatalogService:
    """Load/query/cache the Step 5 Gaia field catalog through one code path."""

    def __init__(
        self,
        params,
        result_dir: Path | str,
        *,
        log_fn: LogFn | None = None,
        stop_fn: StopFn | None = None,
    ):
        self.params = params
        self.result_dir = Path(result_dir)
        self.log_fn = log_fn
        self.stop_fn = stop_fn

    def _log(self, msg: str) -> None:
        if self.log_fn is not None:
            self.log_fn(str(msg))

    def _stopped(self) -> bool:
        try:
            return bool(self.stop_fn and self.stop_fn())
        except Exception:
            return False

    @property
    def P(self):
        return self.params.P

    def _load_gaia_cache_if_ok(self, path: Path) -> pd.DataFrame | None:
        if not path.exists():
            return None
        try:
            tab = Table.read(path, format="ascii.ecsv")
            raw = _normalise_gaia_columns(tab.to_pandas())
            if not {"ra", "dec"} <= set(raw.columns):
                return None
            if "source_id" not in raw.columns and _HAS_GAIA:
                return None
            # Which contract columns this cached catalogue is missing, judged
            # BEFORE the normaliser fills in defaults. `only_used=True` skips
            # ESA-only columns: a VizieR catalogue has no `phot_variable_flag`
            # and never will, so re-querying cannot supply it. Treating that
            # absence as staleness is what made every VizieR-written cache miss
            # on every run — the catalogue was re-fetched each time, and each
            # fetch paid the ESA timeout before falling back.
            pre = pd.DataFrame(columns=[str(c) for c in tab.colnames])
            pre = _normalise_gaia_columns(pre)
            df = raw
            df.attrs["missing_contract_columns"] = [
                c.name for c in missing_columns(pre, only_used=True)]
            return df
        except Exception:
            return None

    def _query_gaia_vizier(self, center: SkyCoord, radius_deg: float, mag_max: float, timeout_s: float) -> pd.DataFrame:
        try:
            from astroquery.utils.tap.core import TapPlus
        except ImportError:
            raise RuntimeError("astroquery.utils.tap not available")
        mag_where = f'AND "I/355/gaiadr3".Gmag <= {mag_max:.4f}' if np.isfinite(mag_max) and mag_max > 0 else ""
        adql = f"""
SELECT
  {vizier_select_clause()}
FROM "I/355/gaiadr3"
WHERE 1=CONTAINS(
    POINT('ICRS', "I/355/gaiadr3".RA_ICRS, "I/355/gaiadr3".DE_ICRS),
    CIRCLE('ICRS', {center.ra.deg:.8f}, {center.dec.deg:.8f}, {radius_deg:.8f})
)
{mag_where}
        """.strip()
        tap = TapPlus(url="https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
        deadline_s = float(getattr(self.P, "gaia_hard_deadline_s", 0.0) or 0.0) or None
        try:
            tab = tap_query_with_deadline(
                tap, adql, timeout_s=timeout_s, deadline_s=deadline_s,
            )
        except Exception as exc:
            raise RuntimeError(f"VizieR fallback query failed: {_exc_brief(exc)}") from exc
        df = _normalise_gaia_columns(tab.to_pandas())
        if "phot_g_mean_mag" in df.columns and np.isfinite(mag_max):
            g = pd.to_numeric(df["phot_g_mean_mag"], errors="coerce")
            df = df[g.notna() & (g <= float(mag_max))]
        return df

    def _query_gaia(self, center: SkyCoord, radius_deg: float, mag_max: float) -> pd.DataFrame:
        _get_gaia()  # lazy import on first actual use
        if not _HAS_GAIA or Gaia is None:
            raise RuntimeError("astroquery.gaia not available")
        if self._stopped():
            raise RuntimeError("stopped")
        timeout_s = float(getattr(self.P, "gaia_timeout_s", 30.0))
        mag_where = f"AND phot_g_mean_mag <= {mag_max:.4f}" if np.isfinite(mag_max) and mag_max > 0 else ""
        # Both SELECT lists are generated from apex.utils.gaia_columns so they
        # cannot drift apart again — see that module for what drifting cost.
        adql = f"""
SELECT {esa_select_clause()}
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', {center.ra.deg:.8f}, {center.dec.deg:.8f}, {radius_deg:.8f})
)
{mag_where}
        """.strip()
        Gaia.ROW_LIMIT = -1
        if self._stopped():
            raise RuntimeError("stopped")
        deadline_s = float(getattr(self.P, "gaia_hard_deadline_s", 0.0) or 0.0) or None
        _t0 = time.time()
        try:
            tab = tap_query_with_deadline(
                Gaia, adql, timeout_s=timeout_s, deadline_s=deadline_s, dump_to_file=False,
            )
        except Exception as exc:
            cause = self._classify_query_failure(exc)
            # Log the exception verbatim, and how long it took to get it. The
            # label alone is not enough to act on: `_classify_query_failure`
            # tests for "timeout" before "connection", so "Connection timed
            # out" — a network fault — is also reported as TIMEOUT. And when
            # the VizieR fallback succeeds (the common case) the original
            # exception was previously discarded entirely, so a run that fell
            # back eight times left no evidence of why.
            #
            # The elapsed time separates the two cases on sight: a hard
            # deadline fires at ~max(2*timeout_s, 60) s, while a refused or
            # unreachable connection returns in well under a second.
            _dt = time.time() - _t0
            self._log(
                f"[Gaia][WARN] ESA TAP failed [{cause}] after {_dt:.1f}s "
                f"(timeout_s={timeout_s:.0f}, hard_deadline_s="
                f"{deadline_s if deadline_s else max(timeout_s * 2.0, 60.0):.0f}): "
                f"{_exc_brief(exc)}"
            )
            if cause in {"IP_BANNED", "SERVER_JOB_LOST", "SERVER_DOWN", "TIMEOUT", "NETWORK_ERROR"}:
                self._log(f"[Gaia][WARN] ESA TAP failed [{cause}], trying VizieR fallback.")
                _t1 = time.time()
                try:
                    df_viz = self._query_gaia_vizier(center, radius_deg, mag_max, timeout_s)
                    self._log(
                        f"[Gaia][WARN] ESA TAP failed [{cause}] -> VizieR fallback used "
                        f"(N={len(df_viz)}, mag_max={mag_max:.2f}, {time.time() - _t1:.1f}s)."
                    )
                    df_viz.attrs["gaia_source"] = "vizier_fallback"
                    return df_viz
                except Exception as exc2:
                    self._log(
                        f"[Gaia][WARN] VizieR fallback failed after "
                        f"{time.time() - _t1:.1f}s [{gaia_failure_reason(exc2)}]: "
                        f"{_exc_brief(exc2)}"
                    )
                    raise RuntimeError(
                        f"Gaia TAP async query failed [{cause}]: {_exc_brief(exc)}; "
                        f"VizieR fallback also failed: {_exc_brief(exc2)}"
                    ) from exc
            raise RuntimeError(f"Gaia TAP async query failed [{cause}]: {_exc_brief(exc)}") from exc
        if "phot_g_mean_mag" in tab.colnames and np.isfinite(mag_max):
            tab = tab[np.isfinite(tab["phot_g_mean_mag"]) & (tab["phot_g_mean_mag"] <= mag_max)]
        return _normalise_gaia_columns(tab.to_pandas())

    def _classify_query_failure(self, exc: Exception) -> str:
        """Label a query failure, most specific cause first.

        Ordering is the whole content of this function, and two of the tests
        it now passes used to fail:

        * SSL was checked last, after the generic "connection" test. An expired
          or untrusted certificate reports through a connection error, so a
          fixable local trust-store problem was filed as NETWORK_ERROR — the
          one label that suggests waiting rather than acting.
        * A bare "timed out" was checked before any connection test, so
          "Connection timed out" — a host that cannot be reached — read as a
          slow server.

        TIMEOUT and NETWORK_ERROR both fall through to the same VizieR
        fallback, so this changes what gets written down, not what runs.
        """
        err_str = str(exc).lower()
        if ssl_error_reason(exc):
            return "SSL_CERTIFICATE_VERIFY_FAILED"
        if "ip" in err_str and any(w in err_str for w in ("disabled", "blocked", "banned", "heavy")):
            return "IP_BANNED"
        if "404" in err_str or "job not found" in err_str:
            return "SERVER_JOB_LOST"
        if any(c in err_str for c in ("503", "502", "500")):
            return "SERVER_DOWN"
        if _OWN_DEADLINE_MARKER in err_str:
            return "TIMEOUT"
        if any(p in err_str for p in _CONNECTION_FAULT_PHRASES):
            return "NETWORK_ERROR"
        if "timeout" in err_str or "timed out" in err_str:
            return "TIMEOUT"
        if any(w in err_str for w in ("connection", "refused")):
            return "NETWORK_ERROR"
        return "UNKNOWN"

    def _cache_mag_max(self, df: pd.DataFrame, meta: dict | None) -> float:
        try:
            if isinstance(meta, dict) and "mag_max" in meta:
                value = float(meta.get("mag_max"))
                if np.isfinite(value):
                    return value
        except Exception:
            pass
        try:
            if "phot_g_mean_mag" in df.columns:
                g = pd.to_numeric(df["phot_g_mean_mag"], errors="coerce")
                if g.notna().any():
                    return float(g.max())
        except Exception:
            pass
        return np.nan

    def _filter_cache_by_mag(self, df: pd.DataFrame, mag_max: float) -> pd.DataFrame:
        if not np.isfinite(mag_max) or "phot_g_mean_mag" not in df.columns:
            return df
        g = pd.to_numeric(df["phot_g_mean_mag"], errors="coerce")
        keep = g.notna() & (g <= float(mag_max))
        return df.loc[keep].copy()

    def _cache_covers_field(self, df: pd.DataFrame, center: SkyCoord, radius_deg: float) -> bool:
        try:
            ra = pd.to_numeric(df.get("ra"), errors="coerce")
            dec = pd.to_numeric(df.get("dec"), errors="coerce")
        except Exception:
            return False
        valid = ra.notna() & dec.notna()
        n_valid = int(valid.sum())
        if n_valid <= 0:
            return False
        if n_valid < 50:
            return True
        ra_values = ra[valid].to_numpy(float)
        dec_values = dec[valid].to_numpy(float)
        cos_dec = float(np.cos(np.deg2rad(float(center.dec.deg))))
        if not np.isfinite(cos_dec) or cos_dec <= 0:
            cos_dec = 1.0
        dx = (ra_values - float(center.ra.deg)) * cos_dec
        dy = dec_values - float(center.dec.deg)
        if dx.size == 0 or dy.size == 0:
            return False
        side_frac = 0.60 if n_valid >= 200 else 0.45
        need = float(radius_deg) * side_frac
        min_x = float(np.nanmin(dx))
        max_x = float(np.nanmax(dx))
        min_y = float(np.nanmin(dy))
        max_y = float(np.nanmax(dy))
        return bool(
            np.isfinite(min_x) and np.isfinite(max_x)
            and np.isfinite(min_y) and np.isfinite(max_y)
            and min_x <= -need and max_x >= need
            and min_y <= -need and max_y >= need
        )

    def _mag_limit(self) -> float:
        mag_max_cfg = float(getattr(self.P, "gaia_mag_max", 18.0))
        wcs_mag_max = float(getattr(self.P, "gaia_wcs_mag_max", 18.0))
        if np.isfinite(mag_max_cfg) and np.isfinite(wcs_mag_max) and wcs_mag_max > 0:
            return float(min(mag_max_cfg, wcs_mag_max))
        return mag_max_cfg

    def load_or_query(self, center: SkyCoord, radius_deg: float) -> tuple[pd.DataFrame, str]:
        # Probe astroquery lazily here too — the cache branch below reads
        # _HAS_GAIA before any _query_gaia() call would set it.
        _get_gaia()
        step5_out = step5_wcs_dir(self.result_dir)
        step5_out.mkdir(parents=True, exist_ok=True)
        cache_path = step5_out / "gaia_fov.ecsv"
        meta_path = step5_out / "gaia_fov_meta.json"
        retry = int(getattr(self.P, "gaia_retry", 2))
        backoff_s = float(getattr(self.P, "gaia_backoff_s", 6.0))
        mag_max = self._mag_limit()
        allow_no_cache = bool(getattr(self.P, "gaia_allow_no_cache", True))

        df_cache = self._load_gaia_cache_if_ok(cache_path)
        meta = None
        if df_cache is not None and _HAS_GAIA:
            stale = list(df_cache.attrs.get("missing_contract_columns") or ())
            if stale:
                # Say why. Dropping a cache costs a network round trip, and a
                # silent one is indistinguishable from a slow server.
                self._log(
                    f"[Gaia] cached catalog lacks {', '.join(stale)} — re-querying "
                    f"so the same quality cuts apply as on a fresh catalog"
                )
                df_cache = None
        if df_cache is not None:
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text(encoding="utf-8"))
                except Exception:
                    meta = None
            cached_mag_max = self._cache_mag_max(df_cache, meta)
            mag_ok = (
                not np.isfinite(mag_max)
                or not np.isfinite(cached_mag_max)
                or cached_mag_max + 1e-6 >= mag_max
            )
            if isinstance(meta, dict):
                try:
                    cached_ra = float(meta.get("center_ra_deg", 0))
                    cached_dec = float(meta.get("center_dec_deg", 0))
                    dist_deg = np.hypot(center.ra.deg - cached_ra, center.dec.deg - cached_dec)
                    # Only require the cache to be centred on the same field.
                    # The cached radius is NOT compared here: different solvers
                    # (ASTAP vs the Internal quad solver) estimate the Gaia query
                    # radius slightly differently (FOV diagonal × fudge), so a
                    # strict ``cached_radius >= radius_deg*0.9`` test made the
                    # Internal solver reject ASTAP's freshly-written cache and
                    # re-query — which could hang on a slow ESA TAP async job.
                    # Whether the cache actually spans the requested field is
                    # verified independently by ``_cache_covers_field`` below,
                    # using the real star distribution rather than the recorded
                    # radius, so dropping the radius gate here loses no safety.
                    same_field = bool(dist_deg < 0.03)
                except Exception:
                    same_field = False
            else:
                same_field = True
            coverage_ok = self._cache_covers_field(df_cache, center, radius_deg)
            if same_field and mag_ok and coverage_ok:
                return self._filter_cache_by_mag(df_cache, mag_max), "cache"

        if not _HAS_GAIA:
            if allow_no_cache:
                return pd.DataFrame(), "no_gaia_module"
            raise RuntimeError("astroquery.gaia not available and no cache")

        last_err = None
        for attempt in range(1, max(1, retry) + 1):
            if self._stopped():
                raise RuntimeError("stopped")
            try:
                df = self._query_gaia(center, radius_deg, mag_max)
                source = df.attrs.get("gaia_source", "esa") if hasattr(df, "attrs") else "esa"
                df = _normalise_gaia_columns(df)
                df.attrs["gaia_source"] = source
                try:
                    Table.from_pandas(df).write(cache_path, format="ascii.ecsv", overwrite=True)
                    meta_path.write_text(
                        json.dumps(
                            {
                                "center_ra_deg": float(center.ra.deg),
                                "center_dec_deg": float(center.dec.deg),
                                "radius_deg": float(radius_deg),
                                "mag_max": float(mag_max),
                                "n_stars": int(len(df)),
                                "gaia_source": source,
                            },
                            indent=2,
                        ),
                        encoding="utf-8",
                    )
                except Exception:
                    pass
                return df, source
            except Exception as exc:
                if self._stopped():
                    raise RuntimeError("stopped")
                last_err = exc
                if attempt < retry:
                    slept = 0.0
                    while slept < backoff_s:
                        if self._stopped():
                            raise RuntimeError("stopped")
                        dt = min(0.25, backoff_s - slept)
                        time.sleep(dt)
                        slept += dt

        df_cache = self._load_gaia_cache_if_ok(cache_path)
        if df_cache is not None:
            cached_mag_max = self._cache_mag_max(df_cache, None)
            mag_ok = (
                not np.isfinite(mag_max)
                or not np.isfinite(cached_mag_max)
                or cached_mag_max + 1e-6 >= mag_max
            )
            if mag_ok:
                return self._filter_cache_by_mag(df_cache, mag_max), "cache(after_fail)"
        if allow_no_cache:
            return pd.DataFrame(), f"fail_no_cache:{gaia_failure_reason(last_err)}"
        raise RuntimeError(f"Gaia query failed: {last_err}")
