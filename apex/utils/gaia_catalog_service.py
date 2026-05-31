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


def gaia_failure_reason(exc: Exception | None) -> str:
    if exc is None:
        return "unknown"
    ssl_reason = ssl_error_reason(exc)
    if ssl_reason:
        return ssl_reason
    text = _exc_brief(exc, limit=180)
    low = text.lower()
    if "timeout" in low or "timed out" in low:
        return "timeout"
    if "connection" in low or "unreachable" in low or "refused" in low:
        return "network_error"
    if "server_down" in low or "503" in low or "502" in low or "500" in low:
        return "server_down"
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
    out = df.copy()
    rename: dict[str, str] = {}
    for col in out.columns:
        key = str(col).strip().lower()
        if key in {"source", "source_id"}:
            rename[col] = "source_id"
        elif key in {"ra", "ra_deg", "ra_icrs"}:
            rename[col] = "ra"
        elif key in {"dec", "dec_deg", "de_icrs"}:
            rename[col] = "dec"
        elif key in {"gmag", "phot_g_mean_mag"}:
            rename[col] = "phot_g_mean_mag"
        elif key in {"bpmag", "phot_bp_mean_mag"}:
            rename[col] = "phot_bp_mean_mag"
        elif key in {"rpmag", "phot_rp_mean_mag"}:
            rename[col] = "phot_rp_mean_mag"
        elif key in {"ruwe"}:
            rename[col] = "ruwe"
        elif key in {"plx", "parallax"}:
            rename[col] = "parallax"
        elif key in {"e_plx", "parallax_error"}:
            rename[col] = "parallax_error"
        elif key in {"pmra"}:
            rename[col] = "pmra"
        elif key in {"e_pmra", "pmra_error"}:
            rename[col] = "pmra_error"
        elif key in {"pmde", "pmdec"}:
            rename[col] = "pmdec"
        elif key in {"e_pmde", "pmdec_error"}:
            rename[col] = "pmdec_error"
        elif key in {"phot_variable_flag"}:
            rename[col] = "phot_variable_flag"
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
            raw_cols = [str(c).strip().lower() for c in tab.colnames]
            missing_var_flag = "phot_variable_flag" not in raw_cols
            df = _normalise_gaia_columns(tab.to_pandas())
            if not {"ra", "dec"} <= set(df.columns):
                return None
            if "source_id" not in df.columns and _HAS_GAIA:
                return None
            df.attrs["missing_phot_variable_flag"] = bool(missing_var_flag)
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
  "I/355/gaiadr3".Source AS source_id,
  "I/355/gaiadr3".RA_ICRS AS ra,
  "I/355/gaiadr3".DE_ICRS AS dec,
  "I/355/gaiadr3".Gmag AS phot_g_mean_mag,
  "I/355/gaiadr3".BPmag AS phot_bp_mean_mag,
  "I/355/gaiadr3".RPmag AS phot_rp_mean_mag,
  "I/355/gaiadr3".RUWE AS ruwe,
  "I/355/gaiadr3".Plx AS parallax,
  "I/355/gaiadr3".e_Plx AS parallax_error,
  "I/355/gaiadr3".pmRA AS pmra,
  "I/355/gaiadr3".e_pmRA AS pmra_error,
  "I/355/gaiadr3".pmDE AS pmdec,
  "I/355/gaiadr3".e_pmDE AS pmdec_error
FROM "I/355/gaiadr3"
WHERE 1=CONTAINS(
    POINT('ICRS', "I/355/gaiadr3".RA_ICRS, "I/355/gaiadr3".DE_ICRS),
    CIRCLE('ICRS', {center.ra.deg:.8f}, {center.dec.deg:.8f}, {radius_deg:.8f})
)
{mag_where}
        """.strip()
        tap = TapPlus(url="https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
        try:
            job = tap_launch_job_async(tap, adql, timeout_s=timeout_s)
            tab = job.get_results()
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
        adql = f"""
SELECT
  source_id, ra, dec,
  phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
  phot_variable_flag,
  pmra, pmdec, pmra_error, pmdec_error,
  parallax, parallax_error
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
        try:
            job = tap_launch_job_async(Gaia, adql, timeout_s=timeout_s, dump_to_file=False)
            tab = job.get_results()
        except Exception as exc:
            cause = self._classify_query_failure(exc)
            if cause in {"IP_BANNED", "SERVER_JOB_LOST", "SERVER_DOWN", "TIMEOUT", "NETWORK_ERROR"}:
                self._log(f"[Gaia][WARN] ESA TAP failed [{cause}], trying VizieR fallback.")
                try:
                    df_viz = self._query_gaia_vizier(center, radius_deg, mag_max, timeout_s)
                    self._log(
                        f"[Gaia][WARN] ESA TAP failed [{cause}] -> VizieR fallback used "
                        f"(N={len(df_viz)}, mag_max={mag_max:.2f})."
                    )
                    df_viz.attrs["gaia_source"] = "vizier_fallback"
                    return df_viz
                except Exception as exc2:
                    self._log(f"[Gaia][WARN] VizieR fallback failed: {gaia_failure_reason(exc2)}")
                    raise RuntimeError(
                        f"Gaia TAP async query failed [{cause}]: {_exc_brief(exc)}; "
                        f"VizieR fallback also failed: {_exc_brief(exc2)}"
                    ) from exc
            raise RuntimeError(f"Gaia TAP async query failed [{cause}]: {_exc_brief(exc)}") from exc
        if "phot_g_mean_mag" in tab.colnames and np.isfinite(mag_max):
            tab = tab[np.isfinite(tab["phot_g_mean_mag"]) & (tab["phot_g_mean_mag"] <= mag_max)]
        return _normalise_gaia_columns(tab.to_pandas())

    def _classify_query_failure(self, exc: Exception) -> str:
        err_str = str(exc).lower()
        if "ip" in err_str and any(w in err_str for w in ("disabled", "blocked", "banned", "heavy")):
            return "IP_BANNED"
        if "404" in err_str or "job not found" in err_str:
            return "SERVER_JOB_LOST"
        if any(c in err_str for c in ("503", "502", "500")):
            return "SERVER_DOWN"
        if "timeout" in err_str or "timed out" in err_str:
            return "TIMEOUT"
        if any(w in err_str for w in ("connection", "refused", "unreachable")):
            return "NETWORK_ERROR"
        if ssl_error_reason(exc):
            return "SSL_CERTIFICATE_VERIFY_FAILED"
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
        if df_cache is not None and bool(df_cache.attrs.get("missing_phot_variable_flag", False)) and _HAS_GAIA:
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
