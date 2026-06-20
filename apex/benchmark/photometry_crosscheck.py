"""Photometric cross-validation harness for APEX Step 7 forced photometry.

This module proves — repeatably — that APEX aperture photometry agrees with
INDEPENDENT photometry engines on the *same stars in the same image*. It is the
measurement-level counterpart to the artificial-star and isochrone-recovery
suites: where those test selection + recovery against an injected truth, this
re-measures APEX's own catalog with a different tool and reports the agreement.

It is wired into ``apex validate --suite crosscheck`` (see
``apex.benchmark.validate.run_crosscheck_suite``).

Engines
-------
* ``sep`` (Bertin SExtractor backend, Barbary's ``sep``): a fully independent
  C re-implementation of aperture photometry. ``sep.sum_circle`` with
  ``bkgann=(r_in, r_out)`` does its OWN local-annulus sky subtraction, exactly
  like APEX, but shares NO code with photutils — this is the strong test.
* ``iraf``: thin wrapper around the existing
  :func:`apex.benchmark.iraf_crosscheck.run_iraf_crosscheck`. Skips cleanly when
  IRAF/PyRAF is unavailable.
* ``photutils``: a re-measurement with photutils ``CircularAperture`` +
  ``CircularAnnulus`` ``ApertureStats``. NOTE: APEX itself uses photutils, so
  this is only *semi*-independent (same engine, different call path); it is a
  consistency reference, not proof of cross-tool agreement.

The comparison is differential and zeropoint-free: APEX and the reference
instrumental magnitudes are ``-2.5 log10(flux)``; their median difference is the
zeropoint offset (which absorbs gain, exposure time, and aperture correction —
all constant per frame), and the residual scatter about that offset is the real
agreement metric.

LAYER RULE: this module is Qt-free. It must NOT import apex.gui or PyQt5. The
heavy science imports (sep, photutils, astropy.io.fits) happen lazily inside the
functions so ``import apex.benchmark.photometry_crosscheck`` stays cheap.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

# ── tunables ──────────────────────────────────────────────────────────────────

# APEX forced-photometry aperture/sky scales (FWHM units) and floors — these are
# the production defaults from apex.config.schema.ApertureScalesConfig /
# ApertureRadiiConfig and apex.analysis.forced_photometry._phot_frame. They are
# overridden at runtime by the actual params.P attributes when available.
DEFAULT_APERTURE_SCALE = 0.8     # forced_r_ap_scale
DEFAULT_REF_AP_SCALE = 2.4       # forced_ref_ap_scale
DEFAULT_ANN_IN_SCALE = 4.0       # fitsky_annulus_scale
DEFAULT_DANN_SCALE = 2.0         # fitsky_dannulus_scale
DEFAULT_MIN_R_AP_PX = 4.0        # min_r_ap_px
DEFAULT_ANN_MIN_GAP_PX = 6.0     # annulus_min_gap_px
DEFAULT_GAIN = 1.5

REFERENCE_CHOICES = ("sep", "iraf", "photutils")


# ── small numeric helpers (self-contained; no cross-module coupling) ───────────

def _finite(arr: Iterable[Any]) -> np.ndarray:
    out = np.asarray(list(arr), dtype=float)
    return out[np.isfinite(out)]


def _mad_sigma(values: Iterable[Any]) -> float:
    arr = _finite(values)
    if arr.size == 0:
        return float("nan")
    med = float(np.median(arr))
    return float(1.4826 * np.median(np.abs(arr - med)))


def _rms(values: Iterable[Any]) -> float:
    arr = _finite(values)
    if arr.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(arr * arr)))


def _robust_sigma(values: Iterable[Any]) -> float:
    """68.27 percentile half-width of |x - median| — robust 1-sigma estimate."""
    arr = _finite(values)
    if arr.size == 0:
        return float("nan")
    med = float(np.median(arr))
    return float(np.percentile(np.abs(arr - med), 68.27))


def _pearson(x: Iterable[Any], y: Iterable[Any]) -> float:
    xa = np.asarray(list(x), dtype=float)
    ya = np.asarray(list(y), dtype=float)
    good = np.isfinite(xa) & np.isfinite(ya)
    if int(good.sum()) < 3:
        return float("nan")
    xa = xa[good]
    ya = ya[good]
    if np.std(xa) == 0 or np.std(ya) == 0:
        return float("nan")
    return float(np.corrcoef(xa, ya)[0, 1])


def _finite_float(value: Any, default: float | None = None) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def _attr(params: Any, name: str, default: float) -> float:
    """Read a forced-photometry scale from params.P (falls back to default)."""
    if params is None:
        return float(default)
    holder = getattr(params, "P", params)
    val = _finite_float(getattr(holder, name, None))
    return float(val) if val is not None and val > 0 else float(default)


def _resolve_scales(params: Any) -> dict[str, float]:
    """Resolve APEX aperture/sky scales from params.P, matching _phot_frame."""
    return {
        "aperture_scale": _attr(params, "forced_r_ap_scale", DEFAULT_APERTURE_SCALE),
        "ref_ap_scale": _attr(params, "forced_ref_ap_scale", DEFAULT_REF_AP_SCALE),
        "ann_in_scale": _attr(params, "fitsky_annulus_scale", DEFAULT_ANN_IN_SCALE),
        "dann_scale": _attr(params, "fitsky_dannulus_scale", DEFAULT_DANN_SCALE),
        "min_r_ap_px": _attr(params, "min_r_ap_px", DEFAULT_MIN_R_AP_PX),
        "ann_min_gap_px": _attr(params, "annulus_min_gap_px", DEFAULT_ANN_MIN_GAP_PX),
        "gain": _attr(params, "gain_e_per_adu", DEFAULT_GAIN),
    }


def _radii(fwhm_px: float, scales: dict[str, float]) -> tuple[float, float, float]:
    """APEX aperture / sky-annulus radii (px), matching forced_photometry._phot_frame.

    r_ap  = max(min_r_ap, aperture_scale * fwhm)
    r_ref = max(r_ap + 2, ref_ap_scale * fwhm)
    r_in  = max(r_ref + ann_gap, ann_in_scale * fwhm)
    r_out = r_in + max(ann_gap, dann_scale * fwhm)
    """
    fwhm = float(fwhm_px)
    r_ap = max(scales["min_r_ap_px"], scales["aperture_scale"] * fwhm)
    r_ref = max(r_ap + 2.0, scales["ref_ap_scale"] * fwhm)
    r_in = max(r_ref + scales["ann_min_gap_px"], scales["ann_in_scale"] * fwhm)
    r_out = r_in + max(scales["ann_min_gap_px"], scales["dann_scale"] * fwhm)
    return float(r_ap), float(r_in), float(r_out)


# ── Step 7 contract loaders ────────────────────────────────────────────────────

def _step7_dir(result_dir: Path) -> Path:
    return Path(result_dir) / "step7_forced_phot"


def _read_index(result_dir: Path) -> pd.DataFrame:
    idx = _step7_dir(result_dir) / "photometry_index.csv"
    if not idx.is_file():
        raise FileNotFoundError(f"Step 7 index missing: {idx}")
    return pd.read_csv(idx)


def load_apex_tsv(tsv_path: str | Path, *, snr_min: float | None = None) -> pd.DataFrame:
    """Load an APEX Step 7 photometry TSV, keeping the columns we compare.

    Requires x, y, flux_net_adu. mag_inst/snr are kept when present. Drops
    saturated / non-linear / bad-phot / off-frame rows so the comparison runs on
    clean photometry only.
    """
    path = Path(tsv_path)
    if not path.is_file():
        raise FileNotFoundError(f"Step 7 TSV missing: {path}")
    df = pd.read_csv(path, sep="\t")
    required = {"x", "y", "flux_net_adu"}
    if not required <= set(df.columns):
        missing = sorted(required - set(df.columns))
        raise ValueError(f"Step 7 TSV missing columns: {', '.join(missing)} ({path.name})")

    out = pd.DataFrame()
    out["x"] = pd.to_numeric(df["x"], errors="coerce")
    out["y"] = pd.to_numeric(df["y"], errors="coerce")
    out["flux_net_adu"] = pd.to_numeric(df["flux_net_adu"], errors="coerce")
    out["mag_inst"] = pd.to_numeric(df["mag_inst"], errors="coerce") if "mag_inst" in df.columns else np.nan
    out["snr"] = pd.to_numeric(df["snr"], errors="coerce") if "snr" in df.columns else np.nan
    out["apcorr"] = pd.to_numeric(df["apcorr"], errors="coerce") if "apcorr" in df.columns else 1.0

    good = (
        np.isfinite(out["x"].to_numpy(float))
        & np.isfinite(out["y"].to_numpy(float))
        & np.isfinite(out["flux_net_adu"].to_numpy(float))
        & (out["flux_net_adu"].to_numpy(float) > 0)
    )
    for flag in ("is_saturated", "is_nonlinear", "bad_phot_flag", "off_frame_flag", "centroid_outlier"):
        if flag in df.columns:
            f = df[flag]
            if f.dtype == bool:
                mask = f.fillna(False).to_numpy(bool)
            else:
                mask = f.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "on"}).to_numpy(bool)
            good &= ~mask
    if snr_min is not None and math.isfinite(float(snr_min)):
        good &= out["snr"].to_numpy(float) >= float(snr_min)
    out = out.loc[good].reset_index(drop=True)
    out.insert(0, "row_id", np.arange(len(out), dtype=int))
    return out


# ── image resolution ───────────────────────────────────────────────────────────

def _data_dir_from_params(result_dir: Path, params: Any) -> Path | None:
    """Best-effort science-frame directory.

    Prefers params.P.data_dir; otherwise falls back to the result-dir parent
    (APEX convention: result_dir = data_dir / 'result').
    """
    if params is not None:
        holder = getattr(params, "P", params)
        data_dir = getattr(holder, "data_dir", None)
        if data_dir:
            cand = Path(data_dir)
            if cand.is_dir():
                return cand
    parent = Path(result_dir).parent
    return parent if parent.is_dir() else None


def resolve_frame_image(result_dir: str | Path, frame_name: str, params: Any = None) -> Path:
    """Resolve the SAME pixel image APEX measured for ``frame_name``.

    The Step 7 (x, y) live in this image's pixel frame, so getting this wrong
    silently produces a meaningless comparison. Resolution order:

    1. If a crop is active (``step2_crop/crop_rect.json`` exists) and the cropped
       frame is on disk, use ``step2_crop/cropped/<frame>``.
    2. Otherwise use the raw science frame in ``params.P.data_dir`` (or the
       result-dir parent, the APEX default layout).

    Raises FileNotFoundError when no candidate exists.
    """
    result_dir = Path(result_dir)
    name = Path(frame_name).name

    crop_rect = result_dir / "step2_crop" / "crop_rect.json"
    if crop_rect.exists():
        cropped = result_dir / "step2_crop" / "cropped" / name
        if cropped.exists():
            return cropped

    data_dir = _data_dir_from_params(result_dir, params)
    candidates: list[Path] = []
    if data_dir is not None:
        candidates.append(data_dir / name)
    # Also try the raw frame next to the result dir even if data_dir resolved
    # elsewhere (covers configs where data_dir points at a stale path).
    candidates.append(result_dir.parent / name)
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError(
        f"Could not resolve science image for '{frame_name}'. Tried: "
        + ", ".join(str(c) for c in candidates)
    )


def _load_image(path: Path) -> np.ndarray:
    from astropy.io import fits

    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                # sep requires native byte order, C-contiguous float
                data = np.ascontiguousarray(hdu.data.astype(np.float64))
                return data
    raise ValueError(f"No 2-D image data in {path}")


def _sanity_check_positions(
    data: np.ndarray, x: np.ndarray, y: np.ndarray
) -> dict[str, Any]:
    """Confirm the (x, y) land on real sources: peak pixel >> background.

    Samples the nearest-pixel value at each star and compares the median of
    those to a global background level (sigma-clipped). If the stars are NOT
    well above sky, the image/coordinate frame is almost certainly mismatched.
    """
    h, w = data.shape
    xi = np.clip(np.round(x).astype(int), 0, w - 1)
    yi = np.clip(np.round(y).astype(int), 0, h - 1)
    on_frame = (x >= 0) & (x < w) & (y >= 0) & (y < h)
    star_vals = data[yi[on_frame], xi[on_frame]] if int(on_frame.sum()) else np.asarray([])

    flat = data[np.isfinite(data)]
    if flat.size > 200000:
        rng = np.random.default_rng(0)
        flat = rng.choice(flat, size=200000, replace=False)
    med = float(np.median(flat)) if flat.size else float("nan")
    mad = float(1.4826 * np.median(np.abs(flat - med))) if flat.size else float("nan")
    star_med = float(np.median(star_vals)) if star_vals.size else float("nan")
    n_sigma = float((star_med - med) / mad) if (math.isfinite(mad) and mad > 0) else float("nan")
    # A genuine star catalog should sit many sigma above sky at the peak pixel.
    plausible = bool(math.isfinite(n_sigma) and n_sigma >= 3.0)
    return {
        "image_shape": [int(h), int(w)],
        "n_on_frame": int(on_frame.sum()),
        "n_total": int(x.size),
        "frac_on_frame": float(on_frame.mean()) if x.size else float("nan"),
        "sky_median_adu": med,
        "sky_mad_adu": mad,
        "star_pixel_median_adu": star_med,
        "star_over_sky_sigma": n_sigma,
        "positions_plausible": plausible,
    }


# ── statistics + plots shared across engines ───────────────────────────────────

def _agreement_stats(
    apex_mag: np.ndarray,
    ref_mag: np.ndarray,
    *,
    snr: np.ndarray | None = None,
    snr_min: float = 20.0,
) -> dict[str, Any]:
    """Δmag agreement statistics on the SNR>snr_min subset.

    Δmag = apex_mag - ref_mag. The median is the zeropoint offset (absorbs the
    constant gain/exptime/apcorr difference between the two flux scales); the
    scatter about it is the real agreement metric.
    """
    apex_mag = np.asarray(apex_mag, dtype=float)
    ref_mag = np.asarray(ref_mag, dtype=float)
    good = np.isfinite(apex_mag) & np.isfinite(ref_mag)
    if snr is not None and math.isfinite(float(snr_min)):
        snr = np.asarray(snr, dtype=float)
        good &= np.isfinite(snr) & (snr >= float(snr_min))
    n = int(good.sum())
    if n == 0:
        return {"n_matched": 0}

    am = apex_mag[good]
    rm = ref_mag[good]
    dmag = am - rm
    offset = float(np.median(dmag))
    resid = dmag - offset  # zeropoint-aligned residual
    abs_resid = np.abs(resid)
    return {
        "n_matched": n,
        "dmag_median": offset,
        "dmag_mad": _mad_sigma(dmag),
        "dmag_rms": _rms(resid),
        "robust_sigma": _robust_sigma(dmag),
        "pearson_r": _pearson(am, rm),
        "frac_within_0p05": float(np.mean(abs_resid <= 0.05)),
        "frac_within_0p1": float(np.mean(abs_resid <= 0.10)),
        "mag_min": float(np.min(am)),
        "mag_max": float(np.max(am)),
    }


def _binned_residual(
    apex_mag: np.ndarray,
    dmag_aligned: np.ndarray,
    *,
    n_bins: int = 8,
) -> list[dict[str, Any]]:
    """Median Δmag (zeropoint-aligned) vs APEX magnitude, in equal-count bins."""
    am = np.asarray(apex_mag, dtype=float)
    dm = np.asarray(dmag_aligned, dtype=float)
    good = np.isfinite(am) & np.isfinite(dm)
    am = am[good]
    dm = dm[good]
    if am.size < n_bins * 2:
        return []
    edges = np.quantile(am, np.linspace(0, 1, n_bins + 1))
    edges = np.unique(edges)
    rows: list[dict[str, Any]] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (am >= lo) & (am <= hi)
        if int(sel.sum()) < 2:
            continue
        rows.append(
            {
                "mag_center": float(0.5 * (lo + hi)),
                "n": int(sel.sum()),
                "dmag_median": float(np.median(dm[sel])),
                "dmag_mad": _mad_sigma(dm[sel]),
            }
        )
    return rows


def _write_scatter_plot(
    apex_mag: np.ndarray,
    ref_mag: np.ndarray,
    *,
    stats: dict[str, Any],
    reference: str,
    out_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    good = np.isfinite(apex_mag) & np.isfinite(ref_mag)
    am = apex_mag[good]
    rm = ref_mag[good] + float(stats.get("dmag_median", 0.0))  # align ZP to APEX
    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    if am.size:
        lo = float(min(am.min(), rm.min()))
        hi = float(max(am.max(), rm.max()))
        pad = max(0.05, 0.04 * (hi - lo))
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="black", linewidth=1)
        ax.scatter(am, rm, s=10, alpha=0.55)
        ax.set(xlim=(lo - pad, hi + pad), ylim=(lo - pad, hi + pad))
        ax.invert_xaxis()
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.text(
            0.04, 0.96,
            f"N={stats.get('n_matched')}\n"
            f"offset={stats.get('dmag_median'):.4f}\n"
            f"MAD={stats.get('dmag_mad'):.4f}\n"
            f"r={stats.get('pearson_r'):.5f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.75"},
        )
    ax.set(
        xlabel="APEX instrumental mag (-2.5 log10 flux_net_adu)",
        ylabel=f"{reference} instrumental mag (ZP-aligned to APEX)",
    )
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def _write_residual_plot(
    apex_mag: np.ndarray,
    dmag_aligned: np.ndarray,
    *,
    reference: str,
    out_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    good = np.isfinite(apex_mag) & np.isfinite(dmag_aligned)
    am = apex_mag[good]
    dm = dmag_aligned[good]
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.scatter(am, dm, s=10, alpha=0.5)
    ax.axhline(0.0, color="black", linewidth=1)
    for level, color in ((0.05, "tab:orange"), (0.10, "tab:red")):
        ax.axhline(level, color=color, linewidth=0.8, linestyle="--", alpha=0.7)
        ax.axhline(-level, color=color, linewidth=0.8, linestyle="--", alpha=0.7)
    ax.set(
        xlabel="APEX instrumental mag",
        ylabel=f"APEX - {reference} (zeropoint removed, mag)",
    )
    if am.size:
        ax.invert_xaxis()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


# ── engine: sep ─────────────────────────────────────────────────────────────────

def crosscheck_sep(
    image_path: str | Path,
    step7_tsv: str | Path,
    *,
    fwhm_px: float,
    aperture_scale: float = DEFAULT_APERTURE_SCALE,
    ref_ap_scale: float = DEFAULT_REF_AP_SCALE,
    ann_in_scale: float = DEFAULT_ANN_IN_SCALE,
    dann_scale: float = DEFAULT_DANN_SCALE,
    min_r_ap_px: float = DEFAULT_MIN_R_AP_PX,
    ann_min_gap_px: float = DEFAULT_ANN_MIN_GAP_PX,
    gain: float = DEFAULT_GAIN,
    snr_min: float = 20.0,
    output_dir: str | Path | None = None,
    plot_prefix: str = "sep",
) -> dict[str, Any]:
    """Independently re-measure APEX Step 7 stars with ``sep`` and compare.

    ``sep.sum_circle(data, x, y, r_ap, bkgann=(r_in, r_out), gain=gain)`` performs
    its OWN local-annulus sky subtraction — sharing no code with photutils — so
    agreement here is a genuine cross-tool result. APEX (x, y) are 0-indexed
    pixel centers, matching sep's convention.

    Returns a stats dict with n_matched, dmag_median (zeropoint offset),
    dmag_mad, dmag_rms, robust_sigma, pearson_r, frac_within_0p05/0p1, binned
    residuals, a sanity-check block, and (when output_dir is set) plot paths.
    """
    import sep

    image_path = Path(image_path)
    data = _load_image(image_path)
    apex = load_apex_tsv(step7_tsv)  # keep all SNR for sanity-check; filter later

    scales = {
        "aperture_scale": float(aperture_scale),
        "ref_ap_scale": float(ref_ap_scale),
        "ann_in_scale": float(ann_in_scale),
        "dann_scale": float(dann_scale),
        "min_r_ap_px": float(min_r_ap_px),
        "ann_min_gap_px": float(ann_min_gap_px),
        "gain": float(gain),
    }
    r_ap, r_in, r_out = _radii(fwhm_px, scales)

    x = apex["x"].to_numpy(float)
    y = apex["y"].to_numpy(float)
    sanity = _sanity_check_positions(data, x, y)

    g = float(gain) if (math.isfinite(float(gain)) and float(gain) > 0) else None
    flux, fluxerr, flag = sep.sum_circle(
        data, x, y, float(r_ap), bkgann=(float(r_in), float(r_out)), gain=g
    )
    flux = np.asarray(flux, dtype=float)
    flag = np.asarray(flag, dtype=int)

    apex_flux = apex["flux_net_adu"].to_numpy(float)
    # Remove APEX's per-frame aperture correction so we compare raw aperture
    # sums; the apcorr is constant per frame, so this only re-centers the ZP
    # offset and does not change the scatter.
    apcorr = apex["apcorr"].to_numpy(float)
    apcorr = np.where(np.isfinite(apcorr) & (apcorr > 0), apcorr, 1.0)
    apex_aper_flux = apex_flux / apcorr

    with np.errstate(divide="ignore", invalid="ignore"):
        apex_mag = -2.5 * np.log10(np.where(apex_aper_flux > 0, apex_aper_flux, np.nan))
        sep_mag = -2.5 * np.log10(np.where(flux > 0, flux, np.nan))

    snr = apex["snr"].to_numpy(float)
    valid = (flag == 0) & np.isfinite(sep_mag) & np.isfinite(apex_mag)
    stats = _agreement_stats(
        np.where(valid, apex_mag, np.nan),
        np.where(valid, sep_mag, np.nan),
        snr=snr,
        snr_min=snr_min,
    )

    result: dict[str, Any] = {
        "status": "ok" if stats.get("n_matched", 0) > 0 else "error",
        "reference": "sep",
        "image": str(image_path),
        "n_input": int(len(apex)),
        "n_sep_flag_ok": int(np.sum(flag == 0)),
        "radii_px": {"r_ap": r_ap, "r_in": r_in, "r_out": r_out},
        "fwhm_px": float(fwhm_px),
        "gain": float(gain),
        "snr_min": float(snr_min),
        "sanity": sanity,
        **stats,
    }
    if stats.get("n_matched", 0) == 0:
        result["reason"] = "no valid matched sources after sep measurement / SNR cut"
        return result

    # Binned residual vs magnitude on the analysis subset.
    sub = valid & np.isfinite(snr) & (snr >= float(snr_min))
    am = apex_mag[sub]
    dm = (apex_mag[sub] - sep_mag[sub]) - float(stats["dmag_median"])
    result["binned_dmag_vs_mag"] = _binned_residual(am, dm)

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        scatter_path = out / f"{plot_prefix}_apex_vs_sep_mag.png"
        resid_path = out / f"{plot_prefix}_dmag_vs_mag.png"
        _write_scatter_plot(am, sep_mag[sub], stats=stats, reference="sep", out_path=scatter_path)
        _write_residual_plot(am, dm, reference="sep", out_path=resid_path)
        result["plots"] = [scatter_path.name, resid_path.name]
    return result


# ── engine: photutils (semi-independent) ────────────────────────────────────────

def crosscheck_photutils(
    image_path: str | Path,
    step7_tsv: str | Path,
    *,
    fwhm_px: float,
    aperture_scale: float = DEFAULT_APERTURE_SCALE,
    ref_ap_scale: float = DEFAULT_REF_AP_SCALE,
    ann_in_scale: float = DEFAULT_ANN_IN_SCALE,
    dann_scale: float = DEFAULT_DANN_SCALE,
    min_r_ap_px: float = DEFAULT_MIN_R_AP_PX,
    ann_min_gap_px: float = DEFAULT_ANN_MIN_GAP_PX,
    gain: float = DEFAULT_GAIN,
    snr_min: float = 20.0,
    output_dir: str | Path | None = None,
    plot_prefix: str = "photutils",
) -> dict[str, Any]:
    """Re-measure with photutils CircularAperture + CircularAnnulus.

    NOTE: APEX itself uses photutils, so this engine shares APEX's measurement
    backend. It is a *consistency* reference (different call path, same library),
    not proof of cross-tool agreement; treat ``sep`` / ``iraf`` as the
    independent tests. Sky is the sigma-clipped median in the annulus,
    subtracted per pixel of the aperture (median-sky aperture photometry).
    """
    from photutils.aperture import (
        CircularAperture,
        CircularAnnulus,
        ApertureStats,
        aperture_photometry,
    )
    from astropy.stats import SigmaClip

    image_path = Path(image_path)
    data = _load_image(image_path)
    apex = load_apex_tsv(step7_tsv)

    scales = {
        "aperture_scale": float(aperture_scale),
        "ref_ap_scale": float(ref_ap_scale),
        "ann_in_scale": float(ann_in_scale),
        "dann_scale": float(dann_scale),
        "min_r_ap_px": float(min_r_ap_px),
        "ann_min_gap_px": float(ann_min_gap_px),
        "gain": float(gain),
    }
    r_ap, r_in, r_out = _radii(fwhm_px, scales)

    x = apex["x"].to_numpy(float)
    y = apex["y"].to_numpy(float)
    sanity = _sanity_check_positions(data, x, y)
    positions = np.column_stack([x, y])

    aperture = CircularAperture(positions, r=float(r_ap))
    annulus = CircularAnnulus(positions, r_in=float(r_in), r_out=float(r_out))
    sigclip = SigmaClip(sigma=3.0, maxiters=5)
    bkg_stats = ApertureStats(data, annulus, sigma_clip=sigclip)
    bkg_median = np.asarray(bkg_stats.median, dtype=float)
    phot = aperture_photometry(data, aperture)
    raw_sum = np.asarray(phot["aperture_sum"], dtype=float)
    net_flux = raw_sum - bkg_median * float(aperture.area)

    apex_flux = apex["flux_net_adu"].to_numpy(float)
    apcorr = apex["apcorr"].to_numpy(float)
    apcorr = np.where(np.isfinite(apcorr) & (apcorr > 0), apcorr, 1.0)
    apex_aper_flux = apex_flux / apcorr

    with np.errstate(divide="ignore", invalid="ignore"):
        apex_mag = -2.5 * np.log10(np.where(apex_aper_flux > 0, apex_aper_flux, np.nan))
        pu_mag = -2.5 * np.log10(np.where(net_flux > 0, net_flux, np.nan))

    snr = apex["snr"].to_numpy(float)
    valid = np.isfinite(pu_mag) & np.isfinite(apex_mag)
    stats = _agreement_stats(
        np.where(valid, apex_mag, np.nan),
        np.where(valid, pu_mag, np.nan),
        snr=snr,
        snr_min=snr_min,
    )

    result: dict[str, Any] = {
        "status": "ok" if stats.get("n_matched", 0) > 0 else "error",
        "reference": "photutils",
        "semi_independent": True,
        "image": str(image_path),
        "n_input": int(len(apex)),
        "radii_px": {"r_ap": r_ap, "r_in": r_in, "r_out": r_out},
        "fwhm_px": float(fwhm_px),
        "gain": float(gain),
        "snr_min": float(snr_min),
        "sanity": sanity,
        **stats,
    }
    if stats.get("n_matched", 0) == 0:
        result["reason"] = "no valid matched sources after photutils measurement / SNR cut"
        return result

    sub = valid & np.isfinite(snr) & (snr >= float(snr_min))
    am = apex_mag[sub]
    dm = (apex_mag[sub] - pu_mag[sub]) - float(stats["dmag_median"])
    result["binned_dmag_vs_mag"] = _binned_residual(am, dm)

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        scatter_path = out / f"{plot_prefix}_apex_vs_photutils_mag.png"
        resid_path = out / f"{plot_prefix}_dmag_vs_mag.png"
        _write_scatter_plot(am, pu_mag[sub], stats=stats, reference="photutils", out_path=scatter_path)
        _write_residual_plot(am, dm, reference="photutils", out_path=resid_path)
        result["plots"] = [scatter_path.name, resid_path.name]
    return result


# ── engine: iraf (wrapper around the existing cross-check) ───────────────────────

def _iraf_runtime_available(runtime_cmd: Sequence[str] | None) -> bool:
    import os

    if runtime_cmd:
        head = runtime_cmd[0]
    else:
        head = "wsl" if os.name == "nt" else "python3"
    if shutil.which(head) is not None:
        return True
    # An explicit path may be given instead of a PATH name.
    return Path(head).exists()


def crosscheck_iraf(
    image_path: str | Path,
    step7_tsv: str | Path,
    *,
    fwhm_px: float | None = None,
    gain: float | None = None,
    snr_min: float = 20.0,
    output_dir: str | Path | None = None,
    runtime_cmd: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Thin wrapper around the existing IRAF cross-check.

    Reuses :func:`apex.benchmark.iraf_crosscheck.run_iraf_crosscheck` verbatim;
    this does NOT reimplement any IRAF logic. Returns ``{"status": "skipped",
    "reason": ...}`` when the IRAF/PyRAF runtime is not found on PATH (or at a
    configured explicit/WSL path). Never raises for a missing runtime.
    """
    image_path = Path(image_path)
    step7_tsv = Path(step7_tsv)
    if not image_path.exists():
        return {"status": "skipped", "reason": f"science image not found: {image_path}", "reference": "iraf"}
    if not step7_tsv.exists():
        return {"status": "skipped", "reason": f"Step 7 TSV not found: {step7_tsv}", "reference": "iraf"}
    if not _iraf_runtime_available(runtime_cmd):
        return {
            "status": "skipped",
            "reason": "IRAF/PyRAF runtime not found on PATH (checked wsl/python3 and any configured path)",
            "reference": "iraf",
        }

    try:
        from apex.benchmark.iraf_crosscheck import IRAFCrosscheckConfig, run_iraf_crosscheck
    except Exception as exc:  # noqa: BLE001
        return {"status": "skipped", "reason": f"could not load IRAF cross-check: {exc}", "reference": "iraf"}

    out_root = Path(output_dir) / "iraf_crosscheck" if output_dir is not None else None
    cc = IRAFCrosscheckConfig(
        input_fits=str(image_path),
        step7_tsv=str(step7_tsv),
        output_root=str(out_root) if out_root is not None else "benchmark/runs/iraf_crosscheck",
        mode="phot_fixed_coords",
        min_snr=float(snr_min),
        overwrite=True,
    )
    if fwhm_px is not None and math.isfinite(float(fwhm_px)):
        cc.fwhm_px = float(fwhm_px)
    if gain is not None and math.isfinite(float(gain)) and float(gain) > 0:
        cc.gain_e_per_adu = float(gain)
    if runtime_cmd:
        cc.runtime_cmd = list(runtime_cmd)

    try:
        produced = run_iraf_crosscheck(cc)
    except Exception as exc:  # noqa: BLE001 - surface as skipped, never crash report
        return {"status": "skipped", "reason": f"IRAF cross-check could not complete: {exc}", "reference": "iraf"}

    summary_path = Path(produced) / "phot_fixed_coords" / "fixed_summary.json"
    result: dict[str, Any] = {"status": "ok", "reference": "iraf", "output_dir": str(produced)}
    if summary_path.exists():
        fs = json.loads(summary_path.read_text(encoding="utf-8"))
        # Map IRAF cross-check fields onto the common agreement schema.
        result.update(
            {
                "n_matched": fs.get("n_matched"),
                "dmag_median": _finite_float(fs.get("median_delta_mag_raw")),
                "dmag_mad": _finite_float(fs.get("mad_sigma_delta_mag_zp_aligned"))
                or _finite_float(fs.get("mad_sigma_delta_mag_centered")),
                "dmag_rms": _finite_float(fs.get("rms_delta_mag_zp_aligned"))
                or _finite_float(fs.get("rms_delta_mag_centered")),
                "robust_sigma": _finite_float(fs.get("mad_sigma_delta_mag_centered")),
                "frac_within_0p05": (
                    1.0 - _finite_float(fs.get("outlier_fraction_abs_centered_gt_0p05"), float("nan"))
                    if _finite_float(fs.get("outlier_fraction_abs_centered_gt_0p05")) is not None
                    else None
                ),
                "position_rms_px": _finite_float(fs.get("position_rms_px")),
            }
        )
    return result


# ── frame picker + orchestration ────────────────────────────────────────────────

def _pick_representative_frame(index: pd.DataFrame, *, frame: str | None) -> dict[str, Any]:
    """Pick a frame: explicit name if given, else status==ok with median FWHM."""
    work = index.copy()
    if "file" not in work.columns:
        raise ValueError("photometry_index.csv has no 'file' column")
    if frame is not None:
        name = Path(frame).name
        match = work[work["file"].astype(str).map(lambda p: Path(p).name) == name]
        if match.empty:
            raise ValueError(f"frame '{frame}' not found in photometry_index.csv")
        return match.iloc[0].to_dict()

    if "status" in work.columns:
        ok = work[work["status"].astype(str).str.strip().str.lower() == "ok"]
        if not ok.empty:
            work = ok
    if "fwhm_px" in work.columns:
        fwhm = pd.to_numeric(work["fwhm_px"], errors="coerce")
        valid = work[np.isfinite(fwhm.to_numpy(float))].copy()
        if not valid.empty:
            valid_fwhm = pd.to_numeric(valid["fwhm_px"], errors="coerce").to_numpy(float)
            target = float(np.median(valid_fwhm))
            order = np.argsort(np.abs(valid_fwhm - target))
            return valid.iloc[int(order[0])].to_dict()
    return work.iloc[0].to_dict()


def _tsv_path_for_row(result_dir: Path, row: dict[str, Any]) -> Path:
    """Resolve the Step 7 TSV path for an index row."""
    step7 = _step7_dir(result_dir)
    raw_path = row.get("path")
    if raw_path:
        cand = Path(str(raw_path))
        if not cand.is_absolute():
            cand = result_dir / cand
        if cand.exists():
            return cand
        # Try by basename inside the step7 dir.
        by_name = step7 / Path(str(raw_path)).name
        if by_name.exists():
            return by_name
    # Fall back to the canonical naming: photometry_<file>.tsv
    name = Path(str(row.get("file", ""))).name
    cand = step7 / f"photometry_{name}.tsv"
    if cand.exists():
        return cand
    raise FileNotFoundError(f"Could not resolve Step 7 TSV for frame '{name}' under {step7}")


def run_crosscheck(
    result_dir: str | Path,
    output_dir: str | Path,
    *,
    frame: str | None = None,
    references: Sequence[str] = ("sep",),
    params: Any = None,
    snr_min: float = 20.0,
) -> dict[str, Any]:
    """Auto-pick a representative frame and cross-check it against engine(s).

    Loads photometry_index.csv, picks a frame (explicit ``frame`` or status==ok
    nearest to the median FWHM), resolves the SAME image APEX measured, then runs
    each requested reference engine and aggregates the agreement stats.
    """
    result_dir = Path(result_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    refs = [r.strip().lower() for r in references]
    if "all" in refs:
        refs = list(REFERENCE_CHOICES)
    unknown = [r for r in refs if r not in REFERENCE_CHOICES]
    if unknown:
        raise ValueError(f"Unknown crosscheck reference(s): {', '.join(unknown)}")

    try:
        index = _read_index(result_dir)
    except FileNotFoundError as exc:
        return {"status": "error", "reason": str(exc), "result_dir": str(result_dir)}

    try:
        row = _pick_representative_frame(index, frame=frame)
    except ValueError as exc:
        return {"status": "error", "reason": str(exc), "result_dir": str(result_dir)}

    frame_name = Path(str(row.get("file", ""))).name
    fwhm_px = _finite_float(row.get("fwhm_px"))
    scales = _resolve_scales(params)
    gain = scales["gain"]

    try:
        tsv_path = _tsv_path_for_row(result_dir, row)
    except FileNotFoundError as exc:
        return {"status": "error", "reason": str(exc), "frame": frame_name}

    try:
        image_path = resolve_frame_image(result_dir, frame_name, params)
    except FileNotFoundError as exc:
        return {"status": "error", "reason": str(exc), "frame": frame_name, "step7_tsv": str(tsv_path)}

    if fwhm_px is None or fwhm_px <= 0:
        return {
            "status": "error",
            "reason": f"no usable fwhm_px in photometry_index.csv for frame '{frame_name}'",
            "frame": frame_name,
        }

    plot_dir = output_dir / "plots"
    aggregate: dict[str, Any] = {
        "status": "ok",
        "result_dir": str(result_dir),
        "frame": frame_name,
        "filter": str(row.get("filter", "")),
        "fwhm_px": float(fwhm_px),
        "gain": float(gain),
        "snr_min": float(snr_min),
        "image": str(image_path),
        "step7_tsv": str(tsv_path),
        "aperture_scales": scales,
        "references": {},
    }

    common = dict(
        fwhm_px=float(fwhm_px),
        aperture_scale=scales["aperture_scale"],
        ref_ap_scale=scales["ref_ap_scale"],
        ann_in_scale=scales["ann_in_scale"],
        dann_scale=scales["dann_scale"],
        min_r_ap_px=scales["min_r_ap_px"],
        ann_min_gap_px=scales["ann_min_gap_px"],
        gain=float(gain),
        snr_min=float(snr_min),
        output_dir=plot_dir,
    )

    for ref in refs:
        try:
            if ref == "sep":
                aggregate["references"]["sep"] = crosscheck_sep(image_path, tsv_path, **common)
            elif ref == "photutils":
                aggregate["references"]["photutils"] = crosscheck_photutils(image_path, tsv_path, **common)
            elif ref == "iraf":
                aggregate["references"]["iraf"] = crosscheck_iraf(
                    image_path, tsv_path,
                    fwhm_px=float(fwhm_px), gain=float(gain), snr_min=float(snr_min),
                    output_dir=output_dir,
                )
        except ImportError as exc:
            aggregate["references"][ref] = {"status": "skipped", "reason": f"{ref} not installed: {exc}", "reference": ref}
        except Exception as exc:  # noqa: BLE001 - never crash the suite on one engine
            aggregate["references"][ref] = {"status": "error", "reason": str(exc), "reference": ref}

    # Sanity verdict from the first engine that ran a position sanity check.
    for ref_result in aggregate["references"].values():
        sanity = ref_result.get("sanity") if isinstance(ref_result, dict) else None
        if sanity is not None:
            aggregate["positions_plausible"] = bool(sanity.get("positions_plausible"))
            aggregate["sanity"] = sanity
            break

    (output_dir / "crosscheck.json").write_text(
        json.dumps(aggregate, indent=2, sort_keys=True, default=str), encoding="utf-8"
    )
    return aggregate
