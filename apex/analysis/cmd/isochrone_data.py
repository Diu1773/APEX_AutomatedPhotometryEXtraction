"""Shared, Qt-free data preparation for CMD isochrone fitting.

These helpers were previously private to ``validation/cmd_step12_realdata.py``.
They are promoted into the ``apex`` package so the validation runner, the headless
pipeline, **and the GUI** can all use one implementation (the GUI must not import
from ``validation/``). Everything here is pure data/IO + numpy — no Qt, no
plotting — so it respects the layer rule (gui → analysis, analysis ↛ gui).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from apex.utils.gaia_transforms import build_color_pairs, filter_bands_from_columns

__all__ = [
    "PARSEC_ISO_COL",
    "EXTINCTION_R",
    "repo_root",
    "resolve_wide_csv",
    "infer_color_mag",
    "default_iso_file",
    "magnitude_column",
    "membership_mask",
    "load_multiband_isochrone",
    "extract_multicolor_arrays",
]

# PARSEC CMD 3.x column indices (0-based) for each photometric band.
PARSEC_ISO_COL = {
    "u": 28, "g": 29, "r": 30, "i": 31, "z": 32,
    "B": 30, "V": 31, "R": 32, "I": 33,
}

# Total-to-selective extinction A_band / E(B-V) (Cardelli/SF2011, R_V=3.1).
EXTINCTION_R = {
    "U": 4.902, "B": 4.035, "V": 3.116, "R": 2.634, "I": 1.903,
    "u": 4.239, "g": 3.303, "r": 2.285, "i": 1.698, "z": 1.263,
}


def repo_root() -> Path:
    """Repository root (…/apex/analysis/cmd/isochrone_data.py → repo)."""
    return Path(__file__).resolve().parents[3]


def resolve_wide_csv(result_dir: Path) -> Path:
    """Locate the Step-10 wide CMD CSV under an APEX result directory."""
    result_dir = Path(result_dir)
    candidates = [
        result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv",
        result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide.csv",
        result_dir / "median_by_ID_filter_wide_cmd.csv",
        result_dir / "median_by_ID_filter_wide.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("Step10 wide CMD CSV not found under result_dir")


def infer_color_mag(df: pd.DataFrame) -> Tuple[Tuple[str, str], str]:
    """Pick a sensible default (color, mag_band) from the available bands."""
    bands = filter_bands_from_columns(df.columns, "mag_std_")
    pairs = build_color_pairs(bands, adjacent_only=True)
    for pair, mag in [("g-r", "g"), ("B-V", "V"), ("V-I", "V"), ("r-i", "r")]:
        b1, b2 = pair.split("-")
        if (b1, b2) in pairs and mag in bands:
            return (b1, b2), mag
    if not pairs:
        raise ValueError("No mag_std_* bands found in Step10 CSV")
    b1, b2 = pairs[0]
    return (b1, b2), b1


def default_iso_file(color: Tuple[str, str]) -> Path:
    """Default bundled PARSEC grid for a color (Johnson vs SDSS auto-pick)."""
    root = repo_root()
    if any(b in {"B", "V", "R", "I"} for b in color):
        return root / "isochrone" / "PARSEC" / "johnson" / "iso_data_cmd39_johnson.dat"
    return root / "isochrone" / "PARSEC" / "sdss" / "iso_data_cmd39.dat"


def magnitude_column(df: pd.DataFrame, band: str) -> Tuple[np.ndarray, np.ndarray, str]:
    """Return (values, errors, source_column) for a band, preferring mag_std."""
    std = f"mag_std_{band}"
    inst = f"mag_inst_{band}"
    if std in df.columns:
        values = pd.to_numeric(df[std], errors="coerce").to_numpy(float)
        source = std
    elif inst in df.columns:
        values = pd.to_numeric(df[inst], errors="coerce").to_numpy(float)
        source = inst
    else:
        raise KeyError(f"Missing magnitude column for band {band}")

    err_candidates = [f"mag_std_err_{band}", f"mag_inst_err_{band}", f"mag_inst_med_err_{band}"]
    err = None
    for col in err_candidates:
        if col in df.columns:
            err = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
            break
    if err is None:
        err = np.full(len(values), 0.03, dtype=float)
    err = np.where(np.isfinite(err) & (err > 0), err, 0.03)
    return values, err, source


def membership_mask(
    df: pd.DataFrame, *, nsig: float = 3.0, plx_nsig: float = 3.0,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Boolean mask of probable cluster members from Gaia proper motion + parallax.

    Auto-detects the cluster clump in (pmra, pmdec) as the densest 2-D histogram
    peak, refines the centre by iterative robust (MAD) clipping, then keeps stars
    within ``nsig`` robust sigma of the PM centroid AND ``plx_nsig`` of the
    members' median parallax. Stars lacking astrometry are treated as non-members
    (conservative). Returns an all-True mask as a graceful no-op when the
    astrometry columns are absent.
    """
    need = {"pmra", "pmdec", "parallax"}
    if not need.issubset(df.columns):
        return np.ones(len(df), dtype=bool), {"applied": False, "reason": "no_astrometry"}
    pmra = pd.to_numeric(df["pmra"], errors="coerce").to_numpy(float)
    pmdec = pd.to_numeric(df["pmdec"], errors="coerce").to_numpy(float)
    plx = pd.to_numeric(df["parallax"], errors="coerce").to_numpy(float)
    fin = np.isfinite(pmra) & np.isfinite(pmdec) & np.isfinite(plx)
    if int(fin.sum()) < 20:
        return np.ones(len(df), dtype=bool), {"applied": False, "reason": "too_few_astrometric"}

    x, y, p = pmra[fin], pmdec[fin], plx[fin]
    lo_x, hi_x = np.percentile(x, [2, 98])
    lo_y, hi_y = np.percentile(y, [2, 98])
    bw = 0.5
    bx = np.arange(lo_x, hi_x + bw, bw)
    by = np.arange(lo_y, hi_y + bw, bw)
    if len(bx) < 2 or len(by) < 2:
        return np.ones(len(df), dtype=bool), {"applied": False, "reason": "degenerate_pm_range"}
    hist, xe, ye = np.histogram2d(x, y, bins=[bx, by])
    ip, jp = np.unravel_index(int(np.argmax(hist)), hist.shape)
    cx = 0.5 * (xe[ip] + xe[ip + 1])
    cy = 0.5 * (ye[jp] + ye[jp + 1])

    sel = np.hypot(x - cx, y - cy) < 1.5
    sx = sy = 0.3
    for _ in range(5):
        if int(sel.sum()) < 10:
            break
        cx = float(np.median(x[sel]))
        cy = float(np.median(y[sel]))
        sx = max(1.4826 * float(np.median(np.abs(x[sel] - cx))), 0.2)
        sy = max(1.4826 * float(np.median(np.abs(y[sel] - cy))), 0.2)
        sel = (np.abs(x - cx) < nsig * sx) & (np.abs(y - cy) < nsig * sy)

    pc = float(np.median(p[sel]))
    ps = max(1.4826 * float(np.median(np.abs(p[sel] - pc))), 0.05)
    member_fin = sel & (np.abs(p - pc) < plx_nsig * ps)

    mask = np.zeros(len(df), dtype=bool)
    mask[np.flatnonzero(fin)[member_fin]] = True
    meta = {
        "applied": True,
        "pm_center": [round(cx, 3), round(cy, 3)],
        "pm_sigma": [round(sx, 3), round(sy, 3)],
        "parallax_center": round(pc, 3),
        "parallax_sigma": round(ps, 3),
        "n_astrometric": int(fin.sum()),
        "n_members": int(mask.sum()),
        "nsig": nsig,
        "plx_nsig": plx_nsig,
    }
    return mask, meta


def load_multiband_isochrone(
    iso_file: Path,
    bands: Sequence[str],
    age_bounds: Tuple[float, float],
    mh_bounds: Tuple[float, float],
    cache_dir: Optional[Path] = None,
) -> Tuple[np.ndarray, Dict[str, int], Dict[str, Any]]:
    """Load a compact PARSEC subset carrying every requested band.

    Returns ``(iso_data, band_col_in_subset, meta)`` where ``iso_data`` columns
    are ``[pad, mh, log_age, band_0, band_1, ..., mass]``.
    """
    iso_file = Path(iso_file)
    bands = list(bands)
    missing = [b for b in bands if b not in PARSEC_ISO_COL]
    if missing:
        raise KeyError(f"No PARSEC column mapping for: {', '.join(missing)}")

    band_src_cols = {b: PARSEC_ISO_COL[b] for b in bands}
    usecols = sorted({1, 2, 5, *band_src_cols.values()})

    cache_path = None
    if cache_dir is not None:
        stat = iso_file.stat()
        key_text = json.dumps(
            {
                "iso_file": str(iso_file.resolve()),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "bands": bands,
                "age_bounds": list(age_bounds),
                "mh_bounds": list(mh_bounds),
                "kind": "multiband",
            },
            sort_keys=True,
        )
        digest = hashlib.sha1(key_text.encode("utf-8")).hexdigest()[:16]
        cache_path = Path(cache_dir) / f"parsec_multiband_{digest}.npy"
        if cache_path.exists():
            data = np.load(cache_path)
            band_col = {b: 3 + k for k, b in enumerate(bands)}
            return data, band_col, {"cache_hit": True, "cache_path": str(cache_path)}

    chunks: List[np.ndarray] = []
    for chunk in pd.read_csv(
        iso_file, sep=r"\s+", comment="#", header=None,
        usecols=usecols, chunksize=200_000, engine="c",
    ):
        age = pd.to_numeric(chunk[2], errors="coerce").to_numpy(float)
        mh = pd.to_numeric(chunk[1], errors="coerce").to_numpy(float)
        mask = (
            np.isfinite(age) & np.isfinite(mh)
            & (age >= age_bounds[0]) & (age <= age_bounds[1])
            & (mh >= mh_bounds[0]) & (mh <= mh_bounds[1])
        )
        if not mask.any():
            continue
        sub = chunk.loc[mask]
        cols = [np.zeros(len(sub)), sub[1].to_numpy(float), sub[2].to_numpy(float)]
        cols += [sub[band_src_cols[b]].to_numpy(float) for b in bands]
        cols.append(sub[5].to_numpy(float))
        chunks.append(np.column_stack(cols))

    if not chunks:
        raise ValueError("No isochrone rows survived the selected age/[M/H] bounds")
    iso_data = np.vstack(chunks)
    band_col = {b: 3 + k for k, b in enumerate(bands)}
    meta = {"cache_hit": False, "cache_path": str(cache_path) if cache_path else None}
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, iso_data)
    return iso_data, band_col, meta


def extract_multicolor_arrays(
    df: pd.DataFrame,
    colors: Sequence[Tuple[str, str]],
    mag_band: str,
    data_snr_min: float,
    max_stars: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """Extract an (N, D) observed array + per-axis errors for the fit.

    Axis order: one column per colour (b1 − b2) then the magnitude band.
    """
    colors = [tuple(c) for c in colors]
    bands_needed = sorted({b for pair in colors for b in pair} | {mag_band})
    vals: Dict[str, np.ndarray] = {}
    errs: Dict[str, np.ndarray] = {}
    sources: Dict[str, str] = {}
    for b in bands_needed:
        v, e, src = magnitude_column(df, b)
        vals[b], errs[b], sources[b] = v, e, src

    axes, axis_errs, labels = [], [], []
    for b1, b2 in colors:
        axes.append(vals[b1] - vals[b2])
        axis_errs.append(np.sqrt(errs[b1] ** 2 + errs[b2] ** 2))
        labels.append(f"{b1}-{b2}")
    axes.append(vals[mag_band])
    axis_errs.append(errs[mag_band])
    labels.append(mag_band)

    obs = np.column_stack(axes)
    err = np.column_stack(axis_errs)

    mask = np.isfinite(obs).all(axis=1) & np.isfinite(err).all(axis=1)
    for b in bands_needed:
        snr_col = f"snr_{b}"
        if snr_col in df.columns and data_snr_min > 0:
            snr = pd.to_numeric(df[snr_col], errors="coerce").to_numpy(float)
            mask &= np.isfinite(snr) & (snr >= data_snr_min)

    idx = np.flatnonzero(mask)
    if max_stars > 0 and len(idx) > max_stars:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(idx, size=max_stars, replace=False))

    meta = {
        "n_total_rows": int(len(df)),
        "n_after_filters": int(len(idx)),
        "axes": labels,
        "columns": sources,
    }
    return obs[idx], err[idx], labels, meta
