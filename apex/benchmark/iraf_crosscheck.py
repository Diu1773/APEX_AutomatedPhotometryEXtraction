"""IRAF/PyRAF cross-checks for APEX CMD photometry outputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from apex.utils.constants import INSTRUMENTAL_ZMAG


FIXED_MODE = "phot_fixed_coords"
DAOFIND_MODE = "daofind_phot"
BOTH_MODE = "both"


@dataclass
class IRAFCrosscheckConfig:
    input_fits: str = ""
    step7_tsv: str = ""
    output_root: str = "benchmark/runs/iraf_crosscheck"
    mode: str = BOTH_MODE
    filter_name: str = ""
    fwhm_px: float | None = None
    aperture_scale_fwhm: float = 1.0
    annulus_scale_fwhm: float = 4.0
    dannulus_scale_fwhm: float = 2.0
    min_aperture_px: float = 4.0
    annulus_min_gap_px: float = 6.0
    max_fixed_coords: int = 500
    min_snr: float = 20.0
    sample_seed: int = 20260608
    zmag: float = 25.0
    gain_e_per_adu: float | None = None
    rdnoise_e: float | None = None
    sigma_adu: float | None = None
    exptime: float | None = None
    recenter: bool = False
    daofind_threshold_grid: list[float] = field(default_factory=lambda: [12.0, 9.0, 7.0, 5.0])
    daofind_max_sources: int = 2500
    daofind_max_ratio_to_apex: float = 2.0
    daofind_match_radius_px: float = 2.0
    overwrite: bool = False
    runtime_cmd: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class IRAFFrameMetadata:
    filter_name: str
    fwhm_px: float
    r_ap: float
    r_in: float
    r_out: float
    gain_e_per_adu: float
    rdnoise_e: float
    sigma_adu: float
    apcorr: float
    n_apex_sources: int
    frame_zeropoint_mag: float
    exptime: float = 1.0


def windows_to_wsl_path(path: str | Path) -> str:
    """Convert a Windows path to the corresponding WSL mount path."""
    text = str(path)
    if text.startswith("/"):
        return text
    resolved = str(Path(text).expanduser().resolve())
    if len(resolved) >= 3 and resolved[1] == ":" and resolved[2] in {"\\", "/"}:
        drive = resolved[0].lower()
        rest = resolved[2:].replace("\\", "/")
        return f"/mnt/{drive}{rest}"
    return resolved.replace("\\", "/")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normal_filter(value: Any) -> str:
    return str(value).strip().lower().rstrip("'")


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def _finite(values: Iterable[Any]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    return arr[np.isfinite(arr)]


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


def _median_or_nan(values: Any) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else float("nan")


def _finite_float(value: Any, default: float | None = None) -> float | None:
    try:
        out = float(value)
    except Exception:
        return default
    return out if math.isfinite(out) else default


def _default_runtime_cmd() -> list[str]:
    return ["wsl", "python3"] if os.name == "nt" else ["python3"]


def _uses_wsl(runtime_cmd: list[str]) -> bool:
    return bool(runtime_cmd) and Path(runtime_cmd[0]).name.lower() == "wsl"


def _runtime_path(path: str | Path, *, use_wsl: bool) -> str:
    resolved = Path(path).expanduser().resolve()
    return windows_to_wsl_path(resolved) if use_wsl else str(resolved)


def _clean_output_dir(path: Path, *, overwrite: bool) -> None:
    if path.exists() and any(path.iterdir()) and not overwrite:
        raise FileExistsError(f"IRAF cross-check output already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _read_step7(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Step 7 photometry table is missing: {path}")
    return pd.read_csv(path, sep="\t")


def _lookup_index_row(step7_path: Path, input_fits: Path) -> dict[str, Any]:
    index_path = step7_path.parent / "photometry_index.csv"
    if not index_path.is_file():
        return {}
    index = pd.read_csv(index_path)
    if "file" not in index.columns:
        return {}
    filename = input_fits.name
    rows = index[index["file"].astype(str).map(Path).map(lambda p: p.name) == filename]
    if rows.empty:
        rows = index[index["path"].astype(str) == str(step7_path)] if "path" in index.columns else rows
    return rows.iloc[0].to_dict() if not rows.empty else {}


def _lookup_frame_zeropoint(step7_path: Path, input_fits: Path, filter_name: str) -> float:
    zp_path = step7_path.parent.parent / "cmd_zeropoint" / "frame_zeropoint.csv"
    if not zp_path.is_file():
        return float("nan")
    try:
        frame_zp = pd.read_csv(zp_path)
    except Exception:
        return float("nan")
    required = {"file", "filter", "zp_frame"}
    if not required <= set(frame_zp.columns):
        return float("nan")
    work = frame_zp.copy()
    work["file_key"] = work["file"].astype(str).map(lambda text: Path(text).name)
    work["filter_key"] = work["filter"].map(_normal_filter)
    rows = work[
        (work["file_key"] == input_fits.name)
        & (work["filter_key"] == _normal_filter(filter_name))
    ]
    if rows.empty:
        rows = work[work["file_key"] == input_fits.name]
    if rows.empty:
        return float("nan")
    return _finite_float(rows.iloc[0].get("zp_frame"), float("nan")) or float("nan")


def _mode_or_empty(values: pd.Series) -> str:
    cleaned = values.dropna().astype(str).str.strip()
    if cleaned.empty:
        return ""
    modes = cleaned.mode()
    return str(modes.iloc[0]) if not modes.empty else str(cleaned.iloc[0])


def _positive_median(table: pd.DataFrame, column: str, default: float) -> float:
    if column not in table.columns:
        return float(default)
    vals = _num(table[column]).to_numpy(float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    return float(np.median(vals)) if vals.size else float(default)


def _resolve_metadata(
    config: IRAFCrosscheckConfig,
    step7: pd.DataFrame,
    index_row: dict[str, Any],
    step7_path: Path,
    input_fits: Path,
) -> IRAFFrameMetadata:
    filter_name = _normal_filter(
        config.filter_name
        or index_row.get("filter", "")
        or (_mode_or_empty(step7["filter"]) if "filter" in step7.columns else "")
        or (_mode_or_empty(step7["FILTER"]) if "FILTER" in step7.columns else "")
        or "unknown"
    )
    fwhm = _finite_float(config.fwhm_px)
    if fwhm is None:
        fwhm = _finite_float(index_row.get("fwhm_px"))
    if fwhm is None or fwhm <= 0:
        raise ValueError("fwhm_px is required when it cannot be read from photometry_index.csv")

    r_ap = max(float(config.min_aperture_px), float(config.aperture_scale_fwhm) * fwhm)
    r_in = max(r_ap + float(config.annulus_min_gap_px), float(config.annulus_scale_fwhm) * fwhm)
    r_out = r_in + max(float(config.annulus_min_gap_px), float(config.dannulus_scale_fwhm) * fwhm)

    gain = _finite_float(config.gain_e_per_adu)
    if gain is None or gain <= 0:
        gain = _finite_float(index_row.get("gain_e_per_adu"))
    if gain is None or gain <= 0:
        gain = _positive_median(step7, "gain_e_per_adu", 1.0)

    rdnoise = _finite_float(config.rdnoise_e)
    if rdnoise is None or rdnoise < 0:
        rdnoise = _finite_float(index_row.get("rdnoise_e"))
    if rdnoise is None or rdnoise < 0:
        rdnoise = _positive_median(step7, "rdnoise_e", 0.0)

    sigma = _finite_float(config.sigma_adu)
    if sigma is None or sigma <= 0:
        sigma = _positive_median(step7, "bkg_std_adu", float("nan"))
    if sigma is None or not math.isfinite(sigma) or sigma <= 0:
        sigma = _positive_median(step7, "sky_std", 1.0)

    apcorr = _finite_float(index_row.get("apcorr"))
    if apcorr is None or apcorr <= 0:
        apcorr = _positive_median(step7, "apcorr", 1.0)

    # Exposure time drives IRAF `datapars.itime` so IRAF reports a count-rate
    # magnitude on the same scale as APEX mag_inst (which is normalized by
    # exptime). Without this, mixed-exposure frames would be compared on
    # different flux scales.
    exptime = _finite_float(config.exptime)
    if exptime is None or exptime <= 0:
        exptime = _finite_float(index_row.get("exptime"))
    if exptime is None or exptime <= 0:
        exptime = _positive_median(step7, "exptime", 1.0)

    n_sources = _finite_float(index_row.get("n_detected"))
    if n_sources is None:
        n_sources = _finite_float(index_row.get("n_master"))
    if n_sources is None:
        n_sources = float(len(step7))

    frame_zp = _lookup_frame_zeropoint(step7_path, input_fits, filter_name)

    return IRAFFrameMetadata(
        filter_name=filter_name,
        fwhm_px=float(fwhm),
        r_ap=float(r_ap),
        r_in=float(r_in),
        r_out=float(r_out),
        gain_e_per_adu=float(gain),
        rdnoise_e=float(rdnoise),
        sigma_adu=float(sigma),
        apcorr=float(apcorr),
        n_apex_sources=int(max(0, round(float(n_sources)))),
        frame_zeropoint_mag=float(frame_zp),
        exptime=float(exptime),
    )


def load_apex_reference_table(
    step7_path: str | Path,
    *,
    min_snr: float = 20.0,
) -> pd.DataFrame:
    """Load finite, non-flagged Step 7 sources for IRAF cross-checks."""
    table = _read_step7(Path(step7_path))
    required = {"x", "y", "mag_inst"}
    if not required <= set(table.columns):
        missing = sorted(required - set(table.columns))
        raise ValueError(f"Step 7 table is missing columns: {', '.join(missing)}")

    out = table.copy()
    out["x"] = _num(out["x"])
    out["y"] = _num(out["y"])
    out["mag_inst"] = _num(out["mag_inst"])
    if "snr" in out.columns:
        out["snr"] = _num(out["snr"])
    else:
        out["snr"] = np.nan

    good = (
        np.isfinite(out["x"].to_numpy(float))
        & np.isfinite(out["y"].to_numpy(float))
        & np.isfinite(out["mag_inst"].to_numpy(float))
    )
    if math.isfinite(float(min_snr)):
        good &= out["snr"].to_numpy(float) >= float(min_snr)
    for flag in ("bad_phot_flag", "is_saturated", "is_nonlinear", "off_frame_flag", "centroid_outlier"):
        if flag in out.columns:
            good &= ~_bool_series(out[flag]).to_numpy(bool)

    if "detected_flag" in out.columns:
        out["detected_flag"] = _bool_series(out["detected_flag"])

    keep_cols = [
        col
        for col in (
            "file",
            "filter",
            "ID",
            "master_id",
            "source_id",
            "det_uid",
            "x",
            "y",
            "mag_inst",
            "mag_err",
            "snr",
            "flux",
            "flux_e",
            "flux_net_adu",
            "gain_e_per_adu",
            "apcorr",
            "detected_flag",
            "forced_flag",
            "recentered_flag",
            "center_error_px",
            "centroid_shift_px",
        )
        if col in out.columns
    ]
    out = out.loc[good, keep_cols].copy()
    for col in ("mag_err", "flux", "flux_e", "flux_net_adu", "gain_e_per_adu", "apcorr"):
        if col in out.columns:
            out[col] = _num(out[col])
    return out.reset_index(drop=True)


def select_fixed_coordinate_sources(
    reference: pd.DataFrame,
    *,
    max_sources: int = 500,
    seed: int = 20260608,
    bins: int = 8,
) -> pd.DataFrame:
    """Select a reproducible magnitude-stratified subset for fixed-coordinate phot."""
    if int(max_sources) <= 0 or len(reference) <= int(max_sources):
        selected = reference.copy()
    else:
        rng = np.random.default_rng(int(seed))
        work = reference.copy()
        q = min(int(bins), max(1, int(work["mag_inst"].nunique())))
        try:
            work["_mag_bin"] = pd.qcut(work["mag_inst"], q=q, duplicates="drop")
        except ValueError:
            work["_mag_bin"] = "all"
        groups = list(work.groupby("_mag_bin", observed=True, sort=True))
        per_bin = max(1, int(math.ceil(int(max_sources) / max(1, len(groups)))))
        parts: list[pd.DataFrame] = []
        for _, group in groups:
            take = min(len(group), per_bin)
            parts.append(group.sample(n=take, random_state=int(rng.integers(0, 2**31 - 1))))
        selected = pd.concat(parts, ignore_index=False)
        if len(selected) > int(max_sources):
            selected = selected.sample(n=int(max_sources), random_state=int(seed))
        elif len(selected) < int(max_sources):
            remaining = work.drop(index=selected.index, errors="ignore")
            take = min(len(remaining), int(max_sources) - len(selected))
            if take > 0:
                selected = pd.concat(
                    [
                        selected,
                        remaining.sample(n=take, random_state=int(seed) + 1),
                    ],
                    ignore_index=False,
                )
        selected = selected.drop(columns=["_mag_bin"], errors="ignore")
    selected = selected.sort_values(["mag_inst", "x", "y"], kind="stable").reset_index(drop=True)
    selected.insert(0, "iraf_id", np.arange(1, len(selected) + 1, dtype=int))
    return selected


# IRAF numbers the first pixel 1; APEX, numpy and photutils number it 0.
IRAF_ORIGIN_OFFSET = 1.0


def write_iraf_coords(reference: pd.DataFrame, path: str | Path) -> Path:
    """Write an IRAF coordinate file in the same order as ``iraf_id``.

    APEX counts the first pixel 0, like numpy and photutils; IRAF counts it 1.
    Writing APEX's numbers unchanged asks IRAF to measure one pixel down and to
    the left of every star — 1.41 px away — which was the state of this file
    until 2026-08-14. Verified against the pixels: centroids of 60 bright stars
    land 0.13 px from APEX's coordinates and 1.40 px from those coordinates
    read as 1-based.

    With `recenter` off (the default) IRAF cannot walk back onto the star, so
    the aperture stays off-centre. In a sparse field that is close to a pure
    zero-point — every star loses the same fraction — which is why the earlier
    cross-checks still agreed to a few mmag; in a crowded field it is not.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii", newline="\n") as handle:
        for _, row in reference.sort_values("iraf_id").iterrows():
            handle.write(f"{float(row['x']) + IRAF_ORIGIN_OFFSET:.6f} "
                         f"{float(row['y']) + IRAF_ORIGIN_OFFSET:.6f}\n")
    return out


def parse_txdump(path: str | Path) -> pd.DataFrame:
    names = [
        "iraf_id",
        "iraf_x",
        "iraf_y",
        "iraf_mag",
        "iraf_merr",
        "iraf_msky",
        "iraf_stdev",
        "iraf_nsky",
    ]
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return pd.DataFrame(columns=names)
    rows: list[list[Any]] = []
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < len(names):
            continue
        rows.append(parts[: len(names)])
    table = pd.DataFrame(rows, columns=names)
    if table.empty:
        return pd.DataFrame(columns=names)
    table["iraf_id"] = pd.to_numeric(table["iraf_id"], errors="coerce").astype("Int64")
    for col in names[1:]:
        table[col] = pd.to_numeric(table[col].replace("INDEF", np.nan), errors="coerce")
    # Back into APEX's 0-based frame, so `iraf_x - x` is a real disagreement
    # about where the star is rather than a constant one-pixel bookkeeping
    # difference. See `write_iraf_coords` for the other half.
    table["iraf_x"] = table["iraf_x"] - IRAF_ORIGIN_OFFSET
    table["iraf_y"] = table["iraf_y"] - IRAF_ORIGIN_OFFSET
    return table[np.isfinite(table["iraf_x"]) & np.isfinite(table["iraf_y"])].reset_index(drop=True)


def _apex_iraf_unit_mag(reference: pd.DataFrame, zmag: float) -> pd.Series:
    """Express APEX mag_inst on IRAF's count-rate magnitude scale.

    APEX mag_inst = INSTRUMENTAL_ZMAG - 2.5 log10(flux_e/exptime) already carries
    a zeropoint and exposure-time normalization. IRAF `phot` (with datapars.itime
    set to exptime) reports MAG = zmag - 2.5 log10(flux_ADU/exptime). The two
    differ only by the zeropoint offset (zmag - INSTRUMENTAL_ZMAG) and the
    gain/apcorr factors; exptime cancels.
    """
    mag = _num(reference["mag_inst"]) + (float(zmag) - INSTRUMENTAL_ZMAG)
    if "apcorr" in reference.columns:
        apcorr = _num(reference["apcorr"])
        good = np.isfinite(apcorr) & (apcorr > 0)
        mag = mag.where(~good, mag + 2.5 * np.log10(apcorr))
    if "gain_e_per_adu" in reference.columns:
        gain = _num(reference["gain_e_per_adu"])
        good = np.isfinite(gain) & (gain > 0)
        mag = mag.where(~good, mag + 2.5 * np.log10(gain))
    return mag


def _positive_factor(table: pd.DataFrame, column: str, default: float = 1.0) -> pd.Series:
    if column not in table.columns:
        return pd.Series(float(default), index=table.index, dtype=float)
    values = _num(table[column])
    good = np.isfinite(values.to_numpy(float)) & (values.to_numpy(float) > 0)
    return values.where(good, float(default)).astype(float)


def add_iraf_calibrated_equivalent_columns(
    table: pd.DataFrame,
    *,
    zmag: float,
    frame_zeropoint: float | None,
) -> pd.DataFrame:
    """Put IRAF `phot` magnitudes onto APEX's apcorr+ZP-equivalent scale.

    IRAF `phot` (with datapars.itime = exptime) reports a count-rate magnitude
    `MAG = zmag - 2.5 log10(flux_ADU/exptime)`. APEX Step 7 stores an
    aperture-corrected, exposure-normalized instrumental magnitude using
    electrons: `mag_inst = INSTRUMENTAL_ZMAG - 2.5 log10(flux_e/exptime * apcorr)`.
    The two scales differ only by the zeropoint offset and the gain/apcorr
    factors (exptime cancels), so when IRAF's zmag equals INSTRUMENTAL_ZMAG the
    converted IRAF magnitude is identical to APEX mag_inst for matched sources.
    """
    out = table.copy()
    if "iraf_mag" not in out.columns:
        return out

    gain = _positive_factor(out, "gain_e_per_adu", 1.0)
    apcorr = _positive_factor(out, "apcorr", 1.0)
    out["iraf_gain_factor"] = gain
    out["iraf_apcorr_factor"] = apcorr
    out["iraf_mag_inst_apcorr_equiv"] = (
        _num(out["iraf_mag"])
        + (INSTRUMENTAL_ZMAG - float(zmag))
        - 2.5 * np.log10(gain)
        - 2.5 * np.log10(apcorr)
    )

    frame_zp = _finite_float(frame_zeropoint, float("nan"))
    if frame_zp is not None and math.isfinite(frame_zp):
        out["iraf_mag_cal_apcorr_zp"] = out["iraf_mag_inst_apcorr_equiv"] + float(frame_zp)
        if "apex_mag_cal_frame" not in out.columns and "mag_inst" in out.columns:
            out["apex_mag_cal_frame"] = _num(out["mag_inst"]) + float(frame_zp)
    else:
        out["iraf_mag_cal_apcorr_zp"] = np.nan

    if {"apex_mag_cal_frame", "iraf_mag_cal_apcorr_zp"} <= set(out.columns):
        out["delta_mag_apcorr_zp"] = _num(out["apex_mag_cal_frame"]) - _num(out["iraf_mag_cal_apcorr_zp"])
        offset = _median_or_nan(out["delta_mag_apcorr_zp"])
        out["iraf_apcorr_zp_to_apex_offset"] = offset
        out["delta_mag_apcorr_zp_centered"] = out["delta_mag_apcorr_zp"] - offset
        out["iraf_mag_apcorr_zp_aligned_to_apex"] = out["iraf_mag_cal_apcorr_zp"] + offset
        out["delta_mag_apcorr_zp_aligned"] = (
            _num(out["apex_mag_cal_frame"]) - _num(out["iraf_mag_apcorr_zp_aligned_to_apex"])
        )
    else:
        out["delta_mag_apcorr_zp"] = np.nan
        out["iraf_apcorr_zp_to_apex_offset"] = np.nan
        out["delta_mag_apcorr_zp_centered"] = np.nan
        out["iraf_mag_apcorr_zp_aligned_to_apex"] = np.nan
        out["delta_mag_apcorr_zp_aligned"] = np.nan
    return out


def _add_alignment_columns(
    table: pd.DataFrame,
    *,
    zmag: float,
    frame_zeropoint: float | None,
) -> pd.DataFrame:
    out = table.copy()
    # mag_inst already carries INSTRUMENTAL_ZMAG (IRAF-style count-rate scale),
    # so to express it on IRAF's photpars.zmag scale we only re-base by the
    # zeropoint difference — not add zmag a second time.
    out["apex_mag_iraf_zmag"] = _num(out["mag_inst"]) + (float(zmag) - INSTRUMENTAL_ZMAG)
    out["apex_mag_iraf_units"] = _apex_iraf_unit_mag(out, float(zmag))
    out["delta_mag_raw"] = out["apex_mag_iraf_zmag"] - out["iraf_mag"]
    out["delta_mag_units"] = out["apex_mag_iraf_units"] - out["iraf_mag"]
    med = float(np.nanmedian(out["delta_mag_raw"]))
    med_units = float(np.nanmedian(out["delta_mag_units"]))
    out["delta_mag_centered"] = out["delta_mag_raw"] - med
    out["delta_mag_units_centered"] = out["delta_mag_units"] - med_units

    frame_zp = _finite_float(frame_zeropoint, float("nan"))
    if frame_zp is not None and math.isfinite(frame_zp):
        out["apex_mag_cal_frame"] = _num(out["mag_inst"]) + float(frame_zp)
        zp_offset = float(np.nanmedian(out["apex_mag_cal_frame"] - out["iraf_mag"]))
        out["iraf_to_apex_zp_offset"] = zp_offset
        out["iraf_mag_aligned_to_apex"] = out["iraf_mag"] + zp_offset
        out["delta_mag_zp_aligned"] = out["apex_mag_cal_frame"] - out["iraf_mag_aligned_to_apex"]
    else:
        out["apex_mag_cal_frame"] = np.nan
        out["iraf_to_apex_zp_offset"] = np.nan
        out["iraf_mag_aligned_to_apex"] = np.nan
        out["delta_mag_zp_aligned"] = out["delta_mag_centered"]
    out = add_iraf_calibrated_equivalent_columns(
        out,
        zmag=float(zmag),
        frame_zeropoint=frame_zeropoint,
    )
    return out


def build_fixed_comparison(
    reference: pd.DataFrame,
    iraf: pd.DataFrame,
    *,
    zmag: float,
    frame_zeropoint: float | None = None,
) -> pd.DataFrame:
    if reference.empty or iraf.empty:
        return pd.DataFrame()
    left = reference.copy()
    left["iraf_id"] = pd.to_numeric(left["iraf_id"], errors="coerce").astype("Int64")
    merged = left.merge(iraf, on="iraf_id", how="inner", validate="one_to_one")
    if merged.empty:
        return merged
    merged = _add_alignment_columns(
        merged,
        zmag=float(zmag),
        frame_zeropoint=frame_zeropoint,
    )
    merged["iraf_dx_px"] = merged["iraf_x"] - merged["x"]
    merged["iraf_dy_px"] = merged["iraf_y"] - merged["y"]
    merged["iraf_dr_px"] = np.hypot(merged["iraf_dx_px"], merged["iraf_dy_px"])
    return merged


def summarize_mag_comparison(table: pd.DataFrame, *, label: str = "") -> dict[str, Any]:
    if table.empty:
        return {
            "label": label,
            "n_matched": 0,
            "n_finite_residual": 0,
            "median_delta_mag_raw": float("nan"),
            "mad_sigma_delta_mag_centered": float("nan"),
            "p95_abs_delta_mag_centered": float("nan"),
            "rms_delta_mag_centered": float("nan"),
            "median_delta_mag_units": float("nan"),
            "mad_sigma_delta_mag_units_centered": float("nan"),
            "frame_zeropoint_mag": float("nan"),
            "iraf_to_apex_zp_offset": float("nan"),
            "mad_sigma_delta_mag_zp_aligned": float("nan"),
            "p95_abs_delta_mag_zp_aligned": float("nan"),
            "rms_delta_mag_zp_aligned": float("nan"),
            "variance_delta_mag_zp_aligned": float("nan"),
            "median_delta_mag_apcorr_zp": float("nan"),
            "mad_sigma_delta_mag_apcorr_zp_centered": float("nan"),
            "p95_abs_delta_mag_apcorr_zp_centered": float("nan"),
            "rms_delta_mag_apcorr_zp_centered": float("nan"),
            "iraf_apcorr_zp_to_apex_offset": float("nan"),
            "mad_sigma_delta_mag_apcorr_zp_aligned": float("nan"),
            "p95_abs_delta_mag_apcorr_zp_aligned": float("nan"),
            "rms_delta_mag_apcorr_zp_aligned": float("nan"),
            "median_apex_mag_err": float("nan"),
            "median_iraf_merr": float("nan"),
            "median_quadrature_formal_err": float("nan"),
            "residual_mad_over_apex_mag_err": float("nan"),
            "residual_mad_over_iraf_merr": float("nan"),
            "residual_mad_over_quadrature_err": float("nan"),
            "position_rms_px": float("nan"),
            "outlier_fraction_abs_centered_gt_0p05": float("nan"),
        }
    centered = table["delta_mag_centered"].to_numpy(float)
    centered_units = table["delta_mag_units_centered"].to_numpy(float)
    if "delta_mag_zp_aligned" in table.columns:
        zp_aligned = table["delta_mag_zp_aligned"].to_numpy(float)
    else:
        zp_aligned = centered
    finite_centered = centered[np.isfinite(centered)]
    finite_zp_aligned = zp_aligned[np.isfinite(zp_aligned)]
    if "delta_mag_apcorr_zp_centered" in table.columns:
        apcorr_zp_centered = table["delta_mag_apcorr_zp_centered"].to_numpy(float)
    else:
        apcorr_zp_centered = np.asarray([], dtype=float)
    if "delta_mag_apcorr_zp_aligned" in table.columns:
        apcorr_zp_aligned = table["delta_mag_apcorr_zp_aligned"].to_numpy(float)
    else:
        apcorr_zp_aligned = np.asarray([], dtype=float)
    finite_apcorr_zp_centered = apcorr_zp_centered[np.isfinite(apcorr_zp_centered)]
    finite_apcorr_zp_aligned = apcorr_zp_aligned[np.isfinite(apcorr_zp_aligned)]
    abs_centered = np.abs(finite_centered)
    abs_zp_aligned = np.abs(finite_zp_aligned)
    abs_apcorr_zp_centered = np.abs(finite_apcorr_zp_centered)
    abs_apcorr_zp_aligned = np.abs(finite_apcorr_zp_aligned)
    apex_err = pd.to_numeric(table.get("mag_err", pd.Series(dtype=float)), errors="coerce").to_numpy(float)
    iraf_err = pd.to_numeric(table.get("iraf_merr", pd.Series(dtype=float)), errors="coerce").to_numpy(float)
    apex_err = apex_err[np.isfinite(apex_err) & (apex_err > 0)]
    iraf_err = iraf_err[np.isfinite(iraf_err) & (iraf_err > 0)]
    median_apex_err = float(np.median(apex_err)) if apex_err.size else float("nan")
    median_iraf_err = float(np.median(iraf_err)) if iraf_err.size else float("nan")
    if apex_err.size and iraf_err.size:
        aligned_errors = table[["mag_err", "iraf_merr"]].apply(pd.to_numeric, errors="coerce")
        quadrature = np.hypot(
            aligned_errors["mag_err"].to_numpy(float),
            aligned_errors["iraf_merr"].to_numpy(float),
        )
        quadrature = quadrature[np.isfinite(quadrature) & (quadrature > 0)]
    else:
        quadrature = np.asarray([], dtype=float)
    median_quadrature = float(np.median(quadrature)) if quadrature.size else float("nan")
    zp_mad = _mad_sigma(finite_zp_aligned)
    return {
        "label": label,
        "n_matched": int(len(table)),
        "n_finite_residual": int(finite_centered.size),
        "median_delta_mag_raw": float(np.nanmedian(table["delta_mag_raw"])),
        "mad_sigma_delta_mag_centered": _mad_sigma(finite_centered),
        "p95_abs_delta_mag_centered": float(np.nanpercentile(abs_centered, 95)) if abs_centered.size else float("nan"),
        "rms_delta_mag_centered": _rms(finite_centered),
        "median_delta_mag_units": float(np.nanmedian(table["delta_mag_units"])),
        "mad_sigma_delta_mag_units_centered": _mad_sigma(centered_units),
        "frame_zeropoint_mag": _median_or_nan(table["apex_mag_cal_frame"] - table["mag_inst"])
        if "apex_mag_cal_frame" in table.columns
        else float("nan"),
        "iraf_to_apex_zp_offset": _median_or_nan(table["iraf_to_apex_zp_offset"])
        if "iraf_to_apex_zp_offset" in table.columns
        else float("nan"),
        "mad_sigma_delta_mag_zp_aligned": zp_mad,
        "p95_abs_delta_mag_zp_aligned": (
            float(np.nanpercentile(abs_zp_aligned, 95)) if abs_zp_aligned.size else float("nan")
        ),
        "rms_delta_mag_zp_aligned": _rms(finite_zp_aligned),
        "variance_delta_mag_zp_aligned": float(np.var(finite_zp_aligned)) if finite_zp_aligned.size else float("nan"),
        "median_delta_mag_apcorr_zp": _median_or_nan(table.get("delta_mag_apcorr_zp", pd.Series(dtype=float))),
        "mad_sigma_delta_mag_apcorr_zp_centered": _mad_sigma(finite_apcorr_zp_centered),
        "p95_abs_delta_mag_apcorr_zp_centered": (
            float(np.nanpercentile(abs_apcorr_zp_centered, 95)) if abs_apcorr_zp_centered.size else float("nan")
        ),
        "rms_delta_mag_apcorr_zp_centered": _rms(finite_apcorr_zp_centered),
        "iraf_apcorr_zp_to_apex_offset": _median_or_nan(
            table.get("iraf_apcorr_zp_to_apex_offset", pd.Series(dtype=float))
        ),
        "mad_sigma_delta_mag_apcorr_zp_aligned": _mad_sigma(finite_apcorr_zp_aligned),
        "p95_abs_delta_mag_apcorr_zp_aligned": (
            float(np.nanpercentile(abs_apcorr_zp_aligned, 95)) if abs_apcorr_zp_aligned.size else float("nan")
        ),
        "rms_delta_mag_apcorr_zp_aligned": _rms(finite_apcorr_zp_aligned),
        "median_apex_mag_err": median_apex_err,
        "median_iraf_merr": median_iraf_err,
        "median_quadrature_formal_err": median_quadrature,
        "residual_mad_over_apex_mag_err": (
            float(zp_mad / median_apex_err) if np.isfinite(zp_mad) and np.isfinite(median_apex_err) else float("nan")
        ),
        "residual_mad_over_iraf_merr": (
            float(zp_mad / median_iraf_err) if np.isfinite(zp_mad) and np.isfinite(median_iraf_err) else float("nan")
        ),
        "residual_mad_over_quadrature_err": (
            float(zp_mad / median_quadrature)
            if np.isfinite(zp_mad) and np.isfinite(median_quadrature)
            else float("nan")
        ),
        "position_rms_px": _rms(table["iraf_dr_px"].to_numpy(float)) if "iraf_dr_px" in table.columns else float("nan"),
        "outlier_fraction_abs_centered_gt_0p05": (
            float(np.mean(abs_centered > 0.05)) if abs_centered.size else float("nan")
        ),
    }


def build_error_model_by_magnitude(
    table: pd.DataFrame,
    *,
    bin_width: float = 0.5,
    min_count: int = 5,
) -> pd.DataFrame:
    required = {"apex_mag_cal_frame", "delta_mag_zp_aligned", "mag_err", "iraf_merr"}
    if table.empty or not required <= set(table.columns):
        return pd.DataFrame()
    work = table.copy()
    for column in required:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work[
        np.isfinite(work["apex_mag_cal_frame"].to_numpy(float))
        & np.isfinite(work["delta_mag_zp_aligned"].to_numpy(float))
    ].copy()
    if work.empty:
        return pd.DataFrame()
    width = float(bin_width)
    work["mag_bin"] = (np.floor(work["apex_mag_cal_frame"] / width) * width + 0.5 * width).round(3)
    rows: list[dict[str, Any]] = []
    for mag_bin, group in work.groupby("mag_bin", sort=True):
        if len(group) < int(min_count):
            continue
        residual = group["delta_mag_zp_aligned"].to_numpy(float)
        residual = residual - _median_or_nan(residual)
        apex_err = group["mag_err"].to_numpy(float)
        iraf_err = group["iraf_merr"].to_numpy(float)
        quadrature = np.hypot(apex_err, iraf_err)
        rows.append(
            {
                "mag_bin": float(mag_bin),
                "n": int(len(group)),
                "mag_median": float(np.nanmedian(group["apex_mag_cal_frame"])),
                "apex_mag_err_median": _median_or_nan(apex_err),
                "iraf_merr_median": _median_or_nan(iraf_err),
                "quadrature_formal_err_median": _median_or_nan(quadrature),
                "residual_mad_sigma": _mad_sigma(residual),
                "residual_rms": _rms(residual),
                "residual_variance": float(np.nanvar(residual)),
            }
        )
    return pd.DataFrame(rows)


def _write_common_pyraf_setup(meta: IRAFFrameMetadata, config: IRAFCrosscheckConfig) -> str:
    calgorithm = "centroid" if config.recenter else "none"
    maxshift = max(0.0, min(5.0, meta.fwhm_px / 2.0)) if config.recenter else 0.0
    aperture = meta.r_ap
    annulus = meta.r_in
    dannulus = meta.r_out - meta.r_in
    return f"""
from pyraf import iraf
iraf.noao(); iraf.digiphot(); iraf.daophot()
for task_name in ["daofind", "phot", "datapars", "findpars", "centerpars", "fitskypars", "photpars"]:
    try:
        iraf.unlearn(task_name)
    except Exception:
        pass
iraf.datapars.noise = "poisson"
iraf.datapars.fwhmpsf = {meta.fwhm_px:.8f}
iraf.datapars.sigma = {meta.sigma_adu:.8f}
iraf.datapars.readnoise = {meta.rdnoise_e:.8f}
iraf.datapars.epadu = {meta.gain_e_per_adu:.8f}
iraf.datapars.itime = {meta.exptime:.8f}
iraf.centerpars.calgorithm = {calgorithm!r}
iraf.centerpars.maxshift = {maxshift:.8f}
iraf.centerpars.cbox = {max(3.0, 2.0 * meta.fwhm_px):.8f}
iraf.centerpars.cmaxiter = 10
iraf.fitskypars.salgorithm = "median"
iraf.fitskypars.annulus = {annulus:.8f}
iraf.fitskypars.dannulus = {dannulus:.8f}
iraf.fitskypars.smaxiter = 10
iraf.fitskypars.sloclip = 3.0
iraf.fitskypars.shiclip = 3.0
iraf.photpars.apertures = "{aperture:.4f}"
iraf.photpars.zmag = {float(config.zmag):.8f}
iraf.photpars.mkapert = "no"
"""


def _fixed_pyraf_script(
    *,
    image: str,
    coords: str,
    mag: str,
    txdump: str,
    log_path: str,
    meta: IRAFFrameMetadata,
    config: IRAFCrosscheckConfig,
) -> str:
    setup = _write_common_pyraf_setup(meta, config)
    return f"""# Auto-generated by APEX IRAF cross-check.
import json
import os
import time

{setup}

IMAGE = {image!r}
COORDS = {coords!r}
MAG = {mag!r}
TXDUMP = {txdump!r}
LOG_PATH = {log_path!r}

def safe_rm(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception:
        pass

def log(payload):
    with open(LOG_PATH, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\\n")

safe_rm(MAG)
safe_rm(TXDUMP)
t0 = time.time()
iraf.phot(IMAGE, coords=COORDS, output=MAG, verify="no", interactive="no", verbose="no", Stdout=os.devnull)
t1 = time.time()
iraf.txdump(MAG, fields="ID,XCENTER,YCENTER,MAG,MERR,MSKY,STDEV,NSKY", expr="yes", headers="no", Stdout=TXDUMP)
t2 = time.time()
log({{"mode": "phot_fixed_coords", "status": "ok", "phot_seconds": t1 - t0, "txdump_seconds": t2 - t1}})
"""


def _parse_threshold_grid(values: Iterable[float]) -> list[float]:
    parsed = sorted({float(v) for v in values if math.isfinite(float(v)) and float(v) > 0}, reverse=True)
    if not parsed:
        raise ValueError("daofind_threshold_grid must contain at least one positive threshold")
    return parsed


def _daofind_pyraf_script(
    *,
    image: str,
    output_dir: str,
    log_path: str,
    meta: IRAFFrameMetadata,
    config: IRAFCrosscheckConfig,
) -> str:
    setup = _write_common_pyraf_setup(meta, config)
    thresholds = _parse_threshold_grid(config.daofind_threshold_grid)
    max_by_ratio = max(1, int(round(meta.n_apex_sources * float(config.daofind_max_ratio_to_apex))))
    max_sources = max(1, min(int(config.daofind_max_sources), max_by_ratio))
    return f"""# Auto-generated by APEX IRAF daofind cross-check.
import json
import os
import time

{setup}

IMAGE = {image!r}
OUTPUT_DIR = {output_dir!r}
LOG_PATH = {log_path!r}
THRESHOLDS = {thresholds!r}
MAX_SOURCES = {max_sources}

def safe_rm(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception:
        pass

def count_coords(path):
    n = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
                if len(parts) < 2:
                    continue
                try:
                    float(parts[0]); float(parts[1])
                except Exception:
                    continue
                n += 1
    except Exception:
        return 0
    return n

def log(payload):
    with open(LOG_PATH, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\\n")

os.makedirs(OUTPUT_DIR, exist_ok=True)
unsafe_lower_thresholds = False
for threshold in THRESHOLDS:
    tag = str(threshold).replace(".", "p")
    coo = os.path.join(OUTPUT_DIR, f"daofind_thr_{{tag}}.coo")
    mag = os.path.join(OUTPUT_DIR, f"daofind_thr_{{tag}}.mag")
    txt = os.path.join(OUTPUT_DIR, f"daofind_thr_{{tag}}.txt")
    if unsafe_lower_thresholds:
        log({{"mode": "daofind_phot", "threshold": threshold, "status": "skipped_after_unsafe_threshold", "max_sources": MAX_SOURCES}})
        continue
    safe_rm(coo); safe_rm(mag); safe_rm(txt)
    iraf.findpars.threshold = float(threshold)
    iraf.findpars.nsigma = 1.5
    iraf.findpars.sharplo = 0.2
    iraf.findpars.sharphi = 1.2
    iraf.findpars.roundlo = -1.0
    iraf.findpars.roundhi = 1.0
    t0 = time.time()
    iraf.daofind(IMAGE, output=coo, verify="no", interactive="no", verbose="no")
    t1 = time.time()
    n_detected = count_coords(coo)
    if n_detected > MAX_SOURCES:
        unsafe_lower_thresholds = True
        log({{
            "mode": "daofind_phot",
            "threshold": threshold,
            "status": "too_many_sources",
            "n_detected": n_detected,
            "max_sources": MAX_SOURCES,
            "daofind_seconds": t1 - t0,
            "coo": coo,
        }})
        continue
    iraf.phot(IMAGE, coords=coo, output=mag, verify="no", interactive="no", verbose="no", Stdout=os.devnull)
    t2 = time.time()
    iraf.txdump(mag, fields="ID,XCENTER,YCENTER,MAG,MERR,MSKY,STDEV,NSKY", expr="yes", headers="no", Stdout=txt)
    t3 = time.time()
    log({{
        "mode": "daofind_phot",
        "threshold": threshold,
        "status": "ok",
        "n_detected": n_detected,
        "max_sources": MAX_SOURCES,
        "daofind_seconds": t1 - t0,
        "phot_seconds": t2 - t1,
        "txdump_seconds": t3 - t2,
        "coo": coo,
        "mag": mag,
        "txdump": txt,
    }})
"""


def _run_script(script_path: Path, runtime_cmd: list[str], *, use_wsl: bool) -> None:
    script_arg = _runtime_path(script_path, use_wsl=use_wsl)
    completed = subprocess.run(
        [*runtime_cmd, script_arg],
        cwd=str(script_path.parent),
        text=True,
        capture_output=True,
        check=False,
    )
    (script_path.parent / f"{script_path.stem}_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (script_path.parent / f"{script_path.stem}_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"PyRAF script failed ({completed.returncode}): {script_path}\n"
            f"stderr: {completed.stderr[-4000:]}"
        )


def _write_fixed_plots(table: pd.DataFrame, output_dir: Path) -> None:
    if table.empty:
        return
    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    ax.scatter(table["apex_mag_iraf_zmag"], table["delta_mag_centered"], s=10, alpha=0.55)
    ax.axhline(0.0, color="black", linewidth=1)
    ax.set(xlabel="APEX instrumental magnitude on IRAF zmag", ylabel="APEX - IRAF, median removed (mag)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "fixed_delta_vs_mag.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    finite = table["delta_mag_centered"].to_numpy(float)
    finite = finite[np.isfinite(finite)]
    if finite.size:
        ax.hist(finite, bins=min(40, max(10, int(np.sqrt(finite.size)))), alpha=0.85)
    ax.axvline(0.0, color="black", linewidth=1)
    ax.set(xlabel="APEX - IRAF, median removed (mag)", ylabel="N")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "fixed_delta_hist.png", dpi=180)
    plt.close(fig)

    if {"apex_mag_cal_frame", "iraf_mag_aligned_to_apex"} <= set(table.columns):
        x = table["apex_mag_cal_frame"].to_numpy(float)
        y = table["iraf_mag_aligned_to_apex"].to_numpy(float)
        keep = np.isfinite(x) & np.isfinite(y)
        if np.any(keep):
            x = x[keep]
            y = y[keep]
            lo = float(min(np.min(x), np.min(y)))
            hi = float(max(np.max(x), np.max(y)))
            pad = max(0.02, 0.04 * (hi - lo))
            residual = table.loc[keep, "delta_mag_zp_aligned"].to_numpy(float)
            mad = _mad_sigma(residual)
            p95 = float(np.nanpercentile(np.abs(residual[np.isfinite(residual)]), 95))
            fig, ax = plt.subplots(figsize=(6.5, 6.0))
            ax.scatter(x, y, s=11, alpha=0.55)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="black", linewidth=1)
            ax.set(
                xlabel="APEX calibrated magnitude (Step 10 ZP)",
                ylabel="IRAF magnitude aligned to APEX ZP",
                xlim=(lo - pad, hi + pad),
                ylim=(lo - pad, hi + pad),
            )
            ax.invert_xaxis()
            ax.invert_yaxis()
            ax.set_aspect("equal", adjustable="box")
            ax.grid(alpha=0.25)
            ax.text(
                0.04,
                0.96,
                f"N={len(x)}\nMAD={mad:.4f} mag\nP95={p95:.4f} mag",
                transform=ax.transAxes,
                va="top",
                ha="left",
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.75"},
                fontsize=9,
            )
            fig.tight_layout()
            fig.savefig(output_dir / "fixed_apex_vs_iraf_zp_aligned.png", dpi=180)
            plt.close(fig)

    if {"apex_mag_cal_frame", "delta_mag_zp_aligned"} <= set(table.columns):
        finite = table[
            np.isfinite(table["apex_mag_cal_frame"].to_numpy(float))
            & np.isfinite(table["delta_mag_zp_aligned"].to_numpy(float))
        ]
        if not finite.empty:
            fig, ax = plt.subplots(figsize=(7.8, 5.0))
            ax.scatter(
                finite["apex_mag_cal_frame"],
                finite["delta_mag_zp_aligned"],
                s=10,
                alpha=0.55,
            )
            ax.axhline(0.0, color="black", linewidth=1)
            ax.set(
                xlabel="APEX calibrated magnitude (Step 10 ZP)",
                ylabel="APEX - aligned IRAF (mag)",
            )
            ax.invert_xaxis()
            ax.grid(alpha=0.25)
            fig.tight_layout()
            fig.savefig(output_dir / "fixed_zp_aligned_residual_vs_mag.png", dpi=180)
            plt.close(fig)

            values = finite["delta_mag_zp_aligned"].to_numpy(float)
            fig, ax = plt.subplots(figsize=(7.0, 4.5))
            ax.hist(values, bins=min(40, max(10, int(np.sqrt(len(values))))), alpha=0.85)
            ax.axvline(0.0, color="black", linewidth=1)
            ax.set(xlabel="APEX - aligned IRAF (mag)", ylabel="N")
            ax.grid(alpha=0.25)
            fig.tight_layout()
            fig.savefig(output_dir / "fixed_zp_aligned_residual_hist.png", dpi=180)
            plt.close(fig)

    if {"snr", "delta_mag_zp_aligned"} <= set(table.columns):
        finite = table[
            np.isfinite(table["snr"].to_numpy(float))
            & (table["snr"].to_numpy(float) > 0)
            & np.isfinite(table["delta_mag_zp_aligned"].to_numpy(float))
        ]
        if not finite.empty:
            fig, ax = plt.subplots(figsize=(7.5, 4.8))
            ax.scatter(finite["snr"], finite["delta_mag_zp_aligned"], s=10, alpha=0.55)
            ax.axhline(0.0, color="black", linewidth=1)
            ax.set_xscale("log")
            ax.set(xlabel="APEX S/N", ylabel="APEX - aligned IRAF (mag)")
            ax.grid(alpha=0.25, which="both")
            fig.tight_layout()
            fig.savefig(output_dir / "fixed_zp_aligned_residual_vs_snr.png", dpi=180)
            plt.close(fig)

    error_bins = build_error_model_by_magnitude(table)
    if not error_bins.empty:
        error_bins.to_csv(output_dir / "fixed_error_model_by_mag.csv", index=False)
        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.plot(
            error_bins["mag_median"],
            error_bins["apex_mag_err_median"],
            marker="o",
            linewidth=1.5,
            label="APEX median formal error",
        )
        ax.plot(
            error_bins["mag_median"],
            error_bins["iraf_merr_median"],
            marker="s",
            linewidth=1.5,
            label="IRAF median MERR",
        )
        ax.plot(
            error_bins["mag_median"],
            error_bins["residual_mad_sigma"],
            marker="^",
            linewidth=1.8,
            label="APEX-IRAF residual MAD",
        )
        ax.set(
            xlabel="APEX calibrated magnitude (Step 10 ZP)",
            ylabel="Magnitude error / scatter",
        )
        ax.invert_xaxis()
        ax.set_yscale("log")
        ax.grid(alpha=0.25, which="both")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / "fixed_error_model_vs_mag.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.plot(
            error_bins["mag_median"],
            error_bins["residual_mad_sigma"] / error_bins["apex_mag_err_median"],
            marker="o",
            linewidth=1.5,
            label="Residual MAD / APEX error",
        )
        ax.plot(
            error_bins["mag_median"],
            error_bins["residual_mad_sigma"] / error_bins["iraf_merr_median"],
            marker="s",
            linewidth=1.5,
            label="Residual MAD / IRAF MERR",
        )
        ax.axhline(1.0, color="black", linewidth=1)
        ax.set(
            xlabel="APEX calibrated magnitude (Step 10 ZP)",
            ylabel="Scatter-to-formal-error ratio",
        )
        ax.invert_xaxis()
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / "fixed_residual_to_formal_error_ratio.png", dpi=180)
        plt.close(fig)


def _jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return rows


def _match_iraf_to_apex(
    iraf: pd.DataFrame,
    apex: pd.DataFrame,
    *,
    radius_px: float,
    zmag: float = 25.0,
) -> pd.DataFrame:
    if iraf.empty or apex.empty:
        return pd.DataFrame()
    iraf_xy = iraf[["iraf_x", "iraf_y"]].to_numpy(float)
    apex_xy = apex[["x", "y"]].to_numpy(float)
    candidates: list[tuple[float, int, int]] = []
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(apex_xy)
        lists = tree.query_ball_point(iraf_xy, r=float(radius_px))
        for iraf_i, apex_indices in enumerate(lists):
            for apex_i in apex_indices:
                dist = float(np.hypot(*(iraf_xy[iraf_i] - apex_xy[apex_i])))
                candidates.append((dist, iraf_i, int(apex_i)))
    except Exception:
        for iraf_i, xy in enumerate(iraf_xy):
            dist = np.hypot(apex_xy[:, 0] - xy[0], apex_xy[:, 1] - xy[1])
            nearby = np.where(dist <= float(radius_px))[0]
            for apex_i in nearby:
                candidates.append((float(dist[apex_i]), iraf_i, int(apex_i)))

    used_iraf: set[int] = set()
    used_apex: set[int] = set()
    pairs: list[dict[str, Any]] = []
    for dist, iraf_i, apex_i in sorted(candidates, key=lambda item: item[0]):
        if iraf_i in used_iraf or apex_i in used_apex:
            continue
        used_iraf.add(iraf_i)
        used_apex.add(apex_i)
        left = apex.iloc[apex_i].to_dict()
        right = iraf.iloc[iraf_i].to_dict()
        row = {**left, **right}
        row["match_distance_px"] = dist
        pairs.append(row)
    if not pairs:
        return pd.DataFrame()
    matched = pd.DataFrame(pairs)
    matched["apex_mag_iraf_zmag"] = float("nan")
    if "mag_inst" in matched.columns:
        matched["apex_mag_iraf_zmag"] = float(zmag) + _num(matched["mag_inst"])
    return matched


def _summarize_daofind_threshold(
    *,
    log_row: dict[str, Any],
    matched: pd.DataFrame,
    iraf_rows: int,
    apex_rows: int,
    zmag: float,
    frame_zeropoint: float | None,
) -> dict[str, Any]:
    out = dict(log_row)
    out["n_iraf_txdump"] = int(iraf_rows)
    out["n_apex_reference"] = int(apex_rows)
    out["n_matched"] = int(len(matched))
    out["recall_vs_apex_reference"] = float(len(matched) / apex_rows) if apex_rows else float("nan")
    out["matched_fraction_of_iraf"] = float(len(matched) / iraf_rows) if iraf_rows else float("nan")
    if not matched.empty and "mag_inst" in matched.columns and "iraf_mag" in matched.columns:
        comparison = _add_alignment_columns(
            matched.copy(),
            zmag=float(zmag),
            frame_zeropoint=frame_zeropoint,
        )
        comparison["iraf_dr_px"] = comparison["match_distance_px"]
        out.update(summarize_mag_comparison(comparison, label=f"threshold={log_row.get('threshold')}"))
    return out


def _write_daofind_plots(summary: pd.DataFrame, output_dir: Path) -> None:
    if summary.empty or "threshold" not in summary.columns:
        return
    work = summary.copy()
    work["threshold"] = pd.to_numeric(work["threshold"], errors="coerce")
    work["n_detected"] = pd.to_numeric(work.get("n_detected", np.nan), errors="coerce")
    work = work[np.isfinite(work["threshold"])]
    if work.empty:
        return
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(work["threshold"], work["n_detected"], marker="o", linewidth=1.5)
    ax.invert_xaxis()
    ax.set(xlabel="IRAF findpars.threshold", ylabel="Detected sources")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "daofind_sources_vs_threshold.png", dpi=180)
    plt.close(fig)

    if "recall_vs_apex_reference" in work.columns:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        ax.plot(work["threshold"], work["recall_vs_apex_reference"], marker="o", label="Recall vs APEX")
        ax.plot(work["threshold"], work["matched_fraction_of_iraf"], marker="s", label="Matched fraction of IRAF")
        ax.invert_xaxis()
        ax.set(xlabel="IRAF findpars.threshold", ylabel="Fraction")
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / "daofind_match_fractions.png", dpi=180)
        plt.close(fig)


def _write_report(
    output_dir: Path,
    *,
    config: IRAFCrosscheckConfig,
    meta: IRAFFrameMetadata,
    fixed_summary: dict[str, Any] | None,
    daofind_summary: pd.DataFrame | None,
) -> None:
    lines = [
        "# IRAF Cross-Check Report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope",
        "",
        "This package compares APEX Step 7 aperture photometry against IRAF/DAOPHOT through PyRAF.",
        "The fixed-coordinate mode tests the photometry implementation at the same star positions.",
        "The daofind mode is a guarded detector sanity check, not a full pipeline replacement.",
        "",
        "## Frame",
        "",
        f"- Filter: {meta.filter_name}",
        f"- FWHM: {meta.fwhm_px:.3f} px",
        f"- Aperture: r={meta.r_ap:.3f} px, annulus={meta.r_in:.3f}-{meta.r_out:.3f} px",
        f"- Gain/RD noise: {meta.gain_e_per_adu:.6g} e-/ADU, {meta.rdnoise_e:.6g} e-",
        f"- Background sigma: {meta.sigma_adu:.6g} ADU",
        f"- Step 10 frame zeropoint: {meta.frame_zeropoint_mag:.4f} mag",
        "",
        "## Fixed-Coordinate Photometry",
        "",
    ]
    if fixed_summary:
        lines.extend(
            [
                f"- Matched sources: {fixed_summary['n_matched']}",
                f"- Finite magnitude residuals: {fixed_summary['n_finite_residual']}",
                f"- Median raw offset, APEX - IRAF: {fixed_summary['median_delta_mag_raw']:.4f} mag",
                f"- IRAF-to-APEX fitted ZP offset: {fixed_summary['iraf_to_apex_zp_offset']:.4f} mag",
                f"- Robust scatter after median removal: {fixed_summary['mad_sigma_delta_mag_centered']:.4f} mag",
                f"- Robust scatter after ZP alignment: {fixed_summary['mad_sigma_delta_mag_zp_aligned']:.4f} mag",
                f"- IRAF apcorr+ZP-equivalent residual MAD: {fixed_summary['mad_sigma_delta_mag_apcorr_zp_centered']:.4f} mag",
                f"- Median formal error: APEX={fixed_summary['median_apex_mag_err']:.4f} mag, IRAF={fixed_summary['median_iraf_merr']:.4f} mag",
                f"- Residual MAD / formal error: APEX={fixed_summary['residual_mad_over_apex_mag_err']:.3f}, IRAF={fixed_summary['residual_mad_over_iraf_merr']:.3f}",
                f"- 95th percentile absolute residual: {fixed_summary['p95_abs_delta_mag_centered']:.4f} mag",
                f"- >0.05 mag centered residual fraction: {fixed_summary['outlier_fraction_abs_centered_gt_0p05']:.3f}",
            ]
        )
    else:
        lines.append("- Not run")

    lines.extend(["", "## IRAF daofind", ""])
    if daofind_summary is not None and not daofind_summary.empty:
        for _, row in daofind_summary.iterrows():
            lines.append(
                f"- threshold={float(row['threshold']):.2f}: status={row['status']}, "
                f"detected={row.get('n_detected', 'n/a')}, matched={int(row.get('n_matched', 0))}"
            )
    else:
        lines.append("- Not run")

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "- phot_fixed_coords/fixed_apex_vs_iraf_zp_aligned.png",
            "- phot_fixed_coords/fixed_zp_aligned_residual_vs_mag.png",
            "- phot_fixed_coords/fixed_zp_aligned_residual_hist.png",
            "- phot_fixed_coords/fixed_zp_aligned_residual_vs_snr.png",
            "- phot_fixed_coords/fixed_error_model_vs_mag.png",
            "- phot_fixed_coords/fixed_residual_to_formal_error_ratio.png",
            "- daofind_phot/daofind_sources_vs_threshold.png",
            "- daofind_phot/daofind_match_fractions.png",
        ]
    )

    lines.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- APEX Step 7 applies aperture correction and stores flux in an electron convention.",
            "- IRAF `phot` magnitudes can therefore have a constant offset from APEX.",
            "- The report reads Step 10 `zp_frame` when available and aligns IRAF to APEX with a fitted median offset.",
            "- The apcorr+ZP-equivalent IRAF columns apply the same gain, aperture correction, and frame ZP scale to IRAF `phot` output.",
            "- The primary cross-check metric is the ZP-aligned magnitude scatter.",
            "- The formal-error comparison is an error-model diagnostic; shared image noise makes APEX-IRAF residuals smaller than either formal error.",
            "- This is a reference agreement test, not evidence that aperture photometry is intrinsically more accurate than IRAF.",
            "- The daofind grid starts at high threshold and stops lower thresholds after an unsafe source count.",
            "",
            "## Configuration",
            "",
            "```json",
            json.dumps(asdict(config), indent=2, sort_keys=True),
            "```",
            "",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run_iraf_crosscheck(config: IRAFCrosscheckConfig) -> Path:
    mode = str(config.mode).strip().lower()
    if mode not in {FIXED_MODE, DAOFIND_MODE, BOTH_MODE}:
        raise ValueError(f"Unsupported IRAF cross-check mode: {config.mode}")

    input_fits = Path(config.input_fits).expanduser().resolve()
    step7_path = Path(config.step7_tsv).expanduser().resolve()
    if not input_fits.is_file():
        raise FileNotFoundError(f"Input FITS is missing: {input_fits}")
    if not step7_path.is_file():
        raise FileNotFoundError(f"Step 7 TSV is missing: {step7_path}")

    output_root = Path(config.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = _repo_root() / output_root
    output_root = output_root.resolve()
    _clean_output_dir(output_root, overwrite=bool(config.overwrite))

    step7 = _read_step7(step7_path)
    index_row = _lookup_index_row(step7_path, input_fits)
    meta = _resolve_metadata(config, step7, index_row, step7_path, input_fits)
    reference = load_apex_reference_table(step7_path, min_snr=float(config.min_snr))
    if reference.empty:
        raise RuntimeError("No usable APEX Step 7 reference sources after quality filtering")
    reference.to_csv(output_root / "apex_reference_all.csv", index=False)

    runtime_cmd = list(config.runtime_cmd or _default_runtime_cmd())
    use_wsl = _uses_wsl(runtime_cmd)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_fits": str(input_fits),
        "step7_tsv": str(step7_path),
        "effective_config": asdict(config),
        "frame_metadata": asdict(meta),
        "runtime_cmd": runtime_cmd,
    }
    (output_root / "iraf_crosscheck_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    fixed_summary: dict[str, Any] | None = None
    if mode in {FIXED_MODE, BOTH_MODE}:
        fixed_dir = output_root / "phot_fixed_coords"
        fixed_dir.mkdir(parents=True, exist_ok=True)
        fixed_reference = select_fixed_coordinate_sources(
            reference,
            max_sources=int(config.max_fixed_coords),
            seed=int(config.sample_seed),
        )
        fixed_reference.to_csv(fixed_dir / "fixed_reference.csv", index=False)
        coords_path = write_iraf_coords(fixed_reference, fixed_dir / "fixed_coords.coo")
        log_path = fixed_dir / "pyraf_fixed_log.jsonl"
        log_path.write_text("", encoding="utf-8")
        script_path = fixed_dir / "run_fixed_phot.py"
        script_path.write_text(
            _fixed_pyraf_script(
                image=_runtime_path(input_fits, use_wsl=use_wsl),
                coords=_runtime_path(coords_path, use_wsl=use_wsl),
                mag=_runtime_path(fixed_dir / "fixed_phot.mag", use_wsl=use_wsl),
                txdump=_runtime_path(fixed_dir / "fixed_phot.txt", use_wsl=use_wsl),
                log_path=_runtime_path(log_path, use_wsl=use_wsl),
                meta=meta,
                config=config,
            ),
            encoding="utf-8",
            newline="\n",
        )
        _run_script(script_path, runtime_cmd, use_wsl=use_wsl)
        iraf = parse_txdump(fixed_dir / "fixed_phot.txt")
        iraf.to_csv(fixed_dir / "iraf_fixed_txdump.csv", index=False)
        fixed_comparison = build_fixed_comparison(
            fixed_reference,
            iraf,
            zmag=float(config.zmag),
            frame_zeropoint=meta.frame_zeropoint_mag,
        )
        fixed_comparison.to_csv(fixed_dir / "fixed_comparison.csv", index=False)
        fixed_summary = summarize_mag_comparison(fixed_comparison, label="phot_fixed_coords")
        fixed_summary.update({"n_reference": int(len(fixed_reference)), "n_iraf_rows": int(len(iraf))})
        (fixed_dir / "fixed_summary.json").write_text(
            json.dumps(fixed_summary, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        _write_fixed_plots(fixed_comparison, fixed_dir)

    daofind_summary: pd.DataFrame | None = None
    if mode in {DAOFIND_MODE, BOTH_MODE}:
        dao_dir = output_root / "daofind_phot"
        dao_dir.mkdir(parents=True, exist_ok=True)
        log_path = dao_dir / "pyraf_daofind_log.jsonl"
        log_path.write_text("", encoding="utf-8")
        script_path = dao_dir / "run_daofind_phot.py"
        script_path.write_text(
            _daofind_pyraf_script(
                image=_runtime_path(input_fits, use_wsl=use_wsl),
                output_dir=_runtime_path(dao_dir, use_wsl=use_wsl),
                log_path=_runtime_path(log_path, use_wsl=use_wsl),
                meta=meta,
                config=config,
            ),
            encoding="utf-8",
            newline="\n",
        )
        _run_script(script_path, runtime_cmd, use_wsl=use_wsl)
        rows = []
        for log_row in _jsonl(log_path):
            threshold = log_row.get("threshold")
            tag = str(threshold).replace(".", "p")
            txt = dao_dir / f"daofind_thr_{tag}.txt"
            iraf = parse_txdump(txt)
            matched = _match_iraf_to_apex(
                iraf,
                reference,
                radius_px=float(config.daofind_match_radius_px),
                zmag=float(config.zmag),
            )
            if not matched.empty:
                matched = _add_alignment_columns(
                    matched.copy(),
                    zmag=float(config.zmag),
                    frame_zeropoint=meta.frame_zeropoint_mag,
                )
                matched["iraf_dr_px"] = matched["match_distance_px"]
                matched.to_csv(dao_dir / f"matches_thr_{tag}.csv", index=False)
            rows.append(
                _summarize_daofind_threshold(
                    log_row=log_row,
                    matched=matched,
                    iraf_rows=len(iraf),
                    apex_rows=len(reference),
                    zmag=float(config.zmag),
                    frame_zeropoint=meta.frame_zeropoint_mag,
                )
            )
        daofind_summary = pd.DataFrame(rows)
        daofind_summary.to_csv(dao_dir / "daofind_threshold_summary.csv", index=False)
        _write_daofind_plots(daofind_summary, dao_dir)

    _write_report(
        output_root,
        config=config,
        meta=meta,
        fixed_summary=fixed_summary,
        daofind_summary=daofind_summary,
    )
    return output_root
