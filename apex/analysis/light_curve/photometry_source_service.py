"""Resolve one non-mixed photometry table source for LC processing."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from astropy.time import Time

from apex.utils.common_helpers import normalize_filter_key
from apex.utils.io_utils import load_headers_table, load_night_assignments
from apex.utils.photometry_loader import load_frame_photometry
from apex.utils.photometry_provenance import build_photometry_provenance
from apex.utils.step_paths import forced_phot_input_dir
from apex.utils.step_paths_cmd import step8_psf_dir


def _step_data(project_state: Any, result_dir: Path) -> dict:
    if project_state is not None and hasattr(project_state, "get_step_data"):
        try:
            value = project_state.get_step_data("psf_photometry")
            if isinstance(value, dict):
                return value
        except Exception:
            pass
    state_path = result_dir / "project_state.json"
    try:
        payload = json.loads(state_path.read_text(encoding="utf-8"))
        value = payload.get("step_data", {}).get("psf_photometry", {})
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def _signature_matches(path: Path, expected: dict | None) -> bool:
    if not isinstance(expected, dict) or not path.exists():
        return False
    try:
        stat = path.stat()
        return (
            int(expected.get("size", -1)) == int(stat.st_size)
            and int(expected.get("mtime_ns", -1)) == int(stat.st_mtime_ns)
        )
    except (OSError, TypeError, ValueError):
        return False


def resolve_lightcurve_photometry_source(
    result_dir: Path | str,
    project_state: Any = None,
) -> dict:
    """Use complete PSF output unless the PSF step was skipped."""
    root = Path(result_dir)
    aperture_dir = forced_phot_input_dir(root)
    aperture_index = aperture_dir / "photometry_index.csv"
    state = _step_data(project_state, root)
    skipped = bool(state.get("skip_psf", False))
    fallback = {
        **build_photometry_provenance("aperture", "mag", "mag_err"),
        "source": "aperture",
        "directory": aperture_dir,
        "index_path": aperture_index,
        "psf_available": False,
        "psf_skipped": skipped,
        "reason": "PSF was skipped" if skipped else "PSF output unavailable",
    }
    if skipped:
        return fallback

    psf_dir = step8_psf_dir(root)
    psf_index = psf_dir / "photometry_index.csv"
    signature_path = psf_dir / "psf_output_signature.json"
    if not psf_index.exists() or not signature_path.exists():
        return fallback
    try:
        signature = json.loads(signature_path.read_text(encoding="utf-8"))
        signature_mtime_ns = signature_path.stat().st_mtime_ns
        index = pd.read_csv(psf_index)
    except Exception as exc:
        fallback["reason"] = f"Invalid PSF metadata: {exc}"
        return fallback

    frames = [Path(str(value)).name for value in signature.get("frames", [])]
    if not frames or not {"file", "filter"} <= set(index.columns):
        fallback["reason"] = "PSF metadata is incomplete"
        return fallback
    indexed_frames = [Path(str(value)).name for value in index["file"].tolist()]
    if len(indexed_frames) != len(frames) or set(indexed_frames) != set(frames):
        fallback["reason"] = "PSF frame set does not match its signature"
        return fallback
    step7_signature = signature.get("inputs", {}).get("step7_index")
    if not _signature_matches(aperture_index, step7_signature):
        fallback["reason"] = "Step 7 changed after the PSF run"
        return fallback

    frame_inputs = signature.get("inputs", {}).get("frames", [])
    frame_input_map = {
        Path(str(item.get("file", ""))).name: item
        for item in frame_inputs
        if isinstance(item, dict)
    }
    if set(frame_input_map) != set(frames):
        fallback["reason"] = "PSF input signatures are incomplete"
        return fallback

    required = {"det_uid", "mag_psf", "mag_psf_err", "flags_psf"}
    for frame in frames:
        aperture_table = aperture_dir / f"photometry_{frame}.tsv"
        if not _signature_matches(
            aperture_table, frame_input_map[frame].get("step7_tsv")
        ):
            fallback["reason"] = f"Step 7 table changed after the PSF run: {frame}"
            return fallback
        table = psf_dir / f"photometry_{frame}.tsv"
        if not table.exists():
            fallback["reason"] = f"Missing PSF table for {frame}"
            return fallback
        try:
            table_stat = table.stat()
        except OSError as exc:
            fallback["reason"] = f"Cannot stat PSF table for {frame}: {exc}"
            return fallback
        if table_stat.st_size <= 0:
            fallback["reason"] = f"Missing PSF table for {frame}"
            return fallback
        # The completion signature is written only after every output schema
        # has been validated. Recheck a table only when it changed afterwards.
        if table_stat.st_mtime_ns > signature_mtime_ns:
            try:
                columns = set(pd.read_csv(table, sep="\t", nrows=0).columns)
            except Exception as exc:
                fallback["reason"] = f"Cannot read PSF table for {frame}: {exc}"
                return fallback
            if not required <= columns:
                fallback["reason"] = f"Incomplete PSF table for {frame}"
                return fallback

    return {
        **build_photometry_provenance("psf", "mag_psf", "mag_psf_err"),
        "source": "psf",
        "directory": psf_dir,
        "index_path": psf_index,
        "psf_available": True,
        "psf_skipped": False,
        "reason": f"Complete PSF output for {len(frames)} frames",
        "frames": frames,
        "experimental": True,
    }


def load_lightcurve_frame_photometry(
    result_dir: Path | str,
    filename: str,
    source_info: dict,
    filter_name: str = "",
) -> pd.DataFrame | None:
    """Load aperture or PSF data and expose the common ID/mag/error schema."""
    source = str(source_info.get("source", "aperture"))
    if source != "psf":
        frame = load_frame_photometry(Path(result_dir), filename, filter_name)
        if frame is not None and not frame.empty:
            frame = frame.copy()
            frame["photometry_source"] = "aperture"
            frame["mag_input_column"] = str(source_info.get("mag_column", "mag"))
            frame["mag_error_input_column"] = str(
                source_info.get("mag_error_column", "mag_err")
            )
        return frame

    aperture = load_frame_photometry(Path(result_dir), filename, filter_name)
    if aperture is None or aperture.empty or "det_uid" not in aperture.columns:
        return None

    table = Path(source_info["directory"]) / f"photometry_{Path(filename).name}.tsv"
    try:
        psf = pd.read_csv(table, sep="\t")
    except Exception:
        return None
    if psf.empty or "det_uid" not in psf.columns:
        return None

    frame = aperture.copy()
    frame_uid = pd.to_numeric(frame["det_uid"], errors="coerce").astype("Int64")
    psf = psf.copy()
    psf["det_uid"] = pd.to_numeric(psf["det_uid"], errors="coerce").astype("Int64")
    psf = psf[psf["det_uid"].notna() & (psf["det_uid"] >= 0)]
    psf = psf.drop_duplicates(subset="det_uid", keep="first").set_index("det_uid")

    def _map_psf(column: str) -> pd.Series:
        if column not in psf.columns:
            return pd.Series(np.nan, index=frame.index, dtype=float)
        return pd.to_numeric(frame_uid.map(psf[column]), errors="coerce")

    psf_mag = _map_psf("mag_psf")
    psf_err = _map_psf("mag_psf_err")
    flags = _map_psf("flags_psf")
    clean = np.isfinite(flags) & (flags == 0) & np.isfinite(psf_mag)
    frame["mag"] = psf_mag.where(clean, np.nan)
    frame["mag_err"] = psf_err.where(clean, np.nan)
    for output_column, psf_column in (("x", "x_fit"), ("y", "y_fit")):
        fitted = _map_psf(psf_column)
        if output_column in frame.columns:
            existing = pd.to_numeric(frame[output_column], errors="coerce")
            frame[output_column] = fitted.where(np.isfinite(fitted), existing)
        else:
            frame[output_column] = fitted
    frame["flags_psf"] = flags
    for column in (
        "qfit",
        "cfit",
        "reduced_chi2",
        "n_pixels_fit",
        "snr_psf",
        "iter_found",
        "psf_core_r_px",
        "psf_core_cut_px",
    ):
        frame[column] = _map_psf(column)
    frame["psf_qc_clean"] = clean
    frame["photometry_source"] = "psf"
    frame["mag_input_column"] = "mag_psf"
    frame["mag_error_input_column"] = "mag_psf_err"
    return frame


def _time_value(value: Any) -> float:
    try:
        numeric = float(value)
        if np.isfinite(numeric):
            return numeric
    except (TypeError, ValueError):
        pass
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none"}:
        return float("nan")
    try:
        return float(Time(text).jd)
    except Exception:
        return float("nan")


def _first_finite(*values: Any) -> float:
    for value in values:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric):
            return numeric
    return float("nan")


def _bad_measurement_mask(frame: pd.DataFrame) -> pd.Series:
    bad = pd.Series(False, index=frame.index)
    for column in ("bad_phot_flag", "is_saturated", "is_nonlinear", "off_frame_flag"):
        if column not in frame.columns:
            continue
        values = frame[column]
        if pd.api.types.is_bool_dtype(values):
            bad |= values.fillna(False)
        else:
            bad |= values.astype(str).str.strip().str.lower().isin(
                {"true", "1", "yes", "y"}
            )
    return bad


def load_filter_photometry_timeseries(
    result_dir: Path | str,
    filter_name: str,
    project_state: Any = None,
    *,
    should_stop: Callable[[], bool] | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Load one filter into a common long frame/star time-series table."""
    root = Path(result_dir)
    source = resolve_lightcurve_photometry_source(root, project_state)
    if should_stop is not None and should_stop():
        return pd.DataFrame(), source
    aperture_index = forced_phot_input_dir(root) / "photometry_index.csv"
    try:
        index = pd.read_csv(aperture_index)
    except Exception:
        return pd.DataFrame(), source
    if "file" not in index.columns:
        return pd.DataFrame(), source

    filter_column = "filter" if "filter" in index.columns else (
        "FILTER" if "FILTER" in index.columns else None
    )
    wanted_filter = normalize_filter_key(filter_name)
    if filter_column is not None and wanted_filter:
        normalized = index[filter_column].astype(str).map(normalize_filter_key)
        index = index[normalized == wanted_filter].copy()
    if source.get("source") == "psf":
        allowed = {Path(str(value)).name for value in source.get("frames", [])}
        basenames = index["file"].astype(str).map(lambda value: Path(value).name)
        index = index[basenames.isin(allowed)].copy()
    if index.empty:
        return pd.DataFrame(), source

    headers = load_headers_table(root)
    header_rows: dict[str, pd.Series] = {}
    if not headers.empty and "Filename" in headers.columns:
        header_rows = {
            Path(str(row["Filename"])).name: row
            for _, row in headers.iterrows()
        }
    night_assignments = {
        Path(str(key)).name: int(value)
        for key, value in load_night_assignments(root).items()
    }

    output_frames: list[pd.DataFrame] = []
    fallback_nights: dict[str, int] = {}
    for frame_order, (_, index_row) in enumerate(index.iterrows()):
        if should_stop is not None and should_stop():
            break
        filename = Path(str(index_row.get("file", ""))).name
        if not filename:
            continue
        frame = load_lightcurve_frame_photometry(
            root, filename, source, str(index_row.get(filter_column, filter_name) if filter_column else filter_name)
        )
        if should_stop is not None and should_stop():
            break
        if frame is None or frame.empty or "mag" not in frame.columns:
            continue
        identity_column = "source_id" if "source_id" in frame.columns else (
            "ID" if "ID" in frame.columns else None
        )
        if identity_column is None:
            continue
        frame = frame.copy()
        frame["star_id"] = pd.to_numeric(frame[identity_column], errors="coerce").astype("Int64")
        frame["mag"] = pd.to_numeric(frame["mag"], errors="coerce")
        if "mag_err" in frame.columns:
            frame["mag_err"] = pd.to_numeric(frame["mag_err"], errors="coerce")
        else:
            frame["mag_err"] = np.nan
        usable = frame["star_id"].notna() & np.isfinite(frame["mag"])
        usable &= ~_bad_measurement_mask(frame)
        frame = frame[usable].copy()
        if frame.empty:
            continue
        # Windows maps NumPy's platform `int` to int32. Gaia DR3 source IDs
        # require all 64 bits, otherwise distinct stars collapse to low-bit IDs.
        frame["star_id"] = frame["star_id"].astype("int64")
        frame = frame.sort_values("mag_err", na_position="last").drop_duplicates(
            subset="star_id", keep="first"
        )

        header_row = header_rows.get(filename, pd.Series(dtype=object))
        time_value = float("nan")
        for column in ("BJD_TDB", "JD", "jd", "DATE-OBS", "date_obs"):
            if column in index_row.index:
                time_value = _time_value(index_row.get(column))
            if not np.isfinite(time_value) and column in header_row.index:
                time_value = _time_value(header_row.get(column))
            if np.isfinite(time_value):
                break
        if not np.isfinite(time_value):
            time_value = float(frame_order)

        night_id = int(night_assignments.get(filename, 0))
        if night_id <= 0:
            try:
                raw_night = int(float(index_row.get("night_id", 0)))
                night_id = raw_night if raw_night > 0 else 0
            except (TypeError, ValueError):
                night_id = 0
        if night_id <= 0:
            date_value = ""
            for column in ("DATE-OBS", "date_obs"):
                value = index_row.get(column, header_row.get(column, ""))
                if str(value or "").strip():
                    date_value = str(value).split("T", 1)[0]
                    break
            if date_value:
                if date_value not in fallback_nights:
                    fallback_nights[date_value] = len(fallback_nights) + 1
                night_id = fallback_nights[date_value]

        airmass = _first_finite(
            index_row.get("airmass"),
            index_row.get("AIRMASS"),
            header_row.get("airmass"),
            header_row.get("AIRMASS"),
        )
        fwhm = _first_finite(
            index_row.get("fwhm_px"),
            index_row.get("FWHM"),
            header_row.get("fwhm_px"),
            header_row.get("FWHM"),
        )
        frame["frame"] = filename
        frame["time"] = time_value
        frame["night_id"] = night_id
        frame["airmass"] = airmass
        frame["fwhm"] = fwhm
        if "sky" in frame.columns:
            frame["sky"] = pd.to_numeric(frame["sky"], errors="coerce")
        elif "bkg_median_adu" in frame.columns:
            frame["sky"] = pd.to_numeric(frame["bkg_median_adu"], errors="coerce")
        else:
            frame["sky"] = np.nan
        for column in ("x", "y"):
            if column in frame.columns:
                frame[column] = pd.to_numeric(frame[column], errors="coerce")
            else:
                frame[column] = np.nan
        frame["photometry_source"] = source.get("source", "aperture")
        frame["mag_input_column"] = source.get("mag_column", "mag")
        frame["mag_error_input_column"] = source.get("mag_error_column", "mag_err")
        output_frames.append(
            frame[
                [
                    "frame",
                    "time",
                    "night_id",
                    "star_id",
                    "mag",
                    "mag_err",
                    "airmass",
                    "fwhm",
                    "sky",
                    "x",
                    "y",
                    "photometry_source",
                    "mag_input_column",
                    "mag_error_input_column",
                ]
            ]
        )
    if not output_frames:
        return pd.DataFrame(), source
    output = pd.concat(output_frames, ignore_index=True)
    output = output.sort_values(["time", "frame", "star_id"], na_position="last")
    return output.reset_index(drop=True), source
