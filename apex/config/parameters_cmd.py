"""
Parameter management for aperture photometry pipeline
APEX CMD parameter handling.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, Iterable
import hashlib
import types
try:  # Python 3.11+
    import tomllib  # type: ignore
except Exception:  # Python 3.10 and earlier
    import tomli as tomllib  # type: ignore

from apex.config.calibration_section import read_calibration_section
from apex.config.parameter_map import (
    build_settings,
    CANONICAL_SCHEMA_VERSION,
    CMD_TOML_KEY_MAP,
    ensure_schema_version,
    read_schema_version,
    toml_value_for_runtime_attr,
)


def _as_bool(v, default=False):
    """Convert value to boolean"""
    if v is None:
        return default
    return str(v).strip().lower() in ("1", "true", "yes", "y", "on")


def _as_float_or_none(v):
    """Convert value to float or None"""
    try:
        s = str(v).strip()
        return float(s) if s != "" else None
    except:
        return None


def _get_path(data: Dict[str, Any], path: Iterable[str]):
    cur: Any = data
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def _set_path(data: Dict[str, Any], path: Iterable[str], value: Any) -> None:
    """Set a value at a nested path in a dictionary, creating intermediate dicts as needed."""
    cur: Any = data
    keys = list(path)
    for key in keys[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]
    cur[keys[-1]] = value


def _delete_path(data: Dict[str, Any], path: Iterable[str]) -> None:
    cur: Any = data
    keys = list(path)
    for key in keys[:-1]:
        if not isinstance(cur, dict) or key not in cur:
            return
        cur = cur[key]
    if isinstance(cur, dict):
        cur.pop(keys[-1], None)


# The live map is `CMD_TOML_KEY_MAP` in `apex.config.parameter_map`. A second
# copy used to sit here and was overwritten on the next line, so editing it did
# nothing; by the time it was removed it was missing 81 of the map's rows
# (2026-08-16).
TOML_KEY_MAP = CMD_TOML_KEY_MAP


def _read_toml(path: Path) -> Dict[str, Any]:
    """Read TOML config into a flat key dict."""
    if not path.exists():
        return {}

    from apex.config.config_io import check_workspace_identity, load_config_data
    data, _cfg_path = load_config_data(path)  # JSON authority; legacy TOML auto-migrated
    for _issue in check_workspace_identity(data, _cfg_path):
        import warnings as _warnings
        _warnings.warn(f"[workspace] {_issue}")

    raw: Dict[str, Any] = {}
    raw["schema_version"] = read_schema_version(data)

    def set_if(key: str, value: Any):
        if value is not None:
            raw[key] = value

    # Use the global TOML_KEY_MAP for mapping
    for row in TOML_KEY_MAP:
        value = _get_path(data, row[0])
        set_if(row[1], value)

    peak_kernel_scales = _get_path(data, ("detection", "peak", "kernel_scales"))
    if isinstance(peak_kernel_scales, list):
        raw["peak_kernel_scales"] = ",".join(str(v) for v in peak_kernel_scales)

    target_ra = _get_path(data, ("target", "ra_deg"))
    target_dec = _get_path(data, ("target", "dec_deg"))
    target_name = _get_path(data, ("target", "name"))
    set_if("target_name", target_name)
    set_if("target_ra_deg", target_ra)
    set_if("target_dec_deg", target_dec)
    set_if("ra_deg", target_ra)
    set_if("dec_deg", target_dec)

    inst = _get_path(data, ("instrument",)) or {}
    if isinstance(inst, dict):
        set_if("telescope_focal_mm", inst.get("telescope_focal_mm"))
        set_if("camera_pixel_um", inst.get("camera_pixel_um"))
        set_if("camera_binning", inst.get("binning"))
        set_if("binning_default", inst.get("binning"))
        set_if("gain_e_per_adu", inst.get("gain_e_per_adu"))
        set_if("rdnoise_e", inst.get("rdnoise_e"))
        set_if("noise_use_fits_header", inst.get("noise_use_fits_header"))
        set_if("noise_reference_binning", inst.get("noise_reference_binning"))
        set_if("noise_scale_by_binning", inst.get("noise_scale_by_binning"))
        set_if("saturation_adu", inst.get("saturation_adu"))
        set_if("zp_initial", inst.get("zp_initial"))
        set_if("datamin_adu", inst.get("datamin_adu"))
        set_if("datamax_adu", inst.get("datamax_adu"))

    hud5x = _get_path(data, ("hud5x",)) or {}
    if isinstance(hud5x, dict):
        set_if("5x.aperture_scale", hud5x.get("aperture_scale"))
        set_if("5x.center_cbox_scale", hud5x.get("center_cbox_scale"))
        set_if("5x.annulus_in_scale", hud5x.get("annulus_in_scale"))
        set_if("5x.annulus_out_scale", hud5x.get("annulus_out_scale"))
        set_if("5x.sigma_clip", hud5x.get("sigma_clip"))
        set_if("5x.neighbor_mask_scale", hud5x.get("neighbor_mask_scale"))
        set_if("5x.mag_flux", str(hud5x.get("mag_flux")) if hud5x.get("mag_flux") is not None else None)
        set_if("5x.use_header_exptime", hud5x.get("use_header_exptime"))
        set_if("5x.min_r_ap_px", hud5x.get("min_r_ap_px"))
        set_if("5x.min_r_in_px", hud5x.get("min_r_in_px"))
        set_if("5x.min_r_out_px", hud5x.get("min_r_out_px"))

    site = _get_path(data, ("site",)) or {}
    if isinstance(site, dict):
        set_if("site_lat_deg", site.get("lat_deg"))
        set_if("site_lon_deg", site.get("lon_deg"))
        set_if("site_alt_m", site.get("alt_m"))
        set_if("site_tz_offset_hours", site.get("tz_offset_hours"))

    # [calibration] (+ [calibration.overscan]) — kept as a nested dict because
    # its keys map 1:1 onto CalibrationOptions, which owns their defaults.
    calibration = read_calibration_section(data)
    if calibration:
        raw["_calibration"] = calibration

    extra = _get_path(data, ("parameters",)) or {}
    if isinstance(extra, dict):
        for key, value in extra.items():
            raw.setdefault(key, value)

    return raw


def _getf(raw, key, default):
    """Get float value from raw dict"""
    s = str(raw.get(key, "")).strip()
    try:
        return default if s == "" else float(s)
    except:
        return default


def _geti(raw, key, default):
    """Get int value from raw dict"""
    s = str(raw.get(key, "")).strip()
    try:
        if s == "":
            return default
        return int(float(s))
    except:
        return default


class Parameters:
    """
    Main parameter container for the photometry pipeline
    All configuration is stored as a SimpleNamespace for easy attribute access
    """

    def __init__(self, param_file: str | Path = "parameters.toml"):
        from apex.config.config_io import migrate_config_path
        param_file = migrate_config_path(param_file)
        """Initialize parameters from file"""
        self.param_file = Path(param_file)
        self.P = self._load_from_file(self.param_file)
        self.param_hash = self._compute_hash(self.param_file)

    @staticmethod
    def _load_from_file(path: Path) -> types.SimpleNamespace:
        """Load parameters from text file"""
        raw = _read_toml(path)
        parallel_workers = _geti(raw, "max_workers", _geti(raw, "parallel_max_workers", 0))

        rdnoise_candidate = (
            _as_float_or_none(raw.get("rdnoise_e", ""))
            or _as_float_or_none(raw.get("datapar.readnoise", ""))
            or _as_float_or_none(raw.get("readnoise_e", ""))
        )

        # Every setting the key map defines in full — path, name, type,
        # default — is built straight from those rows. What follows is only
        # what a row cannot say: values that need normalising, that differ
        # between CMD and LC, that are computed, or that alias another key.
        values = build_settings(raw, TOML_KEY_MAP)
        values.update(
            schema_version=_geti(raw, "schema_version", CANONICAL_SCHEMA_VERSION),

            # Parallel processing
            max_workers=parallel_workers,  # 0 = auto
            # Backward-compatible alias for older code paths.
            parallel_max_workers=parallel_workers,

            # Detection parameters
            detect_engine=raw.get("detect_engine", "segm"),
            detect_mode=str(raw.get("detect_mode", "normal")).strip().lower() or "normal",
            detect_sigma_by_filter=raw.get("detect_sigma_by_filter", {}) or {},
            detect_keep_max=_geti(raw, "detect_keep_max", 6000),
            deblend_mode=str(raw.get("deblend_mode", "exponential")).strip().lower() or "exponential",
            dao_refine_enable=_as_bool(raw.get("dao_refine_enable", "false"), False),

            # Background parameters
            bkg2d_enable=_as_bool(raw.get("bkg2d_enable", "true"), True),
            bkg2d_downsample=_geti(raw, "bkg2d_downsample", 4),

            # Camera/instrument
            rdnoise_e=rdnoise_candidate,
            # InstrumentConfig and overwrote the right value with the default one.
            camera_binning=_geti(raw, "camera_binning", 0) or None,

            # Detector calibration (Step 0); CalibrationOptions owns the defaults
            calibration=dict(raw.get("_calibration") or {}),

            # 5X HUD viewer parameters
            _hud5={
                "5x.aperture_scale": raw.get("5x.aperture_scale", ""),
                "5x.center_cbox_scale": raw.get("5x.center_cbox_scale", ""),
                "5x.annulus_in_scale": raw.get("5x.annulus_in_scale", ""),
                "5x.annulus_out_scale": raw.get("5x.annulus_out_scale", ""),
                "5x.min_r_ap_px": raw.get("5x.min_r_ap_px", ""),
                "5x.min_r_in_px": raw.get("5x.min_r_in_px", ""),
                "5x.min_r_out_px": raw.get("5x.min_r_out_px", ""),
                "5x.sigma_clip": raw.get("5x.sigma_clip", ""),
                "5x.neighbor_mask_scale": raw.get("5x.neighbor_mask_scale", ""),
                "5x.mag_flux": raw.get("5x.mag_flux", "rate_e"),
                "5x.use_header_exptime": raw.get("5x.use_header_exptime", "true"),
            },

            # Photometry execution flags
            bkg_use_segm_mask=_as_bool(raw.get("bkg_use_segm_mask", "true"), True),

            # ID matching
            idmatch_mode=str(raw.get("idmatch_mode", "normal")).strip().lower() or "normal",

            # Aperture/annulus config
            apcorr_scale_min=_getf(raw, "apcorr_scale_min", 0.5),
            apcorr_scale_max=_getf(raw, "apcorr_scale_max", 5.0),
            apcorr_scale_step=_getf(raw, "apcorr_scale_step", 0.25),
            # loaded into `raw` and was dropped on the way into the namespace.
            apcorr_large_scale=_getf(raw, "apcorr_large_scale", _getf(raw, "apcorr_large_ref_scale", 2.4)),
            apcorr_large_ref_scale=_getf(raw, "apcorr_large_scale", _getf(raw, "apcorr_large_ref_scale", 2.4)),

            # 없으면 None → resolve_wcs_engine 이 내장 엔진을 기본으로 선택.
            wcs_engine=(str(raw.get("wcs_engine", "")).strip().lower() or None),
            wcs_qc_match_radius_arcsec=_getf(raw, "wcs_qc_match_radius_arcsec", 2.0),
            wcs_qc_clip_sigma=_getf(raw, "wcs_qc_clip_sigma", 3.0),
            wcs_qc_require_wcs_ok=_as_bool(raw.get("wcs_qc_require_wcs_ok", "true"), True),
            wcs_qc_min_match_n=_geti(raw, "wcs_qc_min_match_n", 20),
            wcs_qc_min_match_rate=_getf(raw, "wcs_qc_min_match_rate", 0.20),
            wcs_qc_max_rms_px=_getf(raw, "wcs_qc_max_rms_px", 3.5),
            wcs_qc_max_p99_px=_getf(raw, "wcs_qc_max_p99_px", 5.0),
            wcs_qc_min_inlier_rate=_getf(raw, "wcs_qc_min_inlier_rate", 0.50),
            wcs_qc_max_edge_ratio=_getf(raw, "wcs_qc_max_edge_ratio", 0.0),
            wcs_qc_max_center_offset_arcsec=_getf(raw, "wcs_qc_max_center_offset_arcsec", 0.0),
            gaia_pmem_method=str(raw.get("gaia_pmem_method", "gmm3d")).strip().lower() or "gmm3d",
            idmatch_use_qc_pass_only=_as_bool(raw.get("idmatch_use_qc_pass_only", "true"), True),
            idmatch_use_wcs_qc_gate=_as_bool(raw.get("idmatch_use_wcs_qc_gate", "true"), True),
            idmatch_wcs_qc_min_match_rate=_getf(raw, "idmatch_wcs_qc_min_match_rate", 0.20),
            idmatch_wcs_qc_min_match_n=_geti(raw, "idmatch_wcs_qc_min_match_n", 20),
            idmatch_wcs_qc_max_rms_px=_getf(raw, "idmatch_wcs_qc_max_rms_px", 3.5),
            idmatch_wcs_qc_min_inlier_rate=_getf(raw, "idmatch_wcs_qc_min_inlier_rate", 0.50),
            idmatch_wcs_qc_max_p99_px=_getf(raw, "idmatch_wcs_qc_max_p99_px", 5.0),

            # Essential when a band's Gaia reference is 'approx' (Johnson U).
            std_anchor_enable=_as_bool(raw.get("std_anchor_enable", "false"), False),
            std_anchor_catalog=str(raw.get("std_anchor_catalog", "") or ""),
            std_anchor_match_radius=_getf(raw, "std_anchor_match_radius", 1.5),
            std_anchor_min_stars=_geti(raw, "std_anchor_min_stars", 20),
            # enabling it changes every globular-cluster zero point.
            gaia_cstar_cut=bool(raw.get("gaia_cstar_cut", False)),

            # PSF photometry (CMD Step 8)
            psf_mode=str(raw.get("psf_mode", "normal")).strip().lower() or "normal",
            psf_model_mode="per_frame",
            psf_interp_order=_geti(raw, "psf_interp_order", 1),
            psf_epsf_contamination_filter=_as_bool(
                raw.get("psf_epsf_contamination_filter", "true"),
                True,
            ),
            psf_flux_scale_correction=_as_bool(
                raw.get("psf_flux_scale_correction", "false"),
                False,
            ),
            psf_flux_scale_min_neighbor_fwhm=_getf(
                raw, "psf_flux_scale_min_neighbor_fwhm", 4.0
            ),
            psf_flux_scale_max_scatter_mag=_getf(
                raw, "psf_flux_scale_max_scatter_mag", 0.10
            ),
            psf_fit_window_mode=str(
                raw.get("psf_fit_window_mode", "auto")
            ).strip().lower() or "auto",
            psf_fit_encircled_energy=_getf(
                raw, "psf_fit_encircled_energy", 0.90
            ),
            psf_fit_mode=str(raw.get("psf_fit_mode", "new")).strip().lower() or "new",
            psf_core_cut_center_mode=str(raw.get("psf_core_cut_center_mode", "auto")).strip().lower() or "auto",
            psf_grouper_budget_frac=_getf(raw, "psf_grouper_budget_frac", 0.10),
            psf_grouper_budget_cap=_geti(raw, "psf_grouper_budget_cap", 200),
            psf_final_pass_max_iter=_geti(raw, "psf_final_pass_max_iter", 2),
            psf_forced_position_lock=str(raw.get("psf_forced_position_lock", "always")),
            psf_profile_error_frac=_getf(raw, "psf_profile_error_frac", 0.0),
            psf_fit_engine=str(raw.get("psf_fit_engine", "apex_iterative")).strip().lower() or "apex_iterative",
            psf_build_mode=str(raw.get("psf_build_mode", "epsf")).strip().lower() or "epsf",

            # Cross-frame pixel matching (Steps 7–8)
            cross_frame_ransac_tol_px=_getf(raw, "cross_frame_ransac_tol_px", 2.0),
            cross_frame_ransac_max_iter=_geti(raw, "cross_frame_ransac_max_iter", 1000),
            cross_frame_ransac_min_inliers=_geti(raw, "cross_frame_ransac_min_inliers", 6),
            cross_frame_match_tol_px=_getf(raw, "cross_frame_match_tol_px", 3.0),

            # Step 2 crop rectangle (config-driven, headless). Default: skip.
            crop_enable=_as_bool(raw.get("crop_enable", "false"), False),
            crop_x0=_as_float_or_none(raw.get("crop_x0", "")),
            crop_y0=_as_float_or_none(raw.get("crop_y0", "")),
            crop_x1=_as_float_or_none(raw.get("crop_x1", "")),
            crop_y1=_as_float_or_none(raw.get("crop_y1", "")),
        )
        P = types.SimpleNamespace(**values)

        # Store raw dict for compatibility
        P._raw = raw

        # Setup directory paths
        P.data_dir = Path(P.data_dir)
        P.result_dir = Path(P.result_dir) if P.result_dir else (P.data_dir / "result")
        P.result_dir.mkdir(parents=True, exist_ok=True)
        P.cache_dir = (P.result_dir / str(P.cache_dir))
        P.cache_dir.mkdir(parents=True, exist_ok=True)

        pix_guess = P.fwhm_pix_guess
        P.fwhm_seed_px = float(pix_guess if pix_guess is not None else 6.0)
        P._fwhm_seed_from = "pixel"

        return P

    @staticmethod
    def _compute_hash(path: Path) -> str:
        """Compute SHA1 hash of parameter file content (excluding comments)"""
        try:
            txt = Path(path).read_text(encoding="utf-8", errors="ignore")
            lines = []
            for ln in txt.splitlines():
                s = ln.strip()
                if (not s) or s.startswith("#"):
                    continue
                if "#" in s:
                    s = s.split("#", 1)[0].strip()
                lines.append(s)
            norm = "\n".join(lines).encode("utf-8")
            return hashlib.sha1(norm).hexdigest()
        except Exception:
            return "NO_PARAM"

    def save(self, path: Path):
        """Save current parameters to file"""
        path = Path(path)
        data = {"parameters": {}}

        for key, value in vars(self.P).items():
            if key.startswith("_"):
                continue
            if isinstance(value, Path):
                value = str(value)
            data["parameters"][key] = value

        from apex.config.config_io import resolve_config_path, save_config_data
        out = resolve_config_path(path)
        if not save_config_data(out, data):
            raise RuntimeError(f"could not write workspace config: {out}")

    def get(self, name: str, default=None):
        """Get parameter value with fallback to raw dict"""
        if hasattr(self.P, name):
            val = getattr(self.P, name)
            if not (val is None or (isinstance(val, str) and val.strip() == "")):
                return val
        if hasattr(self.P, "_raw") and name in self.P._raw:
            rawv = self.P._raw[name]
            if isinstance(default, bool):
                return _as_bool(rawv, default)
            if isinstance(default, int):
                fv = _as_float_or_none(rawv)
                return int(fv) if fv is not None else default
            if isinstance(default, float):
                fv = _as_float_or_none(rawv)
                return fv if fv is not None else default
            return rawv
        return default

    def get_file_path(self, filename: str) -> Path:
        """Resolve a filename to the original FITS path (multi-night safe)."""
        path_map = getattr(self.P, "file_path_map", None)
        if isinstance(path_map, dict):
            mapped = path_map.get(filename)
            if mapped:
                return Path(mapped)
        return Path(self.P.data_dir) / filename

    def save_toml(self, path: Path | str | None = None) -> bool:
        """Persist current parameters to the workspace JSON.

        The name is historical — every persisted config is JSON now; a legacy
        ``.toml`` argument is resolved to its JSON authority first.
        """
        from apex.config.config_io import load_config_data, save_config_data
        try:
            data, param_path = load_config_data(
                path or getattr(self, "param_file", "parameters.toml"))
        except Exception:
            param_path = None
            data = {}
        if param_path is None:
            from apex.config.config_io import resolve_config_path
            param_path = resolve_config_path(
                path or getattr(self, "param_file", "parameters.toml"))
        ensure_schema_version(data)

        for row in TOML_KEY_MAP:
            path_keys, attr = row[0], row[1]
            if not hasattr(self.P, attr):
                continue
            val = getattr(self.P, attr)
            if val is None:
                if attr in {
                    "gain_e_per_adu",
                    "rdnoise_e",
                    "target_ra_deg",
                    "target_dec_deg",
                }:
                    _delete_path(data, path_keys)
                continue
            _set_path(data, path_keys, toml_value_for_runtime_attr(attr, val, path_keys))

        # Keep wcs_refine.enable in sync with the flat runtime flag if present.
        if hasattr(self.P, "wcs_refine_enable"):
            _set_path(data, ("wcs_refine", "enable"), bool(getattr(self.P, "wcs_refine_enable")))

        if not save_config_data(param_path, data):
            return False
        self.param_file = Path(param_path)
        self.param_hash = self._compute_hash(param_path)
        return True

    def print_summary(self):
        """Print parameter summary"""
        P = self.P
        print("\n==================== PARAM SUMMARY ====================")
        print(f"DATA_DIR      : {P.data_dir}")
        print(f"RESULT_DIR    : {P.result_dir}")
        print(f"CACHE_DIR     : {P.cache_dir}")
        print(
            f"resume_mode   : {P.resume_mode} | force_redetect={P.force_redetect} | "
            f"force_rephot={P.force_rephot}"
        )
        print(
            f"idmatch_mode  : {getattr(P, 'idmatch_mode', 'normal')} | "
            f"gaia_only={getattr(P, 'idmatch_use_gaia_refs_only', False)}"
        )
        print(f"parallel_mode : {P.parallel_mode} | max_workers={P.max_workers}")
        print(f"FWHM seed     : {P.fwhm_seed_px:.2f} px (from={getattr(P, '_fwhm_seed_from', '?')})")
        print(f"FWHM range    : {P.fwhm_px_min:.2f} ~ {P.fwhm_px_max:.2f} px | elong_max={P.fwhm_elong_max} | iso_min_sep={P.iso_min_sep_pix}px")
        print(f"bkg2d detect  : {P.bkg2d_in_detect} | box={P.bkg2d_box}")
        print(f"detect_mode   : {getattr(P, 'detect_mode', 'normal')}")
        print(f"detect_sigma  : base={P.detect_sigma} | by_filter={P.detect_sigma_by_filter}")
        print(f"deblend       : enable={P.deblend_enable} nthresh={P.deblend_nthresh} cont={P.deblend_cont} dilate={P.segm_dilate_radius_px}")
        print(f"clip          : sat_adu={P.saturation_adu}")
        print(f"camera        : gain={P.gain_e_per_adu} e-/ADU | rdnoise={P.rdnoise_e} e- | zp_init={P.zp_initial}")
        print("=======================================================\n")


def read_params(path: str | Path = "parameters.toml") -> Parameters:
    """
    Load parameters from file

    Args:
        path: Path to parameter file

    Returns:
        Parameters object
    """
    return Parameters(path)
