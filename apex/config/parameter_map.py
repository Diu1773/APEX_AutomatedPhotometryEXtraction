"""
Shared parameter schema metadata.

This module is the foundation for moving CMD/LC parameter handling toward a
single canonical TOML schema. The current runtime still exposes flat
``params.P.<name>`` attributes for compatibility; later phases can move the
large per-mode maps here without changing the public TOML shape again.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

# Row defaults that are named constants rather than literals live with the rest
# of the measured constants; `apex.utils` never imports `apex.config`, so this
# stays inside the allowed direction.
from apex.utils.constants import DEPTH_QC_TOLERANCE_MAG


CANONICAL_SCHEMA_VERSION = 1
LEGACY_SCHEMA_VERSION = 0


@dataclass(frozen=True)
class ParameterMapEntry:
    """Map one canonical TOML path to one legacy runtime attribute."""

    toml_path: tuple[str, ...]
    attr: str
    modes: frozenset[str] = frozenset({"cmd", "lc"})
    legacy: bool = False
    note: str = ""

    @property
    def dotted_path(self) -> str:
        return ".".join(self.toml_path)


COMMON_FOUNDATION_MAP: tuple[ParameterMapEntry, ...] = (
    ParameterMapEntry(("io", "data_dir"), "data_dir"),
    ParameterMapEntry(("io", "filename_prefix"), "filename_prefix"),
    ParameterMapEntry(("io", "result_dir"), "result_dir"),
    ParameterMapEntry(("io", "cache_dir"), "cache_dir"),
    ParameterMapEntry(("parallel", "mode"), "parallel_mode"),
    ParameterMapEntry(("parallel", "max_workers"), "max_workers"),
    ParameterMapEntry(("ui", "log_tail"), "ui_log_tail"),
    ParameterMapEntry(("fwhm", "guess_arcsec"), "fwhm_guess_arcsec"),
    ParameterMapEntry(("fwhm", "px_min"), "fwhm_px_min"),
    ParameterMapEntry(("fwhm", "px_max"), "fwhm_px_max"),
    ParameterMapEntry(("detection", "engine"), "detect_engine"),
    ParameterMapEntry(("detection", "sigma"), "detect_sigma"),
    ParameterMapEntry(("detection", "minarea_pix"), "minarea_pix"),
    ParameterMapEntry(("background", "in_detect"), "bkg2d_in_detect"),
    ParameterMapEntry(("photometry", "mode"), "aperture_mode"),
    ParameterMapEntry(("instrument", "rdnoise_e"), "rdnoise_e"),
)


LEGACY_ALIAS_MAP: tuple[ParameterMapEntry, ...] = (
    ParameterMapEntry(
        ("fwhm", "iso_min_sep_pix"),
        "iso_min_sep_pix",
        legacy=True,
        note="Use *_px naming in new canonical keys when this is migrated.",
    ),
    ParameterMapEntry(
        ("master", "iso_min_sep_pix"),
        "master_iso_min_sep_pix",
        legacy=True,
        note="Use *_px naming in new canonical keys when this is migrated.",
    ),
)


# Full legacy runtime maps. Keep order stable because duplicate runtime
# attributes can intentionally write/read multiple TOML paths.
#
# A row is either
#
#     (dotted path, attribute)                      -- mapped, not built here
#     (dotted path, attribute, kind, default)       -- mapped and built here
#
# The four-part row is the whole definition of a setting: where it comes from
# in the file, what the runtime calls it, what type it is, and what it falls
# back to. `build_settings` turns those rows straight into the namespace, so a
# new setting is one row and cannot be half-added.
#
# Two-part rows are settings the map knows and the namespace does not build.
# Each is one of three things, all deliberate: handled explicitly by the loader
# because it needs normalising or differs per mode, read out of `P._raw`, or a
# setting that genuinely does not reach the code yet. The third group is listed
# in `docs/audit/CONFIG_REACHABILITY.md` and pinned by
# `tests/test_config_settings_reach_code.py`. Until 2026-08-16 that third group
# was invisible: the key map and a 452-line hand-written `SimpleNamespace(...)`
# call held the same list twice, and 74 settings were only in one of them, so
# they loaded and then vanished with nothing in the log.


def coerce_setting(kind: str, value: Any, default: Any) -> Any:
    """Turn one raw config value into what the runtime expects.

    Mirrors the `_getf` / `_geti` / `_as_bool` helpers the loaders used inline,
    including their habit of accepting strings: legacy TOML wrote numbers and
    booleans as text, and workspaces migrated from it still carry them.
    """
    if kind == "str":
        return default if value is None else value
    if kind == "bool":
        # `_as_bool`: anything unrecognised is False, not the default. Only an
        # absent key falls back, and the loader's `set_if` never stores None.
        if value is None:
            return default
        return str(value).strip().lower() in ("1", "true", "yes", "y", "on")

    # `_getf` / `_geti` go through `str()` first, so a JSON `true` in a numeric
    # field is the string "True" and fails the conversion — it must not become
    # 1.0 via `float(True)`, which is what a bare isinstance check would do.
    text = "" if value is None else str(value).strip()
    if kind == "float_or_none":
        try:
            return float(text) if text != "" else None
        except (TypeError, ValueError):
            return None
    if text == "":
        return default
    try:
        return int(float(text)) if kind == "int" else float(text)
    except (TypeError, ValueError):
        return default


def build_settings(raw: dict, key_map: Iterable[tuple]) -> dict[str, Any]:
    """Every setting the map defines fully, as attribute -> value.

    Rows without a kind are skipped; the loader sets those itself. Later rows
    win, matching the order-dependent behaviour of the map (an attribute may
    legitimately be reachable from more than one path).
    """
    values: dict[str, Any] = {}
    for row in key_map:
        if len(row) < 4:
            continue
        attr, kind, default = row[1], row[2], row[3]
        values[attr] = coerce_setting(kind, raw.get(attr), default)
    return values

COMMON_TOML_KEY_MAP: tuple[tuple, ...] = (
    (('io', 'data_dir'), 'data_dir', 'str', "."),
    (('io', 'filename_prefix'), 'filename_prefix', 'str', ""),
    (('io', 'result_dir'), 'result_dir', 'str', ""),
    (('io', 'cache_dir'), 'cache_dir', 'str', "cache"),
    (('io', 'night_gap_hours'), 'night_gap_hours', 'float', 8.0),
    (('parallel', 'mode'), 'parallel_mode', 'str', "thread"),
    (('parallel', 'max_workers'), 'max_workers'),
    (('parallel', 'resume_mode'), 'resume_mode', 'bool', True),
    (('parallel', 'force_redetect'), 'force_redetect', 'bool', False),
    (('parallel', 'force_rephot'), 'force_rephot', 'bool', False),
    (('parallel', 'detect_cache_strategy'), 'detect_cache_strategy', 'str', "mtime"),
    (('ui', 'log_tail'), 'ui_log_tail', 'int', 300),
    (('ui', 'detect_progress_bar'), 'detect_progress_bar', 'bool', True),
    (('alignment', 'ref_index'), 'align_ref_index', 'int', 0),
    (('alignment', 'global_align'), 'global_align', 'bool', True),
    (('alignment', 'global_ref_filter'), 'global_ref_filter', 'str', "r"),
    (('alignment', 'global_ref_index'), 'global_ref_index', 'int', 0),
    (('fwhm', 'guess_px'), 'fwhm_pix_guess', 'float_or_none', None),
    (('fwhm', 'guess_arcsec'), 'fwhm_guess_arcsec', 'float_or_none', None),
    (('fwhm', 'arcsec_min'), 'fwhm_arcsec_min', 'float_or_none', None),
    (('fwhm', 'arcsec_max'), 'fwhm_arcsec_max', 'float_or_none', None),
    (('fwhm', 'px_min'), 'fwhm_px_min', 'float', 3.5),
    (('fwhm', 'px_max'), 'fwhm_px_max', 'float', 12.0),
    (('fwhm', 'qc_max_sources'), 'fwhm_qc_max_sources', 'int', 40),
    (('fwhm', 'elong_max'), 'fwhm_elong_max', 'float', 1.3),
    (('fwhm', 'iso_min_sep_pix'), 'iso_min_sep_pix', 'float', 18.0),
    (('fwhm', 'measure_max'), 'fwhm_measure_max', 'int', 25),
    (('fwhm', 'min_sources'), 'fwhm_min_sources', 'int', 15),
    (('fwhm', 'candidate_max'), 'fwhm_candidate_max', 'int', 200),
    (('fwhm', 'measure_all_sources'), 'fwhm_measure_all_sources', 'bool', False),
    (('fwhm', 'dr'), 'fwhm_dr', 'float', 0.5),
    (('clip', 'min_adu'), 'clip_min_adu'),
    (('clip', 'max_adu'), 'clip_max_adu'),
    (('detection', 'engine'), 'detect_engine'),
    (('detection', 'sigma'), 'detect_sigma', 'float', 3.2),
    (('detection', 'sigma_g'), 'detect_sigma_g', 'float_or_none', None),
    (('detection', 'sigma_r'), 'detect_sigma_r', 'float_or_none', None),
    (('detection', 'sigma_i'), 'detect_sigma_i', 'float_or_none', None),
    (('detection', 'sigma_by_filter'), 'detect_sigma_by_filter'),
    (('detection', 'minarea_pix'), 'minarea_pix', 'int', 3),
    (('detection', 'keep_max'), 'detect_keep_max', 'int', 6000),
    (('detection', 'dilate_radius_px'), 'segm_dilate_radius_px', 'int', 4),
    (('detection', 'deblend', 'enable'), 'deblend_enable', 'bool', True),
    (('detection', 'deblend', 'nthresh'), 'deblend_nthresh', 'int', 64),
    (('detection', 'deblend', 'contrast'), 'deblend_cont', 'float', 0.0025),
    (('detection', 'deblend', 'max_labels'), 'deblend_max_labels', 'int', 4000),
    (('detection', 'deblend', 'label_hard_max'), 'deblend_label_hard_max', 'int', 7000),
    (('detection', 'deblend', 'nlevels_soft'), 'deblend_nlevels_soft', 'int', 32),
    (('detection', 'deblend', 'contrast_soft'), 'deblend_contrast_soft', 'float', 0.005),
    (('detection', 'peak', 'enable'), 'peak_pass_enable', 'bool', True),
    (('detection', 'peak', 'nsigma'), 'peak_nsigma', 'float', 3.2),
    (('detection', 'peak', 'min_sep_px'), 'peak_min_sep_px', 'float', 4.0),
    (('detection', 'peak', 'max_add'), 'peak_max_add', 'int', 600),
    (('detection', 'peak', 'max_elong'), 'peak_max_elong', 'float', 1.6),
    (('detection', 'peak', 'sharp_lo'), 'peak_sharp_lo', 'float', 0.12),
    (('detection', 'peak', 'skip_if_nsrc_ge'), 'peak_skip_if_nsrc_ge', 'int', 4500),
    (('detection', 'dao', 'enable'), 'dao_refine_enable'),
    (('detection', 'dao', 'fwhm_px'), 'dao_fwhm_px', 'float', 6.0),
    (('detection', 'dao', 'sharp_lo'), 'dao_sharp_lo', 'float', 0.2),
    (('detection', 'dao', 'sharp_hi'), 'dao_sharp_hi', 'float', 1.0),
    (('detection', 'dao', 'round_lo'), 'dao_round_lo', 'float', -0.5),
    (('detection', 'dao', 'round_hi'), 'dao_round_hi', 'float', 0.5),
    (('detection', 'dao', 'match_tol_px'), 'dao_match_tol_px', 'float', 2.0),
    (('background', 'enable'), 'bkg2d_enable'),
    (('background', 'in_detect'), 'bkg2d_in_detect', 'bool', True),
    (('background', 'box'), 'bkg2d_box', 'int', 64),
    (('background', 'filter_size'), 'bkg2d_filter_size', 'int', 3),
    (('background', 'edge_method'), 'bkg2d_edge_method', 'str', "pad"),
    (('background', 'method'), 'bkg2d_method', 'str', "median"),
    (('background', 'downsample'), 'bkg2d_downsample', 'int', 4),
    (('qc', 'gate_enable'), 'gate_enable', 'bool', True),
    (('qc', 'sky_sigma_max_e'), 'gate_sky_sigma_max_e'),
    (('qc', 'nsrc_min'), 'gate_nsrc_min'),
    (('qc', 'keep_positions_if_fail'), 'keep_positions_if_qc_fail'),
    (('photometry', 'mode'), 'aperture_mode', 'str', "apcorr"),
    (('photometry', 'recenter'), 'recenter_aperture', 'bool', True),
    (('photometry', 'use_segm_mask'), 'bkg_use_segm_mask'),
    (('photometry', 'min_snr_for_mag'), 'min_snr_for_mag', 'float', 3.0),
    (('photometry', 'max_recenter_shift'), 'max_recenter_shift', 'float', 2.0),
    (('photometry', 'centroid_outlier_px'), 'centroid_outlier_px', 'float', 1.0),
    (('photometry', 'registration', 'match_radius_px'), 'registration_match_radius_px', 'float', 4.0),
    (('photometry', 'registration', 'min_anchors'), 'registration_min_anchors', 'int', 6),
    (('photometry', 'sky_sigma_mode'), 'sky_sigma_mode', 'str', "local"),
    (('photometry', 'sky_sigma_includes_rn'), 'sky_sigma_includes_rn', 'bool', True),
    (('photometry', 'sky_sigma_min_n_sky'), 'sky_sigma_min_n_sky', 'int', 50),
    (('photometry', 'use_qc_pass_only'), 'phot_use_qc_pass_only', 'bool', False),
    (('photometry', 'ref_require_apcorr_candidate'), 'phot_ref_require_apcorr_candidate', 'bool', True),
    (('photometry', 'ref_apcorr_min_keep'), 'phot_ref_apcorr_min_keep', 'int', 8),
    (('photometry', 'scales', 'aperture_scale'), 'phot_aperture_scale', 'float', 1.0, {'label': 'Aperture Scale (× FWHM)', 'lo': 0.5, 'hi': 10.0, 'step': 0.1, 'decimals': 2}),
    (('photometry', 'scales', 'annulus_scale'), 'fitsky_annulus_scale', 'float', 4.0, {'label': 'Annulus Inner Scale (× FWHM)', 'lo': 1.0, 'hi': 20.0, 'step': 0.1, 'decimals': 2}),
    (('photometry', 'scales', 'dannulus_scale'), 'fitsky_dannulus_scale', 'float', 2.0, {'label': 'Annulus Outer Width (× FWHM)', 'lo': 0.5, 'hi': 20.0, 'step': 0.1, 'decimals': 2}),
    (('photometry', 'scales', 'center_cbox_scale'), 'center_cbox_scale', 'float', 1.5),
    (('photometry', 'scales', 'annulus_min_gap_px'), 'annulus_min_gap_px', 'float', 6.0),
    (('photometry', 'scales', 'annulus_min_width_px'), 'annulus_min_width_px', 'float', 12.0),
    (('photometry', 'radii', 'min_r_ap_px'), 'min_r_ap_px', 'float', 4.0),
    (('photometry', 'radii', 'min_r_in_px'), 'min_r_in_px', 'float', 12.0),
    (('photometry', 'radii', 'min_r_out_px'), 'min_r_out_px', 'float', 20.0),
    (('photometry', 'radii', 'sigma_clip'), 'annulus_sigma_clip', 'float', 3.0, {'label': 'Sigma Clipping (σ)', 'lo': 1.0, 'hi': 10.0, 'step': 0.1, 'decimals': 1}),
    (('photometry', 'radii', 'max_iter'), 'fitsky_max_iter', 'int', 5),
    (('photometry', 'radii', 'neighbor_mask_scale'), 'annulus_neighbor_mask_scale', 'float', 1.3),
    (('photometry', 'apcorr', 'apply'), 'apcorr_apply', 'bool', True),
    (('photometry', 'apcorr', 'small_scale'), 'apcorr_small_scale', 'float', 0.8),
    (('photometry', 'apcorr', 'large_scale'), 'apcorr_large_scale'),
    (('photometry', 'apcorr', 'min_n'), 'apcorr_use_min_n', 'int', 20),
    (('photometry', 'apcorr', 'scatter_max'), 'apcorr_scatter_max', 'float', 0.05),
    (('photometry', 'apcorr', 'optimize_scales'), 'apcorr_optimize_scales'),
    (('photometry', 'apcorr', 'small_scale_min'), 'apcorr_small_scale_min'),
    (('photometry', 'apcorr', 'small_scale_max'), 'apcorr_small_scale_max'),
    (('photometry', 'apcorr', 'large_scale_min'), 'apcorr_large_scale_min'),
    (('photometry', 'apcorr', 'large_scale_max'), 'apcorr_large_scale_max'),
    (('photometry', 'apcorr', 'scale_step'), 'apcorr_scale_step'),
    (('photometry', 'apcorr', 'min_gap_fwhm'), 'apcorr_min_gap_fwhm'),
    (('photometry', 'apcorr', 'max_pairs'), 'apcorr_max_pairs'),
    (('photometry', 'apcorr', 'max_sources'), 'apcorr_max_sources', 'int', 250),
    (('photometry', 'apcorr', 'scale_min'), 'apcorr_scale_min'),
    (('photometry', 'apcorr', 'scale_max'), 'apcorr_scale_max'),
    (('photometry', 'apcorr', 'min_snr'), 'apcorr_min_snr', 'float', 40.0),
    (('photometry', 'apcorr', 'isolation_factor'), 'apcorr_isolation_factor', 'float', 2.5),
    (('photometry', 'depth_qc', 'tolerance_mag'), 'depth_qc_tolerance_mag', 'float', DEPTH_QC_TOLERANCE_MAG),
    (('photometry', 'depth_qc', 'min_snr'), 'depth_qc_min_snr', 'float', 40.0),
    (('source_quality', 'fwhm_ratio_lo'), 'source_quality_fwhm_ratio_lo', 'float', 0.6),
    (('source_quality', 'fwhm_ratio_hi'), 'source_quality_fwhm_ratio_hi', 'float', 1.6),
    (('source_quality', 'anchor_neighbor_fwhm_mult'), 'source_quality_anchor_neighbor_fwhm_mult', 'float', 2.0),
    (('source_quality', 'anchor_flux_pct'), 'source_quality_anchor_flux_pct', 'float', 60.0),
    (('source_quality', 'apcorr_flux_pct'), 'source_quality_apcorr_flux_pct', 'float', 60.0),
    (('source_quality', 'psf_seed_flux_pct'), 'source_quality_psf_seed_flux_pct', 'float', 30.0),
    (('source_quality', 'edge_fwhm_mult'), 'source_quality_edge_fwhm_mult', 'float', 1.0),
    (('wcs', 'engine'), 'wcs_engine'),
    (('wcs', 'astap_exe'), 'astap_exe', 'str', "astap_cli.exe"),
    (('wcs', 'timeout_s'), 'astap_timeout_s', 'float', 120.0),
    (('wcs', 'astap_search_radius_deg'), 'astap_search_radius_deg', 'float', 8.0),
    (('wcs', 'astap_database'), 'astap_database', 'str', "D80"),
    (('wcs', 'astap_fov_fudge'), 'astap_fov_fudge', 'float', 1.0),
    (('wcs', 'astap_downsample'), 'astap_downsample_z', 'int', 2),
    (('wcs', 'astap_max_stars'), 'astap_max_stars_s', 'int', 500),
    (('wcs', 'astnet_local_enable'), 'astnet_local_enable', 'bool', False),
    (('wcs', 'astnet_local_use_wsl'), 'astnet_local_use_wsl', 'bool', True),
    (('wcs', 'astnet_local_command'), 'astnet_local_command', 'str', "solve-field"),
    (('wcs', 'astnet_local_timeout_s'), 'astnet_local_timeout_s', 'float', 300.0),
    (('wcs', 'astnet_local_downsample'), 'astnet_local_downsample', 'int', 2),
    (('wcs', 'astnet_local_scale_low'), 'astnet_local_scale_low', 'float', 0.0),
    (('wcs', 'astnet_local_scale_high'), 'astnet_local_scale_high', 'float', 0.0),
    (('wcs', 'astnet_local_radius_deg'), 'astnet_local_radius_deg', 'float', 8.0),
    (('wcs', 'astnet_local_keep_outputs'), 'astnet_local_keep_outputs', 'bool', True),
    (('wcs', 'astnet_local_use_cache'), 'astnet_local_use_cache', 'bool', True),
    (('wcs', 'astnet_local_max_objs'), 'astnet_local_max_objs', 'int', 2000),
    (('wcs', 'astnet_local_cpulimit_s'), 'astnet_local_cpulimit_s', 'float', 30.0),
    (('wcs', 'astnet_blind_retry_on_fail'), 'astnet_blind_retry_on_fail', 'bool', True),
    (('wcs', 'astnet_blind_cpulimit_s'), 'astnet_blind_cpulimit_s', 'float', 120.0),
    (('wcs', 'platesolve_gaia_radius_scale'), 'platesolve_gaia_radius_scale', 'float', 1.35),
    (('wcs', 'max_workers'), 'wcs_max_workers', 'int', 1),
    (('wcs', 'require_qc_pass'), 'wcs_require_qc_pass', 'bool', True),
    (('wcs', 'refine_enable'), 'wcs_refine_enable', 'bool', True),
    (('wcs_qc', 'match_radius_arcsec'), 'wcs_qc_match_radius_arcsec'),
    (('wcs_qc', 'clip_sigma'), 'wcs_qc_clip_sigma', 'float', 3.0),
    (('wcs_qc', 'require_wcs_ok'), 'wcs_qc_require_wcs_ok', 'bool', True),
    (('wcs_qc', 'min_match_n'), 'wcs_qc_min_match_n'),
    (('wcs_qc', 'min_match_rate'), 'wcs_qc_min_match_rate'),
    (('wcs_qc', 'max_rms_px'), 'wcs_qc_max_rms_px', 'float', 2.5),
    (('wcs_qc', 'max_p99_px'), 'wcs_qc_max_p99_px'),
    (('wcs_qc', 'min_inlier_rate'), 'wcs_qc_min_inlier_rate', 'float', 0.50),
    (('wcs_qc', 'max_edge_ratio'), 'wcs_qc_max_edge_ratio', 'float', 0.0),
    (('wcs_qc', 'max_center_offset_arcsec'), 'wcs_qc_max_center_offset_arcsec', 'float', 0.0),
    (('wcs_refine', 'enable'), 'wcs_refine_enable', 'bool', True),
    (('wcs_refine', 'max_match'), 'wcs_refine_max_match', 'int', 600),
    (('wcs_refine', 'match_r_fwhm'), 'wcs_refine_match_r_fwhm', 'float', 1.6),
    (('wcs_refine', 'min_match'), 'wcs_refine_min_match', 'int', 50),
    (('gaia', 'radius_fudge'), 'gaia_radius_fudge', 'float', 1.35),
    (('gaia', 'mag_max'), 'gaia_mag_max', 'float', 18.0),
    (('gaia', 'wcs_mag_max'), 'gaia_wcs_mag_max', 'float', 18.0),
    (('refbuild', 'wcs_match_radius_arcsec'), 'ref_wcs_match_radius_arcsec', 'float', 2.0, {'label': 'WCS match radius (arcsec)', 'lo': 0.1, 'hi': 30.0, 'step': 0.1, 'decimals': 2}),
    (('gaia', 'snr_calib_min'), 'gaia_snr_calib_min', 'float', 20.0),
    (('gaia', 'gi_min'), 'gaia_gi_min', 'float', -0.5),
    (('gaia', 'gi_max'), 'gaia_gi_max', 'float', 4.5),
    (('gaia', 'retry'), 'gaia_retry', 'int', 2),
    (('gaia', 'timeout_s'), 'gaia_timeout_s', 'float', 30.0),
    (('gaia', 'hard_deadline_s'), 'gaia_hard_deadline_s', 'float', 0.0),
    (('gaia', 'backoff_s'), 'gaia_backoff_s', 'float', 6.0),
    (('gaia', 'allow_no_cache'), 'gaia_allow_no_cache', 'bool', True),
    (('gaia', 'g_limit'), 'idmatch_gaia_g_limit', 'float', 18.0, {'label': 'Gaia G limit (hybrid ID)', 'lo': 10.0, 'hi': 25.0, 'step': 0.5, 'decimals': 2}),
    (('idmatch', 'gaia_g_limit'), 'idmatch_gaia_g_limit', 'float', 18.0, {'label': 'Gaia G limit (hybrid ID)', 'lo': 10.0, 'hi': 25.0, 'step': 0.5, 'decimals': 2}),
    (('idmatch', 'match_r_fwhm'), 'idmatch_match_r_fwhm', 'float', 0.8),
    (('idmatch', 'two_pass_enable'), 'idmatch_two_pass_enable', 'bool', True),
    (('idmatch', 'tight_radius_arcsec'), 'idmatch_tight_radius_arcsec', 'float', 1.0),
    (('idmatch', 'loose_radius_arcsec'), 'idmatch_loose_radius_arcsec', 'float', 3.0),
    (('idmatch', 'adaptive_retry_threshold'), 'idmatch_adaptive_retry_threshold', 'float', 0.5),
    (('idmatch', 'fwhm_adaptive_floor'), 'idmatch_fwhm_adaptive_floor', 'bool', True),
    (('idmatch', 'geom_correction_enable'), 'idmatch_geom_correction_enable', 'bool', True),
    (('idmatch', 'min_correction_pairs'), 'idmatch_min_correction_pairs', 'int', 3),
    (('idmatch', 'min_affine_pairs'), 'idmatch_min_affine_pairs', 'int', 6),
    (('idmatch', 'tol_px'), 'idmatch_tol_px', 'float', 2.0),
    (('idmatch', 'tol_arcsec'), 'idmatch_tol_arcsec', 'float_or_none', None),
    (('idmatch', 'force'), 'force_idmatch', 'bool', False),
    (('idmatch', 'use_qc_pass_only'), 'idmatch_use_qc_pass_only'),
    (('idmatch', 'use_wcs_qc_gate'), 'idmatch_use_wcs_qc_gate'),
    (('idmatch', 'wcs_qc_min_match_rate'), 'idmatch_wcs_qc_min_match_rate'),
    (('idmatch', 'wcs_qc_min_match_n'), 'idmatch_wcs_qc_min_match_n'),
    (('idmatch', 'wcs_qc_max_rms_px'), 'idmatch_wcs_qc_max_rms_px'),
    (('idmatch', 'wcs_qc_min_inlier_rate'), 'idmatch_wcs_qc_min_inlier_rate'),
    (('idmatch', 'wcs_qc_max_p99_px'), 'idmatch_wcs_qc_max_p99_px'),
    (('master', 'n_master'), 'N_master'),
    (('master', 'iso_min_sep_pix'), 'master_iso_min_sep_pix'),
    (('master', 'keep_max'), 'master_keep_max'),
    (('master', 'flux_quantile'), 'master_flux_quantile'),
    (('master', 'filter_keep'), 'master_filter_keep'),
    (('master', 'ref_frame'), 'ref_frame', 'str', None),
    (('master', 'min_frames_xy'), 'master_min_frames_xy', 'int', 1),
    (('master', 'preserve_ids'), 'master_preserve_ids', 'bool', True),
    (('master', 'force_build'), 'force_master_build', 'bool', False),
    (('master_editor', 'search_radius_px'), 'search_radius_px', 'float', 7.0),
    (('master_editor', 'bulk_drop_box_px'), 'bulk_drop_box_px', 'int', 200),
    (('master_editor', 'gaia_add_max_sep_arcsec'), 'gaia_add_max_sep_arcsec', 'float', 2.0),
    (('match', 'tol_px'), 'match_tol_px', 'float', 1.0),
    (('match', 'wcs_radius_arcsec'), 'wcs_match_radius_arcsec'),
    (('match', 'min_gaia_matches'), 'min_master_gaia_matches', 'int', 10),
    (('match', 'pixel_scale_arcsec'), 'pixel_scale_arcsec', 'float_or_none', None),
    (('overlay', 'max_labels'), 'overlay_max_labels', 'int', 2000),
    (('overlay', 'label_fontsize'), 'overlay_label_fontsize', 'float', 6.0),
    (('overlay', 'label_offset_px'), 'overlay_label_offset_px', 'float', 3.0),
    (('overlay', 'show_id_when_no_mag'), 'overlay_show_id_when_no_mag', 'bool', False),
    (('overlay', 'use_phot_centroid'), 'overlay_use_phot_centroid', 'bool', True),
    (('overlay', 'show_ref_pos'), 'overlay_show_ref_pos', 'bool', True),
    (('overlay', 'show_shift_vectors'), 'overlay_show_shift_vectors', 'bool', False),
    (('overlay', 'shift_max_vectors'), 'overlay_shift_max_vectors', 'int', 300),
    (('overlay', 'shift_min_px'), 'overlay_shift_min_px', 'float', 1.5),
    (('overlay', 'inspect_index'), 'inspect_index', 'int', 0),
    (('ui', 'canvas_px'), 'ui_canvas_px'),
    (('transform', 'save_src2ref'), 'save_src2ref_tforms'),
    (('detection', 'peak', 'kernel_scales'), 'peak_kernel_scales', 'str', "0.9,1.3"),
    (('target', 'name'), 'target_name', 'str', ""),
    (('target', 'ra_deg'), 'target_ra_deg', 'float_or_none', None),
    (('target', 'dec_deg'), 'target_dec_deg', 'float_or_none', None),
    (('simbad', 'timeout_s'), 'simbad_timeout_s', 'float', 20.0),
    # --- Extinction (Airmass Fit) tool. Seventeen knobs the window offered and
    # could not persist; the file's `extinction_fit` section named none of them.
    (('extinction_fit', 'snr_min'), 'extinction_snr_min', 'float', 10.0),
    (('extinction_fit', 'min_good_stars'), 'extinction_min_good_stars', 'int', 3),
    (('extinction_fit', 'min_points_color'), 'extinction_min_points_color', 'int', 30),
    (('extinction_fit', 'min_points_quadratic'), 'extinction_min_points_quadratic', 'int', 20),
    (('extinction_fit', 'use_quadratic'), 'extinction_use_quadratic', 'bool', True),
    (('extinction_fit', 'use_color_dependent'), 'extinction_use_color_dependent', 'bool', True),
    (('extinction_fit', 'delta_x_enable'), 'extinction_delta_x_enable', 'bool', True),
    (('extinction_fit', 'delta_x_min'), 'extinction_delta_x_min', 'float', 0.3),
    (('extinction_fit', 'frame_qc_method'), 'extinction_frame_qc_method', 'str', "mad"),
    (('extinction_fit', 'frame_qc_sigma'), 'extinction_frame_qc_sigma', 'float', 3.0),
    (('extinction_fit', 'star_min_frames'), 'extinction_star_min_frames', 'int', 8),
    (('extinction_fit', 'star_rms_max'), 'extinction_star_rms_max', 'float', 0.10),
    (('extinction_fit', 'star_snr_med_min'), 'extinction_star_snr_med_min', 'float', 10.0),
    (('extinction_fit', 'star_use_weights'), 'extinction_star_use_weights', 'bool', True),
    (('extinction_fit', 'varstar_method'), 'extinction_varstar_method', 'str', "mad"),
    (('extinction_fit', 'varstar_min_frames'), 'extinction_varstar_min_frames', 'int', 5),
    (('extinction_fit', 'varstar_sigma'), 'extinction_varstar_sigma', 'float', 3.0),
    # --- Reference build and airmass windows, same shape.
    (('refbuild', 'build_mode'), 'ref_build_mode', 'str', "hybrid"),
    (('refbuild', 'compare_exclude_split'), 'step6_compare_exclude_split', 'bool', True),
    (('airmass', 'update_source'), 'airmass_update_source', 'str', "auto"),
    # --- Step 6 (master catalogue) is a shared step, and these lived in the
    # LC-only map, so the CMD window offered them and could save none of them.
    # Every default here is the literal the code was already using (2026-08-16).
    (('refbuild', 'sat_drop_pct'), 'ref_select_sat_pct', 'float', 20.0, {'label': 'Drop top saturation frames', 'lo': 0.0, 'hi': 100.0, 'step': 1.0, 'decimals': 1, 'suffix': '%'}),
    (('refbuild', 'elong_drop_pct'), 'ref_select_elong_pct', 'float', 20.0, {'label': 'Drop top elongation frames', 'lo': 0.0, 'hi': 100.0, 'step': 1.0, 'decimals': 1, 'suffix': '%'}),
    (('refbuild', 'per_date'), 'ref_per_date', 'bool', True, {'label': 'Per-date reference'}),
    (('refbuild', 'master_union'), 'ref_master_union', 'bool', True, {'label': 'Union master (all frames)'}),
    (('refbuild', 'union_min_frames'), 'ref_union_min_frames', 'int', 1, {'label': 'Union min detections/star', 'lo': 1, 'hi': 1000}),
    (('refbuild', 'ref_cat_max_sources'), 'ref_cat_max_sources', 'int', 0, {'label': 'Ref catalog max sources (0=all)', 'lo': 0, 'hi': 50000}),
    (('refbuild', 'ref_cat_min_sources'), 'ref_cat_min_sources', 'int', 50, {'label': 'Ref catalog min sources', 'lo': 0, 'hi': 50000}),
    (('refbuild', 'ref_cat_max_elong'), 'ref_cat_max_elong', 'float', 1.5, {'label': 'Ref max elongation', 'lo': 0.0, 'hi': 10.0, 'step': 0.1, 'decimals': 2}),
    (('refbuild', 'ref_cat_max_abs_round'), 'ref_cat_max_abs_round', 'float', 0.4, {'label': 'Ref max |roundness|', 'lo': 0.0, 'hi': 5.0, 'step': 0.05, 'decimals': 2}),
    (('refbuild', 'ref_cat_sharp_min'), 'ref_cat_sharp_min', 'float', 0.2, {'label': 'Ref sharpness min', 'lo': -5.0, 'hi': 5.0, 'step': 0.1, 'decimals': 2}),
    (('refbuild', 'ref_cat_sharp_max'), 'ref_cat_sharp_max', 'float', 1.0, {'label': 'Ref sharpness max', 'lo': -5.0, 'hi': 5.0, 'step': 0.1, 'decimals': 2}),
    (('refbuild', 'ref_cat_min_peak_adu'), 'ref_cat_min_peak_adu', 'float', 0.0, {'label': 'Ref min peak/flux (0=off)', 'lo': 0.0, 'hi': 1000000000.0, 'step': 1.0, 'decimals': 1}),
    (('refbuild', 'wcs_min_match_rate'), 'ref_wcs_min_match_rate', 'float', 0.2, {'label': 'WCS min match rate', 'lo': 0.0, 'hi': 1.0, 'step': 0.05, 'decimals': 2}),
    (('refbuild', 'wcs_min_match_n'), 'ref_wcs_min_match_n', 'int', 50, {'label': 'WCS min match count', 'lo': 0, 'hi': 100000}),
    (('refbuild', 'wcs_max_sep_med_arcsec'), 'ref_wcs_max_sep_med_arcsec', 'float', 1.5, {'label': 'WCS max sep median (arcsec)', 'lo': 0.0, 'hi': 30.0, 'step': 0.1, 'decimals': 2}),
    (('refbuild', 'wcs_max_sep_p90_arcsec'), 'ref_wcs_max_sep_p90_arcsec', 'float', 2.5, {'label': 'WCS max sep p90 (arcsec)', 'lo': 0.0, 'hi': 60.0, 'step': 0.1, 'decimals': 2}),
    (('refbuild', 'wcs_max_dup_rate'), 'ref_wcs_max_dup_rate', 'float', 0.1, {'label': 'WCS max duplicate rate', 'lo': 0.0, 'hi': 1.0, 'step': 0.05, 'decimals': 2}),
    # --- The airmass and extinction tools open in both modes, and the night
    # parsing keys sit in every config's `io` section, yet all of these lived
    # in the LC-only map: CMD read none of them and its windows could save
    # none of them. Defaults match the literals CMD was using (2026-08-16).
    (('io', 'night_parse_mode'), 'night_parse_mode', 'str', "regex"),
    (('io', 'night_parse_regex'), 'night_parse_regex', 'str', r".*_(\d{8})"),
    (('io', 'night_parse_split_delim'), 'night_parse_split_delim', 'str', "_"),
    (('io', 'night_parse_split_index'), 'night_parse_split_index', 'int', -1),
    (('io', 'night_parse_last_digits'), 'night_parse_last_digits', 'int', 8),
    (('io', 'airmass_formula'), 'airmass_formula', 'str', "Kasten & Young (1989)"),
    (('io', 'airmass_update_mode'), 'airmass_update_mode', 'str', "overwrite"),
    (('extinction_fit', 'clip_sigma'), 'extfit_clip_sigma', 'float', 3.0),
    (('extinction_fit', 'fit_iters'), 'extfit_fit_iters', 'int', 5),
    # --- Steps 4 and 5 are shared; these were CMD-only, so the LC windows
    # offered them and could not save them (2026-08-16).
    (('detection', 'mode'), 'detect_mode'),
    (('wcs', 'astap_annotate_variables'), 'astap_annotate_variables', 'bool', False),
)

CMD_ONLY_TOML_KEY_MAP: tuple[tuple, ...] = (
    (('detection', 'deblend', 'mode'), 'deblend_mode'),
    (('wcs', 'header_coord_warn_sep_deg'), 'wcs_header_coord_warn_sep_deg', 'float', 0.5),
    (('wcs', 'header_coord_max_sep_deg'), 'wcs_header_coord_max_sep_deg', 'float', 5.0),
    (('gaia', 'match_tol_arcsec'), 'ref_wcs_match_radius_arcsec', 'float', 2.0, {'label': 'WCS match radius (arcsec)', 'lo': 0.1, 'hi': 30.0, 'step': 0.1, 'decimals': 2}),
    (('gaia', 'derived_enable'), 'gaia_derived_enable', 'bool', True),
    (('gaia', 'pmem_method'), 'gaia_pmem_method'),
    (('gaia', 'pmem_ruwe_max'), 'gaia_pmem_ruwe_max', 'float', 2.0),
    (('gaia', 'pmem_min_visibility_periods'), 'gaia_pmem_min_visibility_periods', 'float', 8.0),
    (('gaia', 'pmem_min_valid'), 'gaia_pmem_min_valid', 'int', 30),
    (('gaia', 'pmem_min_fit'), 'gaia_pmem_min_fit', 'int', 25),
    (('idmatch', 'mode'), 'idmatch_mode'),
    (('idmatch', 'use_gaia_refs_only'), 'idmatch_use_gaia_refs_only', 'bool', False),
    (('master_editor', 'membership_overlay_enable'), 'step8_membership_overlay_enable', 'bool', False),
    (('master_editor', 'membership_threshold'), 'step8_membership_threshold', 'float', 0.5),
    (('cmd', 'membership_mode'), 'cmd_membership_mode', 'str', "normal"),
    (('cmd', 'membership_compare'), 'cmd_membership_compare', 'bool', True),
    (('instrument', 'gain_e_per_adu'), 'gain_e_per_adu', 'float_or_none', None),
    (('instrument', 'rdnoise_e'), 'rdnoise_e'),
    (('instrument', 'noise_use_fits_header'), 'noise_use_fits_header', 'bool', False),
    (('instrument', 'noise_reference_binning'), 'noise_reference_binning', 'float_or_none', None),
    (('instrument', 'noise_scale_by_binning'), 'noise_scale_by_binning', 'bool', False),
    (('instrument', 'saturation_adu'), 'saturation_adu', 'float', 60000.0),
    (('instrument', 'datamin_adu'), 'datamin_adu', 'float', 0.1),
    (('instrument', 'datamax_adu'), 'datamax_adu', 'float', 60000.0),
    (('instrument', 'binning'), 'binning_default', 'int', 2),
    (('instrument', 'telescope_focal_mm'), 'telescope_focal_mm', 'float', 3947.0),
    (('instrument', 'camera_pixel_um'), 'camera_pixel_um', 'float', 3.76),
    (('instrument', 'zp_initial'), 'zp_initial', 'float', 25.0, {'label': 'Preview ZP', 'lo': 0.0, 'hi': 50.0, 'step': 0.1, 'decimals': 2, 'tooltip': 'Display-only zero point for Step 3 preview magnitude.'}),
    (('site', 'lat_deg'), 'site_lat_deg', 'float', 0.0),
    (('site', 'lon_deg'), 'site_lon_deg', 'float', 0.0),
    (('site', 'alt_m'), 'site_alt_m', 'float', 0.0),
    (('site', 'tz_offset_hours'), 'site_tz_offset_hours', 'float', 0.0),
    (('cmd', 'psf_match_radius_px'), 'psf_cmd_match_radius_px', 'float', 1.0),
    # --- Isochrone fit (Step 12). `colors` has no default on purpose: the step
    # refuses to run without it, because a default age window cannot reach a
    # globular and a missing reddening prior rails an open cluster at the floor.
    (('isochrone', 'colors'), 'iso_colors', 'str', ""),
    (('isochrone', 'mag_band'), 'iso_mag_band', 'str', "g"),
    (('isochrone', 'age_min'), 'iso_age_min', 'float', 8.5),
    (('isochrone', 'age_max'), 'iso_age_max', 'float', 10.2),
    # The [M/H] box must reach a globular. -1.0 does not: M13 sits near -1.5, and
    # because the fit narrows this box around an [M/H] prior and never widens it,
    # a metal-poor prior inside a (-1.0, 0.5) box leaves lo > hi and the grid
    # comes out empty. -2.2 is the floor of the run that recovered M13's
    # literature metallicity (-1.559).
    (('isochrone', 'mh_min'), 'iso_mh_min', 'float', -2.2),
    (('isochrone', 'mh_max'), 'iso_mh_max', 'float', 0.5),
    (('isochrone', 'dm_min'), 'iso_dm_min', 'float', 5.0),
    (('isochrone', 'dm_max'), 'iso_dm_max', 'float', 18.0),
    (('isochrone', 'ecolor_min'), 'iso_ecolor_min', 'float', 0.0),
    # E(colour), not E(B-V) — for B-V they are nearly the same (R_B - R_V ~ 1.0).
    # 1.0 is what the desktop dialog has always used; 0.5 would have been a
    # silent narrowing of it.
    (('isochrone', 'ecolor_max'), 'iso_ecolor_max', 'float', 1.0),
    (('isochrone', 'mh_prior'), 'iso_mh_prior', 'str', ""),
    (('isochrone', 'ecolor_prior'), 'iso_ecolor_prior', 'str', ""),
    (('isochrone', 'dm_prior'), 'iso_dm_prior', 'str', ""),
    (('isochrone', 'parallax_distance_prior'), 'iso_parallax_distance_prior', 'bool', False),
    (('isochrone', 'parallax_dm_sigma'), 'iso_parallax_dm_sigma', 'float', 0.05),
    (('isochrone', 'parallax_dm_window'), 'iso_parallax_dm_window', 'float', 0.06),
    (('isochrone', 'use_membership'), 'iso_use_membership', 'bool', True),
    (('isochrone', 'data_snr_min'), 'iso_data_snr_min', 'float', 20.0),
    (('isochrone', 'fit_snr_min'), 'iso_fit_snr_min', 'float', 5.0),
    (('isochrone', 'max_stars'), 'iso_max_stars', 'int', 500),
    (('isochrone', 'n_walkers'), 'iso_n_walkers', 'int', 32),
    (('isochrone', 'n_burn'), 'iso_n_burn', 'int', 600),
    (('isochrone', 'n_steps'), 'iso_n_steps', 'int', 2000),
    (('isochrone', 'f_bin'), 'iso_f_bin', 'float', 0.3),
    (('isochrone', 'f_field'), 'iso_f_field', 'float', 0.1),
    (('isochrone', 'err_floor'), 'iso_err_floor', 'float', 0.02),
    (('isochrone', 'seed'), 'iso_seed', 'int', 2024),
)

PSF_TOML_KEY_MAP: tuple[tuple, ...] = (
    (('psf', 'epsf_oversampling'), 'psf_epsf_oversampling', 'int', 2),
    (('psf', 'interp_order'), 'psf_interp_order'),
    (('psf', 'epsf_maxiters'), 'psf_epsf_maxiters', 'int', 5),
    (('psf', 'epsf_size_px'), 'psf_epsf_size_px', 'int', 25),
    (('psf', 'epsf_size_fwhm_mult'), 'psf_epsf_size_fwhm_mult', 'float', 4.0),
    (('psf', 'n_stars_max'), 'psf_n_stars_max', 'int', 0),
    (('psf', 'isolation_fwhm_mult'), 'psf_isolation_fwhm_mult', 'float', 3.0),
    (('psf', 'epsf_contamination_filter'), 'psf_epsf_contamination_filter'),
    (('psf', 'flux_scale_correction'), 'psf_flux_scale_correction'),
    (('psf', 'flux_scale_min_snr'), 'psf_flux_scale_min_snr', 'float', 50.0),
    (('psf', 'flux_scale_min_stars'), 'psf_flux_scale_min_stars', 'int', 8),
    (('psf', 'flux_scale_min_neighbor_fwhm'), 'psf_flux_scale_min_neighbor_fwhm'),
    (('psf', 'flux_scale_max_scatter_mag'), 'psf_flux_scale_max_scatter_mag'),
    (('psf', 'flux_percentile_lo'), 'psf_flux_percentile_lo', 'float', 75.0),
    (('psf', 'flux_percentile_hi'), 'psf_flux_percentile_hi', 'float', 95.0),
    (('psf', 'fit_shape_px'), 'psf_fit_shape_px', 'int', 5),
    (('psf', 'fit_shape_fwhm_mult'), 'psf_fit_shape_fwhm_mult', 'float', 2.4),
    (('psf', 'fit_window_mode'), 'psf_fit_window_mode'),
    (('psf', 'fit_encircled_energy'), 'psf_fit_encircled_energy'),
    (('psf', 'max_iter'), 'psf_max_iter', 'int', 2),
    (('psf', 'fitter_max_iter'), 'psf_fitter_max_iter', 'int', 6),
    (('psf', 'fit_mode'), 'psf_fit_mode'),
    (('psf', 'redetect_sigma'), 'psf_redetect_sigma', 'float', 3.5),
    (('psf', 'redetect_sigma_g'), 'psf_redetect_sigma_g', 'float', float("nan")),
    (('psf', 'redetect_sigma_r'), 'psf_redetect_sigma_r', 'float', float("nan")),
    (('psf', 'redetect_sigma_i'), 'psf_redetect_sigma_i', 'float', float("nan")),
    (('psf', 'epsf_sharp_lo'), 'psf_epsf_sharp_lo', 'float', 0.3),
    (('psf', 'epsf_sharp_hi'), 'psf_epsf_sharp_hi', 'float', 0.8),
    (('psf', 'epsf_round_abs_max'), 'psf_epsf_round_abs_max', 'float', 0.5),
    (('psf', 'epsf_elong_max'), 'psf_epsf_elong_max', 'float', 1.3),
    (('psf', 'model_mode'), 'psf_model_mode'),
    (('psf', 'fit_engine'), 'psf_fit_engine'),
    (('psf', 'build_mode'), 'psf_build_mode'),
    (('psf', 'field_mode'), 'psf_mode'),
    (('psf', 'parallel_workers'), 'psf_parallel_workers', 'int', 0),
    (('psf', 'duplicate_radius_fwhm_mult'), 'psf_duplicate_radius_fwhm_mult', 'float', 0.8),
    (('psf', 'duplicate_radius_px'), 'psf_duplicate_radius_px', 'float', float("nan")),
    (('psf', 'new_sources_cap_per_iter'), 'psf_new_sources_cap_per_iter', 'int', 70),
    (('psf', 'new_sources_cap_frac'), 'psf_new_sources_cap_frac', 'float', 0.02),
    (('psf', 'fit_init_max_sources'), 'psf_fit_init_max_sources', 'int', 0),
    (('psf', 'core_cut_enable'), 'psf_core_cut_enable', 'bool', False),
    (('psf', 'core_cut_center_mode'), 'psf_core_cut_center_mode'),
    (('psf', 'core_cut_x_px'), 'psf_core_cut_x_px', 'float', 0.0),
    (('psf', 'core_cut_y_px'), 'psf_core_cut_y_px', 'float', 0.0),
    (('psf', 'core_cut_radius_px'), 'psf_core_cut_radius_px', 'float', 0.0),
    (('psf', 'core_cut_radius_fwhm_mult'), 'psf_core_cut_radius_fwhm_mult', 'float', 20.0),
    (('psf', 'core_cut_auto_cell_fwhm_mult'), 'psf_core_cut_auto_cell_fwhm_mult', 'float', 8.0),
    (('psf', 'core_cut_auto_min_density_ratio'), 'psf_core_cut_auto_min_density_ratio', 'float', 1.5),
    (('psf', 'core_cut_auto_min_sources'), 'psf_core_cut_auto_min_sources', 'int', 50),
    (('psf', 'core_cut_max_exclude_frac'), 'psf_core_cut_max_exclude_frac', 'float', 0.70),
    (('psf', 'substar_iters'), 'psf_substar_iters', 'int', 1),
    (('psf', 'substar_neighbor_r_fwhm_mult'), 'psf_substar_neighbor_r_fwhm_mult', 'float', 8.0),
    (('psf', 'substar_max_sources'), 'psf_substar_max_sources', 'int', 1500),
    (('psf', 'conv_new_frac'), 'psf_conv_new_frac', 'float', 0.02),
    (('psf', 'postfit_snr_min'), 'psf_postfit_snr_min', 'float', 3.0),
    (('psf', 'postfit_qfit_max'), 'psf_postfit_qfit_max', 'float', 3.0),
    (('psf', 'postfit_reduced_chi2_max'), 'psf_postfit_reduced_chi2_max', 'float', 25.0),
    (('psf', 'blend_residual_ratio'), 'psf_blend_residual_ratio', 'float', 0.3),
    (('psf', 'use_grouper'), 'psf_use_grouper', 'bool', False),
    (('psf', 'grouper_max_size'), 'psf_grouper_max_size', 'int', 3),
    (('psf', 'grouper_radius_fwhm'), 'psf_grouper_radius_fwhm', 'float', 1.5),
    (('psf', 'grouper_budget_frac'), 'psf_grouper_budget_frac'),
    (('psf', 'grouper_budget_cap'), 'psf_grouper_budget_cap'),
    (('psf', 'final_pass_max_iter'), 'psf_final_pass_max_iter'),
    (('psf', 'forced_position_lock'), 'psf_forced_position_lock'),
    (('psf', 'profile_error_frac'), 'psf_profile_error_frac'),
    (('psf', 'forced_match_radius_fwhm'), 'psf_forced_match_radius_fwhm', 'float', 1.25),
    (('psf', 'redetect_sharp_lo'), 'psf_redetect_sharp_lo', 'float', 0.15),
    (('psf', 'redetect_sharp_hi'), 'psf_redetect_sharp_hi', 'float', 0.95),
    (('psf', 'redetect_round_abs_max'), 'psf_redetect_round_abs_max', 'float', 0.8),
    (('psf', 'flux_conv_threshold'), 'psf_flux_conv_threshold', 'float', 0.01),
    (('psf', 'use_error_image'), 'psf_use_error_image', 'bool', False),
    (('psf', 'save_residuals'), 'psf_save_residuals', 'bool', True),
    (('psf', 'save_model_image'), 'psf_save_model_image', 'bool', True),
    (('psf', 'shared_filter_epsf'), 'psf_shared_filter_epsf', 'bool', False),
    (('psf', 'save_all_iter_residuals'), 'psf_save_all_iter_residuals', 'bool', False),
    (('psf', 'min_epsf_stars'), 'psf_min_epsf_stars', 'int', 10),
)

CMD_ONLY_TOML_KEY_MAP += (
    (('cross_frame', 'ransac_tol_px'), 'cross_frame_ransac_tol_px'),
    (('cross_frame', 'ransac_max_iter'), 'cross_frame_ransac_max_iter'),
    (('cross_frame', 'ransac_min_inliers'), 'cross_frame_ransac_min_inliers'),
    (('cross_frame', 'match_tol_px'), 'cross_frame_match_tol_px'),
    # ── CMD-only keys migrated from COMMON_TOML_KEY_MAP ──────────────────
    (('cmd', 'snr_calib_min'), 'cmd_snr_calib_min', 'float', 20.0),
    (('cmd', 'max_sources'), 'cmd_max_sources'),
    (('cmd', 'apply_extinction'), 'cmd_apply_extinction', 'bool', False),
    (('cmd', 'extinction_mode'), 'cmd_extinction_mode', 'str', "absorb"),
    (('cmd', 'frame_zp_min_n'), 'frame_zp_min_n', 'int', 5),
    (('cmd', 'zp', 'clip_sigma'), 'zp_clip_sigma', 'float', 3.0),
    (('cmd', 'zp', 'fit_iters'), 'zp_fit_iters', 'int', 5),
    (('cmd', 'zp', 'slope_absmax'), 'zp_slope_absmax', 'float', 1.0),
    (('cmd', 'color', 'clip_sigma'), 'color_clip_sigma'),
    (('cmd', 'color', 'fit_iters'), 'color_fit_iters'),
    (('cmd', 'color', 'slope_absmax'), 'color_slope_absmax'),
    # ── Step 10 external standard-star anchor (Gaia-independent ZP re-anchor) ──
    (('cmd', 'standard_anchor', 'enable'), 'std_anchor_enable'),
    (('cmd', 'standard_anchor', 'catalog'), 'std_anchor_catalog'),
    (('cmd', 'standard_anchor', 'match_radius_arcsec'), 'std_anchor_match_radius'),
    (('cmd', 'standard_anchor', 'min_stars'), 'std_anchor_min_stars'),
    (('gaia', 'zp_slope_absmax'), 'gaia_zp_slope_absmax', 'float', 1.0),
    (('gaia', 'color_slope_absmax'), 'gaia_color_slope_absmax', 'float', 2.0),
    (('isochrone', 'file_path'), 'iso_file_path', 'str', ""),
    (('isochrone', 'age_init'), 'iso_age_init', 'float', 9.7),
    (('isochrone', 'mh_init'), 'iso_mh_init', 'float', -0.1),
    (('isochrone', 'eg_r_init'), 'iso_eg_r_init', 'float', 0.0033),
    (('isochrone', 'dm_init'), 'iso_dm_init', 'float', 9.46),
)

LC_ONLY_TOML_KEY_MAP: tuple[tuple, ...] = (
    (('io', 'night_parse_include_unmatched'), 'night_parse_include_unmatched', 'bool', False),
    (('io', 'airmass_update_header'), 'airmass_update_header', 'bool', False),
    (('instrument', 'gain_e_per_adu'), 'gain_e_per_adu', 'float_or_none', None),
    (('instrument', 'rdnoise_e'), 'rdnoise_e'),
    (('instrument', 'noise_use_fits_header'), 'noise_use_fits_header', 'bool', False),
    (('instrument', 'noise_reference_binning'), 'noise_reference_binning', 'float_or_none', None),
    (('instrument', 'noise_scale_by_binning'), 'noise_scale_by_binning', 'bool', False),
    (('instrument', 'saturation_adu'), 'saturation_adu', 'float', 60000.0),
    (('instrument', 'binning'), 'binning_default', 'int', 2),
    (('instrument', 'telescope_focal_mm'), 'telescope_focal_mm', 'float', 3947.0),
    (('instrument', 'camera_pixel_um'), 'camera_pixel_um', 'float', 3.76),
    (('instrument', 'zp_initial'), 'zp_initial', 'float', 25.0),
    (('site', 'lat_deg'), 'site_lat_deg', 'float', 0.0),
    (('site', 'lon_deg'), 'site_lon_deg', 'float', 0.0),
    (('site', 'alt_m'), 'site_alt_m', 'float', 0.0),
    (('site', 'tz_offset_hours'), 'site_tz_offset_hours', 'float', 0.0),
    (('photometry', 'use_original_frames'), 'phot_use_original_frames'),
    (('wcs', 'propagate_max_shift_px'), 'wcs_propagate_max_shift_px', 'float', 50.0),
    (('wcs', 'propagate_min_match'), 'wcs_propagate_min_match', 'int', 12),
    (('wcs', 'propagate_sigma_clip'), 'wcs_propagate_sigma_clip', 'float', 3.0),
    (('idmatch', 'init_r_fwhm'), 'idmatch_init_r_fwhm', 'float', 5.0),
    (('idmatch', 'ratio_max'), 'idmatch_ratio_max', 'float', 0.7),
    (('idmatch', 'min_pairs'), 'idmatch_min_pairs', 'int', 15),
    (('idmatch', 'transform_mode'), 'idmatch_transform_mode', 'str', "similarity"),
    (('idmatch', 'mutual_nearest'), 'idmatch_mutual_nearest', 'bool', True),
    (('light_curve', 'color_index_by_filter'), 'lightcurve_color_index_by_filter'),
    (('light_curve', 'color_term_by_filter'), 'lightcurve_color_term_by_filter'),
    (('extinction_fit', 'order'), 'extfit_order', 'int', 1),
    (('extinction_fit', 'min_points'), 'extfit_min_points', 'int', 5),
    (('extinction_fit', 'use_color_terms'), 'extfit_use_color_terms', 'bool', False),
    (('extinction_fit', 'color_index_by_filter'), 'extfit_color_index_by_filter'),
    (('extinction_fit', 'color_c1_by_filter'), 'extfit_color_c1_by_filter'),
    (('extinction_fit', 'color_c2_by_filter'), 'extfit_color_c2_by_filter'),
    # Keys read by the shared Extinction (Airmass Fit) tool — available in LC mode
    # too (registry modes default to cmd+lc). These were over-eagerly moved to
    # CMD_ONLY; restore them for LC so parameters.toml values are honored there.
    (('cmd', 'snr_calib_min'), 'cmd_snr_calib_min', 'float', 20.0),
    (('cmd', 'zp', 'clip_sigma'), 'zp_clip_sigma', 'float', 3.0),
    (('cmd', 'zp', 'fit_iters'), 'zp_fit_iters', 'int', 5),
    (('light_curve', 'comparison_auto_pool_max'), 'comparison_auto_pool_max', 'int', 30),
    (('light_curve', 'comparison_auto_ensemble_max'), 'comparison_auto_ensemble_max', 'int', 12),
)

# ── Canonical composed maps — callers must use these ─────────────────────────
# (The expanded hand-written literals that lived here were removed in favour of
#  the authoritative composition below.  Use git blame to recover them.)
CMD_TOML_KEY_MAP = COMMON_TOML_KEY_MAP + PSF_TOML_KEY_MAP + CMD_ONLY_TOML_KEY_MAP
LC_TOML_KEY_MAP = COMMON_TOML_KEY_MAP + PSF_TOML_KEY_MAP + LC_ONLY_TOML_KEY_MAP


def path_to_dotted(path: Iterable[str]) -> str:
    return ".".join(str(part) for part in path)


def get_toml_path(data: dict[str, Any], path: Iterable[str]) -> Any:
    """Return a nested TOML value, or None when any path segment is missing."""

    cur: Any = data
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def set_toml_path(data: dict[str, Any], path: Iterable[str], value: Any) -> None:
    """Set a nested TOML value, creating intermediate dictionaries as needed."""

    cur: Any = data
    keys = list(path)
    if not keys:
        raise ValueError("TOML path must not be empty")
    for key in keys[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]
    cur[keys[-1]] = value


def toml_value_for_runtime_attr(attr: str, value: Any, path_keys: Iterable[str]) -> Any:
    """Convert a flat runtime parameter value to a TOML-compatible value."""

    path_tuple = tuple(path_keys)
    if isinstance(value, Path):
        if path_tuple == ("io", "cache_dir"):
            return value.name
        return str(value)
    if attr == "peak_kernel_scales":
        if isinstance(value, str):
            items = [v.strip() for v in value.split(",") if v.strip()]
            out: list[Any] = []
            for item in items:
                try:
                    out.append(float(item))
                except Exception:
                    out.append(item)
            return out
        if isinstance(value, (list, tuple)):
            return list(value)
    return value


def duplicate_toml_paths(
    mapping: Iterable[tuple],
) -> dict[str, list[str]]:
    """Return TOML paths mapped more than once."""

    paths: dict[str, list[str]] = {}
    for row in mapping:
        paths.setdefault(path_to_dotted(row[0]), []).append(row[1])
    return {path: attrs for path, attrs in paths.items() if len(attrs) > 1}


def duplicate_runtime_attrs(
    mapping: Iterable[tuple],
) -> dict[str, list[str]]:
    """Return runtime attrs intentionally or accidentally mapped from multiple paths."""

    attrs: dict[str, list[str]] = {}
    for row in mapping:
        attrs.setdefault(row[1], []).append(path_to_dotted(row[0]))
    return {attr: sorted(paths) for attr, paths in attrs.items() if len(paths) > 1}


def toml_key_map_for_mode(mode: str) -> tuple[tuple, ...]:
    """Return the full TOML key map for a workflow mode."""

    mode_key = str(mode).strip().lower()
    if mode_key == "cmd":
        return CMD_TOML_KEY_MAP
    if mode_key == "lc":
        return LC_TOML_KEY_MAP
    raise ValueError(f"Unknown parameter mode: {mode!r}")


def read_schema_version(data: dict[str, Any] | None) -> int:
    """Return top-level schema_version, treating missing/invalid as legacy."""

    if not isinstance(data, dict):
        return LEGACY_SCHEMA_VERSION
    value = data.get("schema_version", LEGACY_SCHEMA_VERSION)
    try:
        return int(value)
    except Exception:
        return LEGACY_SCHEMA_VERSION


def ensure_schema_version(
    data: dict[str, Any],
    version: int = CANONICAL_SCHEMA_VERSION,
) -> dict[str, Any]:
    """Set the canonical top-level schema_version in a TOML data dict."""

    data["schema_version"] = int(version)
    return data


def entries_for_mode(mode: str) -> tuple[ParameterMapEntry, ...]:
    """Return foundation entries relevant to a workflow mode."""

    mode_key = str(mode).strip().lower()
    return tuple(entry for entry in COMMON_FOUNDATION_MAP if mode_key in entry.modes)


def legacy_aliases_for_mode(mode: str) -> tuple[ParameterMapEntry, ...]:
    """Return known legacy aliases relevant to a workflow mode."""

    mode_key = str(mode).strip().lower()
    return tuple(entry for entry in LEGACY_ALIAS_MAP if mode_key in entry.modes)
