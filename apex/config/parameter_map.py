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

COMMON_TOML_KEY_MAP: tuple[tuple[tuple[str, ...], str], ...] = (
    (('io', 'data_dir'), 'data_dir'),
    (('io', 'filename_prefix'), 'filename_prefix'),
    (('io', 'result_dir'), 'result_dir'),
    (('io', 'cache_dir'), 'cache_dir'),
    (('io', 'night_gap_hours'), 'night_gap_hours'),
    (('parallel', 'mode'), 'parallel_mode'),
    (('parallel', 'max_workers'), 'max_workers'),
    (('parallel', 'resume_mode'), 'resume_mode'),
    (('parallel', 'force_redetect'), 'force_redetect'),
    (('parallel', 'force_rephot'), 'force_rephot'),
    (('parallel', 'detect_cache_strategy'), 'detect_cache_strategy'),
    (('ui', 'log_tail'), 'ui_log_tail'),
    (('ui', 'detect_progress_bar'), 'detect_progress_bar'),
    (('alignment', 'ref_index'), 'align_ref_index'),
    (('alignment', 'global_align'), 'global_align'),
    (('alignment', 'global_ref_filter'), 'global_ref_filter'),
    (('alignment', 'global_ref_index'), 'global_ref_index'),
    (('fwhm', 'guess_px'), 'fwhm_pix_guess'),
    (('fwhm', 'guess_arcsec'), 'fwhm_guess_arcsec'),
    (('fwhm', 'arcsec_min'), 'fwhm_arcsec_min'),
    (('fwhm', 'arcsec_max'), 'fwhm_arcsec_max'),
    (('fwhm', 'px_min'), 'fwhm_px_min'),
    (('fwhm', 'px_max'), 'fwhm_px_max'),
    (('fwhm', 'qc_max_sources'), 'fwhm_qc_max_sources'),
    (('fwhm', 'elong_max'), 'fwhm_elong_max'),
    (('fwhm', 'iso_min_sep_pix'), 'iso_min_sep_pix'),
    (('fwhm', 'measure_max'), 'fwhm_measure_max'),
    (('fwhm', 'min_sources'), 'fwhm_min_sources'),
    (('fwhm', 'candidate_max'), 'fwhm_candidate_max'),
    (('fwhm', 'measure_all_sources'), 'fwhm_measure_all_sources'),
    (('fwhm', 'dr'), 'fwhm_dr'),
    (('clip', 'min_adu'), 'clip_min_adu'),
    (('clip', 'max_adu'), 'clip_max_adu'),
    (('detection', 'engine'), 'detect_engine'),
    (('detection', 'sigma'), 'detect_sigma'),
    (('detection', 'sigma_g'), 'detect_sigma_g'),
    (('detection', 'sigma_r'), 'detect_sigma_r'),
    (('detection', 'sigma_i'), 'detect_sigma_i'),
    (('detection', 'sigma_by_filter'), 'detect_sigma_by_filter'),
    (('detection', 'minarea_pix'), 'minarea_pix'),
    (('detection', 'keep_max'), 'detect_keep_max'),
    (('detection', 'dilate_radius_px'), 'segm_dilate_radius_px'),
    (('detection', 'deblend', 'enable'), 'deblend_enable'),
    (('detection', 'deblend', 'nthresh'), 'deblend_nthresh'),
    (('detection', 'deblend', 'contrast'), 'deblend_cont'),
    (('detection', 'deblend', 'max_labels'), 'deblend_max_labels'),
    (('detection', 'deblend', 'label_hard_max'), 'deblend_label_hard_max'),
    (('detection', 'deblend', 'nlevels_soft'), 'deblend_nlevels_soft'),
    (('detection', 'deblend', 'contrast_soft'), 'deblend_contrast_soft'),
    (('detection', 'peak', 'enable'), 'peak_pass_enable'),
    (('detection', 'peak', 'nsigma'), 'peak_nsigma'),
    (('detection', 'peak', 'min_sep_px'), 'peak_min_sep_px'),
    (('detection', 'peak', 'max_add'), 'peak_max_add'),
    (('detection', 'peak', 'max_elong'), 'peak_max_elong'),
    (('detection', 'peak', 'sharp_lo'), 'peak_sharp_lo'),
    (('detection', 'peak', 'skip_if_nsrc_ge'), 'peak_skip_if_nsrc_ge'),
    (('detection', 'dao', 'enable'), 'dao_refine_enable'),
    (('detection', 'dao', 'fwhm_px'), 'dao_fwhm_px'),
    (('detection', 'dao', 'sharp_lo'), 'dao_sharp_lo'),
    (('detection', 'dao', 'sharp_hi'), 'dao_sharp_hi'),
    (('detection', 'dao', 'round_lo'), 'dao_round_lo'),
    (('detection', 'dao', 'round_hi'), 'dao_round_hi'),
    (('detection', 'dao', 'match_tol_px'), 'dao_match_tol_px'),
    (('background', 'enable'), 'bkg2d_enable'),
    (('background', 'in_detect'), 'bkg2d_in_detect'),
    (('background', 'box'), 'bkg2d_box'),
    (('background', 'filter_size'), 'bkg2d_filter_size'),
    (('background', 'edge_method'), 'bkg2d_edge_method'),
    (('background', 'method'), 'bkg2d_method'),
    (('background', 'downsample'), 'bkg2d_downsample'),
    (('qc', 'gate_enable'), 'gate_enable'),
    (('qc', 'sky_sigma_max_e'), 'gate_sky_sigma_max_e'),
    (('qc', 'nsrc_min'), 'gate_nsrc_min'),
    (('qc', 'keep_positions_if_fail'), 'keep_positions_if_qc_fail'),
    (('photometry', 'mode'), 'aperture_mode'),
    (('photometry', 'recenter'), 'recenter_aperture'),
    (('photometry', 'use_segm_mask'), 'bkg_use_segm_mask'),
    (('photometry', 'min_snr_for_mag'), 'min_snr_for_mag'),
    (('photometry', 'max_recenter_shift'), 'max_recenter_shift'),
    (('photometry', 'centroid_outlier_px'), 'centroid_outlier_px'),
    (('photometry', 'registration', 'match_radius_px'), 'registration_match_radius_px'),
    (('photometry', 'registration', 'min_anchors'), 'registration_min_anchors'),
    (('photometry', 'sky_sigma_mode'), 'sky_sigma_mode'),
    (('photometry', 'sky_sigma_includes_rn'), 'sky_sigma_includes_rn'),
    (('photometry', 'sky_sigma_min_n_sky'), 'sky_sigma_min_n_sky'),
    (('photometry', 'use_qc_pass_only'), 'phot_use_qc_pass_only'),
    (('photometry', 'ref_require_apcorr_candidate'), 'phot_ref_require_apcorr_candidate'),
    (('photometry', 'ref_apcorr_min_keep'), 'phot_ref_apcorr_min_keep'),
    (('photometry', 'scales', 'aperture_scale'), 'phot_aperture_scale'),
    (('photometry', 'scales', 'annulus_scale'), 'fitsky_annulus_scale'),
    (('photometry', 'scales', 'dannulus_scale'), 'fitsky_dannulus_scale'),
    (('photometry', 'scales', 'center_cbox_scale'), 'center_cbox_scale'),
    (('photometry', 'scales', 'annulus_min_gap_px'), 'annulus_min_gap_px'),
    (('photometry', 'scales', 'annulus_min_width_px'), 'annulus_min_width_px'),
    (('photometry', 'radii', 'min_r_ap_px'), 'min_r_ap_px'),
    (('photometry', 'radii', 'min_r_in_px'), 'min_r_in_px'),
    (('photometry', 'radii', 'min_r_out_px'), 'min_r_out_px'),
    (('photometry', 'radii', 'sigma_clip'), 'annulus_sigma_clip'),
    (('photometry', 'radii', 'max_iter'), 'fitsky_max_iter'),
    (('photometry', 'radii', 'neighbor_mask_scale'), 'annulus_neighbor_mask_scale'),
    (('photometry', 'apcorr', 'apply'), 'apcorr_apply'),
    (('photometry', 'apcorr', 'small_scale'), 'apcorr_small_scale'),
    (('photometry', 'apcorr', 'large_scale'), 'apcorr_large_scale'),
    (('photometry', 'apcorr', 'min_n'), 'apcorr_use_min_n'),
    (('photometry', 'apcorr', 'scatter_max'), 'apcorr_scatter_max'),
    (('photometry', 'apcorr', 'optimize_scales'), 'apcorr_optimize_scales'),
    (('photometry', 'apcorr', 'small_scale_min'), 'apcorr_small_scale_min'),
    (('photometry', 'apcorr', 'small_scale_max'), 'apcorr_small_scale_max'),
    (('photometry', 'apcorr', 'large_scale_min'), 'apcorr_large_scale_min'),
    (('photometry', 'apcorr', 'large_scale_max'), 'apcorr_large_scale_max'),
    (('photometry', 'apcorr', 'scale_step'), 'apcorr_scale_step'),
    (('photometry', 'apcorr', 'min_gap_fwhm'), 'apcorr_min_gap_fwhm'),
    (('photometry', 'apcorr', 'max_pairs'), 'apcorr_max_pairs'),
    (('photometry', 'apcorr', 'max_sources'), 'apcorr_max_sources'),
    (('photometry', 'apcorr', 'scale_min'), 'apcorr_scale_min'),
    (('photometry', 'apcorr', 'scale_max'), 'apcorr_scale_max'),
    (('photometry', 'apcorr', 'min_snr'), 'apcorr_min_snr'),
    (('photometry', 'apcorr', 'isolation_factor'), 'apcorr_isolation_factor'),
    (('source_quality', 'fwhm_ratio_lo'), 'source_quality_fwhm_ratio_lo'),
    (('source_quality', 'fwhm_ratio_hi'), 'source_quality_fwhm_ratio_hi'),
    (('source_quality', 'anchor_neighbor_fwhm_mult'), 'source_quality_anchor_neighbor_fwhm_mult'),
    (('source_quality', 'anchor_flux_pct'), 'source_quality_anchor_flux_pct'),
    (('source_quality', 'apcorr_flux_pct'), 'source_quality_apcorr_flux_pct'),
    (('source_quality', 'psf_seed_flux_pct'), 'source_quality_psf_seed_flux_pct'),
    (('source_quality', 'edge_fwhm_mult'), 'source_quality_edge_fwhm_mult'),
    (('wcs', 'astap_exe'), 'astap_exe'),
    (('wcs', 'timeout_s'), 'astap_timeout_s'),
    (('wcs', 'astap_search_radius_deg'), 'astap_search_radius_deg'),
    (('wcs', 'astap_database'), 'astap_database'),
    (('wcs', 'astap_fov_fudge'), 'astap_fov_fudge'),
    (('wcs', 'astap_downsample'), 'astap_downsample_z'),
    (('wcs', 'astap_max_stars'), 'astap_max_stars_s'),
    (('wcs', 'astnet_local_enable'), 'astnet_local_enable'),
    (('wcs', 'astnet_local_use_wsl'), 'astnet_local_use_wsl'),
    (('wcs', 'astnet_local_command'), 'astnet_local_command'),
    (('wcs', 'astnet_local_timeout_s'), 'astnet_local_timeout_s'),
    (('wcs', 'astnet_local_downsample'), 'astnet_local_downsample'),
    (('wcs', 'astnet_local_scale_low'), 'astnet_local_scale_low'),
    (('wcs', 'astnet_local_scale_high'), 'astnet_local_scale_high'),
    (('wcs', 'astnet_local_radius_deg'), 'astnet_local_radius_deg'),
    (('wcs', 'astnet_local_keep_outputs'), 'astnet_local_keep_outputs'),
    (('wcs', 'astnet_local_use_cache'), 'astnet_local_use_cache'),
    (('wcs', 'astnet_local_max_objs'), 'astnet_local_max_objs'),
    (('wcs', 'astnet_local_cpulimit_s'), 'astnet_local_cpulimit_s'),
    (('wcs', 'astnet_blind_retry_on_fail'), 'astnet_blind_retry_on_fail'),
    (('wcs', 'astnet_blind_cpulimit_s'), 'astnet_blind_cpulimit_s'),
    (('wcs', 'platesolve_gaia_radius_scale'), 'platesolve_gaia_radius_scale'),
    (('wcs', 'max_workers'), 'wcs_max_workers'),
    (('wcs', 'require_qc_pass'), 'wcs_require_qc_pass'),
    (('wcs', 'refine_enable'), 'wcs_refine_enable'),
    (('wcs_qc', 'match_radius_arcsec'), 'wcs_qc_match_radius_arcsec'),
    (('wcs_qc', 'clip_sigma'), 'wcs_qc_clip_sigma'),
    (('wcs_qc', 'require_wcs_ok'), 'wcs_qc_require_wcs_ok'),
    (('wcs_qc', 'min_match_n'), 'wcs_qc_min_match_n'),
    (('wcs_qc', 'min_match_rate'), 'wcs_qc_min_match_rate'),
    (('wcs_qc', 'max_rms_px'), 'wcs_qc_max_rms_px'),
    (('wcs_qc', 'max_p99_px'), 'wcs_qc_max_p99_px'),
    (('wcs_qc', 'min_inlier_rate'), 'wcs_qc_min_inlier_rate'),
    (('wcs_qc', 'max_edge_ratio'), 'wcs_qc_max_edge_ratio'),
    (('wcs_qc', 'max_center_offset_arcsec'), 'wcs_qc_max_center_offset_arcsec'),
    (('wcs_refine', 'enable'), 'wcs_refine_enable'),
    (('wcs_refine', 'max_match'), 'wcs_refine_max_match'),
    (('wcs_refine', 'match_r_fwhm'), 'wcs_refine_match_r_fwhm'),
    (('wcs_refine', 'min_match'), 'wcs_refine_min_match'),
    (('gaia', 'radius_fudge'), 'gaia_radius_fudge'),
    (('gaia', 'mag_max'), 'gaia_mag_max'),
    (('gaia', 'wcs_mag_max'), 'gaia_wcs_mag_max'),
    (('refbuild', 'wcs_match_radius_arcsec'), 'ref_wcs_match_radius_arcsec'),
    (('gaia', 'snr_calib_min'), 'gaia_snr_calib_min'),
    (('gaia', 'gi_min'), 'gaia_gi_min'),
    (('gaia', 'gi_max'), 'gaia_gi_max'),
    (('gaia', 'retry'), 'gaia_retry'),
    (('gaia', 'timeout_s'), 'gaia_timeout_s'),
    (('gaia', 'hard_deadline_s'), 'gaia_hard_deadline_s'),
    (('gaia', 'backoff_s'), 'gaia_backoff_s'),
    (('gaia', 'allow_no_cache'), 'gaia_allow_no_cache'),
    (('gaia', 'g_limit'), 'idmatch_gaia_g_limit'),
    (('idmatch', 'gaia_g_limit'), 'idmatch_gaia_g_limit'),
    (('idmatch', 'match_r_fwhm'), 'idmatch_match_r_fwhm'),
    (('idmatch', 'two_pass_enable'), 'idmatch_two_pass_enable'),
    (('idmatch', 'tight_radius_arcsec'), 'idmatch_tight_radius_arcsec'),
    (('idmatch', 'loose_radius_arcsec'), 'idmatch_loose_radius_arcsec'),
    (('idmatch', 'adaptive_retry_threshold'), 'idmatch_adaptive_retry_threshold'),
    (('idmatch', 'fwhm_adaptive_floor'), 'idmatch_fwhm_adaptive_floor'),
    (('idmatch', 'geom_correction_enable'), 'idmatch_geom_correction_enable'),
    (('idmatch', 'min_correction_pairs'), 'idmatch_min_correction_pairs'),
    (('idmatch', 'min_affine_pairs'), 'idmatch_min_affine_pairs'),
    (('idmatch', 'tol_px'), 'idmatch_tol_px'),
    (('idmatch', 'tol_arcsec'), 'idmatch_tol_arcsec'),
    (('idmatch', 'force'), 'force_idmatch'),
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
    (('master', 'ref_frame'), 'ref_frame'),
    (('master', 'min_frames_xy'), 'master_min_frames_xy'),
    (('master', 'preserve_ids'), 'master_preserve_ids'),
    (('master', 'force_build'), 'force_master_build'),
    (('master_editor', 'search_radius_px'), 'search_radius_px'),
    (('master_editor', 'bulk_drop_box_px'), 'bulk_drop_box_px'),
    (('master_editor', 'gaia_add_max_sep_arcsec'), 'gaia_add_max_sep_arcsec'),
    (('match', 'tol_px'), 'match_tol_px'),
    (('match', 'wcs_radius_arcsec'), 'wcs_match_radius_arcsec'),
    (('match', 'min_gaia_matches'), 'min_master_gaia_matches'),
    (('match', 'pixel_scale_arcsec'), 'pixel_scale_arcsec'),
    (('overlay', 'max_labels'), 'overlay_max_labels'),
    (('overlay', 'label_fontsize'), 'overlay_label_fontsize'),
    (('overlay', 'label_offset_px'), 'overlay_label_offset_px'),
    (('overlay', 'show_id_when_no_mag'), 'overlay_show_id_when_no_mag'),
    (('overlay', 'use_phot_centroid'), 'overlay_use_phot_centroid'),
    (('overlay', 'show_ref_pos'), 'overlay_show_ref_pos'),
    (('overlay', 'show_shift_vectors'), 'overlay_show_shift_vectors'),
    (('overlay', 'shift_max_vectors'), 'overlay_shift_max_vectors'),
    (('overlay', 'shift_min_px'), 'overlay_shift_min_px'),
    (('overlay', 'inspect_index'), 'inspect_index'),
    (('ui', 'canvas_px'), 'ui_canvas_px'),
    (('transform', 'save_src2ref'), 'save_src2ref_tforms'),
    (('detection', 'peak', 'kernel_scales'), 'peak_kernel_scales'),
    (('target', 'name'), 'target_name'),
    (('target', 'ra_deg'), 'target_ra_deg'),
    (('target', 'dec_deg'), 'target_dec_deg'),
    (('simbad', 'timeout_s'), 'simbad_timeout_s'),
)

CMD_ONLY_TOML_KEY_MAP: tuple[tuple[tuple[str, ...], str], ...] = (
    (('detection', 'mode'), 'detect_mode'),
    (('detection', 'deblend', 'mode'), 'deblend_mode'),
    (('wcs', 'astap_annotate_variables'), 'astap_annotate_variables'),
    (('wcs', 'header_coord_warn_sep_deg'), 'wcs_header_coord_warn_sep_deg'),
    (('wcs', 'header_coord_max_sep_deg'), 'wcs_header_coord_max_sep_deg'),
    (('gaia', 'match_tol_arcsec'), 'ref_wcs_match_radius_arcsec'),
    (('gaia', 'derived_enable'), 'gaia_derived_enable'),
    (('gaia', 'pmem_method'), 'gaia_pmem_method'),
    (('gaia', 'pmem_ruwe_max'), 'gaia_pmem_ruwe_max'),
    (('gaia', 'pmem_min_visibility_periods'), 'gaia_pmem_min_visibility_periods'),
    (('gaia', 'pmem_min_valid'), 'gaia_pmem_min_valid'),
    (('gaia', 'pmem_min_fit'), 'gaia_pmem_min_fit'),
    (('idmatch', 'mode'), 'idmatch_mode'),
    (('idmatch', 'use_gaia_refs_only'), 'idmatch_use_gaia_refs_only'),
    (('master_editor', 'membership_overlay_enable'), 'step8_membership_overlay_enable'),
    (('master_editor', 'membership_threshold'), 'step8_membership_threshold'),
    (('cmd', 'membership_mode'), 'cmd_membership_mode'),
    (('cmd', 'membership_compare'), 'cmd_membership_compare'),
    (('instrument', 'gain_e_per_adu'), 'gain_e_per_adu'),
    (('instrument', 'rdnoise_e'), 'rdnoise_e'),
    (('instrument', 'noise_use_fits_header'), 'noise_use_fits_header'),
    (('instrument', 'noise_reference_binning'), 'noise_reference_binning'),
    (('instrument', 'noise_scale_by_binning'), 'noise_scale_by_binning'),
    (('instrument', 'saturation_adu'), 'saturation_adu'),
    (('instrument', 'datamin_adu'), 'datamin_adu'),
    (('instrument', 'datamax_adu'), 'datamax_adu'),
    (('instrument', 'binning'), 'binning_default'),
    (('instrument', 'telescope_focal_mm'), 'telescope_focal_mm'),
    (('instrument', 'camera_pixel_um'), 'camera_pixel_um'),
    (('instrument', 'zp_initial'), 'zp_initial'),
    (('site', 'lat_deg'), 'site_lat_deg'),
    (('site', 'lon_deg'), 'site_lon_deg'),
    (('site', 'alt_m'), 'site_alt_m'),
    (('site', 'tz_offset_hours'), 'site_tz_offset_hours'),
    (('psf', 'epsf_oversampling'), 'psf_epsf_oversampling'),
    (('psf', 'epsf_size_px'), 'psf_epsf_size_px'),
    (('psf', 'epsf_size_fwhm_mult'), 'psf_epsf_size_fwhm_mult'),
    (('psf', 'n_stars_max'), 'psf_n_stars_max'),
    (('psf', 'isolation_fwhm_mult'), 'psf_isolation_fwhm_mult'),
    (('psf', 'flux_percentile_lo'), 'psf_flux_percentile_lo'),
    (('psf', 'flux_percentile_hi'), 'psf_flux_percentile_hi'),
    (('psf', 'fit_shape_px'), 'psf_fit_shape_px'),
    (('psf', 'fit_shape_fwhm_mult'), 'psf_fit_shape_fwhm_mult'),
    (('psf', 'max_iter'), 'psf_max_iter'),
    (('psf', 'redetect_sigma'), 'psf_redetect_sigma'),
    (('psf', 'redetect_sigma_g'), 'psf_redetect_sigma_g'),
    (('psf', 'redetect_sigma_r'), 'psf_redetect_sigma_r'),
    (('psf', 'redetect_sigma_i'), 'psf_redetect_sigma_i'),
    (('psf', 'model_mode'), 'psf_model_mode'),
    (('psf', 'fit_engine'), 'psf_fit_engine'),
    (('psf', 'build_mode'), 'psf_build_mode'),
    (('psf', 'field_mode'), 'psf_mode'),
    (('psf', 'parallel_workers'), 'psf_parallel_workers'),
    (('psf', 'duplicate_radius_fwhm_mult'), 'psf_duplicate_radius_fwhm_mult'),
    (('psf', 'duplicate_radius_px'), 'psf_duplicate_radius_px'),
    (('psf', 'new_sources_cap_per_iter'), 'psf_new_sources_cap_per_iter'),
    (('psf', 'new_sources_cap_frac'), 'psf_new_sources_cap_frac'),
    (('psf', 'fit_init_max_sources'), 'psf_fit_init_max_sources'),
    (('psf', 'core_cut_enable'), 'psf_core_cut_enable'),
    (('psf', 'core_cut_center_mode'), 'psf_core_cut_center_mode'),
    (('psf', 'core_cut_x_px'), 'psf_core_cut_x_px'),
    (('psf', 'core_cut_y_px'), 'psf_core_cut_y_px'),
    (('psf', 'core_cut_radius_px'), 'psf_core_cut_radius_px'),
    (('psf', 'core_cut_radius_fwhm_mult'), 'psf_core_cut_radius_fwhm_mult'),
    (('psf', 'core_cut_auto_cell_fwhm_mult'), 'psf_core_cut_auto_cell_fwhm_mult'),
    (('psf', 'core_cut_auto_min_density_ratio'), 'psf_core_cut_auto_min_density_ratio'),
    (('psf', 'core_cut_auto_min_sources'), 'psf_core_cut_auto_min_sources'),
    (('psf', 'core_cut_max_exclude_frac'), 'psf_core_cut_max_exclude_frac'),
    (('psf', 'substar_neighbor_r_fwhm_mult'), 'psf_substar_neighbor_r_fwhm_mult'),
    (('psf', 'substar_max_sources'), 'psf_substar_max_sources'),
    (('psf', 'conv_new_frac'), 'psf_conv_new_frac'),
    (('psf', 'use_grouper'), 'psf_use_grouper'),
    (('psf', 'redetect_sharp_lo'), 'psf_redetect_sharp_lo'),
    (('psf', 'redetect_sharp_hi'), 'psf_redetect_sharp_hi'),
    (('psf', 'redetect_round_abs_max'), 'psf_redetect_round_abs_max'),
    (('psf', 'flux_conv_threshold'), 'psf_flux_conv_threshold'),
    (('psf', 'use_error_image'), 'psf_use_error_image'),
    (('psf', 'save_residuals'), 'psf_save_residuals'),
    (('psf', 'save_model_image'), 'psf_save_model_image'),
    (('cross_frame', 'ransac_tol_px'), 'cross_frame_ransac_tol_px'),
    (('cross_frame', 'ransac_max_iter'), 'cross_frame_ransac_max_iter'),
    (('cross_frame', 'ransac_min_inliers'), 'cross_frame_ransac_min_inliers'),
    (('cross_frame', 'match_tol_px'), 'cross_frame_match_tol_px'),
    # ── CMD-only keys migrated from COMMON_TOML_KEY_MAP ──────────────────
    (('cmd', 'snr_calib_min'), 'cmd_snr_calib_min'),
    (('cmd', 'max_sources'), 'cmd_max_sources'),
    (('cmd', 'apply_extinction'), 'cmd_apply_extinction'),
    (('cmd', 'extinction_mode'), 'cmd_extinction_mode'),
    (('cmd', 'frame_zp_min_n'), 'frame_zp_min_n'),
    (('cmd', 'zp', 'clip_sigma'), 'zp_clip_sigma'),
    (('cmd', 'zp', 'fit_iters'), 'zp_fit_iters'),
    (('cmd', 'zp', 'slope_absmax'), 'zp_slope_absmax'),
    (('cmd', 'color', 'clip_sigma'), 'color_clip_sigma'),
    (('cmd', 'color', 'fit_iters'), 'color_fit_iters'),
    (('cmd', 'color', 'slope_absmax'), 'color_slope_absmax'),
    (('gaia', 'zp_slope_absmax'), 'gaia_zp_slope_absmax'),
    (('gaia', 'color_slope_absmax'), 'gaia_color_slope_absmax'),
    (('isochrone', 'file_path'), 'iso_file_path'),
    (('isochrone', 'age_init'), 'iso_age_init'),
    (('isochrone', 'mh_init'), 'iso_mh_init'),
    (('isochrone', 'eg_r_init'), 'iso_eg_r_init'),
    (('isochrone', 'dm_init'), 'iso_dm_init'),
)

LC_ONLY_TOML_KEY_MAP: tuple[tuple[tuple[str, ...], str], ...] = (
    (('io', 'night_parse_mode'), 'night_parse_mode'),
    (('io', 'night_parse_regex'), 'night_parse_regex'),
    (('io', 'night_parse_split_delim'), 'night_parse_split_delim'),
    (('io', 'night_parse_split_index'), 'night_parse_split_index'),
    (('io', 'night_parse_last_digits'), 'night_parse_last_digits'),
    (('io', 'night_parse_include_unmatched'), 'night_parse_include_unmatched'),
    (('io', 'airmass_formula'), 'airmass_formula'),
    (('io', 'airmass_update_header'), 'airmass_update_header'),
    (('io', 'airmass_update_mode'), 'airmass_update_mode'),
    (('instrument', 'gain_e_per_adu'), 'gain_e_per_adu'),
    (('instrument', 'rdnoise_e'), 'rdnoise_e'),
    (('instrument', 'noise_use_fits_header'), 'noise_use_fits_header'),
    (('instrument', 'noise_reference_binning'), 'noise_reference_binning'),
    (('instrument', 'noise_scale_by_binning'), 'noise_scale_by_binning'),
    (('instrument', 'saturation_adu'), 'saturation_adu'),
    (('instrument', 'binning'), 'binning_default'),
    (('instrument', 'telescope_focal_mm'), 'telescope_focal_mm'),
    (('instrument', 'camera_pixel_um'), 'camera_pixel_um'),
    (('instrument', 'zp_initial'), 'zp_initial'),
    (('site', 'lat_deg'), 'site_lat_deg'),
    (('site', 'lon_deg'), 'site_lon_deg'),
    (('site', 'alt_m'), 'site_alt_m'),
    (('site', 'tz_offset_hours'), 'site_tz_offset_hours'),
    (('photometry', 'use_original_frames'), 'phot_use_original_frames'),
    (('wcs', 'propagate_max_shift_px'), 'wcs_propagate_max_shift_px'),
    (('wcs', 'propagate_min_match'), 'wcs_propagate_min_match'),
    (('wcs', 'propagate_sigma_clip'), 'wcs_propagate_sigma_clip'),
    (('idmatch', 'init_r_fwhm'), 'idmatch_init_r_fwhm'),
    (('idmatch', 'ratio_max'), 'idmatch_ratio_max'),
    (('idmatch', 'min_pairs'), 'idmatch_min_pairs'),
    (('idmatch', 'transform_mode'), 'idmatch_transform_mode'),
    (('idmatch', 'mutual_nearest'), 'idmatch_mutual_nearest'),
    (('refbuild', 'sat_drop_pct'), 'ref_select_sat_pct'),
    (('refbuild', 'elong_drop_pct'), 'ref_select_elong_pct'),
    (('refbuild', 'per_date'), 'ref_per_date'),
    (('refbuild', 'master_union'), 'ref_master_union'),
    (('refbuild', 'union_min_frames'), 'ref_union_min_frames'),
    (('refbuild', 'ref_cat_max_sources'), 'ref_cat_max_sources'),
    (('refbuild', 'ref_cat_min_sources'), 'ref_cat_min_sources'),
    (('refbuild', 'ref_cat_max_elong'), 'ref_cat_max_elong'),
    (('refbuild', 'ref_cat_max_abs_round'), 'ref_cat_max_abs_round'),
    (('refbuild', 'ref_cat_sharp_min'), 'ref_cat_sharp_min'),
    (('refbuild', 'ref_cat_sharp_max'), 'ref_cat_sharp_max'),
    (('refbuild', 'ref_cat_min_peak_adu'), 'ref_cat_min_peak_adu'),
    (('refbuild', 'wcs_min_match_rate'), 'ref_wcs_min_match_rate'),
    (('refbuild', 'wcs_min_match_n'), 'ref_wcs_min_match_n'),
    (('refbuild', 'wcs_max_sep_med_arcsec'), 'ref_wcs_max_sep_med_arcsec'),
    (('refbuild', 'wcs_max_sep_p90_arcsec'), 'ref_wcs_max_sep_p90_arcsec'),
    (('refbuild', 'wcs_max_dup_rate'), 'ref_wcs_max_dup_rate'),
    (('light_curve', 'color_index_by_filter'), 'lightcurve_color_index_by_filter'),
    (('light_curve', 'color_term_by_filter'), 'lightcurve_color_term_by_filter'),
    (('extinction_fit', 'order'), 'extfit_order'),
    (('extinction_fit', 'min_points'), 'extfit_min_points'),
    (('extinction_fit', 'clip_sigma'), 'extfit_clip_sigma'),
    (('extinction_fit', 'fit_iters'), 'extfit_fit_iters'),
    (('extinction_fit', 'use_color_terms'), 'extfit_use_color_terms'),
    (('extinction_fit', 'color_index_by_filter'), 'extfit_color_index_by_filter'),
    (('extinction_fit', 'color_c1_by_filter'), 'extfit_color_c1_by_filter'),
    (('extinction_fit', 'color_c2_by_filter'), 'extfit_color_c2_by_filter'),
    # Keys read by the shared Extinction (Airmass Fit) tool — available in LC mode
    # too (registry modes default to cmd+lc). These were over-eagerly moved to
    # CMD_ONLY; restore them for LC so parameters.toml values are honored there.
    (('cmd', 'snr_calib_min'), 'cmd_snr_calib_min'),
    (('cmd', 'zp', 'clip_sigma'), 'zp_clip_sigma'),
    (('cmd', 'zp', 'fit_iters'), 'zp_fit_iters'),
)

# ── Canonical composed maps — callers must use these ─────────────────────────
# (The expanded hand-written literals that lived here were removed in favour of
#  the authoritative composition below.  Use git blame to recover them.)
CMD_TOML_KEY_MAP = COMMON_TOML_KEY_MAP + CMD_ONLY_TOML_KEY_MAP
LC_TOML_KEY_MAP = COMMON_TOML_KEY_MAP + LC_ONLY_TOML_KEY_MAP


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
    mapping: Iterable[tuple[Iterable[str], str]],
) -> dict[str, list[str]]:
    """Return TOML paths mapped more than once."""

    paths: dict[str, list[str]] = {}
    for path, attr in mapping:
        paths.setdefault(path_to_dotted(path), []).append(attr)
    return {path: attrs for path, attrs in paths.items() if len(attrs) > 1}


def duplicate_runtime_attrs(
    mapping: Iterable[tuple[Iterable[str], str]],
) -> dict[str, list[str]]:
    """Return runtime attrs intentionally or accidentally mapped from multiple paths."""

    attrs: dict[str, list[str]] = {}
    for path, attr in mapping:
        attrs.setdefault(attr, []).append(path_to_dotted(path))
    return {attr: sorted(paths) for attr, paths in attrs.items() if len(paths) > 1}


def toml_key_map_for_mode(mode: str) -> tuple[tuple[tuple[str, ...], str], ...]:
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
