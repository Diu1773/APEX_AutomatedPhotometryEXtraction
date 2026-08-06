from __future__ import annotations

import textwrap

import pytest

try:
    import tomllib  # type: ignore
except Exception:
    import tomli as tomllib  # type: ignore

from apex.config.parameter_map import (
    CANONICAL_SCHEMA_VERSION,
    CMD_ONLY_TOML_KEY_MAP,
    CMD_TOML_KEY_MAP,
    COMMON_TOML_KEY_MAP,
    LC_ONLY_TOML_KEY_MAP,
    LC_TOML_KEY_MAP,
    PSF_TOML_KEY_MAP,
    duplicate_runtime_attrs,
    duplicate_toml_paths,
    ensure_schema_version,
    entries_for_mode,
    get_toml_path,
    read_schema_version,
    set_toml_path,
    toml_key_map_for_mode,
    toml_value_for_runtime_attr,
)
from apex.config import parameters_cmd, parameters_lc
from apex.config.parameters_cmd import Parameters as CmdParameters
from apex.config.parameters_lc import Parameters as LcParameters


def _write_minimal_toml(tmp_path, *, include_filename_prefix=True):
    data_dir = (tmp_path / "data").as_posix()
    result_dir = (tmp_path / "result").as_posix()
    path = tmp_path / "parameters.toml"
    prefix_line = 'filename_prefix = "pp_"' if include_filename_prefix else ""
    path.write_text(
        textwrap.dedent(
            f"""
            schema_version = {CANONICAL_SCHEMA_VERSION}

            [io]
            data_dir = "{data_dir}"
            {prefix_line}
            result_dir = "{result_dir}"
            cache_dir = "cache"

            [instrument]
            telescope_focal_mm = 3947.0
            camera_pixel_um = 3.76
            binning = 2
            gain_e_per_adu = 0.1
            rdnoise_e = 1.39
            saturation_adu = 65000.0

            [detection]
            engine = "sep"
            sigma = 3.2
            minarea_pix = 3

            [fwhm]
            guess_arcsec = 2.5
            px_min = 3.0
            px_max = 10.0
            measure_max = 25
            min_sources = 15
            candidate_max = 200
            measure_all_sources = false

            [background]
            in_detect = true
            box = 64

            [photometry.registration]
            match_radius_px = 6.5
            min_anchors = 4

            [photometry.apcorr]
            min_snr = 55.0
            isolation_factor = 2.8

            [source_quality]
            fwhm_ratio_lo = 0.7
            fwhm_ratio_hi = 1.5
            anchor_neighbor_fwhm_mult = 2.1
            anchor_flux_pct = 65.0
            apcorr_flux_pct = 70.0
            psf_seed_flux_pct = 35.0
            edge_fwhm_mult = 1.2

            [psf]
            fit_engine = "apex_iterative"
            n_stars_max = 0
            grouper_max_size = 3
            substar_iters = 1
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )
    return path


def test_schema_version_helpers():
    data = {}

    assert read_schema_version({}) == 0
    assert read_schema_version({"schema_version": "2"}) == 2
    assert read_schema_version({"schema_version": "bad"}) == 0
    assert ensure_schema_version(data)["schema_version"] == CANONICAL_SCHEMA_VERSION


def test_toml_path_helpers_and_value_conversion(tmp_path):
    data = {}

    set_toml_path(data, ("detection", "peak", "kernel_scales"), [0.9, 1.3])
    assert get_toml_path(data, ("detection", "peak", "kernel_scales")) == [0.9, 1.3]
    assert get_toml_path(data, ("missing", "path")) is None

    assert toml_value_for_runtime_attr("peak_kernel_scales", "0.9, 1.3", ("detection", "peak", "kernel_scales")) == [0.9, 1.3]
    assert toml_value_for_runtime_attr("cache_dir", tmp_path / "cache", ("io", "cache_dir")) == "cache"
    assert toml_value_for_runtime_attr("result_dir", tmp_path / "result", ("io", "result_dir")) == str(tmp_path / "result")


def test_foundation_map_has_common_modes():
    cmd_paths = {entry.dotted_path for entry in entries_for_mode("cmd")}
    lc_paths = {entry.dotted_path for entry in entries_for_mode("lc")}

    assert "detection.engine" in cmd_paths
    assert "detection.engine" in lc_paths
    assert "instrument.rdnoise_e" in cmd_paths


def test_mode_toml_maps_are_centralized():
    assert parameters_cmd.TOML_KEY_MAP is CMD_TOML_KEY_MAP
    assert parameters_lc.TOML_KEY_MAP is LC_TOML_KEY_MAP
    assert toml_key_map_for_mode("cmd") is CMD_TOML_KEY_MAP
    assert toml_key_map_for_mode("lc") is LC_TOML_KEY_MAP
    assert len(COMMON_TOML_KEY_MAP) > len(CMD_ONLY_TOML_KEY_MAP)
    assert len(COMMON_TOML_KEY_MAP) > len(LC_ONLY_TOML_KEY_MAP)
    assert PSF_TOML_KEY_MAP
    assert (("psf", "fit_engine"), "psf_fit_engine") in PSF_TOML_KEY_MAP


def test_parameter_maps_have_no_duplicate_toml_paths_and_known_attr_aliases():
    assert duplicate_toml_paths(CMD_TOML_KEY_MAP) == {}
    assert duplicate_toml_paths(LC_TOML_KEY_MAP) == {}

    assert duplicate_runtime_attrs(CMD_TOML_KEY_MAP) == {
        "idmatch_gaia_g_limit": ["gaia.g_limit", "idmatch.gaia_g_limit"],
        "ref_wcs_match_radius_arcsec": ["gaia.match_tol_arcsec", "refbuild.wcs_match_radius_arcsec"],
        "wcs_refine_enable": ["wcs.refine_enable", "wcs_refine.enable"],
    }
    assert duplicate_runtime_attrs(LC_TOML_KEY_MAP) == {
        "idmatch_gaia_g_limit": ["gaia.g_limit", "idmatch.gaia_g_limit"],
        "wcs_refine_enable": ["wcs.refine_enable", "wcs_refine.enable"],
    }


def test_cmd_and_lc_load_top_level_schema_version(tmp_path):
    param_path = _write_minimal_toml(tmp_path)

    cmd = CmdParameters(param_path)
    lc = LcParameters(param_path)

    assert cmd.P.schema_version == CANONICAL_SCHEMA_VERSION
    assert lc.P.schema_version == CANONICAL_SCHEMA_VERSION
    assert cmd.P.detect_engine == "sep"
    assert lc.P.detect_engine == "sep"
    assert cmd.P.cache_dir.name == "cache"
    assert lc.P.cache_dir.name == "cache"
    assert cmd.P.registration_match_radius_px == 6.5
    assert lc.P.registration_match_radius_px == 6.5
    assert cmd.P.registration_min_anchors == 4
    assert lc.P.registration_min_anchors == 4
    assert cmd.P.apcorr_min_snr == 55.0
    assert lc.P.apcorr_min_snr == 55.0
    assert cmd.P.apcorr_isolation_factor == 2.8
    assert lc.P.apcorr_isolation_factor == 2.8
    assert cmd.P.source_quality_anchor_neighbor_fwhm_mult == 2.1
    assert lc.P.source_quality_anchor_neighbor_fwhm_mult == 2.1
    assert cmd.P.source_quality_apcorr_flux_pct == 70.0
    assert lc.P.source_quality_apcorr_flux_pct == 70.0
    assert cmd.P.psf_fit_engine == "apex_iterative"
    assert lc.P.psf_fit_engine == "apex_iterative"
    assert cmd.P.psf_n_stars_max == 0
    assert lc.P.psf_n_stars_max == 0
    assert cmd.P.psf_epsf_contamination_filter is True
    assert lc.P.psf_epsf_contamination_filter is True
    assert cmd.P.psf_flux_scale_correction is False
    assert lc.P.psf_flux_scale_correction is False
    assert cmd.P.psf_flux_scale_min_snr == 50.0
    assert lc.P.psf_flux_scale_min_stars == 8
    assert cmd.P.psf_flux_scale_min_neighbor_fwhm == 4.0
    assert lc.P.psf_flux_scale_max_scatter_mag == 0.10
    assert cmd.P.psf_grouper_max_size == 3
    assert lc.P.psf_grouper_max_size == 3
    assert cmd.P.psf_grouper_radius_fwhm == 1.5
    assert lc.P.psf_grouper_radius_fwhm == 1.5
    assert cmd.P.psf_forced_match_radius_fwhm == 1.25
    assert lc.P.psf_forced_match_radius_fwhm == 1.25
    assert cmd.P.psf_substar_iters == 1
    assert lc.P.psf_substar_iters == 1
    assert cmd.P.psf_fitter_max_iter == 6
    assert lc.P.psf_fitter_max_iter == 6
    assert cmd.P.psf_postfit_snr_min == 3.0
    assert cmd.P.psf_fit_shape_fwhm_mult == 2.4
    assert lc.P.psf_fit_window_mode == "auto"
    assert cmd.P.psf_fit_encircled_energy == 0.90
    assert cmd.P.psf_postfit_qfit_max == 3.0
    assert lc.P.psf_postfit_reduced_chi2_max == 25.0
    assert lc.P.psf_blend_residual_ratio == 0.3


def test_missing_filename_prefix_defaults_to_all_fits(tmp_path):
    param_path = _write_minimal_toml(tmp_path, include_filename_prefix=False)

    cmd = CmdParameters(param_path)
    lc = LcParameters(param_path)

    assert cmd.P.filename_prefix == ""
    assert lc.P.filename_prefix == ""


def test_missing_manual_noise_values_load_as_header_fallback_candidates(tmp_path):
    param_path = _write_minimal_toml(tmp_path)
    text = param_path.read_text(encoding="utf-8")
    text = text.replace("gain_e_per_adu = 0.1\n", "")
    text = text.replace("rdnoise_e = 1.39\n", "")
    param_path.write_text(text, encoding="utf-8")

    cmd = CmdParameters(param_path)
    lc = LcParameters(param_path)

    assert cmd.P.gain_e_per_adu is None
    assert cmd.P.rdnoise_e is None
    assert lc.P.gain_e_per_adu is None
    assert lc.P.rdnoise_e is None


def test_save_toml_writes_canonical_schema_version(tmp_path):
    param_path = _write_minimal_toml(tmp_path)
    params = CmdParameters(param_path)

    assert params.save_toml()

    with param_path.open("rb") as fh:
        data = tomllib.load(fh)
    assert data["schema_version"] == CANONICAL_SCHEMA_VERSION
    assert data["detection"]["engine"] == "sep"


def _read_authority(param_path):
    """Saved settings land in the workspace JSON authority, not the TOML."""
    import json
    from apex.config.config_io import resolve_config_path
    auth = resolve_config_path(param_path)
    return json.loads(auth.read_text(encoding="utf-8"))


def test_save_toml_removes_blank_manual_noise_values(tmp_path):
    param_path = _write_minimal_toml(tmp_path)
    params = CmdParameters(param_path)

    params.P.gain_e_per_adu = None
    params.P.rdnoise_e = None
    assert params.save_toml()

    data = _read_authority(param_path)
    assert "gain_e_per_adu" not in data["instrument"]
    assert "rdnoise_e" not in data["instrument"]


@pytest.mark.parametrize("ParamsCls", [CmdParameters, LcParameters])
def test_save_toml_removes_blank_target_coordinates(tmp_path, ParamsCls):
    param_path = _write_minimal_toml(tmp_path)
    param_path.write_text(
        param_path.read_text(encoding="utf-8")
        + "\n[target]\nname = \"NGC457\"\nra_deg = 298.30625\ndec_deg = 58.27806\n",
        encoding="utf-8",
    )
    params = ParamsCls(param_path)

    assert params.P.target_name == "NGC457"
    params.P.target_name = "M37"
    params.P.target_ra_deg = None
    params.P.target_dec_deg = None
    assert params.save_toml()

    data = _read_authority(param_path)
    assert data["target"]["name"] == "M37"
    assert "ra_deg" not in data["target"]
    assert "dec_deg" not in data["target"]


def test_save_toml_preserves_forced_phot_quality_knobs(tmp_path):
    param_path = _write_minimal_toml(tmp_path)
    params = CmdParameters(param_path)

    assert params.save_toml()

    with param_path.open("rb") as fh:
        data = tomllib.load(fh)
    assert data["photometry"]["registration"]["match_radius_px"] == 6.5
    assert data["photometry"]["registration"]["min_anchors"] == 4
    assert data["photometry"]["apcorr"]["min_snr"] == 55.0
    assert data["photometry"]["apcorr"]["isolation_factor"] == 2.8
    assert data["source_quality"]["anchor_neighbor_fwhm_mult"] == 2.1
    assert data["source_quality"]["apcorr_flux_pct"] == 70.0


def test_pydantic_schema_accepts_schema_version(tmp_path):
    pytest.importorskip("pydantic")
    from apex.config.schema import Parameters

    param_path = _write_minimal_toml(tmp_path)
    params = Parameters.from_toml(param_path)

    assert params.schema_version == CANONICAL_SCHEMA_VERSION
    assert params.detection.engine.value == "sep"
