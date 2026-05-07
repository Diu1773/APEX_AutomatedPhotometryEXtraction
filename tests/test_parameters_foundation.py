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


def _write_minimal_toml(tmp_path):
    data_dir = (tmp_path / "data").as_posix()
    result_dir = (tmp_path / "result").as_posix()
    path = tmp_path / "parameters.toml"
    path.write_text(
        textwrap.dedent(
            f"""
            schema_version = {CANONICAL_SCHEMA_VERSION}

            [io]
            data_dir = "{data_dir}"
            filename_prefix = "pp_"
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


def test_save_toml_writes_canonical_schema_version(tmp_path):
    pytest.importorskip("tomli_w")
    param_path = _write_minimal_toml(tmp_path)
    params = CmdParameters(param_path)

    assert params.save_toml()

    with param_path.open("rb") as fh:
        data = tomllib.load(fh)
    assert data["schema_version"] == CANONICAL_SCHEMA_VERSION
    assert data["detection"]["engine"] == "sep"


def test_pydantic_schema_accepts_schema_version(tmp_path):
    pytest.importorskip("pydantic")
    from apex.config.schema import Parameters

    param_path = _write_minimal_toml(tmp_path)
    params = Parameters.from_toml(param_path)

    assert params.schema_version == CANONICAL_SCHEMA_VERSION
    assert params.detection.engine.value == "sep"
