"""[calibration] TOML wiring (P3).

Before this the section was documented in parameters.example.toml but nothing
parsed it, so editing it had no effect and the GUI panel reset every restart.
"""

from __future__ import annotations

import pytest

from apex.analysis.calibration import CalibrationOptions
from apex.config.calibration_section import (
    calibration_toml_sections,
    read_calibration_section,
)
from apex.utils.param_file import update_param_file

TOML = """
[io]
data_dir = "."

[calibration]
enabled = true
combine_method = "sigmaclip_mean"
maxiters = 7
dark_scale = false
temp_match_tol_c = 0.2
strict_temp = true

[calibration.overscan]
enable = true
edge = "top"
width = 48
trim = false
"""


def _write(tmp_path, text=TOML):
    path = tmp_path / "parameters.toml"
    path.write_text(text, encoding="utf-8")
    return path


def _parse(path):
    """Read back what update_param_file wrote — the JSON authority when it
    exists (post-TOML-removal), else the original file."""
    import json
    from apex.config.config_io import resolve_config_path
    auth = resolve_config_path(path)
    if auth.exists():
        return json.loads(auth.read_text(encoding="utf-8"))
    import tomllib
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def test_read_calibration_section_flattens_overscan(tmp_path):
    section = read_calibration_section(_parse(_write(tmp_path)))
    assert section["combine_method"] == "sigmaclip_mean"
    assert section["maxiters"] == 7
    assert section["dark_scale"] is False
    assert section["overscan_enable"] is True
    assert section["overscan_edge"] == "top"
    assert section["overscan_width"] == 48
    assert section["overscan_trim"] is False


def test_read_calibration_section_absent():
    assert read_calibration_section({}) == {}
    assert read_calibration_section(None) == {}
    assert read_calibration_section({"calibration": "not a table"}) == {}


def test_options_from_mapping_applies_and_defaults(tmp_path):
    opts = CalibrationOptions.from_mapping(read_calibration_section(_parse(_write(tmp_path))))
    assert opts.combine_method == "sigmaclip_mean"
    assert opts.maxiters == 7
    assert opts.dark_scale is False
    assert opts.temp_match_tol_c == pytest.approx(0.2)
    assert opts.strict_temp is True
    assert opts.overscan_edge == "top"
    # untouched keys keep the dataclass defaults; 'enabled' is not a field
    assert opts.flat_min == CalibrationOptions().flat_min
    assert not hasattr(opts, "enabled")


def test_options_from_mapping_survives_bad_values():
    opts = CalibrationOptions.from_mapping(
        {"maxiters": "not a number", "sigma_low": "2.5", "dark_scale": "false"})
    assert opts.maxiters == CalibrationOptions().maxiters    # bad value ignored
    assert opts.sigma_low == pytest.approx(2.5)              # string coerced
    assert opts.dark_scale is False


def test_settings_round_trip_through_the_param_file(tmp_path):
    """What the Step 0 window saves must come back on the next load."""
    path = _write(tmp_path)
    settings = CalibrationOptions(
        combine_method="mean", temp_match_tol_c=0.5, strict_temp=False,
        overscan_enable=True, overscan_width=64,
    ).to_mapping()
    settings.pop("gain")
    settings.pop("readnoise")

    assert update_param_file(path, calibration_toml_sections(settings))

    reloaded = CalibrationOptions.from_mapping(read_calibration_section(_parse(path)))
    assert reloaded.combine_method == "mean"
    assert reloaded.temp_match_tol_c == pytest.approx(0.5)
    assert reloaded.strict_temp is False
    assert reloaded.overscan_enable is True
    assert reloaded.overscan_width == 64
    # unrelated tables must survive the rewrite
    assert _parse(path)["io"]["data_dir"] == "."
    # 'enabled' is not a CalibrationOptions field and must not be dropped
    assert _parse(path)["calibration"]["enabled"] is True


def test_update_param_file_missing_file(tmp_path):
    assert not update_param_file(tmp_path / "nope.toml", {("calibration",): {"maxiters": 3}})


def test_params_expose_the_calibration_table(tmp_path):
    """Both parameter models must surface P.calibration for the GUI/pipeline."""
    from apex.config.parameters_cmd import read_params as read_cmd
    from apex.config.parameters_lc import read_params as read_lc

    path = _write(tmp_path)
    for read_params in (read_cmd, read_lc):
        params = read_params(path)
        opts = CalibrationOptions.from_mapping(params.P.calibration)
        assert opts.combine_method == "sigmaclip_mean"
        assert opts.temp_match_tol_c == pytest.approx(0.2)
        assert opts.overscan_edge == "top"


def test_the_shipped_example_matches_the_dataclass_defaults():
    """The shipped example must not document a value the code no longer uses.

    It was `parameters.example.toml` until 2026-08-18; TOML left the runtime
    that day and the template went with it.
    """
    import json
    from pathlib import Path

    example = Path(__file__).absolute().parents[1] / "parameters.example.json"
    section = read_calibration_section(json.loads(example.read_text(encoding="utf-8")))
    defaults = CalibrationOptions()
    known = set(CalibrationOptions.field_names())
    for key, value in section.items():
        if key not in known:            # e.g. 'enabled'
            continue
        current = getattr(defaults, key)
        if isinstance(current, float):
            assert value == pytest.approx(current), key
        else:
            assert value == current, key
