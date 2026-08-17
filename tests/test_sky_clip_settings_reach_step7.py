"""The sky-annulus clipping settings have to reach the code that clips.

`photometry.radii.sigma_clip` and `.max_iter` are real config keys that load
into `annulus_sigma_clip` / `fitsky_max_iter`. Step 7's forced photometry read
`phot_sigma_clip` / `phot_max_iter` instead — names nothing ever sets — so
`getattr(..., default)` silently returned the literal every time. Step 4's
detection read the right names, so the same config key changed one half of the
pipeline and not the other, with nothing in the log to say so.

This is the third instance of the shape: a value the config accepts, a run that
succeeds, and a setting that does nothing (grouper knobs 2026-08-14, mode
presets 2026-08-15, this one 2026-08-16). The check is cheap — every attribute
a module reads off `P` with a default should be an attribute `P` actually has.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from apex.config.parameters_cmd import read_params

REPO = Path(__file__).absolute().parents[1]
READS_A_DEFAULT = re.compile(r'getattr\(\s*[pP]\s*,\s*"([a-z_][a-z0-9_]*)"\s*,')

# Modules whose settings must exist. Not the whole tree: GUI widgets legitimately
# probe for attributes that only some modes define.
WATCHED = [
    "apex/analysis/forced_photometry.py",
    "apex/analysis/detection.py",
    "apex/benchmark/runner.py",
]

# Names Step 7 reads off `P` that `parameters_cmd` does not define. Each one is
# a real disconnect. Adding to this list must be a deliberate act with the
# reason written down.
KNOWN_ABSENT = {
    # Set by the headless runner at run time, never by the config.
    "force_redetect": "runtime",
    # Defined in parameters_lc but not parameters_cmd — the same cut exists in
    # one mode and silently falls back to a literal in the other.
    "ref_cat_max_elong": "lc-only",
    "ref_cat_max_abs_round": "lc-only",
    "ref_cat_sharp_min": "lc-only",
    "ref_cat_sharp_max": "lc-only",
}


@pytest.fixture(scope="module")
def parameters(tmp_path_factory):
    config = tmp_path_factory.mktemp("ws") / "apex_config.json"
    config.write_text(json.dumps({"io": {"result_dir": ".", "data_dir": "."}}),
                      encoding="utf-8")
    return read_params(config).P


def test_the_sky_annulus_settings_are_the_ones_that_exist(parameters):
    """The pair that was orphaned and is now wired."""
    assert hasattr(parameters, "annulus_sigma_clip")
    assert hasattr(parameters, "fitsky_max_iter")
    assert not hasattr(parameters, "phot_sigma_clip")
    assert not hasattr(parameters, "phot_max_iter")


def test_step7_reads_the_settings_that_exist():
    source = (REPO / "apex/analysis/forced_photometry.py").read_text(encoding="utf-8")
    assert 'getattr(P, "annulus_sigma_clip"' in source
    assert 'getattr(P, "fitsky_max_iter"' in source
    assert "phot_sigma_clip" not in source
    assert "phot_max_iter" not in source


def test_the_same_config_key_reaches_detection_and_forced_photometry(tmp_path):
    """One key, both halves of the pipeline — that was the actual breakage."""
    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({
        "io": {"result_dir": ".", "data_dir": "."},
        "photometry": {"radii": {"sigma_clip": 2.5, "max_iter": 9}},
    }), encoding="utf-8")
    P = read_params(config).P
    assert P.annulus_sigma_clip == pytest.approx(2.5)
    assert P.fitsky_max_iter == 9


def test_the_aperture_scales_come_from_the_config_now(parameters):
    """Was a two-sided break: the code read `forced_*_scale`, which nothing set,
    while the config's `photometry.apcorr.small_scale` / `.large_scale` reached
    `raw` and were dropped by the namespace constructor. Both halves worked and
    were never connected, so fifty config files on disk carried an aperture
    radius no run had used.

    Wiring it was made behaviour-neutral first — 0.8 / 2.4 everywhere, the radii
    every stored product was measured with — and raised to 1.0 / 3.0 on
    2026-08-17, which is what the configs had been asking for all along.

    0.8 was never chosen. It looked deliberate because it landed exactly on
    `small_scale_min`, but nothing reads that key, or `small_scale_max`,
    `large_scale_min/max`, `scale_step`, `scale_min/max`, or `optimize_scales`:
    the config advertises an aperture optimiser that does not exist. What made
    the value cost something is that the enclosed-flux curve is 2.4x steeper at
    0.8 than at 1.0, and an aperture correction is one scalar per frame — it
    removes the mean offset, not the part that varies star to star.
    """
    assert parameters.apcorr_small_scale == pytest.approx(1.0)
    assert parameters.apcorr_large_scale == pytest.approx(3.0)

    source = (REPO / "apex/analysis/forced_photometry.py").read_text(encoding="utf-8")
    assert 'getattr(P, "apcorr_small_scale"' in source
    assert 'getattr(P, "apcorr_large_scale"' in source
    assert "forced_r_ap_scale" not in source
    assert "forced_ref_ap_scale" not in source


def test_both_modes_share_one_aperture(tmp_path):
    """Step 7 is shared, so CMD and LC must not disagree about its radii.

    `parameters_lc` declared 1.0 / 3.0 while `parameters_cmd` had no field at
    all — the moment the read was wired up, the same engine would have measured
    two different apertures depending on which mode opened the workspace.
    """
    from apex.config.parameters_lc import read_params as read_lc

    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({"io": {"result_dir": ".", "data_dir": "."}}),
                      encoding="utf-8")
    cmd_P, lc_P = read_params(config).P, read_lc(config).P
    assert cmd_P.apcorr_small_scale == lc_P.apcorr_small_scale
    assert cmd_P.apcorr_large_scale == lc_P.apcorr_large_scale


def test_a_configured_aperture_actually_arrives(tmp_path):
    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({
        "io": {"result_dir": ".", "data_dir": "."},
        "photometry": {"apcorr": {"small_scale": 1.1, "large_scale": 3.3}},
    }), encoding="utf-8")
    P = read_params(config).P
    assert P.apcorr_small_scale == pytest.approx(1.1)
    assert P.apcorr_large_scale == pytest.approx(3.3)


@pytest.mark.parametrize("module", WATCHED)
def test_no_new_setting_goes_missing(module, parameters):
    """The general guard: a name not on the documented list must exist."""
    source = (REPO / module).read_text(encoding="utf-8")
    missing = sorted({
        name for name in READS_A_DEFAULT.findall(source)
        if name not in KNOWN_ABSENT and not hasattr(parameters, name)
    })
    assert not missing, (f"{module} 이 없는 설정을 읽는다: {missing} — "
                         f"진짜 고아인지 확인하고 KNOWN_ABSENT 에 이유와 함께 넣을 것")


def test_every_declaration_of_the_radii_agrees():
    """One number, several places that each carry their own literal.

    The map row is the authority, but the two loaders set the attribute
    explicitly (overriding the map) and the benchmark harness mirrors the
    pipeline so a cross-check compares like with like. A benchmark left behind
    would quietly measure APEX at one radius against IRAF at another and report
    the difference as a result.
    """
    import re

    sources = {
        "apex/config/parameter_map.py": [
            (r"'apcorr_small_scale', 'float', ([\d.]+)", "1.0"),
            (r"'apcorr_large_scale', 'float', ([\d.]+)", "3.0"),
        ],
        "apex/analysis/forced_photometry.py": [
            (r'getattr\(P, "apcorr_small_scale", ([\d.]+)\)', "1.0"),
            (r'getattr\(P, "apcorr_large_scale", ([\d.]+)\)', "3.0"),
        ],
        "apex/benchmark/photometry_crosscheck.py": [
            (r"DEFAULT_APERTURE_SCALE = ([\d.]+)", "1.0"),
            (r"DEFAULT_REF_AP_SCALE = ([\d.]+)", "3.0"),
        ],
        "apex/benchmark/iraf_crosscheck.py": [
            (r"aperture_scale_fwhm: float = ([\d.]+)", "1.0"),
        ],
    }
    for path, checks in sources.items():
        text = (REPO / path).read_text(encoding="utf-8")
        for pattern, expected in checks:
            found = re.search(pattern, text)
            assert found, f"{path} 에서 {pattern} 를 못 찾았다"
            assert found.group(1) == expected, (
                f"{path}: {found.group(1)} — 나머지와 어긋난다 (기대 {expected})")


def test_the_shipped_template_matches_the_code():
    """Otherwise only newly created workspaces get a different aperture."""
    import json

    example = json.loads((REPO / "parameters.example.json").read_text(encoding="utf-8"))
    apcorr = example["photometry"]["apcorr"]
    assert apcorr["small_scale"] == pytest.approx(1.0)
    assert apcorr["large_scale"] == pytest.approx(3.0)
