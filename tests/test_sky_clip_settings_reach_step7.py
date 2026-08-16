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

    Wiring it was made behaviour-neutral on purpose: the defaults here and the
    values written into every workspace are 0.8 / 2.4, the radii every stored
    product was measured with. Raising them to the 1.0 / 3.0 the template used
    to claim is a real change and needs a reprocess.
    """
    assert parameters.apcorr_small_scale == pytest.approx(0.8)
    assert parameters.apcorr_large_scale == pytest.approx(2.4)

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
