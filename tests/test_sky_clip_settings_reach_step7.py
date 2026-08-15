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
# a real disconnect, but of three different kinds, and only the first was safe
# to fix without changing what the pipeline computes. Adding to this list must
# be a deliberate act with the reason written down.
KNOWN_ABSENT = {
    # Set by the headless runner at run time, never by the config.
    "force_redetect": "runtime",
    # Defined in parameters_lc but not parameters_cmd — the same cut exists in
    # one mode and silently falls back to a literal in the other.
    "ref_cat_max_elong": "lc-only",
    "ref_cat_max_abs_round": "lc-only",
    "ref_cat_sharp_min": "lc-only",
    "ref_cat_sharp_max": "lc-only",
    # A two-sided disconnect: the code reads `forced_*_scale`, which nothing
    # sets, while the config offers `photometry.apcorr.small_scale` /
    # `.large_scale` (-> apcorr_small_scale / apcorr_large_scale), which nothing
    # reads. So the aperture radii are settable in the config and the setting
    # does nothing. Wiring them is not behaviour-neutral — apcorr_large_scale is
    # 3.0 in every workspace while Step 7 uses 2.4 — so it needs a decision, not
    # a quiet fix (2026-08-16).
    "forced_r_ap_scale": "two-sided",
    "forced_ref_ap_scale": "two-sided",
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


def test_the_aperture_scale_disconnect_is_still_there_and_still_two_sided():
    """Pinned, not fixed. If someone wires it, this fails and they must say so.

    The config offers the aperture scales and nothing reads them; the code has
    knobs for the aperture scales and nothing sets them. Connecting the two
    changes r_ref from 2.4 to 3.0 x FWHM, which moves every growth curve.
    """
    source = (REPO / "apex/analysis/forced_photometry.py").read_text(encoding="utf-8")
    assert 'getattr(P, "forced_r_ap_scale"' in source
    config_source = (REPO / "apex/config/parameters_cmd.py").read_text(encoding="utf-8")
    assert "apcorr_small_scale" in config_source
    reading = [path for path in (REPO / "apex").rglob("*.py")
               if "config" not in path.parts
               and "apcorr_small_scale" in path.read_text(encoding="utf-8", errors="replace")]
    assert not reading, f"이제 읽는 곳이 있다 — 연결했으면 문서를 고칠 것: {reading}"


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
