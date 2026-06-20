"""Qt-free tests for the unified validation harness.

These run a deliberately small artificial-star experiment so the whole suite
stays well under a minute. They assert that the report + manifest are produced
and that the completeness / RMS / position metrics are present and finite.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from apex.benchmark.synthetic_frame import make_synthetic_reference_frame
from apex.benchmark.validate import (
    run_artificial_star_suite,
    run_iraf_suite,
    run_known_targets_suite,
    run_validation,
    write_report,
)


# Small but realistic parameters so the run is fast yet non-degenerate.
_FAST_SYNTH = {"size": 384, "n_stars": 90}
_FAST_KW = dict(
    seed=20260620,
    trials=2,
    stars_per_trial=20,
    bootstrap_samples=20,
    synth_overrides=_FAST_SYNTH,
)


def _is_finite_number(value) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


@pytest.fixture(scope="module")
def artificial_star_result(tmp_path_factory):
    out = tmp_path_factory.mktemp("artificial_star")
    return run_artificial_star_suite(out, config=None, **_FAST_KW), out


def test_synthetic_frame_generation(tmp_path):
    from astropy.io import fits

    path = make_synthetic_reference_frame(
        tmp_path / "frame.fits", seed=3, size=320, n_stars=60
    )
    assert path.exists()
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    assert data.shape == (320, 320)
    assert header["FILTER"]
    assert float(header["EGAIN"]) > 0
    assert float(header["EXPTIME"]) > 0


def test_artificial_star_metrics_present_and_finite(artificial_star_result):
    result, _ = artificial_star_result
    assert result["status"] == "ok", result.get("reason")
    # Core accuracy metrics must be finite numbers.
    for key in (
        "completeness_m50",
        "completeness_m90",
        "position_rmse_px",
        "photometric_bias_median_mag",
        "photometric_scatter_mad_mag",
    ):
        assert _is_finite_number(result[key]), f"{key}={result[key]!r}"
    # False-positive rate is a finite fraction in [0, 1].
    fp = result["false_positive_rate"]
    assert _is_finite_number(fp) and 0.0 <= fp <= 1.0
    # At least one of bright/faint RMS is measured and finite.
    rms_values = [result["photometric_rms_bright_mag"], result["photometric_rms_faint_mag"]]
    assert any(_is_finite_number(v) for v in rms_values)
    # The completeness curve must span a real range (non-degenerate).
    assert result["completeness_non_degenerate"] is True
    # Per-magnitude completeness table is populated.
    assert result["completeness_by_magnitude"]


def test_iraf_suite_skips_cleanly_without_inputs(tmp_path):
    # No config at all -> skipped, never raises.
    result = run_iraf_suite(tmp_path, config=None)
    assert result["status"] == "skipped"
    assert result.get("reason")

    # Config with non-existent inputs -> still skipped, not an error.
    import types

    ns = types.SimpleNamespace(
        iraf={"input_fits": str(tmp_path / "nope.fits"), "step7_tsv": str(tmp_path / "nope.tsv")}
    )
    result2 = run_iraf_suite(tmp_path, config=ns)
    assert result2["status"] == "skipped"


def test_known_targets_synthetic_recovery_runs(tmp_path):
    result = run_known_targets_suite(tmp_path, config=None)
    assert result["status"] == "ok"
    syn = result["synthetic_isochrone"]
    assert syn["status"] == "ok"
    # The recovery should pass against the fitter's tolerances.
    assert result["passed"] is True
    assert syn["parameter_checks"]
    # Real clusters are skipped when none are configured.
    assert result["real_clusters"]["status"] == "skipped"


def test_write_report_produces_files(tmp_path):
    results = {
        "timestamp": "2026-06-20T00:00:00+00:00",
        "git": {"commit": "deadbeef" * 5, "branch": "main", "dirty": False},
        "versions": {"python": "3.11.0"},
        "suites": {
            "artificial-star": {
                "status": "ok",
                "completeness_m50": 17.3,
                "completeness_m90": 16.3,
                "completeness_m10": 18.2,
                "completeness_width_mag": 0.3,
                "completeness_m50_ci95": [17.0, 17.6],
                "completeness_non_degenerate": True,
                "photometric_bias_median_mag": 0.001,
                "photometric_scatter_mad_mag": 0.04,
                "photometric_rms_bright_mag": 0.02,
                "photometric_rms_faint_mag": 0.12,
                "position_rmse_px": 0.12,
                "new_detections": 40,
                "new_false_detections": 0,
                "false_positive_rate": 0.0,
                "n_trials": 2,
                "stars_per_trial": 20,
                "n_injected": 40,
                "injection_mag_range": [14.5, 19.0],
                "fwhm_px": 3.5,
                "zeropoint_mag": 25.0,
                "completeness_by_magnitude": [
                    {"magnitude_center": 15.5, "n_eligible": 10, "completeness": 1.0, "mag_bias_median": 0.0}
                ],
                "plots": [],
            },
            "iraf": {"status": "skipped", "reason": "no IRAF"},
            "known-targets": {
                "status": "ok",
                "passed": True,
                "synthetic_isochrone": {"status": "ok", "passed": True, "n_stars": 80, "reduced_chi2": 1.1, "parameter_checks": {}},
                "real_clusters": {"status": "skipped", "reason": "none configured"},
            },
        },
    }
    paths = write_report(results, tmp_path, update_docs=False)
    assert paths["summary"].exists()
    assert paths["manifest"].exists()
    assert paths["docs"].exists()

    text = paths["summary"].read_text(encoding="utf-8")
    assert "# APEX Validation Report" in text
    assert "Interpretation & caveats" in text
    assert "m50" in text

    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    assert "suites" in manifest
    assert manifest["suites"]["artificial-star"]["completeness_m50"] == 17.3


def test_run_validation_rejects_unknown_suite(tmp_path):
    with pytest.raises(ValueError):
        run_validation("not-a-suite", tmp_path, config=None)
