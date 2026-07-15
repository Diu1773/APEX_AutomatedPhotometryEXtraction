import json

import numpy as np
import pandas as pd

from validation.run_psf_artificial_stars import (
    _copy_baseline_detection,
    _read_background_rms,
    _read_fwhm,
)


def test_read_fwhm_prefers_frame_summary_over_source_median(tmp_path):
    detection_path = tmp_path / "detect_frame.fit.csv"
    pd.DataFrame(
        {
            "fwhm_px": [1.1, 1.4, 8.5, 8.7],
            "fwhm_status": ["out_of_range", "out_of_range", "ok", "ok"],
        }
    ).to_csv(detection_path, index=False)
    detection_path.with_suffix(".json").write_text(
        json.dumps({"fwhm_px": 8.6739}), encoding="utf-8"
    )

    assert _read_fwhm(detection_path, None) == 8.6739


def test_read_fwhm_uses_bounded_sources_without_summary(tmp_path):
    detection_path = tmp_path / "detect_frame.fit.csv"
    pd.DataFrame(
        {
            "fwhm_px": [1.1, 1.4, 8.5, 8.7],
            "fwhm_status": ["out_of_range", "out_of_range", "ok", "ok"],
        }
    ).to_csv(detection_path, index=False)

    assert _read_fwhm(detection_path, None) == 8.6


def test_read_background_rms_prefers_step4_summary(tmp_path):
    detection_path = tmp_path / "detect_frame.fit.csv"
    detection_path.write_text("x,y\n1,1\n", encoding="utf-8")
    detection_path.with_suffix(".json").write_text(
        json.dumps({"bkg_rms": 27.48}), encoding="utf-8"
    )

    value = _read_background_rms(
        detection_path,
        None,
        np.zeros((20, 20), dtype=np.float32),
    )

    assert value == 27.48


def test_copy_baseline_detection_preserves_frame_summary(tmp_path):
    baseline = tmp_path / "baseline"
    source_dir = baseline / "step4_detection"
    source_dir.mkdir(parents=True)
    source = source_dir / "detect_frame.fit.csv"
    source.write_text("x,y\n1,2\n", encoding="utf-8")
    source.with_suffix(".json").write_text(
        json.dumps({"fwhm_px": 8.7}), encoding="utf-8"
    )
    trial_result = tmp_path / "trial" / "result"

    target = _copy_baseline_detection(baseline, trial_result, "frame.fit")

    assert target.read_text(encoding="utf-8") == "x,y\n1,2\n"
    assert json.loads(target.with_suffix(".json").read_text(encoding="utf-8")) == {
        "fwhm_px": 8.7
    }
