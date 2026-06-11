from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.utils.step_paths import step7_forced_phot_dir


def _qa_report_module():
    return pytest.importorskip("apex.gui.tools.qa_report")


def test_qa_filter_keys_are_canonical():
    qa_report = _qa_report_module()

    _qa_filter_key = qa_report._qa_filter_key
    assert _qa_filter_key("b") == "B"
    assert _qa_filter_key("sdss-r") == "r"
    assert _qa_filter_key("  V  ") == "V"
    assert _qa_filter_key(np.nan) == ""


def test_frame_quality_absolute_fwhm_uses_arcsec(tmp_path):
    qa_report = _qa_report_module()

    phot_dir = step7_forced_phot_dir(tmp_path)
    phot_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {"file": "frame_a", "fwhm_used": 3.0, "r_ap": 4.0, "r_in": 7.0, "r_out": 10.0},
            {"file": "frame_b", "fwhm_used": 5.0, "r_ap": 4.0, "r_in": 7.0, "r_out": 10.0},
        ]
    ).to_csv(phot_dir / "aperture_by_frame.csv", index=False)

    rows = []
    for frame in ("frame_a", "frame_b"):
        for idx in range(4):
            rows.append(
                {
                    "frame": frame,
                    "ID": idx + 1,
                    "FILTER": "B",
                    "mag": 12.0 + idx,
                    "snr": 50.0,
                    "mag_err": 0.02,
                }
            )
    df = pd.DataFrame(rows)

    worker = qa_report.QAReportWorker(
        tmp_path,
        params={
            "frame_flag_mode": "absolute",
            "frame_fwhm_mode": "absolute",
            "frame_fwhm_abs": 2.0,
            "pixel_scale_arcsec": 0.5,
        },
    )
    result = worker._compute_frame_quality(df)
    out = result["data"].set_index("frame")

    assert out.loc["frame_a", "fwhm_arcsec"] == pytest.approx(1.5)
    assert out.loc["frame_b", "fwhm_arcsec"] == pytest.approx(2.5)
    assert bool(out.loc["frame_a", "fwhm_flag"]) is False
    assert bool(out.loc["frame_b", "fwhm_flag"]) is True
