"""Step 4's QC products must not need a window.

The window has written `frame_quality_summary.csv`, `step4_qc_overview.png` and
`step4_detection_overlay_examples.png` since it existed. A headless run left
that directory with no figures at all — the export was a method on the dialog.
Same shape as the Step 8 comparison figure and the Step 12 isochrone plots.

Two things had to move, and the second only showed up by rendering the result:

  * the drawing, which is now in `apex.analysis.detection_qc`
  * the QC verdict. `evaluate_frame_qc` already lived in `apex.analysis`, but
    only the window called it. Without a status every mask in the overview is
    empty, so the first headless figure came out blank under a title reading
    "PASS=0 REVIEW=0 FAIL=0". It rendered without error and said nothing.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from apex.analysis.detection_qc import (
    auto_qc_thresholds, build_frame_table, export_qc_products,
)


def _workspace(tmp_path: Path, n: int = 6) -> SimpleNamespace:
    det = tmp_path / "step4_detection"
    det.mkdir(parents=True)
    rng = np.random.default_rng(20260818)
    for i in range(n):
        name = f"frame-{i:03d}-{'BVR'[i % 3]}.fit"
        (det / f"detect_{name}.json").write_text(json.dumps({
            "filter": "BVR"[i % 3],
            "bkg_median": 600.0 + 40 * i,
            "bkg_rms": 28.0 + 0.4 * i,
            "fwhm_px": float(6.5 + 0.2 * rng.standard_normal()),
            "n_sources": int(900 + 40 * rng.standard_normal()),
            "n_raw_detections": 1200,
            "n_after_shape_filter": 1000,
            "median_elongation": 1.15,
            "median_roundness": 0.02,
        }), encoding="utf-8")

    scan = tmp_path / "step1_file_selection"
    scan.mkdir(parents=True)
    pd.DataFrame({
        "Filename": [f"frame-{i:03d}-{'BVR'[i % 3]}.fit" for i in range(n)],
        "AIRMASS": np.linspace(1.01, 1.35, n),
        "JD": np.linspace(2460000.0, 2460000.2, n),
    }).to_csv(scan / "headers.csv", index=False)

    return SimpleNamespace(P=SimpleNamespace(
        result_dir=str(tmp_path), cache_dir=str(tmp_path / "cache"),
        fwhm_elong_max=1.3))


def test_the_frame_table_comes_back_from_disk(tmp_path):
    params = _workspace(tmp_path)
    df = build_frame_table(tmp_path, params)
    assert len(df) == 6
    assert df["fwhm_med"].notna().all()
    assert df["n_sources"].gt(0).all()


def test_the_verdict_is_computed_headless(tmp_path):
    """The blank-figure bug: numbers present, status empty, nothing drawn."""
    params = _workspace(tmp_path)
    df = build_frame_table(tmp_path, params)
    assert "qc_status" in df.columns
    verdicts = set(df["qc_status"].astype(str).str.upper())
    assert verdicts <= {"PASS", "REVIEW", "FAIL"}, verdicts
    assert verdicts != {""}, "판정이 비면 그림이 통째로 빈다"


def test_the_airmass_axis_comes_from_the_header_scan(tmp_path):
    """Otherwise the batch plots against a bare index and the window doesn't."""
    params = _workspace(tmp_path)
    df = build_frame_table(tmp_path, params)
    assert df["airmass"].notna().all()
    assert df["airmass"].max() > df["airmass"].min()


def test_the_three_products_are_written(tmp_path):
    params = _workspace(tmp_path)
    written = export_qc_products(tmp_path, params=params)
    names = {p.name for p in written}
    assert "frame_quality_summary.csv" in names
    assert "step4_qc_overview.png" in names
    for path in written:
        assert path.stat().st_size > 0


def test_an_empty_workspace_writes_nothing_rather_than_a_blank_figure(tmp_path):
    (tmp_path / "step4_detection").mkdir(parents=True)
    assert export_qc_products(tmp_path) == []


def test_the_thresholds_follow_the_configured_elongation_cut(tmp_path):
    params = _workspace(tmp_path)
    params.P.fwhm_elong_max = 1.6
    assert auto_qc_thresholds(params).elong_fail == pytest.approx(1.6)


def test_the_window_and_the_batch_call_the_same_export():
    import io
    import tokenize

    repo = Path(__file__).absolute().parents[1]
    for rel in ("apex/gui/workflow/step4_source_detection.py",
                "apex/pipeline/steps/detect.py"):
        source = (repo / rel).read_text(encoding="utf-8")
        code = " ".join(
            tok.string
            for tok in tokenize.generate_tokens(io.StringIO(source).readline)
            if tok.type not in (tokenize.COMMENT, tokenize.STRING)
        )
        assert "detection_qc" in code, f"{rel} 이 공용 모듈을 안 쓴다"
