"""Steps 6 and 7 must leave pictures, not only tables.

Both windows drew into canvases and saved nothing, so a headless run produced
`ref_frame_stats.csv`, `centering_stats.csv` and `apcorr_summary.csv` with no
figure of any of them.

Step 7's growth curve is the awkward one and worth naming: unlike every other
product it never reaches disk. `run_forced_photometry` hands it to `apcorr_cb`
and the window plots it live, so the batch has to catch it during the run. It
also arrived without the frame's name — the window knew that from the row the
user had clicked — which would have left every figure untitled and unfileable.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from apex.analysis.forcedphot_qc import (
    draw_center_shift, draw_growth_curve, export_forcedphot_qc, one_per_filter,
)
from apex.analysis.refbuild_qc import export_refbuild_qc


def _params(tmp_path):
    return SimpleNamespace(P=SimpleNamespace(
        result_dir=str(tmp_path), cache_dir=str(tmp_path / "cache"),
        centroid_outlier_px=1.0))


def test_the_growth_curve_draws_from_the_run_payload():
    gc = {
        "radii_px": list(np.linspace(2.0, 20.0, 12)),
        "enclosed_frac": list(np.linspace(0.25, 1.0, 12)),
        "mag_err": list(np.linspace(0.02, 0.03, 12)),
        "r_ap_px": 6.1, "r_ref_px": 18.4, "r_opt_px": 6.8,
        "fwhm_px": 7.67, "apcorr": 1.5317, "fname": "frame-B.fit",
    }
    fig = Figure(figsize=(6, 6))
    assert draw_growth_curve(fig.add_subplot(211), fig.add_subplot(212), gc)

    fig2 = Figure(figsize=(6, 6))
    assert not draw_growth_curve(fig2.add_subplot(211), fig2.add_subplot(212), {})


def test_one_curve_per_filter_is_kept():
    curves = {f"f{i}.fit": {"filter": "BVR"[i % 3]} for i in range(9)}
    picked = one_per_filter(curves)
    assert len(picked) == 3
    assert {v["filter"] for v in picked.values()} == {"B", "V", "R"}


def test_the_centre_shift_reads_the_photometry_table(tmp_path):
    out = tmp_path / "step7_forced_phot"
    out.mkdir(parents=True)
    rng = np.random.default_rng(20260818)
    pd.DataFrame({
        "center_error_px": np.abs(rng.normal(0.4, 0.3, 400)),
        "mag_err": np.abs(rng.normal(0.03, 0.01, 400)),
        "centroid_outlier": rng.random(400) < 0.1,
    }).to_csv(out / "photometry_frame-B.fit.tsv", sep="\t", index=False)

    fig = Figure(figsize=(10, 4))
    assert draw_center_shift(fig.add_subplot(121), fig.add_subplot(122),
                             tmp_path, "frame-B.fit",
                             {"center_error_p90_px": 1.05,
                              "centroid_outlier_rate": 0.1057})

    fig2 = Figure(figsize=(10, 4))
    assert not draw_center_shift(fig2.add_subplot(121), fig2.add_subplot(122),
                                 tmp_path, "missing.fit")


def test_the_export_writes_both_kinds(tmp_path):
    out = tmp_path / "step7_forced_phot"
    out.mkdir(parents=True)
    pd.DataFrame({
        "center_error_px": [0.2, 0.5, 1.4], "mag_err": [0.01, 0.02, 0.05],
        "centroid_outlier": [False, False, True],
    }).to_csv(out / "photometry_frame-B.fit.tsv", sep="\t", index=False)
    pd.DataFrame([{"file": "frame-B.fit", "center_error_p90_px": 1.0,
                   "centroid_outlier_rate": 0.05}]).to_csv(
        out / "centering_stats.csv", index=False)

    written = export_forcedphot_qc(tmp_path, _params(tmp_path), {
        "frame-B.fit": {"filter": "B", "radii_px": [2.0, 5.0, 9.0],
                        "enclosed_frac": [0.4, 0.8, 1.0],
                        "mag_err": [0.03, 0.02, 0.025],
                        "r_ap_px": 5.0, "r_ref_px": 9.0, "apcorr": 1.4},
    })
    names = {p.name for p in written}
    assert "step7_center_shift.png" in names
    assert any(n.startswith("step7_growth_curve_") for n in names)
    for path in written:
        assert path.stat().st_size > 0


def test_the_run_labels_each_curve_with_its_frame():
    """Without a name the batch cannot title or file the figure."""
    source = (Path(__file__).absolute().parents[1]
              / "apex/analysis/forced_photometry.py").read_text(encoding="utf-8")
    assert 'gc_data["fname"] = fname_r' in source


def test_step6_writes_nothing_without_its_input(tmp_path):
    (tmp_path / "step6_refbuild").mkdir(parents=True)
    assert export_refbuild_qc(_params(tmp_path)) == []


def test_both_steps_call_the_shared_modules():
    import io
    import tokenize

    repo = Path(__file__).absolute().parents[1]
    for rel, needle in (("apex/pipeline/steps/refbuild.py", "refbuild_qc"),
                        ("apex/pipeline/steps/forcedphot.py", "forcedphot_qc"),
                        ("apex/gui/workflow/step6_ref_build.py", "refbuild_qc"),
                        ("apex/gui/workflow/step7_forced_aperture_phot.py", "forcedphot_qc")):
        source = (repo / rel).read_text(encoding="utf-8")
        code = " ".join(
            tok.string
            for tok in tokenize.generate_tokens(io.StringIO(source).readline)
            if tok.type not in (tokenize.COMMENT, tokenize.STRING)
        )
        assert needle in code, f"{rel} 이 공용 모듈을 안 쓴다"
