from pathlib import Path

import pandas as pd
import pytest

from apex.benchmark.cmd_batch import _coarse_m50_seed, _precision_grid, select_cmd_frames


def _write_project(root: Path) -> None:
    index_dir = root / "result" / "step7_forced_phot"
    zp_dir = root / "result" / "cmd_zeropoint"
    index_dir.mkdir(parents=True)
    zp_dir.mkdir(parents=True)
    index_rows = []
    zp_rows = []
    for filter_name in ("g", "r", "i"):
        for number, fwhm in enumerate((8.0, 4.0, 6.0, 5.0, 7.0), start=1):
            filename = f"{filter_name}_{number}.fit"
            (root / filename).touch()
            index_rows.append(
                {
                    "file": filename,
                    "filter": filter_name,
                    "status": "ok",
                    "wcs_ok": True,
                    "fwhm_px": fwhm,
                }
            )
            zp_rows.append(
                {
                    "file": filename,
                    "filter": filter_name,
                    "zp_frame": 25.0 + number / 100,
                    "zp_scatter": 0.02,
                }
            )
    pd.DataFrame(index_rows).to_csv(index_dir / "photometry_index.csv", index=False)
    pd.DataFrame(zp_rows).to_csv(zp_dir / "frame_zeropoint.csv", index=False)


def test_select_cmd_frames_uses_unique_fwhm_quantiles(tmp_path):
    _write_project(tmp_path)

    selected = select_cmd_frames(tmp_path)

    assert len(selected) == 9
    for _, group in selected.groupby("filter_key"):
        by_condition = group.set_index("condition")
        assert by_condition.loc["best", "fwhm_px"] == pytest.approx(4.0)
        assert by_condition.loc["median", "fwhm_px"] == pytest.approx(6.0)
        assert by_condition.loc["worst", "fwhm_px"] == pytest.approx(8.0)
        assert group["file"].is_unique
        assert all(Path(path).is_file() for path in group["input_fits"])


def test_precision_grid_is_centered_and_inclusive():
    grid = _precision_grid(15.0, 0.6, 0.2)

    assert grid == pytest.approx([14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6])


def test_coarse_m50_seed_rejects_unbracketed_bright_grid():
    points = pd.DataFrame(
        {
            "magnitude": [12.0, 12.5, 13.0],
            "completeness": [1.0, 0.96, 0.92],
        }
    )

    with pytest.raises(RuntimeError, match="too bright"):
        _coarse_m50_seed(points, {"m50": 15.0}, "frame.fit")


def test_coarse_m50_seed_falls_back_to_interpolation_when_fit_is_outside_grid():
    points = pd.DataFrame(
        {
            "magnitude": [12.0, 12.5, 13.0],
            "completeness": [0.9, 0.5, 0.1],
        }
    )

    assert _coarse_m50_seed(points, {"m50": 14.0}, "frame.fit") == pytest.approx(12.5)


def test_select_cmd_frames_requires_calibrated_unique_frames(tmp_path):
    _write_project(tmp_path)
    zp_path = tmp_path / "result" / "cmd_zeropoint" / "frame_zeropoint.csv"
    zeropoints = pd.read_csv(zp_path)
    zeropoints = zeropoints[~((zeropoints["filter"] == "g") & (zeropoints.index > 1))]
    zeropoints.to_csv(zp_path, index=False)

    with pytest.raises(RuntimeError, match="Filter g has 2 calibrated frames"):
        select_cmd_frames(tmp_path)
