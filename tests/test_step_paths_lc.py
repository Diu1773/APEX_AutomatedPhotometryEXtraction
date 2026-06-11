from pathlib import Path

from apex.utils.step_paths_lc import (
    find_best_lightcurve_csv,
    list_lightcurve_csvs,
    step9_lc_dir,
    step10_detrend_dir,
)


def _write_csv(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("JD,diff_mag\n1.0,0.01\n", encoding="utf-8")
    return path


def test_find_best_lightcurve_prefers_step10_current(tmp_path):
    raw = _write_csv(step9_lc_dir(tmp_path) / "lightcurve_ID42_raw.csv")
    current = _write_csv(
        step10_detrend_dir(tmp_path) / "lightcurve_ID42_current.csv"
    )

    assert find_best_lightcurve_csv(tmp_path, 42) == current
    assert list_lightcurve_csvs(tmp_path, 42)[:2] == [current, raw]


def test_find_best_lightcurve_falls_back_to_step9_writer_name(tmp_path):
    raw = _write_csv(step9_lc_dir(tmp_path) / "lightcurve_ID7_raw.csv")

    assert find_best_lightcurve_csv(tmp_path, 7) == raw


def test_find_best_lightcurve_skips_corrupt_preferred_output(tmp_path):
    current = step10_detrend_dir(tmp_path) / "lightcurve_ID7_current.csv"
    current.parent.mkdir(parents=True)
    current.write_text("not,a,usable,lightcurve\n", encoding="utf-8")
    raw = _write_csv(step9_lc_dir(tmp_path) / "lightcurve_ID7_raw.csv")

    assert find_best_lightcurve_csv(tmp_path, 7) == raw


def test_list_lightcurve_csvs_does_not_return_missing_candidates(tmp_path):
    assert list_lightcurve_csvs(tmp_path, 99) == []
