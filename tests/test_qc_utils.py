from __future__ import annotations

import pandas as pd

from apex.utils.qc_utils import filter_files_by_qc, filter_frame_df_by_qc, load_frame_excludes


def test_filter_files_by_qc_uses_step4_frame_quality(tmp_path):
    result_dir = tmp_path / "result"
    step4_dir = result_dir / "step4_detection"
    step4_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "file": ["a.fit", "b.fit", "c.fit"],
            "passed": [True, False, True],
        }
    ).to_csv(step4_dir / "frame_quality.csv", index=False)

    files, info = filter_files_by_qc(result_dir, ["a.fit", "b.fit", "c.fit"], require_qc=True)

    assert files == ["a.fit", "c.fit"]
    assert info["applied"] is True
    assert info["kept"] == 2
    assert info["dropped"] == 1


def test_filter_files_by_qc_is_noop_when_disabled(tmp_path):
    files, info = filter_files_by_qc(tmp_path, ["a.fit", "b.fit"], require_qc=False)

    assert files == ["a.fit", "b.fit"]
    assert info["applied"] is False
    assert info["reason"] == "disabled"


def test_filter_files_by_qc_matches_windows_paths_by_basename(tmp_path):
    result_dir = tmp_path / "result"
    step4_dir = result_dir / "step4_detection"
    step4_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "file": ["a.fit", "b.fit"],
            "passed": ["True", "False"],
        }
    ).to_csv(step4_dir / "frame_quality.csv", index=False)

    files, info = filter_files_by_qc(
        result_dir,
        [r"E:\obs\a.fit", r"E:\obs\b.fit"],
        require_qc=True,
    )

    assert files == [r"E:\obs\a.fit"]
    assert info["kept"] == 1


def test_filter_frame_df_by_qc_filters_indexes(tmp_path):
    result_dir = tmp_path / "result"
    step4_dir = result_dir / "step4_detection"
    step4_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "file": ["a.fit", "b.fit", "c.fit"],
            "passed": [True, False, True],
        }
    ).to_csv(step4_dir / "frame_quality.csv", index=False)
    idx = pd.DataFrame({"file": ["a.fit", "b.fit", "c.fit"], "n": [1, 2, 3]})

    filtered, info = filter_frame_df_by_qc(result_dir, idx, require_qc=True)

    assert filtered["file"].tolist() == ["a.fit", "c.fit"]
    assert info["total"] == 3
    assert info["kept"] == 2


def test_load_frame_excludes_parses_string_false(tmp_path):
    pd.DataFrame(
        {
            "file": ["a.fit", "b.fit"],
            "passed": ["False", "True"],
            "exclude_reason": ["manual", ""],
        }
    ).to_csv(tmp_path / "frame_exclude.csv", index=False)

    assert load_frame_excludes(tmp_path) == {"a.fit": {"manual"}}
