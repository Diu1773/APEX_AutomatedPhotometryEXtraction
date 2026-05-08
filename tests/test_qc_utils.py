from __future__ import annotations

import pandas as pd

from apex.utils.qc_utils import filter_files_by_qc


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
