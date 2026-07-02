from types import SimpleNamespace

import pandas as pd

from apex.utils.qc_utils import (
    filter_files_by_qc,
    frame_quality_has_auto_qc,
    should_use_frame_quality_qc,
)
from apex.utils.step_paths import step4_dir


def _write_frame_quality(result_dir, rows):
    out = step4_dir(result_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "frame_quality.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_auto_qc_frame_quality_forces_downstream_qc(tmp_path):
    _write_frame_quality(
        tmp_path,
        [
            {"file": "good.fit", "passed": True, "exclude_reason": ""},
            {"file": "bad.fit", "passed": False, "exclude_reason": "auto_qc_fail,auto:high_fwhm"},
        ],
    )
    params = SimpleNamespace(phot_use_qc_pass_only=False)

    assert frame_quality_has_auto_qc(tmp_path)
    assert should_use_frame_quality_qc(
        tmp_path,
        params,
        "phot_use_qc_pass_only",
        default=False,
    )
    kept, info = filter_files_by_qc(
        tmp_path,
        ["good.fit", "bad.fit"],
        require_qc=should_use_frame_quality_qc(tmp_path, params, "phot_use_qc_pass_only"),
    )

    assert kept == ["good.fit"]
    assert info["applied"] is True
    assert info["dropped"] == 1


def test_manual_frame_quality_does_not_force_disabled_downstream_qc(tmp_path):
    _write_frame_quality(
        tmp_path,
        [
            {"file": "good.fit", "passed": True, "exclude_reason": ""},
            {"file": "manual_bad.fit", "passed": False, "exclude_reason": "manual"},
        ],
    )
    params = SimpleNamespace(phot_use_qc_pass_only=False)

    assert not frame_quality_has_auto_qc(tmp_path)
    assert not should_use_frame_quality_qc(
        tmp_path,
        params,
        "phot_use_qc_pass_only",
        default=False,
    )
