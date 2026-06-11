from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt5")
from apex.gui.workflow.cmd.step8_psf_photometry import (
    _allstar_newton_group,
    _allstar_newton_one,
    _load_detect_positions,
)


def test_step8_detection_loader_prefers_step4_quality_flux(tmp_path):
    result_dir = tmp_path / "result"
    cache_dir = tmp_path / "cache"
    step4_dir = result_dir / "step4_detection"
    cache_dir.mkdir()
    step4_dir.mkdir(parents=True)

    pd.DataFrame(
        {
            "det_uid": [1, 2],
            "x": [float("nan"), 30.0],
            "y": [10.0, 40.0],
            "flux_for_quality": [111.0, 1234.0],
            "dao_flux": [222.0, 999.0],
            "peak_adu": [333.0, 888.0],
            "epsf_candidate": [True, True],
        }
    ).to_csv(step4_dir / "detect_frame.fits.csv", index=False)

    out = _load_detect_positions("frame.fits", cache_dir, result_dir)

    assert out is not None
    assert list(out.index) == [0]
    assert out.loc[0, "det_uid"] == 2
    assert out.loc[0, "flux_init"] == 1234.0


def test_allstar_newton_degenerate_paths_keep_five_value_contract():
    eval_psf = lambda x, y: np.ones_like(x, dtype=float)

    small = _allstar_newton_one(
        np.ones((2, 2)),
        1.0,
        1.0,
        0,
        0,
        10.0,
        eval_psf,
    )
    singular = _allstar_newton_one(
        np.ones((5, 5)),
        2.0,
        2.0,
        0,
        0,
        10.0,
        eval_psf,
    )
    grouped = _allstar_newton_group(
        np.ones((2, 2)),
        [(1.0, 1.0, 10.0), (2.0, 2.0, 8.0)],
        0,
        0,
        eval_psf,
    )

    assert len(small) == 5
    assert len(singular) == 5
    assert all(len(result) == 5 for result in grouped)
