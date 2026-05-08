from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("PyQt5")
from apex.gui.workflow.cmd.step8_psf_photometry import _load_detect_positions


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
