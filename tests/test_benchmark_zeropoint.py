from pathlib import Path
from types import SimpleNamespace

from astropy.io import fits
import pandas as pd
import pytest

from apex.benchmark.runner import resolve_benchmark_zeropoint


def test_resolve_benchmark_zeropoint_prefers_matching_frame(tmp_path):
    zp_dir = tmp_path / "cmd_zeropoint"
    zp_dir.mkdir()
    pd.DataFrame(
        [
            {
                "file": "frame.fit",
                "filter": "r",
                "zp_frame": 25.123,
                "zp_scatter": 0.02,
                "n_ref": 42,
            }
        ]
    ).to_csv(zp_dir / "frame_zeropoint.csv", index=False)
    p = SimpleNamespace(result_dir=tmp_path, zp_initial=25.0)

    value, source, meta = resolve_benchmark_zeropoint(
        Path("frame.fit"),
        fits.Header({"FILTER": "r"}),
        p,
        None,
        allow_initial_fallback=False,
    )

    assert value == pytest.approx(25.123)
    assert source.endswith(":frame")
    assert float(meta["n_ref"]) == pytest.approx(42)


def test_resolve_benchmark_zeropoint_rejects_silent_initial_fallback(tmp_path):
    p = SimpleNamespace(result_dir=tmp_path, zp_initial=25.0)
    with pytest.raises(RuntimeError, match="No calibrated Step 10 zero point"):
        resolve_benchmark_zeropoint(
            Path("missing.fit"),
            fits.Header({"FILTER": "i"}),
            p,
            None,
            allow_initial_fallback=False,
        )
