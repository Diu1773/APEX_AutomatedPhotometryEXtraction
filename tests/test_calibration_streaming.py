"""Streaming master combine must reproduce the in-memory path exactly.

O1 (docs/audit/APEX_PERF_DEV_PLAN.md) replaced "load every calibration frame,
then combine" with a row-band sweep that reads only the rows it needs.  The
whole optimisation rests on the two paths being numerically identical, so that
is what these tests pin.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from apex.analysis import calibration as C


def _write_frame(path, data, **cards):
    hdu = fits.PrimaryHDU(np.asarray(data, dtype=np.float32))
    for key, value in cards.items():
        hdu.header[key] = value
    hdu.writeto(path, overwrite=True)
    return str(path)


@pytest.fixture
def frames(tmp_path):
    """A small stack with structure, outliers and a dead pixel."""
    rng = np.random.default_rng(1234)
    paths = []
    for i in range(7):
        data = rng.normal(500.0, 8.0, (60, 41)).astype(np.float32)
        data[10 + i, 5] += 4000.0          # a per-frame spike for the clipper
        data[0, 0] = np.nan                # dead pixel
        paths.append(_write_frame(tmp_path / f"f{i:02d}.fits", data, EXPTIME=30.0))
    return paths


@pytest.mark.parametrize("method", ["median", "mean", "sigmaclip_mean"])
def test_streaming_matches_in_memory(frames, method):
    opts = C.CalibrationOptions()
    reference = C.combine_frames(
        [C.load_frame(p, opts)[0] for p in frames], method, 3.0, 3.0, 5)
    streamed = C.combine_frames_streaming(frames, method, 3.0, 3.0, 5)
    assert np.array_equal(np.nan_to_num(reference, nan=-1.0),
                          np.nan_to_num(streamed, nan=-1.0))


def test_streaming_is_band_size_invariant(frames):
    """The band size is a memory knob, never a result knob."""
    big = C.combine_frames_streaming(frames, "median", 3.0, 3.0, 5,
                                     chunk_bytes=8 * 1024 * 1024)
    tiny = C.combine_frames_streaming(frames, "median", 3.0, 3.0, 5,
                                      chunk_bytes=1024)  # forces 1-row bands
    assert np.array_equal(np.nan_to_num(big, nan=-1.0),
                          np.nan_to_num(tiny, nan=-1.0))


def test_prepare_receives_matching_band_rows(frames):
    """`prepare` gets the row range so it can slice a master correctly."""
    opts = C.CalibrationOptions()
    reference_master = C.combine_frames(
        [C.load_frame(p, opts)[0] for p in frames], "median", 3.0, 3.0, 5)

    in_memory = C.combine_frames(
        [C.load_frame(p, opts)[0] - reference_master for p in frames],
        "median", 3.0, 3.0, 5)
    streamed = C.combine_frames_streaming(
        frames, "median", 3.0, 3.0, 5,
        prepare=lambda values, _i, y0, y1: values - reference_master[y0:y1])
    assert np.array_equal(np.nan_to_num(in_memory, nan=-1.0),
                          np.nan_to_num(streamed, nan=-1.0))


def test_build_master_bias_uses_streaming_and_matches(frames, monkeypatch):
    opts = C.CalibrationOptions(overscan_enable=False)
    streamed, prov = C.build_master_bias(frames, opts)

    # Force the legacy path and compare.
    monkeypatch.setattr(C, "can_stream_combine", lambda _opts: False)
    in_memory, _ = C.build_master_bias(frames, opts)
    assert np.array_equal(np.nan_to_num(streamed, nan=-1.0),
                          np.nan_to_num(in_memory, nan=-1.0))
    assert prov["n_frames"] == len(frames)


def test_overscan_disables_streaming():
    """Overscan works along columns and may trim, so a row band is not
    self-contained; those runs must keep the original path."""
    assert C.can_stream_combine(C.CalibrationOptions(overscan_enable=False))
    assert not C.can_stream_combine(C.CalibrationOptions(overscan_enable=True))


def test_streaming_rejects_empty_input():
    with pytest.raises(ValueError):
        C.combine_frames_streaming([], "median")
