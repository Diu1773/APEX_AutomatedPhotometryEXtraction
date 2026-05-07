"""
STF (Screen Transfer Function) utilities — ported from AstralImage.

Pure-numpy implementation; no Qt dependency.
"""

from __future__ import annotations

import numpy as np


def _solve_mtf_balance(norm_med: float, target: float = 0.25) -> float:
    """MTF(norm_med, m) = target 가 되는 m을 역산."""
    x = float(norm_med)
    t = float(target)
    denom = x - 2.0 * t * x + t
    if abs(denom) < 1e-10:
        return 0.5
    m = x * (1.0 - t) / denom
    return float(np.clip(m, 0.001, 0.999))


def compute_stf_params(
    data: np.ndarray,
    target_median: float = 0.25,
) -> tuple[float, float, float]:
    """STF 파라미터 (shadow, highlight, mid) 계산.

    Iterative sigma-clipping으로 배경 노이즈 σ 추정 → PI/Siril 스타일.

    Parameters
    ----------
    data : array-like
        1D flat 배열 또는 임의 shape의 numpy 배열.
    target_median : float
        STF 적용 후 목표 중앙값 (기본 0.25).

    Returns
    -------
    (shadow, highlight, mid) : tuple[float, float, float]
    """
    arr = np.asarray(data, dtype=np.float32).ravel()
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0, 1.0, 0.5

    step = max(1, arr.size // 65536)
    sample = arr[::step][:65536]

    med = float(np.nanmedian(sample))
    sigma = 1.4826 * float(np.nanmedian(np.abs(sample - med)))
    for _ in range(5):
        if sigma <= 0:
            break
        clipped = sample[np.abs(sample - med) <= 3.0 * sigma]
        if clipped.size < 64:
            break
        nm = float(np.nanmedian(clipped))
        ns = 1.4826 * float(np.nanmedian(np.abs(clipped - nm)))
        if abs(nm - med) < 1e-10:
            med, sigma = nm, ns
            break
        med, sigma = nm, ns

    shadow    = max(0.0, med - 2.8 * sigma)
    highlight = float(np.percentile(sample, 99.99))
    if highlight <= shadow:
        highlight = shadow + 1.0

    scale    = highlight - shadow
    norm_med = float(np.clip((med - shadow) / scale, 1e-6, 1.0 - 1e-6))
    mid      = _solve_mtf_balance(norm_med, target_median)
    return shadow, highlight, mid


def apply_stf_params(
    data: np.ndarray,
    shadow: float,
    highlight: float,
    mid: float,
) -> np.ndarray:
    """사전 계산된 (shadow, highlight, mid)로 STF 적용 → [0, 1] float32 반환."""
    arr = np.asarray(data, dtype=np.float32)
    finite_mask = np.isfinite(arr)
    if not finite_mask.any():
        return np.zeros_like(arr)
    scale = float(highlight - shadow)
    if scale <= 0:
        return np.zeros_like(arr)
    safe = np.where(finite_mask, arr, shadow).astype(np.float32, copy=False)
    norm = (safe - shadow) / scale
    np.clip(norm, 0.0, 1.0, out=norm)
    a = float(mid - 1.0)
    b = float(2.0 * mid - 1.0)
    c = float(-mid)
    denom = norm * b + c
    denom[np.abs(denom) < 1e-7] = 1e-7
    result = norm * a / denom
    np.clip(result, 0.0, 1.0, out=result)
    return result


def apply_stretch(
    data: np.ndarray,
    mode: str,
    shadow: float = 0.0,
    highlight: float = 1.0,
    mid: float = 0.5,
) -> np.ndarray:
    """mode 별 스트레치 적용 → [0, 1] float32 반환.

    Parameters
    ----------
    mode : str
        'stf' | 'linear' | 'log' | 'asinh' | 'sqrt'
    shadow, highlight, mid : float
        STF 파라미터 (stf 모드) 또는 linear 범위 클리핑 (나머지 모드).
    """
    arr = np.asarray(data, dtype=np.float32)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros_like(arr)

    mode = mode.lower()

    if mode == "stf":
        return apply_stf_params(arr, shadow, highlight, mid)

    # 나머지 모드: shadow/highlight를 black/white point로 사용
    scale = float(highlight - shadow)
    if scale <= 0:
        scale = 1.0
    norm = np.clip((arr.astype(np.float32) - shadow) / scale, 0.0, 1.0)

    if mode == "log":
        s = 10.0
        return np.log1p(norm * s) / np.log1p(s)
    if mode == "asinh":
        s = 10.0
        return np.arcsinh(norm * s) / np.arcsinh(s)
    if mode == "sqrt":
        return np.sqrt(norm)
    # linear
    return norm
