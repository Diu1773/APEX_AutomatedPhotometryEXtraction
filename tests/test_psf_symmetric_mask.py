"""ePSF 참조별에서 점광원 잡음을 걸러내는 PSF 대칭 판정 회귀 검사.

별은 등방 PSF 라 좌우 또는 상하 **양쪽** 이웃이 피크의 일정 비율 이상이지만,
우주선·핫픽셀은 그렇지 않다. 이 판정이 없으면 ePSF 참조별 선택의 「밝고 고립」
기준이 정확히 우주선을 최적해로 뽑는다 (2026-07-29 M67/QHY600: ePSF 가 실제
별보다 2.75배 좁아지고 PSF 플럭스가 구경의 32% 로 떨어졌다).
"""

import numpy as np

from apex.analysis.psf_policy import psf_symmetric_mask


def _gaussian(img, x, y, amp, sigma):
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx]
    img += amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * sigma ** 2))
    return img


def test_star_passes_spike_patterns_rejected():
    img = np.zeros((60, 60), float)
    _gaussian(img, 10, 10, 200.0, 2.0)      # FWHM ~4.7 px
    img[30, 30] = 500.0                      # 고립 스파이크
    img[40, 40] = img[40, 41] = 500.0        # 2픽셀 수평쌍
    img[50, 50] = img[51, 51] = 500.0        # 대각쌍

    xy = np.array([[10, 10], [30, 30], [40, 40], [50, 50]], float)
    keep = psf_symmetric_mask(img, xy, background=0.0, neighbor_frac=0.3)

    assert keep[0], "실제 별이 걸리면 안 된다"
    assert not keep[1], "고립 스파이크는 걸러야 한다"
    assert not keep[2], "2픽셀 쌍은 걸러야 한다"
    assert not keep[3], "대각쌍은 걸러야 한다"


def test_scale_invariant_faint_star_survives():
    """밝기와 무관해야 한다 — 어두운 별도 통과, 아주 밝은 스파이크도 차단."""
    img = np.zeros((40, 40), float)
    _gaussian(img, 10, 10, 5.0, 2.0)     # 아주 어두운 별
    img[25, 25] = 1.0e5                   # 아주 밝은 스파이크

    keep = psf_symmetric_mask(
        img, np.array([[10, 10], [25, 25]], float), background=0.0
    )
    assert keep[0]
    assert not keep[1]


def test_background_is_subtracted_before_ratio():
    """하늘이 높으면 이웃/피크 비가 1 에 가까워져 잡음도 통과한다.

    배경을 빼야 판정이 성립한다는 것을 고정한다.
    """
    sky = 1000.0
    img = np.full((40, 40), sky, float)
    img[20, 20] += 50.0                   # 하늘 위 고립 스파이크

    xy = np.array([[20, 20]], float)
    assert psf_symmetric_mask(img, xy, background=0.0)[0], (
        "배경을 안 빼면 잡음이 통과한다(이 동작을 알고 있어야 한다)"
    )
    assert not psf_symmetric_mask(img, xy, background=sky)[0], (
        "배경을 빼면 걸러진다"
    )


def test_border_and_empty_inputs():
    img = np.zeros((10, 10), float)
    img[0, 0] = 100.0
    # 가장자리는 이웃을 못 보므로 통과시킨다(다른 컷이 처리)
    assert psf_symmetric_mask(img, np.array([[0, 0]], float))[0]
    # 빈 입력
    assert psf_symmetric_mask(img, np.zeros((0, 2))).shape == (0,)


def test_flat_zero_peak_is_kept():
    """피크가 0 이하면 모양 정보가 없으므로 다른 컷에 맡긴다."""
    img = np.zeros((20, 20), float)
    keep = psf_symmetric_mask(img, np.array([[10, 10]], float), background=0.0)
    assert keep[0]
