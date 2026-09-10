"""계수를 냈으면 그 계수가 얼마나 흔들리는지도 같이 내야 한다.

Step 10 은 색항을 적합해 적용하면서 그 불확실도를 내지 않았다. 그래서
`ct = -0.204` 와 `ct = -0.204 ± 0.311` 이 산출물에서 똑같이 보였다.

2026-09-09 에 M67 을 세 기기로 재서 확인한 것 (`validation/ERROR_BUDGET.md`):

    MuSCAT3 r    ct = -0.204   sigma(ct) = 0.311   <- 못 재는 값
    kb26 B       ct = -0.244   sigma(ct) = 0.278   <- 못 재는 값
    나머지 26 조합                sigma(ct) = 0.005~0.057

못 재는 값을 적용하면 색 범위 끝에서 0.12~0.30 등급이 들어간다.

**적합 공분산으로는 이것을 못 낸다.** 그 값은 별을 그대로 두고 잡음만 다시 뽑는
폭이라 절단이 만드는 흔들림을 안 담는다. 그래서 별을 복원추출로 다시 뽑아 매번
처음부터(절단 포함) 다시 맞춘다.
"""

from __future__ import annotations

import numpy as np
import pytest

from apex.analysis.cmd.zeropoint_runner import (
    ZeropointCalibrationRunner,
    bootstrap_zp_ct,
)


class _Fitter:
    """`_robust_linfit` 은 러너의 메서드인데 self 에서 쓰는 것이 로그뿐이다."""

    def _log(self, *a, **k):
        pass

    fit = ZeropointCalibrationRunner._robust_linfit


FIT = _Fitter()


def _fit(x, y, w=None, **kw):
    return FIT.fit(x, y, w=w, **kw)


def _clean(n=400, ct=0.15, zp=-3.2, noise=0.01, seed=7):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.2, 1.6, n)
    y = zp + ct * x + rng.normal(0.0, noise, n)
    return x, y, np.full(n, 1.0 / noise**2)


def test_a_well_measured_colour_term_gets_a_small_uncertainty():
    x, y, w = _clean()
    zp_sig, ct_sig, n_ok = bootstrap_zp_ct(_fit, x, y, w, draws=120)

    assert n_ok >= 100
    assert 0.0 < ct_sig < 0.01, f"깨끗한 자료인데 sigma(ct)={ct_sig}"
    assert 0.0 < zp_sig < 0.01


def _many_fits(n, lo, hi, noise, true_ct, realisations=40, seed0=100):
    """같은 조건으로 자료를 여러 벌 만들어 매번 적합한다 — 반복 측정."""
    out = []
    for k in range(realisations):
        rng = np.random.default_rng(seed0 + k)
        x = rng.uniform(lo, hi, n)
        y = -3.2 + true_ct * x + rng.normal(0.0, noise, n)
        w = np.full(n, 1.0 / noise**2)
        ct = _fit(x, y, w=w)[1]
        if np.isfinite(ct):
            out.append((float(ct), x, y, w))
    return out


@pytest.mark.parametrize(
    "label, n, lo, hi, noise",
    [("색 지렛대가 넉넉할 때", 400, 0.2, 1.6, 0.01),
     ("색 지렛대가 0.2 등급뿐일 때", 60, 0.55, 0.75, 0.08)],
)
def test_the_uncertainty_predicts_the_spread_of_repeated_measurements(
        label, n, lo, hi, noise):
    """σ 가 가져야 할 성질은 하나다 — **다시 재면 얼마나 달라지는지를 맞히는 것.**

    한 벌의 자료로 「이 값은 못 잰다」를 단정하면 운에 기댄다. 참값이 0 인
    자료에서도 적합이 2σ 를 낼 수 있고 실제로 그렇게 나왔다(어떤 씨앗에서
    ct = −0.36, σ = 0.16). 그래서 여기서는 같은 조건의 자료를 마흔 벌 만들어
    **적합값이 실제로 흩어지는 폭**을 재고, 부트스트랩이 그것을 맞히는지 본다.
    """
    fits = _many_fits(n, lo, hi, noise, true_ct=0.0)
    assert len(fits) >= 30

    empirical = float(np.std([f[0] for f in fits], ddof=1))
    boots = [bootstrap_zp_ct(_fit, x, y, w, draws=80)[1] for _ct, x, y, w in fits[:12]]
    boots = [b for b in boots if np.isfinite(b)]
    assert len(boots) >= 10
    predicted = float(np.median(boots))

    ratio = predicted / empirical
    assert 0.5 < ratio < 2.0, (
        f"{label}: 반복 측정의 산포 {empirical:.4f} 인데 부트스트랩은 "
        f"{predicted:.4f} 를 낸다 (비 {ratio:.2f})")


def test_a_short_colour_lever_gives_a_much_larger_uncertainty():
    """지렛대가 짧으면 σ 가 커야 한다 — 상수를 내고 있지 않다는 확인."""
    wide = bootstrap_zp_ct(_fit, *_clean(n=400), draws=120)[1]

    rng = np.random.default_rng(11)
    x = rng.uniform(0.55, 0.75, 60)
    y = -3.2 + rng.normal(0.0, 0.08, 60)
    narrow = bootstrap_zp_ct(_fit, x, y, np.full(60, 1.0 / 0.08**2), draws=120)[1]

    assert narrow > 10 * wide, f"넓을 때 {wide:.4f} · 좁을 때 {narrow:.4f}"


def test_the_same_data_gives_the_same_uncertainty_twice():
    """씨앗이 고정되어야 산출물의 수를 되짚을 수 있다."""
    x, y, w = _clean()
    a = bootstrap_zp_ct(_fit, x, y, w, draws=80)
    b = bootstrap_zp_ct(_fit, x, y, w, draws=80)
    assert a == b


def test_too_few_stars_is_unknown_rather_than_zero():
    x, y, w = _clean(n=12)
    zp_sig, ct_sig, n_ok = bootstrap_zp_ct(_fit, x, y, w, draws=80)

    assert n_ok == 0
    assert not np.isfinite(zp_sig) and not np.isfinite(ct_sig), (
        "잴 수 없을 때 0 을 돌려주면 '완벽히 안다'는 뜻이 된다")


def test_more_stars_shrink_the_uncertainty():
    """부트스트랩이 자료량에 반응하는지 — 상수를 내고 있지 않다는 확인."""
    small = bootstrap_zp_ct(_fit, *_clean(n=40, seed=3), draws=120)[1]
    large = bootstrap_zp_ct(_fit, *_clean(n=800, seed=3), draws=120)[1]

    assert large < small / 2, f"별이 20 배인데 sigma(ct) {small:.4f} -> {large:.4f}"
