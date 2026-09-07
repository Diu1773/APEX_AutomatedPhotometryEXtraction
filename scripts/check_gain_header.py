"""헤더의 GAIN 을 믿을지 PTC 실측을 쓸지 판정한다.

**왜 규칙이 필요한가.** 이 프로젝트는 헤더가 크게 틀린 기기(Moravian C3-61000,
`EGAIN` 이 실제 변환이득의 1/14)와 잘 맞는 기기(LCO MuSCAT3, 1 % 안)를 둘 다
겪었다. 그러니 「헤더를 믿지 마라」도 「믿어라」도 규칙이 못 된다. 문턱이 있어야
하고, 그 문턱은 지어낸 수가 아니라 측정에서 나와야 한다.

## 문턱 둘

**첫째, 차이가 잴 수 있는 차이인가.** 실측에는 오차가 있다. 5 % 차이는 PTC 가
±1 % 면 「헤더가 틀렸다」이고 PTC 가 ±10 % 면 「아무 말도 못 한다」이다. 그래서
차이 자체가 아니라 **차이를 PTC 오차로 나눈 값**으로 판정한다. 3 시그마를 쓴다.

주의: PTC 가 신호 범위 경고를 내면 그때 인용하는 오차는 참값보다 작다. 도구가
경고를 달았으면 이 스크립트는 오차를 세 배로 부풀려 본다 — MuSCAT3 에서 세
카메라의 어긋남(0.85 %)이 인용 오차(0.25 %)의 3.4 배였다.

**둘째, 그 차이가 결과를 움직이는가.** gain 은 등급에 안 들어간다. CCD 잡음식에서
ADU 를 전자로 바꾸는 데만 쓰이므로 **오차막대만 움직인다.** 실제 식으로 재 보면
gain 을 k 배 틀리게 넣었을 때 보고 σ 가 이만큼 움직인다.

    무엇이 지배하나       gain +1%   +3%     +10%    +1400%
    별(밝은 별)           -0.50%   -1.47%   -4.66%   -74%
    하늘                  -0.58%   -1.72%   -5.43%   -76%
    읽기잡음(어두운 별)   -0.91%   -2.69%   -8.40%   -88%   <- 최악

**최악의 경우가 거의 1:1 이다.** 그리고 APEX 의 오차 모형은 자기 정확도가 이미
측정되어 있다 — 보고 σ 대 실제 산포의 비가 **0.976**, pull 표준편차 1.014
(`validation/paper/captions/fig2_error_model.md`, 합성 5,089 회수).

그러므로 **모형 자신이 2.4 % 어긋나 있는데 그보다 작은 gain 오차는 안 보인다.**
2.4 % ÷ 0.9(최악 감도) ≈ 2.7 %, 반올림해서 **3 %** 를 문턱으로 둔다.

## 판정표

    차이가 3σ 안   차이가 3 % 안   판정
    ------------   -------------   --------------------------------------------
    예             예              헤더를 쓴다. 실측은 「확인함」으로 기록.
    예             아니오          어느 쪽이든 되지만 실측이 이만큼 벌어졌다는
                                   것 자체가 이상하다 — PTC 를 다시 본다.
    아니오         예              헤더가 다른 값을 적고 있을 수 있다. 다만
                                   결과는 안 움직이므로 급하지 않다. 기록만.
    아니오         아니오          **헤더가 실측과 다른 것을 뜻한다.** 값을 바꾸기
                                   전에 무엇을 적은 것인지 밝힌다 (Moravian 은
                                   센서의 명목 레지스터 값이었다).

**「헤더가 틀렸으니 실측값을 넣는다」로 끝내지 않는다.** 크게 어긋나면 그것은
측정 오차가 아니라 **헤더가 다른 물리량을 적고 있다는 뜻**이고, 무엇인지 모르는
채로 숫자만 바꾸면 다음에 같은 자리에서 또 걸린다.

실행:
    python -X utf8 scripts/check_gain_header.py --bias <폴더> --flat <폴더>
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parents[1]))

from apex.analysis.detector_ptc import characterize_detector  # noqa: E402
from apex.utils.io_utils import read_fits_header  # noqa: E402

#: 차이가 실측 오차의 몇 배를 넘으면 「우연이 아니다」로 볼 것인가.
SIGMA_K = 3.0

#: 도구가 신호 범위 경고를 냈을 때 인용 오차를 몇 배로 부풀려 볼 것인가.
#: MuSCAT3 세 카메라의 실제 어긋남이 인용 오차의 3.4 배였다.
WARN_INFLATE = 3.0

#: 이 아래면 오차막대가 오차 모형 자신의 정확도(2.4 %) 안에서 움직인다.
#: 2.4 % ÷ 0.9(읽기잡음 지배 상황의 감도) ≈ 2.7 % -> 3 %.
MATTERS_FRAC = 0.03


def _paths(spec: str) -> list[str]:
    if os.path.isdir(spec):
        out: list[str] = []
        for ext in ("*.fits", "*.fit", "*.fits.fz", "*.fts"):
            out += glob.glob(os.path.join(spec, ext))
        return sorted(out)
    return sorted(glob.glob(spec))


def verdict(header_gain: float, ptc_gain: float, ptc_err: float,
            warned: bool) -> tuple[str, str, float, float]:
    """(판정, 설명, 차이 비율, 쓴 오차) 를 돌려준다."""
    diff = ptc_gain / header_gain - 1.0
    err = ptc_err / header_gain
    if warned:
        err *= WARN_INFLATE
    significant = abs(diff) > SIGMA_K * err if err > 0 else True
    matters = abs(diff) > MATTERS_FRAC

    if not significant and not matters:
        return ("헤더를 쓴다",
                "실측과 구별되지 않고, 구별되더라도 오차막대가 안 움직인다.",
                diff, err)
    if not significant and matters:
        return ("PTC 를 다시 본다",
                "실측 오차 안이지만 차이가 3 % 를 넘는다 — 실측이 너무 헐겁다.",
                diff, err)
    if significant and not matters:
        return ("기록만 하고 헤더를 쓴다",
                "차이가 우연은 아니지만 결과를 안 움직인다.", diff, err)
    return ("헤더가 무엇을 적은 것인지 밝힌다",
            "우연이 아니고 결과도 움직인다. 숫자만 바꾸지 말 것.", diff, err)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="헤더 GAIN 과 PTC 실측을 견준다")
    ap.add_argument("--bias", required=True, help="바이어스 폴더 또는 glob")
    ap.add_argument("--flat", required=True, help="플랫 폴더 또는 glob")
    ap.add_argument("--header-gain", type=float, default=None,
                    help="헤더값을 직접 준다 (기본: 플랫 첫 장에서 읽는다)")
    ap.add_argument("--tolerance", type=float, default=0.06,
                    help="플랫 짝 허용 폭 (하늘 플랫은 밝기가 흘러 넓게 잡는다)")
    ap.add_argument("--box-radius", type=int, default=300)
    a = ap.parse_args(argv)

    bias, flat = _paths(a.bias), _paths(a.flat)
    if len(bias) < 2 or len(flat) < 6:
        print(f"자료 부족: 바이어스 {len(bias)} 장, 플랫 {len(flat)} 장 "
              "(각각 2 장, 6 장 이상 필요)")
        return 2

    hg = a.header_gain
    if hg is None:
        h = read_fits_header(flat[0])
        for key in ("EGAIN", "GAIN"):
            try:
                val = float(h.get(key))
            except (TypeError, ValueError):
                continue
            if val > 0:
                hg = val
                break
    if not hg:
        print("헤더에 GAIN/EGAIN 이 없다. --header-gain 으로 준다.")
        return 2

    res = characterize_detector(bias, flat, signal_floor=0.0,
                                box_radius=a.box_radius, tolerance=a.tolerance)
    warned = any("WARNING" in line for line in res.log)
    v, why, diff, err = verdict(hg, res.gain_eff, res.gain_eff_err, warned)

    print(f"헤더 GAIN      : {hg:.4f} e-/ADU")
    print(f"PTC 실측       : {res.gain_eff:.4f} +- {res.gain_eff_err:.4f}"
          f"  ({res.n_pairs} 쌍, R^2 = {res.r_squared:.4f})")
    print(f"차이           : {diff*100:+.2f} %"
          f"   (쓴 오차 {err*100:.2f} %{' — 경고가 있어 세 배로 부풀림' if warned else ''})")
    print(f"우연이 아닌가  : {'예' if abs(diff) > SIGMA_K*err else '아니오'}"
          f"   ({SIGMA_K:.0f}σ = {SIGMA_K*err*100:.2f} %)")
    print(f"결과를 움직이나: {'예' if abs(diff) > MATTERS_FRAC else '아니오'}"
          f"   (문턱 {MATTERS_FRAC*100:.0f} %)")
    print()
    print(f"판정 → **{v}**")
    print(f"       {why}")
    if warned:
        for line in res.log:
            if "WARNING" in line:
                print(f"       PTC 경고: {line.split('WARNING:')[-1].strip()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
