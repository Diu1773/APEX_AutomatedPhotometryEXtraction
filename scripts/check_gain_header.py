"""헤더의 GAIN 을 믿을지 PTC 실측을 쓸지 판정한다.

**왜 규칙이 필요한가.** 이 프로젝트는 헤더가 크게 틀린 기기(Moravian C3-61000,
`EGAIN` 이 실제 변환이득의 1/14)와 잘 맞는 기기(LCO MuSCAT3, 1 % 안)를 둘 다
겪었다. 그러니 「헤더를 믿지 마라」도 「믿어라」도 규칙이 못 된다.

## 어디까지가 측정이고 어디부터가 관례인가

이 구분을 흐리지 않는 것이 이 파일의 요점이다 (`Main/OPERATOR.md` C-198).

**측정인 것**

* **PTC 로 잰 gain 값** — Janesick 의 광자전달곡선. 헤더를 안 쓴다.
* **그 측정의 재현성** — 바이어스를 반으로, 플랫을 짝·홀로, 짝 허용 폭을 바꿔
  가며 여러 번 다시 재고 흩어짐을 본다. MuSCAT3 에서 0.26 % · 0.35 % · 0.48 %
  가 나왔다. **인용 오차를 임의 배수로 부풀리지 않는다** — 그렇게 했다가
  근거가 「카메라 세 대의 어긋남」이라는 잘못된 추론이었던 적이 있다.
* **gain 오차가 오차막대를 얼마나 움직이는가** — 그 검출기의 RN 과 하늘 밝기로
  실제 CCD 식에서 계산한다. 최악(읽기잡음 지배)이 거의 1:1 이다.

    무엇이 지배하나       gain +1%   +3%     +10%    +1400%
    별(밝은 별)           -0.50%   -1.47%   -4.66%   -74%
    하늘                  -0.58%   -1.72%   -5.43%   -76%
    읽기잡음(어두운 별)   -0.91%   -2.69%   -8.40%   -88%

**관례인 것 (2026-09-07 사용자 결정: 관례를 쓰되 관례라고 밝힌다)**

* **3σ** — 「우연이 아니다」의 관습적 컷이다. **유도된 값이 아니다.**
* **3 %** — 재료는 측정이지만 조합과 반올림은 선택이다. APEX 오차 모형의
  보고 σ 대 실제 산포 비가 0.976 이므로 모형 자신이 2.4 % 어긋나 있고
  (`validation/paper/captions/fig2_error_model.md`), 최악 감도 0.9 로 나누면
  2.7 %, 반올림해서 3 %. **다만 fig2 의 밴드별 비는 0.83 ~ 1.04 이고 나는
  중앙값을 골랐다.** 어두운 끝을 쓰면 문턱이 17 % 가 된다.

**재현성은 정밀도이지 정확도가 아니다.** 같은 좁은 신호 구간의 플랫을 다시
갈라 봐야 그 구간에 공통으로 걸린 계통은 안 보인다. 도구가 신호 범위 경고를
냈다면 **진짜 불확실도는 여기서 재는 것보다 크다** — 판정문에 그렇게 적는다.

## 판정표

    차이가 3σ 안   차이가 3 % 안   판정
    ------------   -------------   --------------------------------------------
    예             예              헤더를 쓴다. 실측은 「확인함」으로 기록.
    예             아니오          PTC 를 다시 본다 — 실측이 너무 헐겁다.
    아니오         예              기록만 하고 헤더를 쓴다.
    아니오         아니오          **헤더가 무엇을 적은 것인지 밝힌다.** 값을
                                   바꾸기 전에 무엇인지 알아낸다 (Moravian 은
                                   센서의 명목 레지스터 값이었다).

**자릿수가 다르면 문턱이 필요 없다.** Moravian 의 1292 % 는 어떤 불확실도
모형으로도 「좀 다르네」가 되지 않는다. 그건 헤더가 **다른 물리량**을 적고 있다는
뜻이다.

실행:
    python -X utf8 scripts/check_gain_header.py --bias <폴더> --flat <폴더>
    python -X utf8 scripts/check_gain_header.py ... --fast   (재현성 측정 생략)
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).absolute().parents[1]))

from apex.analysis.detector_ptc import characterize_detector  # noqa: E402
from apex.utils.io_utils import read_fits_header  # noqa: E402

#: 관례. 「우연이 아니다」의 관습적 컷이며 유도된 값이 아니다.
SIGMA_K = 3.0

#: 관례 + 측정. 위 설명의 「관례인 것」 참조. 어두운 끝을 기준으로 하면 17 % 다.
MATTERS_FRAC = 0.03

#: 이 배수를 넘으면 측정 불일치가 아니라 헤더가 다른 물리량을 적은 것이다.
ORDER_OF_MAGNITUDE = 2.0


def _paths(spec: str) -> list[str]:
    if os.path.isdir(spec):
        out: list[str] = []
        for ext in ("*.fits", "*.fit", "*.fits.fz", "*.fts"):
            out += glob.glob(os.path.join(spec, ext))
        return sorted(out)
    return sorted(glob.glob(spec))


def _gain(bias, flat, tol, box):
    try:
        r = characterize_detector(bias, flat, signal_floor=0.0,
                                  box_radius=box, tolerance=tol)
        return r.gain_eff
    except Exception:  # noqa: BLE001 — 한 갈래가 실패해도 나머지로 흩어짐을 본다
        return None


def repeatability(bias, flat, tol, box) -> tuple[float, int]:
    """자료를 여러 방식으로 갈라 다시 재고 흩어짐을 돌려준다.

    이것은 **정밀도**다. 좁은 신호 구간에 공통으로 걸린 계통은 안 잡힌다.
    """
    ests = [g for g in (
        _gain(bias[: len(bias) // 2], flat, tol, box),
        _gain(bias[len(bias) // 2:], flat, tol, box),
        _gain(bias, flat[0::2], min(tol * 2, 0.2), box),
        _gain(bias, flat[1::2], min(tol * 2, 0.2), box),
        _gain(bias, flat, tol * 0.7, box),
        _gain(bias, flat, tol * 1.7, box),
        _gain(bias, flat, tol * 2.5, box),
    ) if g]
    if len(ests) < 3:
        return float("nan"), len(ests)
    a = np.array(ests)
    return float(a.std(ddof=1)), len(a)


def verdict(diff: float, err: float,
            sigma_k: float = SIGMA_K,
            matters_frac: float = MATTERS_FRAC) -> tuple[str, str]:
    """diff·err 은 헤더값에 대한 비율. (판정, 설명)"""
    significant = abs(diff) > sigma_k * err if err > 0 and np.isfinite(err) else True
    matters = abs(diff) > matters_frac
    if not significant and not matters:
        return ("헤더를 쓴다",
                "실측과 구별되지 않고, 구별되더라도 오차막대가 안 움직인다.")
    if not significant and matters:
        return ("PTC 를 다시 본다",
                "실측 오차 안이지만 차이가 문턱을 넘는다 — 실측이 너무 헐겁다.")
    if significant and not matters:
        return ("기록만 하고 헤더를 쓴다",
                "차이가 우연은 아니지만 결과를 안 움직인다.")
    return ("헤더가 무엇을 적은 것인지 밝힌다",
            "우연이 아니고 결과도 움직인다. 숫자만 바꾸지 말 것.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="헤더 GAIN 과 PTC 실측을 견준다")
    ap.add_argument("--bias", required=True, help="바이어스 폴더 또는 glob")
    ap.add_argument("--flat", required=True, help="플랫 폴더 또는 glob")
    ap.add_argument("--header-gain", type=float, default=None,
                    help="헤더값을 직접 준다 (기본: 플랫 첫 장에서 읽는다)")
    ap.add_argument("--tolerance", type=float, default=0.06,
                    help="플랫 짝 허용 폭 (하늘 플랫은 밝기가 흘러 넓게 잡는다)")
    ap.add_argument("--box-radius", type=int, default=300)
    ap.add_argument("--fast", action="store_true",
                    help="재현성 측정을 건너뛰고 인용 오차를 쓴다 (덜 믿음직하다)")
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

    if a.fast:
        sd, n_rep, how = res.gain_eff_err, 0, "인용 오차 (--fast)"
    else:
        sd, n_rep = repeatability(bias, flat, a.tolerance, a.box_radius)
        how = f"재현성 실측 ({n_rep} 회 다시 잼)"
        if not np.isfinite(sd):
            sd, how = res.gain_eff_err, "인용 오차 (재현성 측정 실패)"

    diff, err = res.gain_eff / hg - 1.0, sd / hg
    ratio = res.gain_eff / hg
    v, why = verdict(diff, err)

    print(f"헤더 GAIN        : {hg:.4f} e-/ADU")
    print(f"PTC 실측         : {res.gain_eff:.4f}  ({res.n_pairs} 쌍, "
          f"R^2 = {res.r_squared:.4f})")
    print(f"불확실도         : ±{sd:.4f}  ({err*100:.2f} %)   ← {how}")
    print(f"차이             : {diff*100:+.2f} %   ({abs(diff)/err:.1f}σ)"
          if err > 0 else f"차이             : {diff*100:+.2f} %")
    print()
    if ratio > ORDER_OF_MAGNITUDE or ratio < 1 / ORDER_OF_MAGNITUDE:
        print(f"판정 → **헤더가 무엇을 적은 것인지 밝힌다**")
        print(f"       {ratio:.1f} 배 차이는 측정 불일치가 아니다. 헤더가 다른")
        print(f"       물리량을 적고 있다 (Moravian 은 명목 레지스터 값이었다).")
        return 0

    print(f"판정 → **{v}**")
    print(f"       {why}")
    print()
    print("문턱이 어디서 왔나 —")
    print(f"  {SIGMA_K:.0f}σ  : **관례다.** 유도된 값이 아니다.")
    print(f"  {MATTERS_FRAC*100:.0f} % : 재료는 측정(오차 모형 정확도 2.4 %, 감도 0.9)"
          "이고 조합·반올림은 선택.")
    print(f"        fig2 의 어두운 끝(비 0.83)을 기준하면 17 % 가 된다.")
    print()
    print("문턱을 바꾸면 판정이 이렇게 바뀐다 —")
    for k in (2.0, 3.0, 5.0):
        row = []
        for m in (0.01, 0.03, 0.17):
            mark = " *" if (k, m) == (SIGMA_K, MATTERS_FRAC) else "  "
            row.append(f"{verdict(diff, err, k, m)[0]}{mark}")
        print(f"  {k:.0f}σ : " + " | ".join(f"{r:<24}" for r in row))
    print(f"        (열: 1 % · 3 % · 17 %,  * 가 지금 쓰는 문턱)")

    if warned:
        print()
        for line in res.log:
            if "WARNING" in line:
                print(f"  PTC 경고: {line.split('WARNING:')[-1].strip()}")
        print("  → 신호 범위가 좁으면 재현성은 계통을 못 잡는다. "
              "진짜 불확실도는 위 값보다 크다.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
