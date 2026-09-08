"""2 차 색항이 시야마다 달라지는 이유를 찾는다 — 첫 시도. **답은 여기 없다.**

> **이 스크립트는 뒤에 나온 둘로 대체됐다** (2026-09-08).
> 흔들림과 함께 움직이는 것을 상관계수로 찾으려 했는데, 조합이 여섯뿐이라
> 아무것도 못 갈랐다. **답을 낸 것은 아래 둘이다.**
>
>     validation/color_term_weighting.py   적합 방식을 넷으로 갈라, 흔들림이
>                                          적합기에서 오는지 자료에서 오는지
>     validation/color_term_transfer.py    적합에 안 쓴 별로 시험 — 2 차항이
>                                          무엇이든 사기는 하는지
>
> 결론은 `validation/COLOR_TERM_QUADRATIC.md` 에 있다. 이 파일은 「저장된 ct2
> 열을 쓰면 안 된다」는 것을 밝힌 기록으로 남긴다.

**D-012 가 답하지 않은 것이 이것이다** (2026-08-21). 그 결정은 「같은 카메라의
다섯 시야에서 R 의 ct2 가 −0.243 · −0.003 · +0.309 로 부호까지 뒤집힌다」를
보이고 「그러니 시야마다 적합해 쓰지 말자」로 닫았다. **왜 뒤집히는지는 안 밝혔다.**

같은 결정문에 B 는 반대라고 적혀 있다 — 시야를 넘어 공유한 −0.11 이 그 시야
자기 적합보다 나았다(산포 0.0373 대 0.0408). 그러니 **어떤 것은 진짜이고 어떤
것은 아닌데 무엇이 그 둘을 가르는지 모른다.** 사장님이 물은 것이 그것이다.

## 저장된 ct2 열은 쓸 수 없다

`zp_fit_coefficients.csv` 의 `ct2` 는 **적용된 값**이지 **적합된 값**이 아니다.
스물넷 중 열이 `+0.000` 인데 「0 으로 쟀다」가 아니다 — D-012 이후 실행은 적용을
안 해서 0 이고, 그 전 실행은 임의 상한 `|ct2| ≤ 0.25` 에서 걸러져 0 이다.
**처음에 이 열로 상관을 냈다가 무의미한 수를 냈다.**

그래서 여기서는 **적합을 다시 한다.** 보정성 표(`gaia_sdss_calibrator_by_ID.csv`)에
기준 등급과 기기 등급의 차이(`delta_*`)와 색(`color_*`)이 다 있으므로 그 자리에서
`delta = zp + ct·색 + ct2·색²` 을 시그마 절단으로 적합한다.
**기존 산출물은 하나도 안 건드린다 — 읽기만 한다.**

절단은 설정과 같게 맞췄다(3σ, 5 회). APEX 안의 적합기와 같은 코드는 아니므로
절대값은 소수점 아래에서 다를 수 있고, **여기서 보는 것은 시야 사이의 흔들림**이다.

## 무엇이 흔들림을 만드는지 후보

    색의 지렛대   적합에 쓴 별들의 색 범위. 짧으면 곡률이 안 잡힌다.
    별 수         적을수록 흔들린다.
    2차의 이득    2차를 넣어 잔차가 실제로 줄었는가. 안 줄면 곡률이 없다는 뜻.
    성단 종류     구상(금속 결핍) 대 산개(태양 조성). 기준이 Gaia 를 옮긴 값이고
                  그 변환은 항성 종족에 따라 다르게 틀리므로, 2 차항이 그 차이를
                  흡수할 수 있다.

**결론을 내려는 스크립트가 아니다.** 무엇이 함께 움직이는지 보이고 다음에 무엇을
재야 할지 좁히는 것이 목적이다. 시야가 일곱뿐이라 상관계수 하나로 못 정한다.

실행:
    python -X utf8 validation/color_term_why_it_varies.py
"""
from __future__ import annotations

import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOTS = ("E:/APEX_validation/reprocess/*/result",
         "E:/observed_Analysis/*/*/result")

#: 성단 종류 — 기준 변환이 항성 종족에 따라 다르게 틀리는지 보려고.
KIND = {
    "M13": "구상", "M3": "구상", "M5": "구상",
    "M67": "산개", "NGC6811": "산개", "M37": "산개", "NGC457": "산개",
}

#: 밴드마다 영점 적합이 쓰는 색지수 (APEX 규약).
BAND_COLOUR = {"g": "g_r", "r": "g_r", "i": "r_i",
               "B": "B_V", "V": "B_V", "R": "V_R", "I": "V_R"}


def field_of(path: str) -> str:
    for name in KIND:
        if re.search(rf"[\\/]{name}[\\/]", path):
            return name
    return "?"


def robust_poly(x, y, deg, clip_sigma=3.0, iters=5):
    """시그마 절단 다항 적합. 설정의 clip_sigma=3 · fit_iters=5 에 맞췄다."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x, float)[m], np.asarray(y, float)[m]
    if len(x) < deg + 5:
        return None
    keep = np.ones(len(x), bool)
    for _ in range(iters):
        c = np.polyfit(x[keep], y[keep], deg)
        r = y - np.polyval(c, x)
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        if not np.isfinite(s) or s <= 0:
            break
        new = np.abs(r - np.median(r[keep])) <= clip_sigma * s
        if new.sum() < deg + 5 or bool((new == keep).all()):
            break
        keep = new
    c = np.polyfit(x[keep], y[keep], deg)
    resid = y[keep] - np.polyval(c, x[keep])
    rms = float(np.sqrt(np.mean(resid ** 2)))
    sigma_ct2 = float("nan")
    if deg == 2:
        V = np.vander(x[keep], 3)
        try:
            cov = np.linalg.inv(V.T @ V) * rms ** 2
            sigma_ct2 = float(np.sqrt(cov[0, 0]))
        except np.linalg.LinAlgError:
            pass
    return dict(coef=c, n=int(keep.sum()), rms=rms, sigma_ct2=sigma_ct2,
                span=float(np.percentile(x[keep], 95) - np.percentile(x[keep], 5)))


def refit(result_dir: Path) -> list[dict]:
    p = result_dir / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv"
    if not p.exists():
        return []
    try:
        d = pd.read_csv(p)
    except Exception:  # noqa: BLE001 — 못 읽는 워크스페이스는 건너뛴다
        return []
    out = []
    for band, colour in BAND_COLOUR.items():
        dcol, ccol = f"delta_{band}", f"color_{colour}"
        if dcol not in d.columns or ccol not in d.columns:
            continue
        x = pd.to_numeric(d[ccol], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d[dcol], errors="coerce").to_numpy(float)
        lin, qua = robust_poly(x, y, 1), robust_poly(x, y, 2)
        if lin is None or qua is None:
            continue
        out.append(dict(band=band, colour=colour,
                        ct2=float(qua["coef"][0]), ct2_sigma=qua["sigma_ct2"],
                        ct=float(lin["coef"][0]), zp=float(lin["coef"][1]),
                        n=qua["n"], colour_span=qua["span"],
                        rms_linear=lin["rms"], rms_quad=qua["rms"]))
    return out


def main() -> int:
    rows: list[dict] = []
    for pat in ROOTS:
        for rd in sorted(glob.glob(pat)):
            fld = field_of(str(Path(rd) / "x"))
            for r in refit(Path(rd)):
                r.update(field=fld, kind=KIND.get(fld, "?"),
                         workspace=Path(rd).parent.name)
                rows.append(r)
    d = pd.DataFrame(rows)
    if d.empty:
        print("적합할 자료를 못 찾았다.")
        return 1
    d = d.drop_duplicates(subset=["field", "band", "colour"], keep="first")
    print(f"다시 적합한 것 {len(d)} 개 · 시야 {d.field.nunique()} 개")

    print()
    print("=== 같은 (밴드, 색지수)가 시야마다 얼마나 흔들리나 ===")
    head = ("밴드", "색지수", "시야", "ct2 값들", "폭", "부호뒤집힘", "지렛대", "유의도")
    print("{:<5}{:<7}{:>5}{:>34}{:>9}{:>11}{:>9}{:>9}".format(*head))
    print("-" * 89)
    groups: list[dict] = []
    for (band, col), g in d.groupby(["band", "colour"]):
        v = g["ct2"].to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v) < 2:
            continue
        sig = (g["ct2"].abs() / g["ct2_sigma"]).replace([np.inf, -np.inf], np.nan)
        row = dict(band=band, colour=col, n_fields=int(len(v)),
                   ptp=float(v.max() - v.min()),
                   sign_flip=bool(v.max() > 0 and v.min() < 0),
                   colour_span=float(np.nanmedian(g["colour_span"])),
                   median_sig=float(np.nanmedian(sig)),
                   median_n=float(np.nanmedian(g["n"])),
                   gain_mmag=float(np.nanmedian(g["rms_linear"] - g["rms_quad"]) * 1000))
        groups.append(row)
        vals = " ".join(f"{x:+.3f}" for x in sorted(v))
        flip = "예" if row["sign_flip"] else "아니오"
        print("{:<5}{:<7}{:>5}{:>34}{:>9.3f}{:>11}{:>9.2f}{:>9.1f}".format(
            band, col, len(v), vals, row["ptp"], flip,
            row["colour_span"], row["median_sig"]))

    o = pd.DataFrame(groups)
    print()
    print("=== 흔들림(폭)이 무엇과 함께 움직이나 ===")
    for key, label in (("colour_span", "색 지렛대"), ("median_sig", "유의도"),
                       ("median_n", "별 수"), ("gain_mmag", "2차로 준 이득")):
        m = np.isfinite(o["ptp"]) & np.isfinite(o[key])
        if m.sum() >= 3 and float(o.loc[m, key].std()) > 0:
            c = float(np.corrcoef(o.loc[m, "ptp"], o.loc[m, key])[0, 1])
            print(f"  ct2 폭  vs  {label:<12} 상관 {c:+.3f}  (조합 {int(m.sum())} 개)")
    print("  ※ 조합이 대여섯 개뿐이다. 상관계수 하나로 결론 내지 않는다.")

    print()
    print("=== 시야별 원값 ===")
    head2 = ("시야", "종류", "밴드", "색지수", "ct2", "±", "유의도", "N",
             "지렛대", "2차 이득(mmag)")
    print("{:<9}{:<5}{:<5}{:<7}{:>9}{:>8}{:>8}{:>6}{:>8}{:>16}".format(*head2))
    print("-" * 82)
    for _, r in d.sort_values(["band", "field"]).iterrows():
        s = (abs(r["ct2"]) / r["ct2_sigma"]
             if r["ct2_sigma"] and np.isfinite(r["ct2_sigma"]) else float("nan"))
        print("{:<9}{:<5}{:<5}{:<7}{:>+9.3f}{:>8.3f}{:>8.1f}{:>6}{:>8.2f}{:>16.2f}".format(
            r["field"], r["kind"], r["band"], r["colour"], r["ct2"], r["ct2_sigma"],
            s, int(r["n"]), r["colour_span"], (r["rms_linear"] - r["rms_quad"]) * 1000))

    outp = Path("validation/color_term_why_it_varies.json")
    outp.write_text(json.dumps({"refits": d.to_dict("records"), "groups": groups},
                               ensure_ascii=False, indent=2, default=str),
                    encoding="utf-8")
    print()
    print(f"썼다: {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
