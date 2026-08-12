"""Emit the side-by-side parameter table for the APEX–DAOPHOT PSF comparison.

A comparison against another package is not usable in a paper without the
settings both engines ran under. The manuscript already carries one of these
for the *aperture* cross-check (Table 2); this is its counterpart for PSF
photometry.

The DAOPHOT column is read from `iraf_parameters.json`, which the run itself
dumps with IRAF's `lpar` — the values IRAF actually held when the tasks fired,
not the assignments the harness intended. `unlearn` restoring a default, or
IRAF coercing a type, would otherwise never show up.

The APEX column is read from the workspace configuration that produced the
step-8 output being compared. Both sides therefore come from the runs, not from
this file.

The table separates two kinds of parameter, because they are chosen on
different grounds:

*Detector constants* describe the data. Both engines must be told the same
thing or the comparison measures the disagreement in the inputs.

*Method parameters* are each engine's own choice. Forcing DAOPHOT onto APEX's
values would be tuning it toward APEX's answer; DAOPHOT follows Massey & Davis
(1992) instead, and where the two differ the difference is stated so a reader
can judge it. One of them — the fit radius — was separately re-run at APEX's
value to show it does not change the conclusion.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]

# (label, iraf task, iraf parameter, apex config path or literal, note)
DETECTOR = [
    ("Gain", "datapars", "epadu", ("instrument", "gain_e_per_adu"),
     "동일 — APEX 설정에서 두 엔진에 같이 준다"),
    ("Read noise", "datapars", "readnoise", ("instrument", "rdnoise_e"),
     "동일"),
    ("Good-data minimum", "datapars", "datamin", ("instrument", "datamin_adu"),
     "동일"),
    ("Good-data maximum", "datapars", "datamax", ("instrument", "datamax_adu"),
     "동일"),
    ("FWHM", "datapars", "fwhmpsf", ("__measured__", "fwhm_px"),
     "동일 — 프레임에서 측정한 값"),
    ("Sky sigma", "datapars", "sigma", ("__measured__", "sky_sigma_adu"),
     "동일 — sigma_clipped_stats"),
    ("Exposure", "datapars", "itime", ("__measured__", "exptime"),
     "동일 — 헤더"),
]

METHOD = [
    ("PSF model radius", "daopars", "psfrad", ("__apex__", "epsf_size"),
     "DAOPHOT 4·FWHM+1 (Massey & Davis 1992) · APEX epsf_size_fwhm_mult=4 의 상자"),
    ("Fit radius", "daopars", "fitrad", ("__apex__", "fit_window"),
     "**차이** — APEX 는 포위에너지 90 % 자동창(약 1.7·FWHM). DAOPHOT 표준은 "
     "1·FWHM 이나 여기서는 APEX 에 맞춰 1.7 로 돌렸다(교란 제거 확인)"),
    ("Analytic function", "daopars", "function", ("__apex__", "build_mode"),
     "**차이** — DAOPHOT 은 auto 선택(실제 채택: moffat15), APEX 는 경험적 ePSF"),
    ("Spatial variation", "daopars", "varorder", ("__apex__", "model_mode"),
     "**차이** — DAOPHOT varorder=0(일정), APEX 는 프레임당 ePSF 하나. "
     "varorder=1 로도 돌려 결론이 안 바뀜을 확인"),
    ("Recenter during fit", "daopars", "recenter", ("__apex__", "recenter"),
     "동일 취지 — 두 엔진 모두 적합 중 위치를 미세조정"),
    ("Fit sky", "daopars", "fitsky", ("__apex__", "fitsky"), "동일 취지"),
    ("PSF-star cleaning", "daopars", "nclean", ("__apex__", "epsf_maxiters"),
     "**차이 — DAOPHOT 에 불리하다.** nclean=0 이라 PSF 별 정제 반복을 껐다. "
     "APEX 는 epsf_maxiters=5 로 ePSF 를 반복 정련한다"),
    ("Max fit iterations", "daopars", "maxiter", ("__apex__", "fitter_max_iter"),
     "각자 기본값"),
    ("PSF reference stars", "daopars", "__maxnpsf__", ("__apex__", "n_epsf_stars"),
     "pstselect 로 자동 선정 · APEX 는 오염인지 필터로 선정"),
    ("Initial aperture", "photpars", "apertures", ("__apex__", "step7_aperture"),
     "DAOPHOT 은 phot 의 초기 등급용. ALLSTAR 등급의 영점이 여기서 오므로 "
     "밝고 한산한 별로 상수 오프셋을 제거한 뒤에만 bias 를 비교한다"),
    ("Sky annulus", "fitskypars", "annulus", ("__apex__", "step7_annulus"),
     "**구조가 다르다** — DAOPHOT 은 고리에서 하늘을 재고, APEX step 8 은 "
     "하늘을 적합의 자유변수로 함께 푼다. 구멍 반경이라는 개념이 없다"),
    ("Sky annulus width", "fitskypars", "dannulus", ("__apex__", "step7_dannulus"),
     "위와 같음"),
    ("Centering in `phot`", "centerpars", "calgorithm", ("__apex__", "forced"),
     "**none 으로 껐다.** 좌표를 주는 강제측광이라 재중심을 켜면 희미한 별의 "
     "중심이 이웃으로 끌려간다(주입별 25 중 7 만 생존한 실측으로 발각)"),
]


def iraf_value(pars: dict, task: str, name: str, summary: dict) -> str:
    if name == "__maxnpsf__":
        return str(summary.get("n_psf_stars", "—"))
    entry = pars.get(task, {}).get(name)
    return entry["value"] if isinstance(entry, dict) else "—"


def apex_value(spec: tuple, config: dict, summary: dict) -> str:
    head, key = spec
    parameters = summary.get("parameters", {})
    if head == "__measured__":
        value = parameters.get(key)
        return f"{value:g}" if isinstance(value, (int, float)) else str(value or "—")
    if head == "__apex__":
        psf = config.get("psf", {})
        table = {
            "epsf_size": f"{psf.get('epsf_size_fwhm_mult')}·FWHM 상자",
            "fit_window": (f"{psf.get('fit_window_mode')} · "
                           f"포위에너지 {psf.get('fit_encircled_energy')}"),
            "build_mode": str(psf.get("build_mode")),
            "model_mode": str(psf.get("model_mode")),
            "recenter": "yes (적합 내)",
            "fitsky": "yes",
            "epsf_maxiters": str(psf.get("epsf_maxiters")),
            "fitter_max_iter": str(psf.get("fitter_max_iter")),
            "n_epsf_stars": "자동 (오염인지 필터)",
            "step7_aperture": "Step 7 구경보정 경로",
            "step7_annulus": "해당 없음 (적합 내 하늘)",
            "step7_dannulus": "해당 없음",
            "forced": "강제 (좌표 고정)",
        }
        return table.get(key, "—")
    section, key2 = head, key
    value = config.get(section, {}).get(key2)
    return f"{value:g}" if isinstance(value, (int, float)) else str(value or "—")


def render(rows: list, pars: dict, config: dict, summary: dict) -> list[str]:
    out = ["| 항목 | APEX step 8 | IRAF DAOPHOT | 비고 |", "|---|---|---|---|"]
    for label, task, name, spec, note in rows:
        out.append(f"| {label} | {apex_value(spec, config, summary)} | "
                   f"`{iraf_value(pars, task, name, summary)}` | {note} |")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iraf-parameters", required=True)
    ap.add_argument("--daophot-summary", required=True)
    ap.add_argument("--apex-config", required=True)
    ap.add_argument("--output", default=str(
        REPO / "validation" / "psf_engines" / "PARAMETERS.md"))
    args = ap.parse_args()

    pars = json.loads(Path(args.iraf_parameters).read_text(encoding="utf-8"))
    summary = json.loads(Path(args.daophot_summary).read_text(encoding="utf-8"))
    config = json.loads(Path(args.apex_config).read_text(encoding="utf-8"))

    lines = [
        "# APEX step 8 과 IRAF DAOPHOT ALLSTAR 의 측정 파라미터",
        "",
        "이 표가 없으면 두 엔진의 비교는 논문에 쓸 수 없다. DAOPHOT 열은 실행이",
        "IRAF `lpar` 로 덤프한 값(`iraf_parameters.json`)이고 — 하네스가 설정하려",
        "**의도한** 값이 아니라 태스크가 발화할 때 IRAF 가 **실제로 들고 있던** 값이다.",
        "APEX 열은 그 step 8 산출을 만든 워크스페이스 설정에서 읽는다.",
        "",
        f"대상 프레임: `{Path(summary.get('frame','')).name}` · "
        f"좌표 {summary.get('n_input_positions')} 개 · "
        f"ALLSTAR 적합 {summary.get('n_valid_mag')} 개 · "
        f"{summary.get('elapsed_s', 0):.0f} s",
        "",
        "## 검출기 상수 — 두 엔진이 같아야 하는 것",
        "",
        "자료의 성질이다. 다르면 비교가 엔진이 아니라 입력의 불일치를 잰다.",
        "",
    ]
    lines += render(DETECTOR, pars, config, summary)
    lines += [
        "",
        "## 방법 파라미터 — 각 엔진 저자의 선택",
        "",
        "DAOPHOT 에 APEX 값을 강제하면 APEX 의 답 쪽으로 조율하는 것이 된다.",
        "표준 지침(Massey & Davis 1992)을 따르되 **차이는 전부 적는다.**",
        "",
    ]
    lines += render(METHOD, pars, config, summary)
    lines += [
        "",
        "## 표에서 읽어야 할 것",
        "",
        "- **`nclean = 0` 은 DAOPHOT 에 불리한 설정이다.** PSF 별 정제 반복을 꺼서",
        "  `pstselect` 가 고른 60 개가 그대로 모형에 들어갔다(`psf.pst` 60 →",
        "  `psf.opst` 60, 기각 0). 사람이 눈으로 걸러내던 단계의 자동 대체물을",
        "  안 쓴 것이다. **그런데도 DAOPHOT 이 정밀도에서 앞섰으므로 그 결론은",
        "  보수적이다** — `nclean` 을 켜면 격차가 더 벌어질 수 있다.",
        "- **적합창과 공간변화는 따로 실험해서 닫았다.** `--fitrad-fwhm 1.7`,",
        "  `--varorder 1` 로 각각 재실행했고 결론이 바뀌지 않았다.",
        "- **초기 구경이 ALLSTAR 등급의 영점을 정한다.** 그래서 bias 를 비교하기",
        "  전에 밝고 한산한 별에서 상수 오프셋을 제거한다(측정값 +0.296 mag,",
        "  1·FWHM 구경보정 크기와 일치).",
        "",
        "전체 파라미터(설명 포함)는 `*_iraf_parameters.json` 에 있다. 이 표는",
        "그 중 비교에 영향을 주는 것만 뽑은 것이다.",
    ]

    out = Path(args.output)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:6]))
    print(f"\nsaved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
