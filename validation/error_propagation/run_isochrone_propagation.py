"""등급 오차가 성단 파라미터에 얼마나 남는가.

RESEARCH_FRAME.md 의 「그 단계의 오차가 최종 결과에 얼마나 남는가」를 CMD 갈래에
대해 재는 실험이다. 파이프라인 전체를 돌리지 않고 **이소크론 적합 층 하나만**
본다 — 알려진 참 파라미터로 합성 성단을 만들고, 크기를 정한 측광 오차를 얹고,
다시 적합해서 참값에서 얼마나 벗어나는지 잰다.

오차값은 APEX 가 실제로 낸 흩어짐 범위에 맞췄다. 인공별 10,993 개를 되찾아 잰
등급별 MAD 가 11~12 등급에서 0.005, 17~18 등급에서 0.076 이었다
(`validation/paper/논문작업/MAG_ACCURACY_20260903.md`).

**무엇을 재고 무엇을 안 재는가.** MCMC 를 참값 근처에서 출발시킨다(기존 시험과
같은 설정). 그러므로 이 표는 「옳은 골짜기 안에서 잡음이 답을 얼마나 흔드는가」를
재고, 「격자 훑기가 엉뚱한 골짜기에 앉을 위험」은 재지 않는다. 후자는 따로 물어야
한다.

생성 모형과 적합기는 tests/test_isochrone_mcmc.py 의 것을 그대로 쓴다. 여기서
다시 만들면 두 곳이 갈라진다.

실행:  .venv-deploy/Scripts/python.exe -X utf8 validation/error_propagation/run_isochrone_propagation.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from test_isochrone_mcmc import _make_synthetic_grid, _simulate_cluster  # noqa: E402
from apex.analysis.cmd.isochrone_fitter_v2 import FitBounds, IsochroneFitterV2  # noqa: E402
from apex.analysis.cmd.isochrone_mcmc import fit_isochrone_mcmc  # noqa: E402
from apex.utils import run_journal  # noqa: E402

OUT = Path(__file__).absolute().parent
TRUTH = (9.20, 0.0, 9.60, 0.06)          # log_age, [M/H], (m-M), E(g-r)
PHOT_ERRS = (0.005, 0.010, 0.020, 0.030, 0.050, 0.080)
SEEDS = tuple(range(11, 21))             # 10 draws per error level
NAMES = ("log_age", "mh", "dm", "e_color")

WHY = (
    "RESEARCH_FRAME.md 의 「그 단계의 오차가 최종 결과에 얼마나 남는가」를 CMD "
    "갈래에 대해 잰다. 단계마다의 참값 회수 시험은 이미 있으나(주기·SYSREM·"
    "등시선) 통과와 실패만 내고, 그 오차가 최종 성단 파라미터에 얼마나 남는지를 "
    "잰 것이 없다. 파이프라인 전체가 아니라 이소크론 적합 층 하나만 돈다."
)

# 이 실험이 실제로 쓰는 값 전부. 저널에 이대로 들어간다.
SETTINGS = {
    "truth": dict(zip(NAMES, TRUTH)),
    "phot_errs": list(PHOT_ERRS),
    "seeds": list(SEEDS),
    "phot_err_basis": (
        "인공별 10,993 개를 되찾아 잰 등급별 MAD 범위 0.005~0.076 에 맞췄다 "
        "(validation/paper/논문작업/MAG_ACCURACY_20260903.md)"
    ),
    "simulate": {"n_members": 350, "n_binaries": 40, "n_field": 25,
                 "source": "tests/test_isochrone_mcmc.py::_simulate_cluster"},
    "mcmc": {"n_walkers": 24, "n_burn": 150, "n_steps": 450,
             "f_bin": 0.3, "f_field": 0.1, "seed": 2024,
             "init_from_gridscan": "truth"},
    "bounds": {"log_age": [8.6, 9.7], "metallicity": [-0.3, 0.3],
               "distance_mod": [9.0, 10.2], "extinction_gr": [0.0, 0.20]},
    "grid": "tests/test_isochrone_mcmc.py::_make_synthetic_grid (합성 격자)",
    "measures": "옳은 골짜기 안의 흔들림. 격자 훑기가 엉뚱한 골짜기에 앉을 위험은 안 잰다.",
}


def _provenance() -> dict:
    """어느 판의 코드가, 어떤 환경에서 돌았는지."""
    import platform
    import subprocess

    def _git(*args: str) -> str:
        try:
            return subprocess.run(["git", *args], cwd=str(REPO), capture_output=True,
                                  text=True, timeout=20).stdout.strip()
        except Exception:
            return ""

    versions = {}
    for name in ("numpy", "scipy", "emcee", "astropy"):
        try:
            versions[name] = __import__(name).__version__
        except Exception:
            versions[name] = "없음"
    return {
        "git_commit": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": versions,
        "script": "validation/error_propagation/run_isochrone_propagation.py",
    }


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="등급 오차 → 성단 파라미터 전파 측정")
    ap.add_argument("--f-bin", type=float, default=0.3,
                    help="적합이 가정하는 쌍성 비율 (기본 0.3). 합성의 실제 값은 40/390=0.103")
    ap.add_argument("--f-field", type=float, default=0.1,
                    help="적합이 가정하는 배경별 비율 (기본 0.1). 합성의 실제 값은 25/415=0.060")
    ap.add_argument("--phot-errs", type=float, nargs="+", default=list(PHOT_ERRS),
                    help="잴 등급 오차들")
    ap.add_argument("--tag", default="",
                    help="산출 파일 이름 꼬리표. 다른 설정의 실행이 섞이지 않게 한다")
    a = ap.parse_args(argv)

    phot_errs = tuple(a.phot_errs)
    tag = ("_" + a.tag) if a.tag else ""
    settings = {**SETTINGS,
                "phot_errs": list(phot_errs),
                "mcmc": {**SETTINGS["mcmc"], "f_bin": a.f_bin, "f_field": a.f_field},
                "fraction_note": (
                    "합성의 실제 비율은 쌍성 40/390=0.103, 배경별 25/415=0.060 이다. "
                    "적합에 주는 값이 이와 다르면 그 자체가 계통 오차의 원인이 된다."
                )}

    grid = _make_synthetic_grid()
    fitter = IsochroneFitterV2(
        "synthetic.dat", col_mh=1, col_age=2, col_g=3, col_r=4,
        col_mag=5, col_mass=6, iso_data=grid, use_bilinear=True,
    )
    bounds = FitBounds(
        log_age=(8.6, 9.7), metallicity=(-0.3, 0.3),
        distance_mod=(9.0, 10.2), extinction_gr=(0.0, 0.20),
    )

    jsonl = OUT / f"runs{tag}.jsonl"
    done = set()
    if jsonl.exists():                    # resume an interrupted sweep
        for line in jsonl.open(encoding="utf-8"):
            try:
                r = json.loads(line)
                done.add((r["phot_err"], r["seed"]))
            except Exception:
                pass
        print(f"이미 끝난 것 {len(done)} 개 — 건너뛴다", flush=True)

    total = len(phot_errs) * len(SEEDS)
    n = 0
    t_start = time.perf_counter()

    run_id = run_journal.new_run_id()
    run_journal.record_note(
        OUT, WHY + f"  [이번 실행] f_bin={a.f_bin} f_field={a.f_field} "
        f"phot_errs={list(phot_errs)} 산출=runs{tag}.jsonl",
        run_id=run_id, author="isochrone-propagation")
    run_journal.record_run_start(
        OUT, run_id, mode="experiment", source="validation",
        plan=[f"phot_err={pe}" for pe in phot_errs],
        environment=_provenance(),
    )
    print(f"저널 {run_journal.journal_path(OUT)}  run={run_id}", flush=True)
    if done:
        # 이 실행이 다시 계산하지 않는 것들. 저널이 60 회 전부를 설명하게 하려고 남긴다.
        run_journal.record_note(
            OUT,
            "이어 돌린다. 아래 조합은 이 실행이 다시 계산하지 않고 runs.jsonl 의 "
            "기존 줄을 그대로 쓴다. 같은 스크립트가 만든 것이고 파라미터는 이 파일의 "
            "상수라 같으나, 저널 기능을 붙이기 전에 돌아서 그 회차의 step 기록은 "
            "없다: " + ", ".join(f"err={pe} seed={sd}" for pe, sd in sorted(done)),
            run_id=run_id, author="isochrone-propagation")
    with jsonl.open("a", encoding="utf-8") as fh:
        for pe in phot_errs:
            for seed in SEEDS:
                n += 1
                if (pe, seed) in done:
                    continue
                obs_c, obs_m, err_c, err_m = _simulate_cluster(
                    fitter, TRUTH, phot_err=pe, seed=seed)
                t0 = time.perf_counter()
                res = fit_isochrone_mcmc(
                    obs_c, obs_m, err_c, err_m, bounds=bounds,
                    interp_fn=fitter._interpolate_isochrone,
                    imf_fn=fitter._imf_weight,
                    n_walkers=24, n_burn=150, n_steps=450,
                    init_from_gridscan=list(TRUTH),
                    f_bin=a.f_bin, f_field=a.f_field, seed=2024,
                )
                med = (res.log_age_med, res.mh_med, res.dm_med, res.e_color_med)
                row = {
                    "phot_err": pe,
                    "seed": seed,
                    "n_stars": int(obs_c.size),
                    "fit_seconds": round(time.perf_counter() - t0, 2),
                    "acceptance": round(float(res.acceptance_fraction), 4),
                }
                for name, truth, m in zip(NAMES, TRUTH, med):
                    row[f"{name}_med"] = round(float(m), 5)
                    row[f"{name}_err"] = round(float(m) - truth, 5)
                fh.write(json.dumps(row, ensure_ascii=False) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
                run_journal.record_step(
                    OUT, run_id, index=n, key=f"fit_err{pe}_seed{seed}",
                    title=f"phot_err={pe} seed={seed}", status="ok",
                    source="validation", duration_s=row["fit_seconds"],
                    outputs=[jsonl], settings={**settings, "this_fit": row},
                    settings_scope="read_by_step",
                )
                el = time.perf_counter() - t_start
                print(f"[{n:3d}/{total}] err={pe:.3f} seed={seed} "
                      f"Δage={row['log_age_med'] - TRUTH[0]:+.4f} "
                      f"Δdm={row['dm_med'] - TRUTH[2]:+.4f} "
                      f"({el/60:.1f}분 경과)", flush=True)

    summarize(jsonl, tag)
    run_journal.record_run_end(
        OUT, run_id, success=True,
        duration_s=time.perf_counter() - t_start,
    )
    return 0


def summarize(jsonl: Path, tag: str = "") -> None:
    rows = [json.loads(l) for l in jsonl.open(encoding="utf-8") if l.strip()]
    if not rows:
        print("결과 없음")
        return
    out = {"truth": dict(zip(NAMES, TRUTH)), "levels": []}
    print()
    print(f"{'등급오차':>8} {'n':>3} " + " ".join(f"{nm+' 편차':>16}" for nm in NAMES))
    print("-" * 80)
    for pe in sorted({r["phot_err"] for r in rows}):
        sel = [r for r in rows if r["phot_err"] == pe]
        entry = {"phot_err": pe, "n": len(sel)}
        cells = []
        for nm in NAMES:
            e = np.array([r[f"{nm}_err"] for r in sel], dtype=float)
            entry[nm] = {"median": round(float(np.median(e)), 5),
                         "rms": round(float(np.sqrt(np.mean(e ** 2))), 5),
                         "max_abs": round(float(np.max(np.abs(e))), 5)}
            cells.append(f"{np.median(e):+.4f}±{np.std(e):.4f}")
        out["levels"].append(entry)
        print(f"{pe:8.3f} {len(sel):3d} " + " ".join(f"{c:>16}" for c in cells))
    (OUT / f"summary{tag}.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n요약 저장: {OUT / ('summary' + tag + '.json')}")


if __name__ == "__main__":
    raise SystemExit(main())
