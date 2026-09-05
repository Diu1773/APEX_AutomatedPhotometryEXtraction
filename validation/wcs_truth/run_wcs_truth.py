"""알려진 좌표계를 되찾는가 — Step 5 내장 solver 의 참값 시험.

공용 체인에서 참값 시험이 없던 마지막 자리다. 지금까지 WCS 는 ASTAP·
astrometry.net 과 맞춰 보는 것뿐이었는데, 셋 다 Gaia 에 맞추므로 일치가
정확성의 증거가 되지 못한다. 여기서는 **좌표계를 우리가 정한다.**

**무엇을 하나.** CRVAL·CRPIX·CD 를 정해 좌표계를 하나 만들고, 그 좌표계로 하늘의
별을 화소 자리로 옮긴다. 그 화소 자리를 검출 목록으로, 원래 하늘 좌표를 참조
목록으로 solver 에 준다. solver 가 우리가 쓴 좌표계를 되찾아야 한다.

**정답이 규격에 있다.** Greisen & Calabretta (2002) 의 FITS WCS 와 TAN 투영이
정의를 정하므로, 어느 모형을 고르느냐 하는 문제가 없다. 별이 그 화소에 있거나
없다.

**두 가지를 쓸어 본다.**
  A. 중심 잡기 오차 — 검출 위치가 흔들릴 때 어디까지 버티는가
  B. 겨눔 힌트 오차 — 헤더의 겨눔이 얼마나 틀려도 푸는가

실행:
    .venv-deploy/Scripts/python.exe -X utf8 validation/wcs_truth/run_wcs_truth.py
    .venv-deploy/Scripts/python.exe -X utf8 validation/wcs_truth/run_wcs_truth.py \
        --sweep scatter --seeds 2 --tag quick
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.astrometry.solver import solve  # noqa: E402
from apex.utils import run_journal  # noqa: E402

OUT = Path(__file__).absolute().parent

RA0, DEC0 = 200.0, -30.0
TRUE_SCALE = 0.5           # 초각/화소
TRUE_ROLL = 23.7           # 도. 둥근 값을 피해서 0 을 가정한 버그가 드러나게 한다
NX = NY = 1024
N_STARS = 300

SCATTERS_PX = (0.05, 0.10, 0.25, 0.50, 1.00, 2.00)
HINT_ARCMIN = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0)
N_SEEDS = 5

WHY = (
    "공용 체인에서 참값 시험이 없던 마지막 자리가 Step 5 다. 지금까지 WCS 는 "
    "ASTAP·astrometry.net 과 맞춰 보는 것뿐이었는데 셋 다 Gaia 에 맞추므로 일치가 "
    "정확성의 증거가 되지 못한다. 여기서는 좌표계를 우리가 정하고 solver 가 그것을 "
    "되찾는지 본다. 정답이 FITS WCS 규격(Greisen & Calabretta 2002)에 있으므로 "
    "모형 선택이 끼어들지 않는다."
)


def _provenance() -> dict:
    import platform
    import subprocess

    def _git(*a: str) -> str:
        try:
            return subprocess.run(["git", *a], cwd=str(REPO), capture_output=True,
                                  text=True, timeout=20).stdout.strip()
        except Exception:
            return ""

    vers = {}
    for m in ("numpy", "astropy", "scipy"):
        try:
            vers[m] = __import__(m).__version__
        except Exception:
            vers[m] = "없음"
    return {"git_commit": _git("rev-parse", "HEAD"),
            "git_dirty": bool(_git("status", "--porcelain")),
            "python": sys.version.split()[0], "platform": platform.platform(),
            "packages": vers,
            "script": "validation/wcs_truth/run_wcs_truth.py"}


def true_wcs() -> WCS:
    """우리가 정한 좌표계. solver 가 이것을 되찾아야 한다."""
    w = WCS(naxis=2)
    w.wcs.crpix = [NX / 2.0 + 0.5, NY / 2.0 - 0.5]   # 정중앙을 피한다
    w.wcs.crval = [RA0, DEC0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    d = TRUE_SCALE / 3600.0
    r = np.deg2rad(TRUE_ROLL)
    w.wcs.cd = np.array([[-d * np.cos(r), d * np.sin(r)],
                         [d * np.sin(r), d * np.cos(r)]])
    return w


def one_run(scatter_px: float, hint_arcmin: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    w = true_wcs()

    # 시야 안쪽에 별을 뿌린다. 화소 자리를 먼저 잡고 하늘로 옮기면
    # 시야를 벗어나는 별이 안 생긴다.
    px = rng.uniform(30, NX - 30, N_STARS)
    py = rng.uniform(30, NY - 30, N_STARS)
    sky = w.all_pix2world(np.column_stack([px, py]), 0)

    # 검출 목록: 참 화소 자리에 중심 잡기 오차를 얹는다.
    det = np.column_stack([px, py]) + rng.normal(0.0, scatter_px, (N_STARS, 2))
    flux = rng.uniform(1000.0, 60000.0, N_STARS)

    # 겨눔 힌트를 일부러 틀리게 준다.
    ang = rng.uniform(0, 2 * np.pi)
    off = hint_arcmin / 60.0
    hint_ra = RA0 + off * np.cos(ang) / np.cos(np.deg2rad(DEC0))
    hint_dec = DEC0 + off * np.sin(ang)

    t0 = time.perf_counter()
    res = solve(
        sources_xy=det, source_flux=flux,
        gaia_ra=sky[:, 0], gaia_dec=sky[:, 1], gaia_flux=flux,
        approx_ra=hint_ra, approx_dec=hint_dec,
        approx_scale_arcsec=TRUE_SCALE,
        img_shape=(NY, NX),
    )
    dt = time.perf_counter() - t0

    row = {"scatter_px": scatter_px, "hint_arcmin": hint_arcmin, "seed": seed,
           "seconds": round(dt, 2), "converged": bool(res.converged),
           "n_matches": int(res.n_matches),
           "rms_arcsec_reported": round(float(res.rms_arcsec), 5)}

    if not res.converged or res.wcs is None:
        row.update({"scale_err_pct": None, "roll_err_deg": None,
                    "pos_med_arcsec": None, "pos_max_arcsec": None,
                    "sip_order": int(getattr(res, "sip_order", 0))})
        return row

    # 되찾은 좌표계를 참 좌표계와 견준다. 참 화소 자리를 되찾은 계로 하늘에
    # 옮겨서, 참 하늘 좌표와 얼마나 떨어지는지가 이 시험의 답이다.
    got = res.wcs.all_pix2world(np.column_stack([px, py]), 0)
    cosd = np.cos(np.deg2rad(DEC0))
    sep = np.hypot((got[:, 0] - sky[:, 0]) * cosd, got[:, 1] - sky[:, 1]) * 3600.0

    roll_err = abs(float(res.rotation_deg) - TRUE_ROLL)
    roll_err = min(roll_err, abs(roll_err - 360.0))
    row.update({
        "scale_err_pct": round(
            100.0 * (float(res.scale_arcsec_per_px) - TRUE_SCALE) / TRUE_SCALE, 5),
        "roll_err_deg": round(roll_err, 5),
        "pos_med_arcsec": round(float(np.median(sep)), 5),
        "pos_max_arcsec": round(float(np.max(sep)), 5),
        "sip_order": int(getattr(res, "sip_order", 0)),
    })
    return row


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="내장 WCS solver 의 참값 시험")
    ap.add_argument("--sweep", choices=("scatter", "hint", "both"), default="both")
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--tag", default="")
    a = ap.parse_args(argv)

    tag = ("_" + a.tag) if a.tag else ""
    jsonl = OUT / f"runs{tag}.jsonl"

    cases = []
    if a.sweep in ("scatter", "both"):
        cases += [("scatter", sc, 0.0) for sc in SCATTERS_PX]
    if a.sweep in ("hint", "both"):
        cases += [("hint", 0.10, h) for h in HINT_ARCMIN]

    settings = {
        "true_wcs": {"crval": [RA0, DEC0], "scale_arcsec": TRUE_SCALE,
                     "roll_deg": TRUE_ROLL, "crpix_offset": "정중앙에서 ±0.5 화소",
                     "ctype": "RA---TAN / DEC--TAN"},
        "field": {"nx": NX, "ny": NY, "n_stars": N_STARS},
        "cases": [{"sweep": s, "scatter_px": sc, "hint_arcmin": h}
                  for s, sc, h in cases],
        "n_seeds": a.seeds,
        "catalog": "참조 목록은 우리가 만든 참 하늘 좌표다. Gaia 를 안 쓴다.",
        "scoring": ("참 화소 자리를 되찾은 좌표계로 하늘에 옮겨 참 하늘 좌표와의 "
                    "거리를 잰다. 규격이 정답이므로 대조 상대가 필요 없다."),
        "solver_defaults": "apex.analysis.astrometry.solver.solve 의 기본값 그대로",
    }

    done = set()
    if jsonl.exists():
        for line in jsonl.open(encoding="utf-8"):
            try:
                r = json.loads(line)
                done.add((r["scatter_px"], r["hint_arcmin"], r["seed"]))
            except Exception:
                pass

    run_id = run_journal.new_run_id()
    run_journal.record_note(
        OUT, WHY + f"  [이번 실행] sweep={a.sweep} seeds={a.seeds} "
        f"산출=runs{tag}.jsonl", run_id=run_id, author="wcs-truth")
    run_journal.record_run_start(
        OUT, run_id, mode="experiment", source="validation",
        plan=[f"{s}:scatter={sc},hint={h}" for s, sc, h in cases],
        environment=_provenance())
    print(f"저널 {run_journal.journal_path(OUT)}  run={run_id}", flush=True)

    total, n = len(cases) * a.seeds, 0
    t0_all = time.perf_counter()
    with jsonl.open("a", encoding="utf-8") as fh:
        for sweep, sc, hint in cases:
            for seed in range(201, 201 + a.seeds):
                n += 1
                if (sc, hint, seed) in done:
                    continue
                row = {"sweep": sweep, **one_run(sc, hint, seed)}
                fh.write(json.dumps(row, ensure_ascii=False) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
                run_journal.record_step(
                    OUT, run_id, index=n,
                    key=f"wcs_sc{sc}_hint{hint}_seed{seed}",
                    title=f"{sweep} scatter={sc}px hint={hint}'",
                    status="ok" if row["converged"] else "failed",
                    source="validation", duration_s=row["seconds"],
                    outputs=[jsonl], settings={**settings, "this_run": row},
                    settings_scope="read_by_step")
                pos = row["pos_med_arcsec"]
                print(f"[{n:3d}/{total}] {sweep:7s} sc={sc:.2f} hint={hint:.1f}' "
                      f"seed={seed}  conv={row['converged']} "
                      f"match={row['n_matches']:3d} "
                      f"pos_med={pos if pos is not None else '-'}", flush=True)

    summarize(jsonl, tag)
    run_journal.record_run_end(OUT, run_id, success=True,
                               duration_s=time.perf_counter() - t0_all)
    return 0


def summarize(jsonl: Path, tag: str = "") -> None:
    rows = [json.loads(l) for l in jsonl.open(encoding="utf-8") if l.strip()]
    if not rows:
        print("결과 없음")
        return
    out = {"true_wcs": {"scale_arcsec": TRUE_SCALE, "roll_deg": TRUE_ROLL},
           "scatter": [], "hint": []}

    for sweep, key, label in (("scatter", "scatter_px", "중심오차(px)"),
                              ("hint", "hint_arcmin", "겨눔힌트(분각)")):
        sel_all = [r for r in rows if r.get("sweep") == sweep]
        if not sel_all:
            continue
        print()
        print("{:>14} {:>3} {:>6} {:>7} {:>11} {:>11} {:>11}".format(
            label, "n", "수렴", "매칭", "위치중앙(\")", "위치최대(\")", "배율오차(%)"))
        print("-" * 70)
        for v in sorted({r[key] for r in sel_all}):
            sel = [r for r in sel_all if r[key] == v]
            conv = [r for r in sel if r["converged"]]
            e = {key: v, "n": len(sel), "n_converged": len(conv)}
            if conv:
                e["pos_med_arcsec"] = float(np.median(
                    [r["pos_med_arcsec"] for r in conv]))
                e["pos_max_arcsec"] = float(np.max(
                    [r["pos_max_arcsec"] for r in conv]))
                e["scale_err_pct"] = float(np.median(
                    [abs(r["scale_err_pct"]) for r in conv]))
                e["roll_err_deg"] = float(np.median(
                    [r["roll_err_deg"] for r in conv]))
                e["n_matches_med"] = float(np.median(
                    [r["n_matches"] for r in conv]))
                print("{:14.2f} {:3d} {:>6} {:7.0f} {:11.4f} {:11.4f} {:11.5f}".format(
                    v, len(sel), f"{len(conv)}/{len(sel)}", e["n_matches_med"],
                    e["pos_med_arcsec"], e["pos_max_arcsec"], e["scale_err_pct"]))
            else:
                print("{:14.2f} {:3d} {:>6} {:>7} {:>11} {:>11} {:>11}".format(
                    v, len(sel), f"0/{len(sel)}", "-", "-", "-", "-"))
            out[sweep].append(e)

    (OUT / f"summary{tag}.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print()
    print("요약 저장:", OUT / ("summary" + tag + ".json"))


if __name__ == "__main__":
    raise SystemExit(main())
