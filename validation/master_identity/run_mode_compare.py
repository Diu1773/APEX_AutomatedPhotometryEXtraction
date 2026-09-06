"""local 과 hybrid 는 무엇이 다른가 — 번호를 붙이는 방식만 다른지 확인한다.

`ref_build_mode` 는 두 값을 받는다. 그런데 `apex/analysis/refbuild.py` 를 읽으면
**갈라지는 자리가 한 곳뿐**이다 — 마스터 목록을 다 만든 뒤(1396 줄) source_id 를
Gaia 번호로 갈아 끼우느냐 마느냐. 병합 알고리즘 자체는 두 모드가 같은 코드다.

그러면 세 가지를 확인해야 한다.

1. **정체성 품질이 정말 같은가.** 같은 씨앗·같은 오차로 두 모드를 돌려 마스터 줄
   수·분할·미청구가 일치하는지 본다.
2. **hybrid 는 번호를 어떻게 바꾸는가.** Gaia 목록을 채워 놓고 돌려서 매칭된 별은
   Gaia 번호를, 나머지는 음수를 받는지 본다.
3. **Gaia 목록이 없으면 hybrid 는 어떻게 되는가.** 이것이 중요하다. 뒤에서
   다중 관측일 병합기(`apex/analysis/merge/id_match.py`)가 **source_id 가 양수면
   Gaia 번호라고 해석**하기 때문이다.

실행:
    .venv-deploy/Scripts/python.exe -X utf8 validation/master_identity/run_mode_compare.py
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.refbuild import run_refbuild  # noqa: E402
from apex.utils import run_journal  # noqa: E402
from apex.utils.step_paths import step5_wcs_dir, step6_refbuild_dir  # noqa: E402

from run_identity import _build, _kwargs, _provenance, _score  # noqa: E402

OUT = Path(__file__).absolute().parent

WHY = (
    "local 과 hybrid 가 정말 번호 붙이는 방식만 다른지, 그리고 Gaia 목록이 없을 때 "
    "hybrid 가 무슨 번호를 남기는지 확인한다. 뒤의 다중 관측일 병합기가 source_id "
    "부호로 Gaia 여부를 판단하므로 남는 번호의 부호가 실제로 중요하다."
)


def _write_gaia(result_dir: Path, truth: np.ndarray, frac: float, seed: int) -> int:
    """참 별 중 frac 만큼을 Gaia 목록에 넣는다. 번호는 실제 DR3 처럼 큰 정수."""
    rng = np.random.default_rng(seed)
    n = len(truth)
    take = rng.permutation(n)[: int(round(frac * n))]
    take.sort()
    g = Table()
    g["source_id"] = (4295806720000000000 + np.arange(len(take), dtype=np.int64) * 977)
    g["ra"] = truth[take, 0]
    g["dec"] = truth[take, 1]
    g["phot_g_mean_mag"] = rng.uniform(12.0, 17.0, len(take))
    g["phot_bp_mean_mag"] = g["phot_g_mean_mag"] + 0.4
    g["phot_rp_mean_mag"] = g["phot_g_mean_mag"] - 0.4
    g.write(str(step5_wcs_dir(result_dir) / "gaia_fov.ecsv"),
            format="ascii.ecsv", overwrite=True)
    return len(take)


def _sid_stats(master: pd.DataFrame) -> dict:
    sid = pd.to_numeric(master["source_id"], errors="coerce").dropna().astype("int64")
    return {
        "n_rows": int(len(sid)),
        "n_positive": int((sid > 0).sum()),
        "n_negative": int((sid < 0).sum()),
        "max_positive": int(sid[sid > 0].max()) if (sid > 0).any() else 0,
        "looks_like_gaia_id": bool((sid > 10**12).any()),
        "is_sequential_1_to_n": bool(
            (sid > 0).all() and sorted(sid.tolist()) == list(range(1, len(sid) + 1))
        ),
    }


def _one(mode: str, gaia_frac: float, scatter_px: float, seed: int) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix="apex_mode_"))
    try:
        params, names, truth, result_dir = _build(tmp, scatter_px, seed)
        n_gaia = 0
        if gaia_frac > 0:
            n_gaia = _write_gaia(result_dir, truth, gaia_frac, seed)
        kw = _kwargs(params, names)
        kw["ref_build_mode"] = mode
        t0 = time.perf_counter()
        summary = run_refbuild(**kw)
        dt = time.perf_counter() - t0
        master = pd.read_csv(step6_refbuild_dir(result_dir) / "ref_catalog.tsv", sep="\t")
        return {
            "mode": mode, "gaia_frac": gaia_frac, "n_gaia_in_fixture": n_gaia,
            "scatter_px": scatter_px, "seed": seed, "seconds": round(dt, 2),
            "n_sources_reported": int(summary.get("n_sources", -1)),
            **_score(master, truth), **_sid_stats(master),
        }
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


CASES = [
    ("local", 0.0, "local · Gaia 목록 비어 있음"),
    ("hybrid", 0.0, "hybrid · Gaia 목록 비어 있음"),
    ("hybrid", 0.8, "hybrid · 참 별의 80 % 가 Gaia 에 있음"),
    ("local", 0.8, "local · 참 별의 80 % 가 Gaia 에 있음"),
]


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="ref_build_mode 두 값의 차이를 잰다")
    ap.add_argument("--scatter", type=float, default=0.10, help="측성 오차 (픽셀)")
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args(argv)

    jsonl = OUT / "mode_compare.jsonl"
    run_id = run_journal.new_run_id()
    run_journal.record_note(
        OUT, WHY + f"  [이번 실행] scatter={a.scatter}px seeds={a.seeds}",
        run_id=run_id, author="mode-compare")
    run_journal.record_run_start(
        OUT, run_id, mode="experiment", source="validation",
        plan=[f"{m}/gaia={f}" for m, f, _ in CASES], environment=_provenance())
    print(f"저널 {run_journal.journal_path(OUT)}  run={run_id}", flush=True)

    rows: list[dict] = []
    t0_all = time.perf_counter()
    with jsonl.open("a", encoding="utf-8") as fh:
        n = 0
        for mode, frac, label in CASES:
            for seed in range(101, 101 + a.seeds):
                n += 1
                r = _one(mode, frac, a.scatter, seed)
                r["label"] = label
                rows.append(r)
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
                fh.flush()
                run_journal.record_step(
                    OUT, run_id, index=n, key=f"{mode}_gaia{frac}_seed{seed}",
                    title=label + f" seed={seed}", status="ok",
                    source="validation", duration_s=r["seconds"], outputs=[jsonl],
                    settings={"case": label, "this_run": r},
                    settings_scope="read_by_step")
                print(f"[{n:2d}] {label:<34} seed={seed}  "
                      f"master={r['n_master']}/{r['n_true']}  "
                      f"분할={r['n_split']}  양수={r['n_positive']}  "
                      f"음수={r['n_negative']}  1..N={r['is_sequential_1_to_n']}",
                      flush=True)

    print()
    hdr = ("경우", "마스터", "분할", "미청구", "중앙편차", "양수 sid", "음수 sid",
           "Gaia 번호꼴", "1..N 순번")
    print("{:<34} {:>7} {:>5} {:>7} {:>9} {:>9} {:>9} {:>11} {:>10}".format(*hdr))
    print("-" * 112)
    for mode, frac, label in CASES:
        sel = [r for r in rows if r["mode"] == mode and r["gaia_frac"] == frac]
        if not sel:
            continue

        def med(k):
            return float(np.median([r[k] for r in sel]))

        print("{:<34} {:>7.0f} {:>5.0f} {:>7.0f} {:>9.4f} {:>9.0f} {:>9.0f} "
              "{:>11} {:>10}".format(
                  label, med("n_master"), med("n_split"), med("n_unclaimed_true"),
                  med("sep_med_arcsec"), med("n_positive"), med("n_negative"),
                  str(sel[0]["looks_like_gaia_id"]),
                  str(sel[0]["is_sequential_1_to_n"])))

    run_journal.record_run_end(OUT, run_id, success=True,
                               duration_s=time.perf_counter() - t0_all)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
