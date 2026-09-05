"""같은 별이 프레임을 건너 같은 번호를 받는가 — Step 6 마스터 목록의 참값 시험.

서론이 APEX 의 심장이라고 적은 것이 이것이다. 「한 프레임에서 생성한 검출목록보다
여러 프레임을 지나면서 유지되는 source identity 가 분석의 기준이 된다.」
그런데 그 주장을 참값으로 확인한 시험이 없었다. 기존 Step 6 시험 넷은 좌표가
유한한가, 중복이 제거되는가, 한 프레임짜리가 걸러지는가 같은 **작동**만 본다.

**무엇을 재나.** 하늘 좌표를 아는 별을 정해 놓고, 서로 다른 곳을 겨눈 프레임 여러
장에 그 별들을 넣고, 측성 오차를 얹은 뒤, 마스터 목록이 별 하나에 번호 하나를
주는지 본다. 참값이 우리 손에 있으므로 대조 상대가 필요 없다.

**왜 이소크론이 아니라 여기인가.** 이소크론 적합은 어느 등시선 모형을 쓰느냐에
답이 딸려 온다. 여기는 그런 선택이 없다 — 별이 거기 있거나 없다.

**모드.** ref_build_mode="local" 로 돈다. hybrid 는 Gaia 목록과 맞춰 번호를 주므로
프레임 사이 병합이 아니라 Gaia 매칭을 재게 된다.

실행:
    .venv-deploy/Scripts/python.exe -X utf8 validation/master_identity/run_identity.py
    .venv-deploy/Scripts/python.exe -X utf8 validation/master_identity/run_identity.py \
        --scatters 0.05 0.5 --seeds 2 --tag quick
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.refbuild import run_refbuild  # noqa: E402
from apex.config.parameters_cmd import read_params  # noqa: E402
from apex.utils import run_journal  # noqa: E402
from apex.utils.step_paths import (  # noqa: E402
    step4_dir, step5_wcs_dir, step6_refbuild_dir,
)

OUT = Path(__file__).absolute().parent
EXAMPLE_TOML = REPO / "parameters.example.toml"

RA0, DEC0 = 200.0, -30.0
PIX_ARCSEC = 0.5          # 실제 자료에 가까운 값 (M13 CDK 0.393, LCO QHY600 0.744)
PIX_DEG = PIX_ARCSEC / 3600.0
NX = NY = 1024            # 8.5 분각 시야
N_TRUE = 400              # 산개성단 정도의 밀도
N_FRAMES = 8
DITHER_PX = 30.0          # 프레임마다 겨눔이 이만큼 어긋난다
SCATTERS_PX = (0.05, 0.10, 0.25, 0.50, 1.00, 2.00)
N_SEEDS = 5

WHY = (
    "서론이 APEX 의 심장이라고 적은 것 — 같은 별이 프레임을 건너 같은 번호를 "
    "받는가 — 을 참값으로 확인한 시험이 없었다. 기존 Step 6 시험 넷은 작동만 본다. "
    "하늘 좌표를 아는 별을 서로 다른 곳을 겨눈 프레임에 넣고 측성 오차를 얹어, "
    "마스터 목록이 별 하나에 번호 하나를 주는지 잰다. 이소크론과 달리 모형 선택이 "
    "끼어들지 않는다."
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
    for m in ("numpy", "pandas", "astropy", "scipy"):
        try:
            vers[m] = __import__(m).__version__
        except Exception:
            vers[m] = "없음"
    return {"git_commit": _git("rev-parse", "HEAD"),
            "git_dirty": bool(_git("status", "--porcelain")),
            "python": sys.version.split()[0], "platform": platform.platform(),
            "packages": vers,
            "script": "validation/master_identity/run_identity.py"}


def _frame_wcs(dx_px: float, dy_px: float, roll_deg: float) -> WCS:
    """프레임 하나의 좌표계. 겨눔을 흔들고 살짝 돌린다."""
    w = WCS(naxis=2)
    w.wcs.crpix = [NX / 2.0 + dx_px, NY / 2.0 + dy_px]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    r = np.deg2rad(roll_deg)
    w.wcs.cd = np.array([[-PIX_DEG * np.cos(r), PIX_DEG * np.sin(r)],
                         [PIX_DEG * np.sin(r), PIX_DEG * np.cos(r)]])
    return w


def _build(tmp: Path, scatter_px: float, seed: int):
    """참 별 목록과, 그 별들을 담은 합성 Step 4·5 산출물을 만든다."""
    rng = np.random.default_rng(seed)
    data_dir, result_dir = tmp / "data", tmp / "result"
    cache_dir = result_dir / "cache"
    s4, s5 = step4_dir(result_dir), step5_wcs_dir(result_dir)
    for d in (data_dir, result_dir, cache_dir, s4, s5):
        d.mkdir(parents=True, exist_ok=True)

    # 참값: 시야 안쪽에 고르게 뿌린 별. 이 좌표가 정답이다.
    half = 0.40 * NX * PIX_DEG        # 가장자리는 비워 둔다 (겨눔이 흔들리므로)
    t_ra = RA0 + rng.uniform(-half, half, N_TRUE) / np.cos(np.deg2rad(DEC0))
    t_dec = DEC0 + rng.uniform(-half, half, N_TRUE)
    truth = np.column_stack([t_ra, t_dec])

    # Gaia 목록은 비운다 — local 모드라 안 쓰지만 파일은 있어야 한다.
    g = Table()
    g["source_id"] = np.array([], dtype=np.int64)
    for c in ("ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"):
        g[c] = np.array([], dtype=float)
    g.write(str(s5 / "gaia_fov.ecsv"), format="ascii.ecsv", overwrite=True)

    names, file_path_map, rows = [], {}, []
    for fi in range(N_FRAMES):
        w = _frame_wcs(rng.uniform(-DITHER_PX, DITHER_PX),
                       rng.uniform(-DITHER_PX, DITHER_PX),
                       rng.uniform(-0.5, 0.5))
        xy = w.all_world2pix(truth, 0)
        on = ((xy[:, 0] >= 8) & (xy[:, 0] < NX - 8)
              & (xy[:, 1] >= 8) & (xy[:, 1] < NY - 8))
        dxy = xy[on] + rng.normal(0.0, scatter_px, (int(on.sum()), 2))
        n = len(dxy)

        name = f"synthetic_V_20250115_{fi:03d}.fits"
        names.append(name)
        hdr = w.to_header()
        hdr["EXPTIME"], hdr["FILTER"] = 30.0, "V"
        hdr["GAIN"], hdr["RDNOISE"] = 1.5, 8.0
        img = (rng.standard_normal((NY, NX)) * 5.0 + 100.0).astype(np.float32)
        p = data_dir / name
        fits.PrimaryHDU(img, header=hdr).writeto(p, overwrite=True)
        file_path_map[name] = str(p)

        pd.DataFrame({
            "det_uid": np.arange(n, dtype=int),
            "x": dxy[:, 0], "y": dxy[:, 1],
            "elongation": np.full(n, 1.1), "roundness": np.full(n, 0.05),
            "sharpness": np.full(n, 0.5),
            "dao_flux": rng.uniform(1000.0, 50000.0, n),
            "dao_peak": rng.uniform(500.0, 30000.0, n),
            "peak_adu": rng.uniform(500.0, 30000.0, n),
            "fwhm_px": np.full(n, 3.2),
            "anchor_candidate": np.ones(n, dtype=bool),
            "apcorr_candidate": np.ones(n, dtype=bool),
        }).to_csv(s4 / f"detect_{name}.csv", index=False)
        (s4 / f"detect_{name}.json").write_text(json.dumps({
            "filter": "V", "fwhm_px": 3.2, "n_sources": n, "sat_star_count": 0,
            "median_elongation": 1.1, "median_roundness": 0.05,
            "sky_med": 100.0, "sky_sigma": 5.0}), encoding="utf-8")
        rows.append({"file": name, "wcs_ok": True, "match_n": n, "n_match": n,
                     "n_catalog_in_fov": N_TRUE, "match_rate": 0.95,
                     "match_rate_cat": 0.9, "match_rate_eff": 0.95,
                     "resid_med": 0.3, "resid_max": 0.8, "rms_px": 0.2,
                     "wcs_qc_pass": True, "wcs_qc_reason": "",
                     "gaia_source": "fixture"})

    pd.DataFrame(rows).to_csv(s5 / "wcs_solve_summary.csv", index=False)
    params = read_params(EXAMPLE_TOML)
    params.P.data_dir, params.P.result_dir = data_dir, result_dir
    params.P.cache_dir, params.P.file_path_map = cache_dir, file_path_map
    return params, names, truth, result_dir


def _kwargs(params, names):
    return dict(
        params=params, data_dir=params.P.data_dir, result_dir=params.P.result_dir,
        cache_dir=params.P.cache_dir, file_list=names, ref_filter="V",
        sat_drop_pct=20.0, elong_drop_pct=20.0,
        ref_cat_max_sources=0, ref_cat_min_sources=20,
        ref_cat_max_elong=1.5, ref_cat_max_abs_round=0.4,
        ref_cat_sharp_min=0.2, ref_cat_sharp_max=1.0, ref_cat_min_peak_adu=0.0,
        wcs_match_radius_arcsec=2.0, wcs_min_match_rate=0.2, wcs_min_match_n=10,
        wcs_max_sep_med_arcsec=1.5, wcs_max_sep_p90_arcsec=2.5,
        wcs_max_dup_rate=0.5,
        ref_per_date=True, ref_master_union=True, ref_union_min_frames=1,
        ref_build_mode="local", gaia_mag_limit=18.0,
    )


def _score(master: pd.DataFrame, truth: np.ndarray) -> dict:
    """마스터 목록의 각 줄을 가장 가까운 참 별에 붙이고 정체성을 센다."""
    m_ra = pd.to_numeric(master["ra_deg"], errors="coerce").to_numpy(float)
    m_dec = pd.to_numeric(master["dec_deg"], errors="coerce").to_numpy(float)
    ok = np.isfinite(m_ra) & np.isfinite(m_dec)
    m_ra, m_dec = m_ra[ok], m_dec[ok]

    cosd = np.cos(np.deg2rad(DEC0))
    d_ra = (m_ra[:, None] - truth[None, :, 0]) * cosd
    d_dec = m_dec[:, None] - truth[None, :, 1]
    dist = np.hypot(d_ra, d_dec) * 3600.0

    nearest = np.argmin(dist, axis=1)
    sep = dist[np.arange(len(nearest)), nearest]
    claims = np.bincount(nearest, minlength=len(truth))

    return {
        "n_master": int(len(m_ra)),
        "n_true": int(len(truth)),
        "n_split": int((claims >= 2).sum()),           # 참 별 하나가 여러 줄로 갈렸다
        "n_extra_rows": int((claims[claims >= 2] - 1).sum()),
        "n_unclaimed_true": int((claims == 0).sum()),  # 아무 줄도 이 별을 안 가리킨다
        "sep_med_arcsec": round(float(np.median(sep)), 4),
        "sep_p90_arcsec": round(float(np.percentile(sep, 90)), 4),
        "sep_max_arcsec": round(float(np.max(sep)), 4),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="마스터 목록 정체성의 참값 시험")
    ap.add_argument("--scatters", type=float, nargs="+", default=list(SCATTERS_PX),
                    help="얹을 측성 오차 (픽셀)")
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--n-stars", type=int, default=N_TRUE,
                    help="참 별 개수. 밀도 축을 넓힐 때 쓴다 (기본 400)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args(argv)

    global N_TRUE
    N_TRUE = int(a.n_stars)

    tag = ("_" + a.tag) if a.tag else ""
    jsonl = OUT / f"runs{tag}.jsonl"
    settings = {
        "field": {"nx": NX, "ny": NY, "pix_arcsec": PIX_ARCSEC,
                  "ra0": RA0, "dec0": DEC0},
        "truth": {"n_stars": N_TRUE, "layout": "시야 안쪽 80 % 에 균일",
                  "mean_separation_arcsec": round(
                      (0.80 * NX * PIX_ARCSEC) / max(np.sqrt(N_TRUE), 1.0), 2)},
        "frames": {"n": N_FRAMES, "dither_px": DITHER_PX, "roll_deg": "+-0.5"},
        "scatters_px": list(a.scatters),
        "n_seeds": a.seeds,
        "refbuild": {"mode": "local", "union": True, "min_frames": 1,
                     "match_radius_arcsec": 2.0},
        "scoring": ("마스터 각 줄을 가장 가까운 참 별에 붙인다. 한 참 별을 두 줄 "
                    "이상이 가리키면 분할이고, 아무 줄도 안 가리키면 그 별은 "
                    "목록에 없거나 다른 별에 병합된 것이다."),
    }

    done = set()
    if jsonl.exists():
        for line in jsonl.open(encoding="utf-8"):
            try:
                r = json.loads(line)
                done.add((r["scatter_px"], r["seed"]))
            except Exception:
                pass

    run_id = run_journal.new_run_id()
    run_journal.record_note(
        OUT, WHY + f"  [이번 실행] scatters={list(a.scatters)} seeds={a.seeds} "
        f"산출=runs{tag}.jsonl", run_id=run_id, author="master-identity")
    run_journal.record_run_start(
        OUT, run_id, mode="experiment", source="validation",
        plan=[f"scatter_px={s}" for s in a.scatters], environment=_provenance())
    print(f"저널 {run_journal.journal_path(OUT)}  run={run_id}", flush=True)
    if done:
        run_journal.record_note(
            OUT, "이어 돌린다. 다시 계산하지 않는 조합: "
            + ", ".join(f"scatter={s} seed={d}" for s, d in sorted(done)),
            run_id=run_id, author="master-identity")

    total, n = len(a.scatters) * a.seeds, 0
    t0_all = time.perf_counter()
    with jsonl.open("a", encoding="utf-8") as fh:
        for sc in a.scatters:
            for seed in range(101, 101 + a.seeds):
                n += 1
                if (sc, seed) in done:
                    continue
                tmp = Path(tempfile.mkdtemp(prefix="apex_ident_"))
                try:
                    params, names, truth, result_dir = _build(tmp, sc, seed)
                    t0 = time.perf_counter()
                    summary = run_refbuild(**_kwargs(params, names))
                    dt = time.perf_counter() - t0
                    cat = step6_refbuild_dir(result_dir) / "ref_catalog.tsv"
                    master = pd.read_csv(cat, sep="\t")
                    row = {"scatter_px": sc, "seed": seed,
                           "seconds": round(dt, 2),
                           "n_sources_reported": int(summary.get("n_sources", -1)),
                           **_score(master, truth)}
                finally:
                    shutil.rmtree(tmp, ignore_errors=True)
                fh.write(json.dumps(row, ensure_ascii=False) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
                run_journal.record_step(
                    OUT, run_id, index=n, key=f"identity_sc{sc}_seed{seed}",
                    title=f"scatter={sc}px seed={seed}", status="ok",
                    source="validation", duration_s=row["seconds"],
                    outputs=[jsonl], settings={**settings, "this_run": row},
                    settings_scope="read_by_step")
                print(f"[{n:3d}/{total}] scatter={sc:.2f}px seed={seed}  "
                      f"master={row['n_master']}/{row['n_true']}  "
                      f"split={row['n_split']}  unclaimed={row['n_unclaimed_true']}  "
                      f"sep_med={row['sep_med_arcsec']}", flush=True)

    summarize(jsonl, tag)
    run_journal.record_run_end(OUT, run_id, success=True,
                               duration_s=time.perf_counter() - t0_all)
    return 0


def summarize(jsonl: Path, tag: str = "") -> None:
    rows = [json.loads(l) for l in jsonl.open(encoding="utf-8") if l.strip()]
    if not rows:
        print("결과 없음")
        return
    out = {"levels": []}
    header = ("측성오차(px)", "n", "마스터", "분할", "미청구", "중앙편차", "최대편차")
    print()
    print("{:>12} {:>3} {:>8} {:>7} {:>8} {:>11} {:>11}".format(*header))
    print("-" * 66)
    for sc in sorted({r["scatter_px"] for r in rows}):
        sel = [r for r in rows if r["scatter_px"] == sc]

        def col(key):
            return np.array([r[key] for r in sel], dtype=float)

        e = {"scatter_px": sc, "n": len(sel),
             "n_master_med": float(np.median(col("n_master"))),
             "n_split_med": float(np.median(col("n_split"))),
             "n_unclaimed_med": float(np.median(col("n_unclaimed_true"))),
             "sep_med_arcsec": float(np.median(col("sep_med_arcsec"))),
             "sep_max_arcsec": float(np.max(col("sep_max_arcsec")))}
        out["levels"].append(e)
        print("{:12.2f} {:3d} {:8.0f} {:7.0f} {:8.0f} {:11.4f} {:11.4f}".format(
            sc, len(sel), e["n_master_med"], e["n_split_med"],
            e["n_unclaimed_med"], e["sep_med_arcsec"], e["sep_max_arcsec"]))
    (OUT / f"summary{tag}.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print()
    print("요약 저장:", OUT / ("summary" + tag + ".json"))


if __name__ == "__main__":
    raise SystemExit(main())
