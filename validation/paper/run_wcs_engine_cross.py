"""Solve the same frames with all three WCS engines and compare the solutions.

Table 1 lists the built-in quad plate solver as self-implemented, to be checked
against ASTAP and astrometry.net. This runs the production entry point
(`apex.analysis.wcs_solve.run_wcs_solve`) three times over the same frames, once
per engine, into separate result trees so nothing is shared between them.

What is compared is the *solution*, not the runtime: for each frame we map a
grid of pixel positions through each engine's WCS and measure the angular
separation between engines, plus each engine's own Gaia residual as recorded by
the acceptance gate.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\run_wcs_engine_cross.py
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

import numpy as np
import tomllib
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import warnings
warnings.filterwarnings("ignore")

SRC = Path(r"E:\APEX_validation\reprocess\M13")
OUT = Path(r"E:\APEX_validation\wcs_engines")
ENGINES = ["internal", "astap", "astnet"]
N_FRAMES = 8
GRID = 5          # 5x5 pixel grid per frame for the engine-to-engine comparison


def load_params():
    """Same loader the headless step runners use. Paths are not taken from the
    config here — `run_wcs_solve` receives data/result/cache directories as
    arguments, so each engine writes into its own tree with identical settings.
    """
    from apex.config.parameters_cmd import read_params
    return read_params(SRC / "parameters.toml")


def solve(engine: str, files: list[str], P) -> dict:
    """Each engine gets its own copy of the frames.

    The solution is written into the FITS header of the file in *data_dir*, not
    into result_dir, so engines sharing a data_dir overwrite each other. Copying
    is the whole fix.
    """
    from apex.analysis.wcs_solve import run_wcs_solve
    rdir = OUT / engine
    if rdir.exists():
        shutil.rmtree(rdir)
    ddir = rdir / "sci"
    ddir.mkdir(parents=True, exist_ok=True)
    for f in files:
        shutil.copy2(SRC / "sci" / f, ddir / f)
    # The solver reads detections from step4_dir(result_dir); without them it
    # falls back to the configured target coordinate and never actually solves.
    # Sharing one detection set is also what we want: same sources, three engines.
    det_src = SRC / "result" / "step4_detection"
    det_dst = rdir / "step4_detection"
    det_dst.mkdir(parents=True, exist_ok=True)
    n_det = 0
    for f in files:
        for ext in (".csv", ".json"):
            p = det_src / f"detect_{f}{ext}"
            if p.exists():
                shutil.copy2(p, det_dst / p.name)
                n_det += ext == ".csv"
    print(f"  [{engine}] 프레임 {len(files)}장 + 검출 {n_det}건 복사", flush=True)
    if n_det < len(files):
        raise RuntimeError(f"검출 결과가 부족하다 ({n_det}/{len(files)}) — 풀 수 없다")
    t0 = time.perf_counter()
    res = run_wcs_solve(
        files, P, str(ddir), str(rdir), str(rdir / "cache"),
        engine=engine, use_cropped=False,
    )
    dt = time.perf_counter() - t0
    print(f"  [{engine}] {dt:.1f}s  {res.get('n_solved', '?')}/{len(files)} solved",
          flush=True)
    return {"engine": engine, "seconds": dt, "result": {
        k: v for k, v in res.items() if isinstance(v, (int, float, str, bool))}}


def compare(files: list[str]) -> dict:
    """Angular separation between engines at a grid of pixels, per frame."""
    rows = []
    for fname in files:
        wcs = {}
        for e in ENGINES:
            p = OUT / e / "sci" / fname     # the engine's own copy carries its solution
            if not p.exists():
                continue
            try:
                w = WCS(fits.getheader(p)).celestial
                if w.has_celestial:
                    wcs[e] = (w, fits.getdata(p).shape)
            except Exception:
                pass
        if len(wcs) < 2:
            rows.append({"file": fname, "engines": list(wcs), "note": "too few"})
            continue
        ny, nx = next(iter(wcs.values()))[1]
        gx, gy = np.meshgrid(np.linspace(nx * 0.1, nx * 0.9, GRID),
                             np.linspace(ny * 0.1, ny * 0.9, GRID))
        sky = {e: SkyCoord(*w.wcs_pix2world(gx.ravel(), gy.ravel(), 0), unit="deg")
               for e, (w, _) in wcs.items()}
        pair = {}
        for i, a in enumerate(ENGINES):
            for b in ENGINES[i + 1:]:
                if a in sky and b in sky:
                    d = sky[a].separation(sky[b]).to(u.arcsec).value
                    pair[f"{a}|{b}"] = {"median": float(np.median(d)),
                                        "max": float(np.max(d))}
        rows.append({"file": fname, "engines": list(wcs), "pairs": pair})
    return {"frames": rows}


def collect_qc() -> dict:
    """Each engine's own Gaia residual from the acceptance gate CSV."""
    import csv
    out = {}
    for e in ENGINES:
        f = OUT / e / "frame_wcs_qc.csv"
        if not f.exists():
            f = OUT / e / "step5_wcs" / "frame_wcs_qc.csv"
        if not f.exists():
            out[e] = None
            continue
        rows = list(csv.DictReader(f.open(encoding="utf-8")))
        def col(name):
            vals = []
            for r in rows:
                try:
                    v = float(r.get(name, "nan"))
                    if np.isfinite(v):
                        vals.append(v)
                except (TypeError, ValueError):
                    pass
            return vals
        rms = col("rms_arcsec") or col("resid_rms_px") or col("rms_px")
        nm = col("n_match") or col("n_matched")
        out[e] = {"n_rows": len(rows),
                  "rms_median": float(np.median(rms)) if rms else None,
                  "rms_max": float(np.max(rms)) if rms else None,
                  "n_match_median": float(np.median(nm)) if nm else None,
                  "columns": list(rows[0].keys()) if rows else []}
    return out


def main() -> int:
    files = sorted(p.name for p in (SRC / "sci").glob("*.fit"))[:N_FRAMES]
    print(f"프레임 {len(files)}장: {files[0]} … {files[-1]}")
    OUT.mkdir(parents=True, exist_ok=True)

    P = load_params()
    timing = {}
    for e in ENGINES:
        try:
            timing[e] = solve(e, files, P)
        except Exception as exc:
            print(f"  [{e}] 실패: {type(exc).__name__}: {exc}", flush=True)
            timing[e] = {"engine": e, "error": f"{type(exc).__name__}: {exc}"}

    cmp_ = compare(files)
    qc = collect_qc()
    payload = {"files": files, "timing": timing, "compare": cmp_, "qc": qc}
    dst = REPO / "validation" / "paper" / "data_wcs_engines"
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "engine_cross.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")

    print("\n=== 엔진 간 좌표 차이 (arcsec, 프레임 중앙값) ===")
    for r in cmp_["frames"]:
        if "pairs" not in r:
            print(f"  {r['file']}: {r.get('note')}  engines={r['engines']}")
            continue
        s = "  ".join(f"{k} {v['median']:.3f}" for k, v in r["pairs"].items())
        print(f"  {r['file'][:26]:<26} {s}")
    print("\n=== 엔진별 Gaia 잔차 ===")
    for e, v in qc.items():
        print(f"  {e:<9} {v}")
    print(f"\n[saved] {dst / 'engine_cross.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
