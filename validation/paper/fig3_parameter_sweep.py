"""Figure — pipeline-parameter and observing-condition sensitivity sweeps.

Full redesign (2026-08-02) of the trials=4 original. Every point is now an
artificial-star benchmark of 30 trials x 40 injected stars run through the
production code path (empirical-PSF injection -> Step-4 detection -> forced
aperture photometry, `apex.benchmark.runner.run_benchmark`), with 95% CIs from
a cluster bootstrap over trials. Twin y-axes are gone: one metric per panel,
2x3 grid with a shared x per column.

* Column 1 — aperture radius on ONE fixed frame (FWHM 6.0 px). The production
  pipeline clamps the forced-aperture radius at `min_r_ap_px` = 4 px, so on the
  old 3.5 px-FWHM frame every requested scale below ~1.1xFWHM collapsed onto
  the same 4 px radius; the old panel drew those duplicate measurements as a
  "flat" segment. The 6.0 px frame puts the whole 0.7-2.5xFWHM range above the
  floor, so each point is a distinct aperture. One deliberately clamped request
  (0.5xFWHM) is kept, plotted at its effective radius with an open marker.
  Top: photometric scatter (MAD). Bottom: position RMSE (aperture-independent).
* Column 2 — sky background: a controlled frame regenerated per level (only
  the sky changes; star field, seed and all other knobs fixed).
  Top: 50% completeness depth m50 (cluster-bootstrap CI from the production
  completeness fit). Bottom: photometric scatter.
* Column 3 — seeing (PSF FWHM): frame regenerated per value, same design.
* Companion (reported in the caption, not drawn): detection threshold
  `detect_sigma` swept 2.0-6.0 on the baseline frame.

Distilled per-point metrics land in data_parameter_sweep/ so the figure can be
rebuilt without the ~45 min of benchmark runs.

Run with the deploy venv:
    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig3_parameter_sweep.py
        [--sweep aperture|background|seeing|threshold]   # run one sweep's points
        [--figure-only]                                  # draw from distilled CSV
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

from apex.benchmark.synthetic_frame import make_synthetic_reference_frame
from apex.benchmark.runner import BenchmarkConfig, run_benchmark


# ── paths ────────────────────────────────────────────────────────────────────
PAPER_DIR = REPO / "validation" / "paper"
FIG_DIR = PAPER_DIR / "figures"
CAP_DIR = PAPER_DIR / "captions"
DATA_DIR = PAPER_DIR / "data_parameter_sweep"
# Short benchmark run root to stay under the Windows MAX_PATH (260) limit;
# NOT under a deep AppData\Local\Temp\claude\... scratchpad path.
SWEEP_ROOT = Path(r"C:\Users\bmffr\AppData\Local\Temp\apx_sweep")

# ── fixed synthetic-frame + benchmark constants ──────────────────────────────
SEED = 20260702          # frame seed AND benchmark base seed (trial t -> SEED+t)
BOOT_SEED = 20260802     # cluster-bootstrap RNG for MAD / RMSE CIs
TRIALS = 30
STARS_PER_TRIAL = 40
N_BOOT = 1000            # bootstrap resamples for MAD / RMSE
M50_BOOT = 500           # completeness-fit cluster bootstrap (production code)
PARAM_FILE = str(REPO / "parameters.toml")

# Baseline synthetic-frame parameters (held fixed except the swept condition).
# 1024 px / 435 stars keeps the star density of the old 640 px / 170-star frame
# (4.15e-4 px^-2) while giving the injection sampler room: at FWHM 6 px the
# min-separation constraint reaches ~62 px and 40 injections hit the random-
# sequential jamming limit of a 640 px frame (placement failures on some
# trial seeds — observed 2026-08-02).
FRAME_BASE = dict(
    seed=SEED,
    n_stars=435,
    size=1024,
    fwhm_px=3.5,
    zeropoint=25.0,
    background=150.0,
    read_noise=5.0,
    gain=1.5,
    mag_min=13.5,
    mag_max=19.0,
    filter_name="r",
    exptime=1.0,
)

# Aperture sweep runs on a FWHM 6.0 px frame so every scale below maps to a
# radius above the 4 px production floor (see module docstring).
APERTURE_FRAME_FWHM = 6.0
DEFAULT_APERTURE_SCALE = 0.8          # forced_r_ap_scale default
MIN_R_AP_PX = 4.0                     # production min_r_ap_px (parameters_cmd)

# Sweep grids. The seeing grid stays inside the production FWHM-QC window
# ([fwhm] px_min = 3.0): below it the frame-FWHM estimator is out of its
# envelope and returns garbage (true 2.5 px measured as 8.5 px from 5 stars,
# observed 2026-08-02), so a 2.5 px point would not be a pipeline measurement.
APERTURE_SCALES = [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 2.5]
BACKGROUNDS = [50.0, 100.0, 150.0, 225.0, 300.0, 450.0, 600.0, 800.0, 1000.0]
FWHMS = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
DETECT_SIGMAS = [2.0, 3.2, 4.0, 5.0, 6.0]     # 3.2 = production default

STARS_COLS = ["trial", "baseline_confounded", "recovered",
              "forced_mag_error", "position_error_px"]


# ── benchmark plumbing ───────────────────────────────────────────────────────
def _make_frame(path: Path, **overrides) -> Path:
    """Generate a frame, or reuse the on-disk one if its parameters match.

    Reuse matters for caching: the generator stamps a DATE header, so
    regenerating an identical frame changes the file bytes and would
    invalidate every cached benchmark point through the sha check.
    """
    params = dict(FRAME_BASE)
    params.update(overrides)
    meta_path = Path(str(path) + ".params.json")
    if path.exists() and meta_path.exists():
        try:
            if json.loads(meta_path.read_text(encoding="utf-8")) == json.loads(
                    json.dumps(params)):
                return path
        except Exception:
            pass
    out = make_synthetic_reference_frame(path, **params)
    meta_path.write_text(json.dumps(params, sort_keys=True), encoding="utf-8")
    return out


def _n_baked(frame: Path) -> int:
    """Stars actually placed in the frame (the generator can fall short of
    n_stars near the packing limit; the header records the real count)."""
    from astropy.io import fits
    with fits.open(frame, memmap=False) as hdul:
        return int(hdul[0].header.get("APEXNST", -1))


def _config(frame: Path, *, aperture_scale, detection_overrides,
            magnitude_min, magnitude_max, m50_boot) -> BenchmarkConfig:
    return BenchmarkConfig(
        input_fits=str(frame),
        parameter_file=PARAM_FILE,
        seed=SEED,
        trials=TRIALS,
        stars_per_trial=STARS_PER_TRIAL,
        magnitude_min=magnitude_min,
        magnitude_max=magnitude_max,
        magnitude_bin_width=1.0,
        completeness_bootstrap_samples=m50_boot,
        save_injected_fits=False,
        zeropoint_mag=25.0,
        isolated_fraction=0.6,
        psf_min_stars=3,
        aperture_scale_fwhm=aperture_scale,
        detection_overrides=detection_overrides,
    )


# Cache key: the config fields that change the *measurements*. The bootstrap
# sample count only changes CI availability and is deliberately excluded.
_KEY_FIELDS = ("seed", "trials", "stars_per_trial", "magnitude_min",
               "magnitude_max", "aperture_scale_fwhm", "detection_overrides",
               "isolated_fraction")


def _sha256(path: Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _cache_valid(out_dir: Path, cfg: BenchmarkConfig, frame: Path) -> bool:
    manifest = out_dir / "manifest.json"
    summary = out_dir / "summary.json"
    stars = out_dir / "stars.csv"
    if not (manifest.exists() and summary.exists() and stars.exists()):
        return False
    try:
        mani = json.loads(manifest.read_text(encoding="utf-8"))
        eff = mani["effective_config"]
    except Exception:
        return False
    for field in _KEY_FIELDS:
        if eff.get(field) != getattr(cfg, field):
            return False
    # the frame is regenerated deterministically every run; a changed frame
    # definition (size, star count, ...) must invalidate the cached point
    if mani.get("input_sha256") != _sha256(frame):
        return False
    return True


def _run_point(frame: Path, out_dir: Path, cfg: BenchmarkConfig) -> dict:
    """Run one benchmark point (or reuse a valid cached run)."""
    if not _cache_valid(out_dir, cfg, frame):
        shutil.rmtree(out_dir, ignore_errors=True)
        run_benchmark(cfg, input_override=frame, output_override=out_dir)
    summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
    summary["_run_dir"] = str(out_dir)
    return summary


# ── cluster bootstrap (over trials) for MAD / position RMSE ──────────────────
# Statistics replicate apex.benchmark.metrics.summarize_results exactly:
#   MAD  : 1.4826 * median|x - median(x)| of forced_mag_error,
#          eligible (= not baseline-confounded) & finite
#   RMSE : sqrt(mean(position_error_px^2)), eligible & recovered & finite
def _point_stats(run_dir: Path) -> dict:
    stars = pd.read_csv(Path(run_dir) / "stars.csv", usecols=STARS_COLS)
    eligible = ~stars["baseline_confounded"].astype(bool)
    recovered = stars["recovered"].astype(bool)
    mag = pd.to_numeric(stars["forced_mag_error"], errors="coerce")
    pos = pd.to_numeric(stars["position_error_px"], errors="coerce")

    mad_by_trial, rmse_by_trial = [], []
    for _t, idx in stars.groupby("trial").groups.items():
        sel = stars.index.isin(idx)
        mad_by_trial.append(mag[sel & eligible & np.isfinite(mag)].to_numpy(float))
        rmse_by_trial.append(
            pos[sel & eligible & recovered & np.isfinite(pos)].to_numpy(float))

    def mad_of(values: np.ndarray) -> float:
        if values.size == 0:
            return np.nan
        return float(1.4826 * np.median(np.abs(values - np.median(values))))

    def rmse_of(values: np.ndarray) -> float:
        if values.size == 0:
            return np.nan
        return float(np.sqrt(np.mean(values ** 2)))

    rng = np.random.default_rng(BOOT_SEED)
    n_tr = len(mad_by_trial)
    boots_mad, boots_rmse = [], []
    for _ in range(N_BOOT):
        pick = rng.integers(0, n_tr, size=n_tr)
        boots_mad.append(mad_of(np.concatenate([mad_by_trial[i] for i in pick])))
        boots_rmse.append(rmse_of(np.concatenate([rmse_by_trial[i] for i in pick])))
    boots_mad = np.asarray(boots_mad, float)
    boots_rmse = np.asarray(boots_rmse, float)

    out = {
        "mad": mad_of(np.concatenate(mad_by_trial)),
        "rmse": rmse_of(np.concatenate(rmse_by_trial)),
        "n_mag": int(sum(v.size for v in mad_by_trial)),
        "n_pos": int(sum(v.size for v in rmse_by_trial)),
    }
    for name, boots in (("mad", boots_mad), ("rmse", boots_rmse)):
        good = boots[np.isfinite(boots)]
        out[f"{name}_lo"] = float(np.percentile(good, 2.5)) if good.size else np.nan
        out[f"{name}_hi"] = float(np.percentile(good, 97.5)) if good.size else np.nan
    return out


# ── sweeps ───────────────────────────────────────────────────────────────────
def run_aperture_sweep() -> pd.DataFrame:
    frame = _make_frame(SWEEP_ROOT / "ap_ref.fits", fwhm_px=APERTURE_FRAME_FWHM)
    rows = []
    for scale in APERTURE_SCALES:
        tag = f"ap_{int(round(scale * 100)):04d}"
        cfg = _config(frame, aperture_scale=scale, detection_overrides={},
                      magnitude_min=14.5, magnitude_max=19.0, m50_boot=0)
        summary = _run_point(frame, SWEEP_ROOT / tag, cfg)
        stats = _point_stats(summary["_run_dir"])
        r_ap = float(summary["aperture_radii_px"]["aperture"])
        fwhm_det = float(summary["fwhm_px"])
        rows.append({
            "sweep": "aperture", "requested_scale": scale,
            "r_ap_px": r_ap, "fwhm_det_px": fwhm_det,
            "clamped": bool(r_ap > scale * fwhm_det + 1e-9),
            "n_baked": _n_baked(frame), **stats,
        })
        print(f"[aperture] scale={scale:.2f} r_ap={r_ap:.2f}px "
              f"clamped={rows[-1]['clamped']} MAD={stats['mad']:.4f} "
              f"[{stats['mad_lo']:.4f},{stats['mad_hi']:.4f}] "
              f"RMSE={stats['rmse']:.3f}px", flush=True)
    return pd.DataFrame(rows)


def _condition_row(summary: dict, stats: dict, frame: Path) -> dict:
    fit = summary.get("completeness_fit") or {}
    return {
        "m50": fit.get("m50"),
        "m50_lo": fit.get("m50_ci95_low"),
        "m50_hi": fit.get("m50_ci95_high"),
        "false_detections": summary.get("new_false_detections"),
        "r_ap_px": float(summary["aperture_radii_px"]["aperture"]),
        "fwhm_det_px": float(summary["fwhm_px"]),
        "n_baked": _n_baked(frame),
        **stats,
    }


def run_background_sweep() -> pd.DataFrame:
    rows = []
    for bkg in BACKGROUNDS:
        tag = f"bkg_{int(round(bkg)):05d}"
        frame = _make_frame(SWEEP_ROOT / f"{tag}_ref.fits", background=bkg)
        cfg = _config(frame, aperture_scale=None, detection_overrides={},
                      magnitude_min=14.0, magnitude_max=20.0, m50_boot=M50_BOOT)
        summary = _run_point(frame, SWEEP_ROOT / tag, cfg)
        stats = _point_stats(summary["_run_dir"])
        rows.append({"sweep": "background", "background": bkg,
                     **_condition_row(summary, stats, frame)})
        print(f"[background] sky={bkg:.0f} ADU m50={rows[-1]['m50']:.3f} "
              f"[{rows[-1]['m50_lo']:.3f},{rows[-1]['m50_hi']:.3f}] "
              f"MAD={stats['mad']:.4f}", flush=True)
    return pd.DataFrame(rows)


def run_seeing_sweep() -> pd.DataFrame:
    rows = []
    for fwhm in FWHMS:
        tag = f"fwhm_{int(round(fwhm * 10)):03d}"
        frame = _make_frame(SWEEP_ROOT / f"{tag}_ref.fits", fwhm_px=fwhm)
        cfg = _config(frame, aperture_scale=None, detection_overrides={},
                      magnitude_min=14.0, magnitude_max=20.0, m50_boot=M50_BOOT)
        summary = _run_point(frame, SWEEP_ROOT / tag, cfg)
        stats = _point_stats(summary["_run_dir"])
        rows.append({"sweep": "seeing", "fwhm_px": fwhm,
                     **_condition_row(summary, stats, frame)})
        print(f"[seeing] fwhm={fwhm:.1f}px m50={rows[-1]['m50']:.3f} "
              f"[{rows[-1]['m50_lo']:.3f},{rows[-1]['m50_hi']:.3f}] "
              f"MAD={stats['mad']:.4f}", flush=True)
    return pd.DataFrame(rows)


def _sigma_actually_used(run_dir: Path) -> float:
    """The sigma the detector really applied (from the baseline detect json).

    Guard against silent no-op overrides: the live parameters.toml defines
    per-filter sigmas ([detection.sigma_by_filter] / detect_sigma_r) which mask
    a bare `detect_sigma` override — observed 2026-08-02: all five "swept"
    thresholds ran at 3.2 and produced bit-identical results.
    """
    js = sorted(Path(run_dir).glob("baseline/step4_detection/detect_*.json"))
    return float(json.loads(js[0].read_text(encoding="utf-8"))["sigma_used"])


def run_threshold_sweep() -> pd.DataFrame:
    frame = _make_frame(SWEEP_ROOT / "thr_ref.fits")
    rows = []
    for sig in DETECT_SIGMAS:
        tag = f"thr_{int(round(sig * 10)):03d}"
        # override the per-filter key too (the synthetic frame is filter "r");
        # otherwise detect_sigma_r from parameters.toml masks the sweep value
        cfg = _config(frame, aperture_scale=None,
                      detection_overrides={"detect_sigma": sig,
                                           "detect_sigma_r": sig},
                      magnitude_min=14.0, magnitude_max=20.0, m50_boot=M50_BOOT)
        summary = _run_point(frame, SWEEP_ROOT / tag, cfg)
        stats = _point_stats(summary["_run_dir"])
        used = _sigma_actually_used(summary["_run_dir"])
        if abs(used - sig) > 1e-6:
            raise RuntimeError(
                f"threshold sweep is a no-op: requested {sig}, detector used {used}")
        rows.append({"sweep": "threshold", "detect_sigma": sig,
                     "sigma_used": used,
                     **_condition_row(summary, stats, frame)})
        print(f"[threshold] sigma={sig:.1f} (used {used:.1f}) "
              f"m50={rows[-1]['m50']:.3f} "
              f"false={rows[-1]['false_detections']}", flush=True)
    return pd.DataFrame(rows)


SWEEPS = {
    "aperture": run_aperture_sweep,
    "background": run_background_sweep,
    "seeing": run_seeing_sweep,
    "threshold": run_threshold_sweep,
}


# ── figure ───────────────────────────────────────────────────────────────────
def _err(lo, hi, y):
    y = np.asarray(y, float)
    return np.vstack([y - np.asarray(lo, float), np.asarray(hi, float) - y])


def make_figure_legacy(ap: pd.DataFrame, bg: pd.DataFrame, se: pd.DataFrame,
                       thr: pd.DataFrame) -> dict:
    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE_COL, 4.9), sharex="col")
    (ax_a, ax_b, ax_c), (ax_d, ax_e, ax_f) = axes
    TITLE_KW = dict(loc="left", fontsize=8.4)

    # ── column 1: aperture radius (fixed frame, FWHM 6 px) ──────────────────
    fwhm_det = float(ap["fwhm_det_px"].iloc[0])
    solid = ap[~ap["clamped"]]
    clamp = ap[ap["clamped"]]
    x_s = solid["r_ap_px"] / fwhm_det
    x_c = clamp["r_ap_px"] / fwhm_det

    opt = solid.loc[solid["mad"].idxmin()]
    x_opt = float(opt["r_ap_px"]) / fwhm_det

    for ax in (ax_a, ax_d):
        ax.axvspan(0.5, MIN_R_AP_PX / fwhm_det, color=PALETTE["grey"], alpha=0.15,
                   lw=0)
        ax.axvline(DEFAULT_APERTURE_SCALE, color=PALETTE["grey"], ls="--", lw=0.9)
    ax_a.text(MIN_R_AP_PX / fwhm_det - 0.02, 0.185, "4 px floor",
              fontsize=6.0, color=PALETTE["grey"], ha="right", va="top",
              rotation=90)
    ax_a.text(DEFAULT_APERTURE_SCALE + 0.02, 0.185,
              f"default {DEFAULT_APERTURE_SCALE:g}", fontsize=6.0,
              color=PALETTE["grey"], ha="left", va="top", rotation=90)

    ax_a.errorbar(x_s, solid["mad"], yerr=_err(solid["mad_lo"], solid["mad_hi"],
                                               solid["mad"]),
                  fmt="o-", ms=3.4, lw=1.1, capsize=2, color=C["data"])
    ax_a.errorbar(x_c, clamp["mad"], yerr=_err(clamp["mad_lo"], clamp["mad_hi"],
                                               clamp["mad"]),
                  fmt="o", ms=4.2, mfc="none", color=C["data"])
    for _i, row in clamp.iterrows():
        ax_a.annotate(f"{row['requested_scale']:g}$\\times$ requested\n"
                      "$\\to$ clamped to 4 px",
                      (row["r_ap_px"] / fwhm_det, row["mad"]),
                      textcoords="offset points", xytext=(-6, 30), fontsize=6.0,
                      color=C["data"], ha="left",
                      arrowprops=dict(arrowstyle="-", lw=0.6, color=C["data"]))
    ax_a.annotate(f"min {opt['mad']:.3f} mag at {x_opt:.2f}$\\times$FWHM",
                  (x_opt, float(opt["mad"])), textcoords="offset points",
                  xytext=(24, -6), fontsize=6.4, color=C["model"],
                  arrowprops=dict(arrowstyle="-", lw=0.7, color=C["model"]))
    ax_a.set_ylim(0.03, 0.22)
    ax_a.set_ylabel("scatter MAD (mag)")
    ax_a.set_title(f"(a) scatter min at {x_opt:.2f}$\\times$FWHM", **TITLE_KW)

    rmse_med = float(np.median(solid["rmse"]))
    ax_d.errorbar(x_s, solid["rmse"], yerr=_err(solid["rmse_lo"],
                                                solid["rmse_hi"], solid["rmse"]),
                  fmt="s-", ms=3.0, lw=1.1, capsize=2, color=C["accent"])
    ax_d.errorbar(x_c, clamp["rmse"], yerr=_err(clamp["rmse_lo"],
                                                clamp["rmse_hi"], clamp["rmse"]),
                  fmt="s", ms=3.8, mfc="none", color=C["accent"])
    ax_d.axhline(rmse_med, color=PALETTE["grey"], lw=0.7, ls=":")
    ax_d.text(0.98, 0.90, f"median {rmse_med:.2f} px",
              transform=ax_d.transAxes, fontsize=6.4, ha="right",
              color=PALETTE["grey"])
    ax_d.set_xlabel("aperture radius / FWHM")
    ax_d.text(0.03, 0.06, f"frame FWHM {fwhm_det:.1f} px",
              transform=ax_d.transAxes, fontsize=6.0, color=PALETTE["grey"])
    ax_d.set_ylabel("position RMSE (px)")
    ax_d.set_ylim(0, max(0.6, 1.6 * float(solid["rmse"].max())))
    ax_d.set_title("(b) astrometry: aperture-independent", **TITLE_KW)

    # ── columns 2-3: observing conditions ───────────────────────────────────
    def depth_panel(ax, x, d, title, logx=False, end_off=(-8, -15)):
        ax.errorbar(x, d["m50"], yerr=_err(d["m50_lo"], d["m50_hi"], d["m50"]),
                    fmt="o-", ms=3.4, lw=1.1, capsize=2, color=C["data"])
        if logx:
            ax.set_xscale("log")
        # normal axis: larger m50 (fainter limit, deeper) plots upward
        ax.set_ylabel(r"depth $m_{50}$ (mag)")
        ax.set_title(title, **TITLE_KW)
        ax.annotate(f"{d['m50'].iloc[0]:.2f}", (x.iloc[0], d["m50"].iloc[0]),
                    textcoords="offset points", xytext=(6, 2), fontsize=6.4,
                    color=C["data"])
        ax.annotate(f"{d['m50'].iloc[-1]:.2f}", (x.iloc[-1], d["m50"].iloc[-1]),
                    textcoords="offset points", xytext=end_off, fontsize=6.4,
                    color=C["data"])

    def scatter_panel(ax, x, d, xlabel, title, logx=False):
        ax.errorbar(x, d["mad"], yerr=_err(d["mad_lo"], d["mad_hi"], d["mad"]),
                    fmt="s-", ms=3.0, lw=1.1, capsize=2, color=C["accent"])
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("scatter MAD (mag)")
        ax.set_title(title, **TITLE_KW)

    # panel letters run down each column: aperture (a,b) / sky (c,d) / seeing (e,f)
    d_bg = float(bg["m50"].iloc[0] - bg["m50"].iloc[-1])
    depth_panel(ax_b, bg["background"], bg,
                f"(c) sky 50$\\to$1000 ADU: $-${d_bg:.2f} mag",
                logx=True, end_off=(-40, 2))
    ax_b.set_xlim(38, 1300)
    scatter_panel(ax_e, bg["background"], bg, "sky background (ADU)",
                  f"(d) sky: scatter {bg['mad'].iloc[0] * 1000:.0f}$\\to$"
                  f"{bg['mad'].iloc[-1] * 1000:.0f} mmag", logx=True)

    d_se = float(se["m50"].iloc[0] - se["m50"].iloc[-1])
    se_lo, se_hi = float(se["fwhm_px"].iloc[0]), float(se["fwhm_px"].iloc[-1])
    depth_panel(ax_c, se["fwhm_px"], se,
                f"(e) seeing {se_lo:g}$\\to${se_hi:g} px: "
                f"$-${d_se:.2f} mag", end_off=(-38, 0))
    scatter_panel(ax_f, se["fwhm_px"], se, "seeing FWHM (px)",
                  f"(f) seeing: scatter {se['mad'].iloc[0] * 1000:.0f}$\\to$"
                  f"{se['mad'].iloc[-1] * 1000:.0f} mmag")

    # ── in-figure data provenance (mandatory since FIGURE_REBUILD_PLAN) ─────
    baked = sorted(set(int(n) for d in (ap, bg, se, thr) for n in d["n_baked"]))
    baked_txt = f"{baked[0]}" if len(baked) == 1 else f"{min(baked)}-{max(baked)}"
    prov = (
        "Synthetic frames: APEX generator (apex.benchmark.synthetic_frame), seed "
        f"{SEED} \u00b7 {FRAME_BASE['size']}$^2$ px \u00b7 {baked_txt} field stars "
        "\u00b7 gain 1.5 e$^-$/ADU \u00b7 "
        "RN 5 e$^-$ \u00b7 ZP 25 \u00b7 r. (a,b): one fixed frame, FWHM "
        f"{APERTURE_FRAME_FWHM:g} px, sky 150 ADU; "
        "(c-f): frame regenerated per condition from the FWHM 3.5 px / sky 150 ADU baseline.\n"
        f"Each point: {TRIALS} trials $\\times$ {STARS_PER_TRIAL} injected stars "
        "(empirical-PSF injection $\\to$ production Step-4 detection $\\to$ forced aperture "
        "photometry). Bars: 95% CI, cluster bootstrap over trials; depth CI from the "
        "production completeness fit."
    )
    fig.tight_layout(pad=0.5, h_pad=1.0, w_pad=1.2, rect=(0, 0.075, 1, 1))
    fig.text(0.005, 0.005, prov, fontsize=5.6, color=PALETTE["grey"], va="bottom")

    paths = save_fig(fig, "fig3_parameter_sweep", FIG_DIR)
    plt.close(fig)
    return paths


# ── caption ──────────────────────────────────────────────────────────────────
def write_caption_legacy(ap, bg, se, thr) -> Path:
    CAP_DIR.mkdir(parents=True, exist_ok=True)
    cap_path = CAP_DIR / "fig3_parameter_sweep.md"

    fwhm_det = float(ap["fwhm_det_px"].iloc[0])
    solid = ap[~ap["clamped"]]
    clamp = ap[ap["clamped"]]
    opt = solid.loc[solid["mad"].idxmin()]
    x_opt = float(opt["r_ap_px"]) / fwhm_det
    rmse_med = float(np.median(solid["rmse"]))

    def tbl(df, cols, fmts):
        head = "| " + " | ".join(cols) + " |"
        sep = "|" + "|".join(["---"] * len(cols)) + "|"
        lines = [head, sep]
        for _i, r in df.iterrows():
            lines.append("| " + " | ".join(f.format(r[c])
                                           for c, f in zip(cols, fmts)) + " |")
        return "\n".join(lines)

    # monotonicity verdict for the scatter trends (the old trials=4 sawtooth)
    def verdict(d):
        diffs = np.diff(d["mad"].to_numpy(float))
        half = 0.5 * (d["mad_hi"] - d["mad_lo"]).to_numpy(float)
        wobble = float(np.max(np.abs(diffs[diffs < 0]))) if (diffs < 0).any() else 0.0
        return wobble, float(np.median(half))

    bg_wob, bg_ci = verdict(bg)
    se_wob, se_ci = verdict(se)

    thr_lo = thr.loc[thr["detect_sigma"].idxmin()]
    thr_hi = thr.loc[thr["detect_sigma"].idxmax()]
    thr_false_max = int(thr["false_detections"].max())

    text = f"""# Figure — parameter & observing-conditions sensitivity (rebuilt 2026-08-02)

**Data: synthetic only.** APEX synthetic generator
(`apex.benchmark.synthetic_frame`), seed {SEED};
{FRAME_BASE['size']}x{FRAME_BASE['size']} px,
{int(ap['n_baked'].iloc[0])} field stars (star density matches the old 640 px
/ 170-star frame), gain 1.5 e-/ADU, read noise 5 e-, zero point 25.0, r band.
Every point is an
artificial-star benchmark through the production code path (empirical-PSF
injection -> Step-4 detection -> forced aperture photometry;
`apex.benchmark.runner.run_benchmark`) with **{TRIALS} trials x
{STARS_PER_TRIAL} injected stars = {TRIALS * STARS_PER_TRIAL} injections per
point** (the old figure used 4 trials; its sawtooth was Monte-Carlo noise, see
verdict below). Error bars are 95% CIs from a cluster bootstrap over trials
({N_BOOT} resamples for MAD/RMSE; the depth CI comes from the production
completeness fit's own {M50_BOOT}-sample cluster bootstrap). Reproduce:
`validation/paper/fig3_parameter_sweep.py` (distilled numbers in
`validation/paper/data_parameter_sweep/`).

**(a) Aperture radius — scatter minimum at {x_opt:.2f}xFWHM.** One fixed frame
with FWHM {APERTURE_FRAME_FWHM:.1f} px (detected {fwhm_det:.2f} px), sky 150
ADU, injections 14.5-19.0 mag. The forced-aperture radius is
`max(min_r_ap_px = 4 px, scale x FWHM)`. **Why FWHM 6 px:** on the old 3.5 px
frame every requested scale <= 1.1xFWHM collapsed onto the same 4 px floor
radius, and the old panel drew those duplicate measurements as a "flat"
segment — a clamp artifact, not physics. On this frame the floor sits at
{MIN_R_AP_PX / fwhm_det:.2f}xFWHM, so 0.7-2.5xFWHM are all distinct apertures.
The one deliberately clamped request (0.5xFWHM -> 4 px, open marker) lands on
the floor. Scatter (MAD of the forced-photometry magnitude error) is minimized
at {x_opt:.2f}xFWHM with MAD {opt['mad']:.4f} mag
[{opt['mad_lo']:.4f}, {opt['mad_hi']:.4f}]; larger apertures admit sky noise
({solid['mad'].iloc[-1]:.3f} mag at 2.5xFWHM). The grey band marks radii the
production pipeline cannot reach (below the 4 px floor); the dashed line is
the default scale {DEFAULT_APERTURE_SCALE:g}.

**(b) Position RMSE is aperture-independent.** Median {rmse_med:.2f} px across
all unclamped apertures; the centroid comes from Step-4 detection, which does
not use the photometry aperture.

**(c,d) Sky background** (frame regenerated per level, injections 14.0-20.0
mag, log x): sky 50 -> 1000 ADU makes the 50% completeness depth
{bg['m50'].iloc[0]:.2f} -> {bg['m50'].iloc[-1]:.2f} mag
({float(bg['m50'].iloc[0] - bg['m50'].iloc[-1]):.2f} mag shallower) and the
scatter {bg['mad'].iloc[0]:.3f} -> {bg['mad'].iloc[-1]:.3f} mag: a brighter sky
raises the shot-noise floor.

**(e,f) Seeing** (frame regenerated per FWHM, sky 150 ADU): FWHM
{se['fwhm_px'].iloc[0]:g} -> {se['fwhm_px'].iloc[-1]:g} px makes the depth
{se['m50'].iloc[0]:.2f} -> {se['m50'].iloc[-1]:.2f} mag
({float(se['m50'].iloc[0] - se['m50'].iloc[-1]):.2f} mag shallower) and the
scatter {se['mad'].iloc[0]:.3f} -> {se['mad'].iloc[-1]:.3f} mag: the same flux
spreads over more pixels, lowering peak SNR. The grid stays inside the
production FWHM-QC window (`[fwhm] px_min = 3.0`): below it the frame-FWHM
estimator is outside its envelope (a true 2.5 px frame was measured as 8.5 px
from 5 stars), so a 2.5 px point would not be a pipeline measurement.

**Sawtooth verdict (old R5 defect).** With 4 trials the scatter curves wobbled
non-monotonically. At 30 trials the largest remaining downward step between
adjacent scatter points is {bg_wob * 1000:.1f} mmag (sky) / {se_wob * 1000:.1f}
mmag (seeing) against median CI half-widths of {bg_ci * 1000:.1f} /
{se_ci * 1000:.1f} mmag — the old sawtooth was Monte-Carlo noise, not
structure. What survives at 30 trials is consistent with monotone trends
within the CIs.

**Detection-threshold companion sweep (not drawn) — the old "threshold-robust"
claim was a no-op artifact.** The old run overrode only `detect_sigma`, which
the per-filter `detect_sigma_r` in parameters.toml silently masks, so all five
"swept" thresholds ran at 3.2-sigma and returned bit-identical results; the
"depth constant to < 0.01 mag over 2-6 sigma" sentence measured nothing. With
the per-filter key overridden too (each point's `sigma_used` is verified equal
to the requested value), the threshold sets the depth directly: m50
{thr_lo['m50']:.2f} at {thr_lo['detect_sigma']:g} sigma ->
{thr_hi['m50']:.2f} at {thr_hi['detect_sigma']:g} sigma
({float(thr_lo['m50'] - thr_hi['m50']):.2f} mag), while false detections stay
at <= {thr_false_max} per {TRIALS}-trial point on this clean low-crowding
field. Lowering the threshold is therefore free depth *on this field*; the
real-frame cost side (false-detection contamination vs threshold) is measured
on real frames in the detection-threshold section (its own figure).
Production default: 3.2 sigma.

## Sweep 1 — aperture (fixed frame, FWHM {APERTURE_FRAME_FWHM:.1f} px)

{tbl(ap, ["requested_scale", "r_ap_px", "mad", "mad_lo", "mad_hi", "rmse"],
     ["{:.2f}", "{:.2f}", "{:.4f}", "{:.4f}", "{:.4f}", "{:.3f}"])}

## Sweep 2 — sky background (frame per value)

{tbl(bg, ["background", "m50", "m50_lo", "m50_hi", "mad", "mad_lo", "mad_hi"],
     ["{:.0f}", "{:.3f}", "{:.3f}", "{:.3f}", "{:.4f}", "{:.4f}", "{:.4f}"])}

## Sweep 3 — seeing / PSF FWHM (frame per value)

{tbl(se, ["fwhm_px", "m50", "m50_lo", "m50_hi", "mad", "mad_lo", "mad_hi"],
     ["{:.1f}", "{:.3f}", "{:.3f}", "{:.3f}", "{:.4f}", "{:.4f}", "{:.4f}"])}

## Sweep 4 — detection threshold (companion, not drawn)

{tbl(thr, ["detect_sigma", "sigma_used", "m50", "m50_lo", "m50_hi",
           "false_detections"],
     ["{:.1f}", "{:.1f}", "{:.3f}", "{:.3f}", "{:.3f}", "{:.0f}"])}
"""
    cap_path.write_text(text, encoding="utf-8")
    return cap_path


def make_figure(ap: pd.DataFrame, bg: pd.DataFrame, se: pd.DataFrame,
                thr: pd.DataFrame) -> dict:
    """Compact three-panel summary; the threshold sweep has its own figure."""
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, 2.75))
    axa, axb, axc = axes
    title_kw = dict(loc="left", fontsize=8.5)

    fwhm_det = float(ap["fwhm_det_px"].iloc[0])
    solid = ap[~ap["clamped"]]
    clamp = ap[ap["clamped"]]
    x_s = solid["r_ap_px"] / fwhm_det
    x_c = clamp["r_ap_px"] / fwhm_det
    optimum = solid.loc[solid["mad"].idxmin()]
    x_opt = float(optimum["r_ap_px"]) / fwhm_det
    rmse_med = float(np.median(solid["rmse"]))

    axa.axvspan(0.5, MIN_R_AP_PX / fwhm_det, color=PALETTE["grey"],
                alpha=0.15, lw=0)
    axa.axvline(DEFAULT_APERTURE_SCALE, color=PALETTE["grey"], ls="--", lw=0.9)
    axa.errorbar(
        x_s, solid["mad"],
        yerr=_err(solid["mad_lo"], solid["mad_hi"], solid["mad"]),
        fmt="o-", ms=3.5, lw=1.1, capsize=2, color=C["data"],
        label="unclamped",
    )
    axa.errorbar(
        x_c, clamp["mad"],
        yerr=_err(clamp["mad_lo"], clamp["mad_hi"], clamp["mad"]),
        fmt="o", ms=4.2, mfc="white", capsize=2, color=C["data"],
        label="4 px floor",
    )
    axa.annotate(
        f"minimum {optimum['mad']:.3f} mag",
        (x_opt, float(optimum["mad"])), xytext=(15, 12),
        textcoords="offset points", fontsize=6.4,
        arrowprops=dict(arrowstyle="-", lw=0.7, color=C["model"]),
    )
    axa.text(
        0.98, 0.95, f"position RMSE {rmse_med:.3f} px\n(aperture-independent)",
        transform=axa.transAxes, ha="right", va="top", fontsize=6.1,
        color=PALETTE["grey"],
    )
    axa.set_xlabel("aperture radius / FWHM")
    axa.set_ylabel("magnitude-error MAD (mag)")
    axa.set_ylim(0.03, 0.22)
    axa.set_title(f"(a) optimum at {x_opt:.2f}$\\times$FWHM", **title_kw)

    def condition_panel(ax, x, data, xlabel, title, logx=False):
        ax.errorbar(
            x, data["m50"],
            yerr=_err(data["m50_lo"], data["m50_hi"], data["m50"]),
            fmt="o-", ms=3.5, lw=1.1, capsize=2, color=C["data"],
        )
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"50% completeness $m_{50}$ (mag)")
        ax.set_title(title, **title_kw)
        ax.annotate(
            f"{data['m50'].iloc[0]:.2f}", (x.iloc[0], data["m50"].iloc[0]),
            xytext=(4, 3), textcoords="offset points", fontsize=6.4,
        )
        ax.annotate(
            f"{data['m50'].iloc[-1]:.2f}", (x.iloc[-1], data["m50"].iloc[-1]),
            xytext=(-4, 3), textcoords="offset points", ha="right", fontsize=6.4,
        )
        ax.text(
            0.03, 0.06,
            f"scatter {data['mad'].iloc[0] * 1000:.0f} → "
            f"{data['mad'].iloc[-1] * 1000:.0f} mmag",
            transform=ax.transAxes, fontsize=6.2, color=PALETTE["grey"],
        )

    d_bg = float(bg["m50"].iloc[0] - bg["m50"].iloc[-1])
    condition_panel(
        axb, bg["background"], bg, "sky background (ADU)",
        f"(b) sky increase: −{d_bg:.2f} mag", logx=True,
    )
    axb.set_xlim(38, 1300)

    d_se = float(se["m50"].iloc[0] - se["m50"].iloc[-1])
    condition_panel(
        axc, se["fwhm_px"], se, "seeing FWHM (px)",
        f"(c) broader PSF: −{d_se:.2f} mag",
    )

    baked = sorted(set(int(n) for d in (ap, bg, se, thr) for n in d["n_baked"]))
    baked_txt = f"{baked[0]}" if len(baked) == 1 else f"{min(baked)}–{max(baked)}"
    provenance = (
        "Synthetic APEX frames · seed "
        f"{SEED} · {FRAME_BASE['size']}$^2$ px · {baked_txt} stars · "
        "gain 1.5 e$^-$/ADU · RN 5 e$^-$ · $r$ band. "
        f"Each point: {TRIALS} trials $\\times$ {STARS_PER_TRIAL} injections; "
        "bars: 95% cluster-bootstrap CI."
    )
    fig.tight_layout(rect=(0, 0.13, 1, 1), w_pad=1.0)
    fig.text(0.005, 0.015, provenance, fontsize=5.6,
             color=PALETTE["grey"], va="bottom")
    paths = save_fig(fig, "fig3_parameter_sweep", FIG_DIR)
    plt.close(fig)
    return paths


def write_caption(ap, bg, se, thr) -> Path:
    CAP_DIR.mkdir(parents=True, exist_ok=True)
    path = CAP_DIR / "fig3_parameter_sweep.md"
    fwhm_det = float(ap["fwhm_det_px"].iloc[0])
    solid = ap[~ap["clamped"]]
    optimum = solid.loc[solid["mad"].idxmin()]
    x_opt = float(optimum["r_ap_px"]) / fwhm_det
    rmse_med = float(np.median(solid["rmse"]))
    text = f"""# Figure — parameter and observing-condition sensitivity

Synthetic APEX frames ({FRAME_BASE['size']}×{FRAME_BASE['size']} px; seed
{SEED}; gain 1.5 e-/ADU; read noise 5 e-). Each point contains {TRIALS}
trials of {STARS_PER_TRIAL} injected stars and follows the production detection
and forced-aperture path; bars are 95% cluster-bootstrap intervals.
**(a)** The magnitude-error MAD is lowest at {x_opt:.2f}×FWHM
({optimum['mad']:.3f} mag). The open marker is the request constrained by the
4 px minimum radius; the median position RMSE is {rmse_med:.3f} px and does not
depend on the photometric aperture. **(b)** Raising the sky from
{bg['background'].iloc[0]:.0f} to {bg['background'].iloc[-1]:.0f} ADU changes
$m_{{50}}$ from {bg['m50'].iloc[0]:.2f} to {bg['m50'].iloc[-1]:.2f} mag and
the MAD from {bg['mad'].iloc[0] * 1000:.0f} to
{bg['mad'].iloc[-1] * 1000:.0f} mmag. **(c)** Broadening the PSF from
{se['fwhm_px'].iloc[0]:.1f} to {se['fwhm_px'].iloc[-1]:.1f} px changes
$m_{{50}}$ from {se['m50'].iloc[0]:.2f} to {se['m50'].iloc[-1]:.2f} mag and
the MAD from {se['mad'].iloc[0] * 1000:.0f} to
{se['mad'].iloc[-1] * 1000:.0f} mmag. The detection-threshold sweep is shown
separately in the threshold-validation figure.
"""
    path.write_text(text, encoding="utf-8")
    return path


# ── main ─────────────────────────────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep", choices=sorted(SWEEPS), default=None,
                        help="run only this sweep's benchmark points (no figure)")
    parser.add_argument("--figure-only", action="store_true",
                        help="draw from the distilled CSV without benchmarks")
    args = parser.parse_args()

    SWEEP_ROOT.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = DATA_DIR / "sweep_points.csv"

    if args.sweep:
        SWEEPS[args.sweep]()
        print(f"[done] sweep '{args.sweep}' points cached under {SWEEP_ROOT}")
        return 0

    if args.figure_only and csv_path.exists():
        allpts = pd.read_csv(csv_path)
        # bool column survives the CSV round-trip as object ("True"/nan)
        allpts["clamped"] = allpts["clamped"].map(
            lambda v: str(v).strip().lower() == "true")
        ap = allpts[allpts["sweep"] == "aperture"].reset_index(drop=True)
        bg = allpts[allpts["sweep"] == "background"].reset_index(drop=True)
        se = allpts[allpts["sweep"] == "seeing"].reset_index(drop=True)
        thr = allpts[allpts["sweep"] == "threshold"].reset_index(drop=True)
    else:
        print("== Sweep 1: aperture radius (fixed frame, FWHM 6 px) ==")
        ap = run_aperture_sweep()
        print("== Sweep 2: sky background (frame per value) ==")
        bg = run_background_sweep()
        print("== Sweep 3: seeing / FWHM (frame per value) ==")
        se = run_seeing_sweep()
        print("== Sweep 4: detection threshold (companion) ==")
        thr = run_threshold_sweep()
        allpts = pd.concat([ap, bg, se, thr], ignore_index=True)
        allpts.to_csv(csv_path, index=False)
        (DATA_DIR / "meta.json").write_text(json.dumps({
            "seed": SEED, "boot_seed": BOOT_SEED, "trials": TRIALS,
            "stars_per_trial": STARS_PER_TRIAL, "n_boot": N_BOOT,
            "m50_boot": M50_BOOT, "frame_base": FRAME_BASE,
            "aperture_frame_fwhm": APERTURE_FRAME_FWHM,
            "grids": {"aperture": APERTURE_SCALES, "background": BACKGROUNDS,
                      "seeing": FWHMS, "detect_sigma": DETECT_SIGMAS},
        }, indent=1), encoding="utf-8")

    paths = make_figure(ap, bg, se, thr)
    cap_path = write_caption(ap, bg, se, thr)
    print("Wrote:")
    for ext, path in paths.items():
        print(f"  {ext}: {path}")
    print(f"  caption: {cap_path}")
    print(f"  distilled: {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
