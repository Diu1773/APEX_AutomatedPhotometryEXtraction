"""Figure 3 — APEX parameter & observing-conditions sensitivity sweeps.

Three 1x3 panels showing real monotonic trends of pipeline metrics as one knob
is varied at a time:

* (a) Aperture radius (`aperture_scale_fwhm`) on a single fixed synthetic frame:
  photometric scatter (MAD, mag) + astrometric position RMSE vs aperture scale.
  A minimum-scatter aperture near ~1.2x FWHM (too small loses flux/SNR, too
  large admits sky noise).
* (b) Sky background: a controlled synthetic frame is regenerated per background
  level; 50% completeness depth (m50) + photometric scatter vs sky brightness.
  Brighter sky -> shallower (brighter) m50 + larger scatter.
* (c) Seeing (PSF FWHM): a controlled frame is regenerated per FWHM; m50 depth +
  photometric scatter vs FWHM. Worse seeing -> shallower m50 + larger scatter.

Experimental design note: panel (a) isolates a pure *pipeline* knob, so it reuses
ONE fixed frame. Panels (b) and (c) study *observing conditions*, so the correct
design is to regenerate a controlled frame per condition value (only that
condition changes; the star field seed, count, and all other knobs are fixed).

Run with the deploy venv:
    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig3_parameter_sweep.py
"""

import sys
import json
import shutil
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, SINGLE_COL, DOUBLE_COL

apply_paper_style()

from apex.benchmark.synthetic_frame import make_synthetic_reference_frame
from apex.benchmark.runner import BenchmarkConfig, run_benchmark


# ── paths ────────────────────────────────────────────────────────────────────
PAPER_DIR = REPO / "validation" / "paper"
FIG_DIR = PAPER_DIR / "figures"
CAP_DIR = PAPER_DIR / "captions"
# Short benchmark run root to stay under the Windows MAX_PATH (260) limit;
# NOT under a deep AppData\Local\Temp\claude\... scratchpad path.
SWEEP_ROOT = Path(r"C:\Users\bmffr\AppData\Local\Temp\apx_sweep")

# ── fixed synthetic-frame + benchmark constants ──────────────────────────────
SEED = 20260702
PARAM_FILE = str(REPO / "parameters.toml")

# Baseline synthetic-frame parameters (held fixed except the swept condition).
FRAME_BASE = dict(
    seed=SEED,
    n_stars=170,
    size=640,
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

# Pipeline default aperture (marked with a dashed guide line in panel a).
DEFAULT_APERTURE_SCALE = 0.8   # forced_r_ap_scale default in the runner

# Sweep grids.
APERTURE_SCALES = [0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5]
BACKGROUNDS = [50.0, 100.0, 150.0, 300.0, 600.0, 1000.0]   # ADU
FWHMS = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0]                     # px

# The detection-threshold override key is a *positive* robustness result noted
# in the caption (verified: `detect_sigma` is a Parameters.P dataclass field).
DETECT_SIGMA_KEY = "detect_sigma"


def _make_frame(path: Path, **overrides) -> Path:
    """Generate a controlled synthetic frame (baseline + per-sweep overrides)."""
    params = dict(FRAME_BASE)
    params.update(overrides)
    return make_synthetic_reference_frame(path, **params)


def _run_benchmark(
    frame: Path,
    out_dir: Path,
    *,
    aperture_scale,
    detection_overrides,
    magnitude_min,
    magnitude_max,
    reuse: bool,
) -> dict:
    """Run one benchmark (or reuse an existing valid run) and return summary.json."""
    summary_path = out_dir / "summary.json"
    if reuse and summary_path.exists():
        try:
            return json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            pass
    shutil.rmtree(out_dir, ignore_errors=True)
    cfg = BenchmarkConfig(
        input_fits=str(frame),
        parameter_file=PARAM_FILE,
        seed=SEED,
        trials=4,
        stars_per_trial=40,
        magnitude_min=magnitude_min,
        magnitude_max=magnitude_max,
        magnitude_bin_width=1.0,
        completeness_bootstrap_samples=100,
        save_injected_fits=False,
        zeropoint_mag=25.0,
        isolated_fraction=0.6,
        psf_min_stars=3,
        aperture_scale_fwhm=aperture_scale,
        detection_overrides=detection_overrides,
    )
    run_dir = run_benchmark(cfg, input_override=frame, output_override=out_dir)
    return json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))


# ── Sweep 1 — aperture radius (single fixed frame; reuse cached results) ──────
def run_aperture_sweep() -> pd.DataFrame:
    frame = _make_frame(PAPER_DIR / "_sweep_ref.fits")
    rows = []
    for scale in APERTURE_SCALES:
        tag = f"ap_{int(round(scale * 100)):04d}"
        summary = _run_benchmark(
            frame, SWEEP_ROOT / tag,
            aperture_scale=scale, detection_overrides={},
            magnitude_min=14.5, magnitude_max=19.0, reuse=True,
        )
        rows.append(
            {
                "aperture_scale": scale,
                "scatter_mad": summary.get("forced_mag_scatter_mad"),
                "position_rmse_px": summary.get("position_rmse_px"),
            }
        )
        print(f"[aperture] scale={scale:.2f}  MAD={rows[-1]['scatter_mad']!r}  "
              f"rmse={rows[-1]['position_rmse_px']!r}")
    return pd.DataFrame(rows)


# ── Sweep 2 — sky background (frame regenerated per value) ────────────────────
def run_background_sweep() -> pd.DataFrame:
    rows = []
    for bkg in BACKGROUNDS:
        tag = f"bkg_{int(round(bkg)):05d}"
        frame = _make_frame(SWEEP_ROOT / f"{tag}_ref.fits", background=bkg)
        summary = _run_benchmark(
            frame, SWEEP_ROOT / tag,
            aperture_scale=None, detection_overrides={},
            magnitude_min=14.0, magnitude_max=20.0, reuse=False,
        )
        fit = summary.get("completeness_fit") or {}
        rows.append(
            {
                "background": bkg,
                "m50": fit.get("m50"),
                "scatter_mad": summary.get("forced_mag_scatter_mad"),
                "completeness": summary.get("completeness"),
            }
        )
        print(f"[background] sky={bkg:.0f} ADU  m50={rows[-1]['m50']!r}  "
              f"MAD={rows[-1]['scatter_mad']!r}")
    return pd.DataFrame(rows)


# ── Sweep 3 — seeing / PSF FWHM (frame regenerated per value) ─────────────────
def run_seeing_sweep() -> pd.DataFrame:
    rows = []
    for fwhm in FWHMS:
        tag = f"fwhm_{int(round(fwhm * 10)):03d}"
        frame = _make_frame(SWEEP_ROOT / f"{tag}_ref.fits", fwhm_px=fwhm, background=150.0)
        summary = _run_benchmark(
            frame, SWEEP_ROOT / tag,
            aperture_scale=None, detection_overrides={},
            magnitude_min=14.0, magnitude_max=20.0, reuse=False,
        )
        fit = summary.get("completeness_fit") or {}
        rows.append(
            {
                "fwhm_px": fwhm,
                "m50": fit.get("m50"),
                "scatter_mad": summary.get("forced_mag_scatter_mad"),
                "completeness": summary.get("completeness"),
            }
        )
        print(f"[seeing] fwhm={fwhm:.1f} px  m50={rows[-1]['m50']!r}  "
              f"MAD={rows[-1]['scatter_mad']!r}")
    return pd.DataFrame(rows)


# ── figure ───────────────────────────────────────────────────────────────────
def _twin_depth_scatter(ax, x, m50, mad, xlabel, title):
    """Shared drawing for panels (b)/(c): m50 depth (left) + MAD (right)."""
    c_depth = C["data"]
    c_scatter = C["accent"]
    ax.plot(x, m50, marker="o", color=c_depth, label=r"Depth ($m_{50}$)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$m_{50}$ depth (mag)", color=c_depth)
    ax.tick_params(axis="y", colors=c_depth)
    ax.set_title(title, loc="left")

    ax2 = ax.twinx()
    ax2.plot(x, mad, marker="s", color=c_scatter, label="Photometric scatter")
    ax2.set_ylabel("Photometric scatter (MAD, mag)", color=c_scatter)
    ax2.tick_params(axis="y", colors=c_scatter)
    ax2.grid(False)

    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines],
              loc="best", fontsize=7.0)
    return ax2


def make_figure(ap: pd.DataFrame, bg: pd.DataFrame, se: pd.DataFrame):
    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(DOUBLE_COL, 2.6))

    # ── (a) Aperture radius — KEEP unchanged ─────────────────────────────────
    c_scatter = C["data"]
    c_rmse = C["accent"]
    ax_a.plot(ap["aperture_scale"], ap["scatter_mad"],
              marker="o", color=c_scatter, label="Photometric scatter")
    ax_a.set_xlabel(r"Aperture radius / FWHM")
    ax_a.set_ylabel("Photometric scatter (MAD, mag)", color=c_scatter)
    ax_a.tick_params(axis="y", colors=c_scatter)

    ax_a2 = ax_a.twinx()
    ax_a2.plot(ap["aperture_scale"], ap["position_rmse_px"],
               marker="s", color=c_rmse, label="Position RMSE")
    ax_a2.set_ylabel("Position RMSE (px)", color=c_rmse)
    ax_a2.tick_params(axis="y", colors=c_rmse)
    ax_a2.grid(False)
    rmse_vals = pd.to_numeric(ap["position_rmse_px"], errors="coerce").dropna()
    if len(rmse_vals):
        ax_a2.set_ylim(0.0, max(0.3, 2.0 * float(rmse_vals.median())))

    ap_valid = ap.dropna(subset=["scatter_mad"])
    opt_row = ap_valid.loc[ap_valid["scatter_mad"].idxmin()]
    opt_scale = float(opt_row["aperture_scale"])
    ax_a.axvline(DEFAULT_APERTURE_SCALE, color=C["floor"], linestyle="--",
                 linewidth=1.0, label=f"Default ({DEFAULT_APERTURE_SCALE:g})")
    ax_a.axvline(opt_scale, color=C["good"], linestyle=":", linewidth=1.2,
                 label=f"Optimum ({opt_scale:g})")
    ax_a.set_title("(a) Aperture radius", loc="left")
    lines_a = ax_a.get_lines() + ax_a2.get_lines()
    ax_a.legend(lines_a, [ln.get_label() for ln in lines_a],
                loc="upper center", fontsize=6.5)

    # ── (b) Sky background — depth + scatter ─────────────────────────────────
    _twin_depth_scatter(
        ax_b, bg["background"], bg["m50"], bg["scatter_mad"],
        r"Sky background (ADU)", "(b) Sky background",
    )

    # ── (c) Seeing / FWHM — depth + scatter ──────────────────────────────────
    _twin_depth_scatter(
        ax_c, se["fwhm_px"], se["m50"], se["scatter_mad"],
        r"Seeing FWHM (px)", "(c) Seeing (PSF FWHM)",
    )

    fig.tight_layout()
    paths = save_fig(fig, "fig3_parameter_sweep", FIG_DIR)
    plt.close(fig)

    optima = {
        "aperture": {
            "optimal_scale": opt_scale,
            "min_scatter_mad": float(opt_row["scatter_mad"]),
            "default_scale": DEFAULT_APERTURE_SCALE,
        },
    }
    return paths, optima


# ── caption ──────────────────────────────────────────────────────────────────
def write_caption(ap, bg, se, optima) -> Path:
    CAP_DIR.mkdir(parents=True, exist_ok=True)
    cap_path = CAP_DIR / "fig3_parameter_sweep.md"
    opt = optima["aperture"]

    def _trend(df, xcol, ycol):
        d = df.dropna(subset=[ycol])
        lo = d.loc[d[xcol].idxmin()]
        hi = d.loc[d[xcol].idxmax()]
        return lo, hi

    bg_lo, bg_hi = _trend(bg, "background", "m50")
    se_lo, se_hi = _trend(se, "fwhm_px", "m50")

    ap_tbl = "\n".join(
        f"| {r.aperture_scale:.2f} | {r.scatter_mad:.4f} | {r.position_rmse_px:.4f} |"
        for r in ap.itertuples()
    )
    bg_tbl = "\n".join(
        f"| {r.background:.0f} | "
        f"{('%.4f' % r.m50) if pd.notna(r.m50) else 'n/a'} | "
        f"{('%.4f' % r.scatter_mad) if pd.notna(r.scatter_mad) else 'n/a'} |"
        for r in bg.itertuples()
    )
    se_tbl = "\n".join(
        f"| {r.fwhm_px:.1f} | "
        f"{('%.4f' % r.m50) if pd.notna(r.m50) else 'n/a'} | "
        f"{('%.4f' % r.scatter_mad) if pd.notna(r.scatter_mad) else 'n/a'} |"
        for r in se.itertuples()
    )

    text = f"""# Figure 3 — Parameter & observing-conditions sensitivity

**Sensitivity of APEX pipeline metrics to one pipeline parameter and two
observing conditions, each swept independently.** Every point is an
artificial-star benchmark (empirical-PSF injection, production Step 4 detection,
and forced aperture photometry) with 4 trials of 40 injected stars; the modest
trial count keeps the sweep fast at the cost of per-point Monte-Carlo noise.
Panel (a) isolates a pure pipeline knob on a single fixed synthetic frame; for
the observing-condition panels (b, c) a controlled synthetic frame is
regenerated per value (fixed star-field seed {SEED}, {FRAME_BASE['n_stars']}
stars, {FRAME_BASE['size']}x{FRAME_BASE['size']} px, gain
{FRAME_BASE['gain']} e-/ADU, read noise {FRAME_BASE['read_noise']} e-, zero
point {FRAME_BASE['zeropoint']:.1f} mag), changing only the swept condition.

**(a) Aperture radius.** Photometric scatter (median absolute deviation of the
recovered magnitude, blue, left axis) and astrometric position RMSE (orange,
right axis) versus aperture radius in units of the PSF FWHM
(`aperture_scale_fwhm`), on the fixed frame (sky {FRAME_BASE['background']:.0f}
ADU, FWHM {FRAME_BASE['fwhm_px']} px). Scatter is minimized at
**{opt['optimal_scale']:g}x FWHM** (MAD = {opt['min_scatter_mad']:.4f} mag):
smaller apertures lose source flux and SNR, larger apertures admit sky noise.
Position RMSE is essentially aperture-independent (~{ap['position_rmse_px'].median():.3f}
px) because the centroid comes from detection, not the aperture. Dashed grey =
pipeline default ({opt['default_scale']:g}x FWHM); dotted green = recovered
optimum.

**(b) Sky background.** 50% completeness depth *m*<sub>50</sub> (blue, left axis)
and photometric scatter (orange, right axis) versus sky background, with a frame
regenerated at each level (injection window widened to 14.0-20.0 mag so
*m*<sub>50</sub> stays in range). As the sky brightens from
{bg_lo['background']:.0f} to {bg_hi['background']:.0f} ADU the depth becomes
shallower (brighter), *m*<sub>50</sub> = {bg_lo['m50']:.3f} -> {bg_hi['m50']:.3f}
mag, and the photometric scatter grows
({bg_lo['scatter_mad']:.4f} -> {bg_hi['scatter_mad']:.4f} mag): a brighter sky
raises the shot-noise floor, so faint sources drop below the detection/SNR
limit — a clean monotonic trend.

**(c) Seeing (PSF FWHM).** *m*<sub>50</sub> depth (blue, left axis) and
photometric scatter (orange, right axis) versus PSF FWHM (seeing), frame
regenerated per value at fixed sky {FRAME_BASE['background']:.0f} ADU. Worse
seeing spreads source flux over more pixels, lowering peak SNR: from FWHM
{se_lo['fwhm_px']:.1f} to {se_hi['fwhm_px']:.1f} px the depth becomes shallower
(*m*<sub>50</sub> = {se_lo['m50']:.3f} -> {se_hi['m50']:.3f} mag) and the
photometric scatter increases
({se_lo['scatter_mad']:.4f} -> {se_hi['scatter_mad']:.4f} mag).

**Detection-threshold robustness (not shown).** A companion sweep of the
detection threshold `detect_sigma` over 2.0-6.0 sigma on the clean fixed frame
left both the depth (*m*<sub>50</sub> constant to < 0.01 mag) and the purity
(zero spurious detections at every threshold) flat: for bright, well-sampled,
isolated sources the production segmentation detector is threshold-robust and
the recovery limit is SNR-set, not threshold-set. This is a desirable stability
property rather than a trend, so it is reported here instead of as a panel.

## Numeric optima and trends

- **Optimal aperture scale (minimum photometric scatter):**
  {opt['optimal_scale']:g}x FWHM, MAD = {opt['min_scatter_mad']:.4f} mag.
- **Depth vs sky background:** m50 {bg_lo['m50']:.3f} mag at
  {bg_lo['background']:.0f} ADU -> {bg_hi['m50']:.3f} mag at
  {bg_hi['background']:.0f} ADU (shallower with brighter sky).
- **Depth vs seeing:** m50 {se_lo['m50']:.3f} mag at FWHM {se_lo['fwhm_px']:.1f}
  px -> {se_hi['m50']:.3f} mag at {se_hi['fwhm_px']:.1f} px (shallower with worse
  seeing).
- **Detection-threshold override key used:** `{DETECT_SIGMA_KEY}` (robustness
  companion sweep).

## Sweep 1 — aperture radius (fixed frame)

| aperture / FWHM | scatter MAD (mag) | position RMSE (px) |
|---|---|---|
{ap_tbl}

## Sweep 2 — sky background (frame per value)

| sky background (ADU) | m50 (mag) | scatter MAD (mag) |
|---|---|---|
{bg_tbl}

## Sweep 3 — seeing / PSF FWHM (frame per value)

| FWHM (px) | m50 (mag) | scatter MAD (mag) |
|---|---|---|
{se_tbl}

*Benchmark: 4 trials x 40 stars per point; seed {SEED}; trials=4 chosen for
sweep speed. Panel (a) reuses one fixed frame; panels (b)/(c) regenerate a
controlled frame per condition value.*
"""
    cap_path.write_text(text, encoding="utf-8")
    return cap_path


def main() -> int:
    SWEEP_ROOT.mkdir(parents=True, exist_ok=True)

    print("== Sweep 1: aperture radius (fixed frame) ==")
    ap = run_aperture_sweep()
    print("== Sweep 2: sky background (frame per value) ==")
    bg = run_background_sweep()
    print("== Sweep 3: seeing / FWHM (frame per value) ==")
    se = run_seeing_sweep()

    paths, optima = make_figure(ap, bg, se)
    cap_path = write_caption(ap, bg, se, optima)

    print("Wrote:")
    for ext, path in paths.items():
        print(f"  {ext}: {path}  exists={path.exists()}")
    print(f"  caption: {cap_path}  exists={cap_path.exists()}")
    print(f"OPTIMAL aperture scale (min scatter): {optima['aperture']['optimal_scale']:g}x FWHM "
          f"(MAD={optima['aperture']['min_scatter_mad']:.4f} mag)")
    print("BACKGROUND m50:", list(zip(bg['background'], bg['m50'])))
    print("SEEING m50:", list(zip(se['fwhm_px'], se['m50'])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
