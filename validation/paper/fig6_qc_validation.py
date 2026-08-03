"""Figure 6 — validation of the automatic frame-QC module (`evaluate_frame_qc`).

Question: does APEX's automatic frame QC (`apex/analysis/frame_qc.py`) correctly
flag frames with known defects while passing clean frames?

Design — the artificial-star methodology lifted one level up, to whole frames:
we synthesize a "night" of 44 frames in which the TRUTH (which frame is
defective, and how) is known by construction, run the REAL production Step-4
detection on every frame (`apex.benchmark.runner.run_apex_detection`, i.e. the
same Qt-free `run_detection` code path the GUI's DetectionWorker executes, with
the repository `parameters.toml`), aggregate the per-frame metrics EXACTLY the
way the production Step-4 window does before calling `evaluate_frame_qc`
(see `FrameQCPanel.load_frames` in `apex/gui/workflow/step4_source_detection.py`:
``sky_med``/``sky_sigma`` are the detector-reported global background level and
RMS, ``fwhm_med`` is the detector's median radial-profile FWHM of the brightest
sources), and compare the PASS/REVIEW/FAIL decisions against the truth.

Frame classes (total 44, every frame a fresh noise + star-field realization):
  * clean x24          — baseline synthetic frame
  * cloud x5           — zeropoint 25.0 -> 24.3 (all stars 0.7 mag fainter;
                          grey transparency loss, sky unchanged)
  * bad_seeing x5      — FWHM 3.5 -> 6.3 px (1.8x bloat)
  * bright_sky x5      — background 150 -> 750 ADU (5x)
  * noisy_readout x5   — frames REALIZED with read noise 25 e-, but the QC
                          input (the "header") still claims 5 e-. This mimics
                          undocumented electronics noise: the CCD-equation
                          check (sky_noise_ratio) is the one that should react.

Honesty rule: thresholds are the shipped `FrameQCThresholds()` defaults; the
per-class outcome is reported exactly as measured (no tuning). In particular a
grey 0.7 mag transparency loss is NOT expected to be reliably detectable at the
detection-metric level — if cloud frames pass, that is a measured limitation
motivating the planned post-photometry transparency QC stage.

Run (full pipeline, ~2-5 min for 44 production detections):
    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig6_qc_validation.py
Re-style without re-running detection (reuses cached frame metrics):
    ... fig6_qc_validation.py --reuse
"""

from __future__ import annotations

import shutil
import sys
import time
from collections import Counter
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

# Windows consoles default to cp949 here; keep unicode prints from crashing.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

from apex.analysis.frame_qc import (  # noqa: E402
    FAIL,
    PASS,
    REVIEW,
    FrameQCThresholds,
    evaluate_frame_qc,
)
from apex.benchmark.runner import run_apex_detection  # noqa: E402
from apex.benchmark.synthetic_frame import make_synthetic_reference_frame  # noqa: E402
from apex.config.parameters_cmd import Parameters  # noqa: E402


# ── configuration ─────────────────────────────────────────────────────────────
# SHORT work root: Windows MAX_PATH safety (never under AppData\...\claude\...).
WORK = Path(r"C:\Users\bmffr\AppData\Local\Temp\apx_qc")
PARAM_FILE = REPO / "parameters.toml"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"
DATADIR = REPO / "validation" / "paper" / "data" / "frame_qc"

SEED0 = 20260702000

# Baseline (clean) synthetic-frame parameters — identical to the tuned
# self-contained artificial-star defaults (apex/benchmark/validate.py).
# min_separation_px is pinned to the clean default (6 x 3.5 = 21 px) for ALL
# classes so every frame carries the same 170-star geometric truth; otherwise
# the bad-seeing frames would silently bake in FEWER stars (the generator's
# separation default scales with FWHM) and confound the source-count metric.
BASE = dict(
    n_stars=170,
    size=640,
    fwhm_px=3.5,
    zeropoint=25.0,
    background=150.0,
    gain=1.5,
    read_noise=5.0,
    mag_min=13.5,
    mag_max=19.0,
    filter_name="r",
    exptime=1.0,
    min_separation_px=21.0,
)

# What the "header" (and thus the QC input) claims for EVERY frame. For
# noisy_readout this is deliberately wrong (frames realized with RN = 25 e-):
# the header lie the sky-noise consistency check is designed to expose.
HEADER_GAIN_E_PER_ADU = 1.5
HEADER_RDNOISE_E = 5.0

TRUE_RDNOISE_NOISY_E = 25.0

CLASSES: list[tuple[str, int, dict]] = [
    ("clean", 24, {}),
    ("cloud", 5, {"zeropoint": 24.3}),
    ("bad_seeing", 5, {"fwhm_px": 6.3}),
    ("bright_sky", 5, {"background": 750.0}),
    ("noisy_readout", 5, {"read_noise": TRUE_RDNOISE_NOISY_E}),
]
TRUTH_ORDER = [name for name, _, _ in CLASSES]
DECISION_ORDER = [PASS, REVIEW, FAIL]

CLASS_STYLE = {
    "clean": dict(color=PALETTE["grey"], label="clean"),
    "cloud": dict(color=PALETTE["purple"], label=r"cloud (ZP $-0.7$ mag)"),
    "bad_seeing": dict(color=PALETTE["vermillion"], label=r"bad seeing (FWHM $\times 1.8$)"),
    "bright_sky": dict(color=PALETTE["orange"], label=r"bright sky ($\times 5$)"),
    "noisy_readout": dict(color=PALETTE["blue"], label=r"noisy readout (RN 25, hdr 5 e$^-$)"),
}
MATRIX_ROW_LABELS = {
    "clean": "clean",
    "cloud": "cloud\n(ZP $-0.7$)",
    "bad_seeing": "bad seeing\n($\\times 1.8$)",
    "bright_sky": "bright sky\n($\\times 5$)",
    "noisy_readout": "noisy readout\n(RN 25 | hdr 5)",
}
DECISION_MARKER = {PASS: "o", REVIEW: "^", FAIL: "X"}


# ── stage 1: synthesize the night + run the production detector ──────────────

def build_night_metrics() -> pd.DataFrame:
    """Generate 44 frames, run REAL Step-4 detection on each, aggregate metrics.

    Returns one row per frame with exactly the columns production Step 4 feeds
    `evaluate_frame_qc` (minus the per-source-quality / candidate-count columns,
    which `evaluate_frame_qc` skips NaN-safely when absent), plus `truth`.
    """
    shutil.rmtree(WORK, ignore_errors=True)
    frames_dir = WORK / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    base_params = Parameters(PARAM_FILE)

    rows: list[dict] = []
    idx = 0
    t_detect_total = 0.0
    for cls_name, n_frames, overrides in CLASSES:
        for _ in range(n_frames):
            seed = SEED0 + idx
            synth = dict(BASE)
            synth.update(overrides)
            fname = f"{idx:02d}_{cls_name}.fits"
            frame_path = frames_dir / fname
            make_synthetic_reference_frame(frame_path, seed=seed, **synth)

            det_dir = WORK / "det" / f"{idx:02d}"
            t0 = time.perf_counter()
            det_df, meta, _ = run_apex_detection(frame_path, base_params, det_dir, {})
            dt = time.perf_counter() - t0
            t_detect_total += dt

            row = {
                "file": fname,
                "truth": cls_name,
                "seed": seed,
                # ── the QC input row, mirroring FrameQCPanel.load_frames ──
                "filter": str(meta.get("filter", "r") or "r").strip(),
                "airmass": 1.0,
                "n_sources": int(meta.get("n_sources", 0) or 0),
                "fwhm_med": float(meta.get("fwhm_px", np.nan)),
                "sky_med": float(meta.get("bkg_median", np.nan)),
                "sky_sigma": float(meta.get("bkg_rms", np.nan)),
                "elong_med": float(meta.get("median_elongation", np.nan)),
                # The "header" values the QC believes — the lie for noisy_readout.
                "gain_e_per_adu": HEADER_GAIN_E_PER_ADU,
                "rdnoise_e": HEADER_RDNOISE_E,
                # bookkeeping (not consumed by evaluate_frame_qc)
                "n_detections_csv": int(len(det_df)),
                "detect_seconds": dt,
            }
            rows.append(row)
            print(
                f"[fig6] {idx + 1:2d}/44 {fname:22s} n_src={row['n_sources']:3d} "
                f"fwhm={row['fwhm_med']:.2f}px sky={row['sky_med']:.1f} "
                f"sky_sig={row['sky_sigma']:.2f} elong={row['elong_med']:.3f} "
                f"({dt:.1f}s)"
            )
            idx += 1

    df = pd.DataFrame(rows)
    print(f"[fig6] total production-detection time: {t_detect_total:.1f} s for {idx} frames")
    return df


# ── stage 2: run the QC under test ────────────────────────────────────────────

def run_qc(metrics: pd.DataFrame) -> tuple[pd.DataFrame, FrameQCThresholds]:
    thr = FrameQCThresholds()
    out = evaluate_frame_qc(metrics, params=None, thresholds=thr)
    return out, thr


def _reason_counter(group: pd.DataFrame) -> Counter:
    tokens: Counter = Counter()
    for reason in group["qc_reasons"].fillna("").astype(str):
        for tok in reason.split(","):
            tok = tok.strip()
            if tok:
                tokens[tok] += 1
    return tokens


def _subcheck_table(out: pd.DataFrame, thr: FrameQCThresholds) -> pd.DataFrame:
    """Attribute WHICH underlying metric crossed its review threshold, per frame.

    Reproduces the individual sub-conditions inside evaluate_frame_qc so the
    caption can say precisely what caught each class.
    """
    sub = pd.DataFrame(index=out.index)
    sub["truth"] = out["truth"]
    sub["fwhm_z_hit"] = out["fwhm_z"] > thr.fwhm_z_review
    sub["fwhm_ratio_hit"] = out["fwhm_model_ratio"] > thr.fwhm_model_ratio_review
    sub["sky_z_hit"] = out["sky_z"] > thr.sky_z_review
    sub["sky_sigma_z_hit"] = out["sky_sigma_z"] > thr.sky_z_review
    sub["noise_ratio_hit"] = out["sky_noise_ratio"] > thr.sky_noise_ratio_review
    sub["nsrc_z_hit"] = out["nsrc_z"] < -thr.nsrc_z_review
    sub["elong_hit"] = out["elong_med"] > thr.elong_review
    sub["depth_hit"] = out["depth_cost_mag"] > thr.depth_cost_review
    return sub.fillna(False)


# ── stage 3: figure ───────────────────────────────────────────────────────────

def draw_decision_matrix(ax, counts: pd.DataFrame) -> None:
    n_rows, n_cols = counts.shape
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.invert_yaxis()
    ax.set_aspect("auto")
    ax.grid(False)
    ax.minorticks_off()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    for i, truth in enumerate(counts.index):
        row_total = int(counts.loc[truth].sum())
        for j, dec in enumerate(counts.columns):
            n = int(counts.loc[truth, dec])
            correct = (dec == PASS) if truth == "clean" else (dec in (REVIEW, FAIL))
            if n > 0:
                base = C["good"] if correct else C["bad"]
                alpha = 0.16 + 0.55 * (n / max(row_total, 1))
            else:
                base, alpha = "#FFFFFF", 1.0
            ax.add_patch(
                Rectangle((j, i), 1, 1, facecolor=base, alpha=alpha,
                          edgecolor="#C8C8C8", linewidth=0.6)
            )
            ax.text(
                j + 0.5, i + 0.5, str(n),
                ha="center", va="center",
                fontsize=9.0 if n else 7.5,
                fontweight="bold" if n else "normal",
                color="#1A1A1A" if n else "#9A9A9A",
            )

    ax.set_xticks([j + 0.5 for j in range(n_cols)])
    ax.set_xticklabels(list(counts.columns), fontsize=8.0)
    ax.xaxis.set_ticks_position("top")
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels(
        [f"{MATRIX_ROW_LABELS[t]}\n$n={int(counts.loc[t].sum())}$" for t in counts.index],
        fontsize=6.8, linespacing=1.15,
    )
    for i, truth in enumerate(counts.index):
        ax.get_yticklabels()[i].set_color(CLASS_STYLE[truth]["color"])
    ax.set_title("(a) QC decision vs injected truth", loc="left", pad=14)


def draw_diagnostic_plane(ax, out: pd.DataFrame, thr: FrameQCThresholds) -> None:
    # Threshold guides from the shipped FrameQCThresholds defaults.
    guide_kw = dict(color=C["floor"], zorder=1)
    ax.axvline(thr.fwhm_model_ratio_review, linestyle="--", linewidth=0.9, **guide_kw)
    ax.axvline(thr.fwhm_model_ratio_fail, linestyle="-", linewidth=1.1, **guide_kw)
    ax.axhline(thr.sky_noise_ratio_review, linestyle="--", linewidth=0.9, **guide_kw)
    ax.axhline(thr.sky_noise_ratio_fail, linestyle="-", linewidth=1.1, **guide_kw)
    label_kw = dict(color=PALETTE["grey"], fontsize=6.3, zorder=2)
    ax.text(thr.fwhm_model_ratio_review, 0.745, " FWHM review 1.30",
            rotation=90, ha="right", va="bottom", **label_kw)
    ax.text(thr.fwhm_model_ratio_fail, 0.745, " FWHM fail 1.60",
            rotation=90, ha="right", va="bottom", **label_kw)
    ax.text(2.13, thr.sky_noise_ratio_review, "noise review 2.5 ",
            ha="right", va="bottom", **label_kw)
    ax.text(2.13, thr.sky_noise_ratio_fail, "noise fail 6.0 ",
            ha="right", va="bottom", **label_kw)

    for truth in TRUTH_ORDER:
        for dec in DECISION_ORDER:
            sel = (out["truth"] == truth) & (out["qc_status"] == dec)
            if not sel.any():
                continue
            marker = DECISION_MARKER[dec]
            color = CLASS_STYLE[truth]["color"]
            kw = dict(s=30, marker=marker, zorder=4, alpha=0.85)
            if marker == "X":
                kw.update(facecolor=color, edgecolor="none")
            else:
                kw.update(facecolor=color, edgecolor="white", linewidth=0.5)
            ax.scatter(
                out.loc[sel, "fwhm_model_ratio"],
                out.loc[sel, "sky_noise_ratio"],
                **kw,
            )

    ax.set_xlim(0.85, 2.15)
    ax.set_yscale("log")
    ax.set_ylim(0.72, 8.5)
    ax.yaxis.set_major_locator(
        mticker.FixedLocator([0.8, 1.0, 1.5, 2.5, 4.0, 6.0, 8.0])
    )
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlabel("FWHM / seeing model  (fwhm_model_ratio)")
    ax.set_ylabel("measured / expected sky noise\n(sky_noise_ratio)")

    class_handles = [
        Line2D([], [], linestyle="", marker="s", markersize=5.5,
               markerfacecolor=CLASS_STYLE[t]["color"], markeredgecolor="none",
               label=CLASS_STYLE[t]["label"])
        for t in TRUTH_ORDER
    ]
    dec_handles = [
        Line2D([], [], linestyle="", marker=DECISION_MARKER[d], markersize=5.5,
               markerfacecolor="#4D4D4D", markeredgecolor="#4D4D4D", label=d)
        for d in DECISION_ORDER
    ]
    leg1 = ax.legend(handles=class_handles, loc="upper left", fontsize=6.4,
                     handletextpad=0.35, borderaxespad=0.25, labelspacing=0.35)
    ax.add_artist(leg1)
    ax.legend(handles=dec_handles, loc="upper right", fontsize=6.4, ncol=1,
              handletextpad=0.35, borderaxespad=0.25, labelspacing=0.35,
              bbox_to_anchor=(1.0, 0.80))
    ax.set_title("(b) QC diagnostic plane", loc="left")


def make_figure(out: pd.DataFrame, counts: pd.DataFrame, thr: FrameQCThresholds):
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 3.0), gridspec_kw={"width_ratios": [1.0, 1.3]}
    )
    draw_decision_matrix(ax_a, counts)
    draw_diagnostic_plane(ax_b, out, thr)
    fig.tight_layout(w_pad=2.2)
    return fig


# ── stage 4: caption ──────────────────────────────────────────────────────────

def _fmt_reasons(counter: Counter) -> str:
    if not counter:
        return "—"
    return ", ".join(f"`{k}` ({v})" for k, v in counter.most_common())


def write_caption(
    out: pd.DataFrame,
    counts: pd.DataFrame,
    thr: FrameQCThresholds,
    sub: pd.DataFrame,
) -> Path:
    med = out.groupby("truth").median(numeric_only=True).reindex(TRUTH_ORDER)
    reasons = {t: _reason_counter(out[out["truth"] == t]) for t in TRUTH_ORDER}
    det_s = float(pd.to_numeric(out.get("detect_seconds"), errors="coerce").sum())

    n_clean = int(counts.loc["clean"].sum())
    n_clean_pass = int(counts.loc["clean", PASS])
    defect_rows = [t for t in TRUTH_ORDER if t != "clean"]
    n_defect = int(counts.loc[defect_rows].sum().sum())
    n_defect_flagged = int(counts.loc[defect_rows, [REVIEW, FAIL]].sum().sum())

    cloud_flagged = int(counts.loc["cloud", [REVIEW, FAIL]].sum())
    noisy_ratio_med = float(med.loc["noisy_readout", "sky_noise_ratio"])
    noisy_ratio_hits = int(sub.loc[sub["truth"] == "noisy_readout", "noise_ratio_hit"].sum())
    noisy_sigz_hits = int(sub.loc[sub["truth"] == "noisy_readout", "sky_sigma_z_hit"].sum())

    # Honest cloud narrative, decided by the measured outcome.
    # How the cloud frames would score if the night-relative count statistics
    # were referenced to the clean frames only (contamination-free scale).
    clean_nsrc = out.loc[out["truth"] == "clean", "n_sources"].to_numpy(float)
    c_med = float(np.median(clean_nsrc))
    c_scale = float(np.median(np.abs(clean_nsrc - c_med)) / 0.6745)
    cloud_nsrc = out.loc[out["truth"] == "cloud", "n_sources"].to_numpy(float)
    cloud_z_cleanref = (cloud_nsrc - c_med) / c_scale if c_scale > 0 else np.full(len(cloud_nsrc), np.nan)
    night_nsrc = out["n_sources"].to_numpy(float)
    n_med_night = float(np.median(night_nsrc))
    n_scale_night = float(np.median(np.abs(night_nsrc - n_med_night)) / 0.6745)
    print(
        f"[fig6] cloud n_src z clean-referenced: {np.round(cloud_z_cleanref, 2).tolist()} "
        f"(clean-only scale {c_scale:.2f} vs whole-night {n_scale_night:.2f})"
    )

    if cloud_flagged == 0:
        cloud_text = (
            "All five cloud frames **PASS** — the expected, and measured, "
            "blind spot. A grey 0.7 mag transparency loss leaves FWHM, sky "
            "level, sky noise, and star shapes untouched (median metrics are "
            "indistinguishable from clean, see table), so the shape, sky, "
            "noise-consistency, and depth checks are *structurally* blind to "
            "it: even the depth-cost proxy is unchanged because extinction "
            "costs *source* flux, not sky noise. The only metric that "
            "responds at all is the detected source count "
            f"(median $N_{{\\rm src}}$ {med.loc['cloud', 'n_sources']:.0f} vs "
            f"{med.loc['clean', 'n_sources']:.0f} clean — an LF-dependent "
            "deficit from the ~0.7 mag slice of stars pushed below the "
            f"detection limit), and its night-relative robust $z$ (median "
            f"{med.loc['cloud', 'nsrc_z']:.1f}) stayed far above the "
            f"$-{thr.nsrc_z_review}$ review cut for two honest reasons: "
            "(i) the deficit itself depends on the luminosity function near "
            "the limit, and (ii) with 20/44 frames defective, the MAD-based "
            "night scale is inflated by the defect population itself "
            f"(clean-only scale {c_scale:.1f} counts vs {n_scale_night:.1f} "
            f"whole-night; referenced to clean frames alone the same cloud "
            f"frames would sit at $z = {cloud_z_cleanref.min():.1f}$ to "
            f"${cloud_z_cleanref.max():.1f}$, i.e. mostly flaggable). Grey "
            "transparency loss is thus only robustly observable relative to a "
            "photometric reference; this measured limitation is precisely the "
            "motivation for the planned post-photometry transparency-QC stage "
            "(frame zeropoint / comparison-star flux monitoring)."
        )
    else:
        cloud_text = (
            f"{cloud_flagged}/5 cloud frames are flagged "
            f"({_fmt_reasons(reasons['cloud'])}): the flag rides on the "
            "source-count deficit produced by the 0.7 mag dimming pushing "
            "stars below the detection limit. This route is luminosity-"
            "function dependent (a field with few stars near the limit would "
            "evade it), so the planned post-photometry transparency QC remains "
            "necessary for robust cloud rejection."
        )

    if noisy_ratio_hits > 0:
        noisy_text = (
            f"the CCD-equation consistency check fired directly: "
            f"`sky_noise_ratio` (median {noisy_ratio_med:.2f}) exceeded the "
            f"{thr.sky_noise_ratio_review} review threshold on "
            f"{noisy_ratio_hits}/5 frames"
        )
    else:
        noisy_text = (
            f"the measured-vs-expected sky-noise ratio rose to a median of "
            f"{noisy_ratio_med:.2f}$\\times$ (clean frames: "
            f"{float(med.loc['clean', 'sky_noise_ratio']):.2f}) but stayed "
            f"below the conservative {thr.sky_noise_ratio_review} review "
            f"threshold at this sky level; the frames were flagged anyway, by "
            f"the *night-relative* sky-noise outlier check (`sky_sigma_z` $>$ "
            f"{thr.sky_z_review} on {noisy_sigz_hits}/5 frames — the same "
            "measured-noise-anomaly family, referenced to the night's own "
            "frames rather than to the CCD equation) and by the depth-cost "
            f"check (median {float(med.loc['noisy_readout', 'depth_cost_mag']):.2f} "
            "mag of estimated depth loss)"
        )

    def row(t: str) -> str:
        return (
            f"| {CLASS_STYLE[t]['label']} | {int(counts.loc[t].sum())} "
            f"| {int(counts.loc[t, PASS])} | {int(counts.loc[t, REVIEW])} "
            f"| {int(counts.loc[t, FAIL])} | {_fmt_reasons(reasons[t])} |"
        )

    def mrow(t: str) -> str:
        m = med.loc[t]
        depth = float(m["depth_cost_mag"])
        depth = 0.0 if abs(depth) < 0.005 else depth  # avoid a "-0.00" cell
        return (
            f"| {CLASS_STYLE[t]['label']} | {m['n_sources']:.0f} | {m['fwhm_med']:.2f} "
            f"| {m['sky_med']:.0f} | {m['sky_sigma']:.2f} | {m['fwhm_model_ratio']:.2f} "
            f"| {m['sky_noise_ratio']:.2f} | {depth:+.2f} |"
        )

    caption = f"""# Figure 6 — Frame-QC validation: injected frame defects vs automatic decisions

**Figure 6.** Validation of the automatic frame-quality-control module
(`apex/analysis/frame_qc.py`, `evaluate_frame_qc`) on a synthetic "night" of 44
frames with defects injected by construction. Each frame is an independent
realization of the self-contained synthetic reference field (170 Gaussian-PSF
stars, uniform true magnitudes in $[13.5, 19.0]$, $640^2$ px, Poisson photon +
Gaussian read noise in the electron domain, gain $g = {HEADER_GAIN_E_PER_ADU}$
e$^-$/ADU, zeropoint 25.0, FWHM 3.5 px, sky 150 ADU, RN 5 e$^-$; star
separation fixed at 21 px for *all* classes so every frame carries the same
170-star geometric truth). Twenty defective frames alter exactly one physical
property: **cloud** (zeropoint $-0.7$ mag; grey transparency loss, sky
unchanged), **bad seeing** (FWHM $\\times 1.8 \\to 6.3$ px), **bright sky**
(background $\\times 5 \\to 750$ ADU), and **noisy readout** (realized with RN
$= {TRUE_RDNOISE_NOISY_E:.0f}$ e$^-$ while the header — and therefore the QC
input — still claims {HEADER_RDNOISE_E:.0f} e$^-$: an undocumented-electronics
"header lie" that only a measured-vs-expected noise check can expose). Every
frame was processed by the **production Step-4 detection service** (the same
Qt-free `run_detection` code path the GUI executes; SEP engine, $\\sigma =
3.2$, repository `parameters.toml`), and the per-frame metrics were assembled
exactly as the production Step-4 window assembles them before calling
`evaluate_frame_qc`: `fwhm_med` is the detector's median radial-profile FWHM of
the brightest $\\leq 25$ sources, and `sky_med` / `sky_sigma` are the
detector-reported global background level and RMS (SEP `Background`
`globalback` / `globalrms`) — no auxiliary pixel-level estimate was
substituted. All frames enter QC with `filter = "r"`, `airmass = 1.0`,
`gain_e_per_adu = {HEADER_GAIN_E_PER_ADU}`, `rdnoise_e = {HEADER_RDNOISE_E}`;
the per-source-quality and candidate-count columns are omitted
(`evaluate_frame_qc` skips absent metrics NaN-safely). Because the airmass is
constant, the Kolmogorov seeing model falls back to the night-median FWHM, so
`fwhm_model_ratio` is the frame FWHM over the night median. Decisions use the
shipped conservative `FrameQCThresholds()` defaults (review/fail): FWHM $z$
{thr.fwhm_z_review}/{thr.fwhm_z_fail}, FWHM model ratio
{thr.fwhm_model_ratio_review}/{thr.fwhm_model_ratio_fail}, sky $z$
{thr.sky_z_review}/{thr.sky_z_fail}, sky-noise ratio
{thr.sky_noise_ratio_review}/{thr.sky_noise_ratio_fail}, $N_{{\\rm src}}$ $z$
$-{thr.nsrc_z_review}$/$-{thr.nsrc_z_fail}$, elongation
{thr.elong_review}/{thr.elong_fail}, depth cost
{thr.depth_cost_review}/{thr.depth_cost_fail} mag (estimated limiting-depth
loss vs the night median from $F_{{\\rm lim}} \\propto \\sigma_{{\\rm sky,e}}
\\cdot \\mathrm{{FWHM}}$). **(a)** Decision matrix, truth class vs
QC decision; green shading marks correct outcomes (clean $\\to$ PASS, defect
$\\to$ REVIEW or FAIL — the policy is deliberately conservative, reserving FAIL
for frames likely to poison downstream steps), red shading marks wrong ones.
**(b)** The two headline diagnostics added by `evaluate_frame_qc`: the seeing
ratio (`fwhm_model_ratio`, abscissa) against the measured-to-expected sky-noise
ratio (`sky_noise_ratio` $= \\sigma_{{\\rm sky}} g \\, / \\sqrt{{B g +
\\mathrm{{RN}}^2}}$, log ordinate), colored by truth class with marker shape
encoding the decision (circle PASS, triangle REVIEW, cross FAIL); grey lines
are the REVIEW (dashed) and FAIL (solid) thresholds of the two ratios.

## Per-class outcome (measured, defaults untouched)

| Truth class | N | PASS | REVIEW | FAIL | QC reasons fired (frames) |
|---|---|---|---|---|---|
{chr(10).join(row(t) for t in TRUTH_ORDER)}

Headline: **{n_clean_pass}/{n_clean} clean frames PASS** and
**{n_defect_flagged}/{n_defect} defective frames are flagged** (REVIEW or
FAIL) at the shipped default thresholds.

## Median per-frame metrics by class

| Truth class | $N_{{\\rm src}}$ | FWHM (px) | sky (ADU) | $\\sigma_{{\\rm sky}}$ (ADU) | FWHM ratio | noise ratio | depth cost (mag) |
|---|---|---|---|---|---|---|---|
{chr(10).join(mrow(t) for t in TRUTH_ORDER)}

## Which check caught what

* **bad seeing** — the only class driven to **FAIL**, by the FWHM checks:
  night-relative robust $z$ (`fwhm_z`, median
  {med.loc['bad_seeing', 'fwhm_z']:.1f} $\\gg$ {thr.fwhm_z_fail}) and the
  seeing-model ratio (median {med.loc['bad_seeing', 'fwhm_model_ratio']:.2f}
  $>$ {thr.fwhm_model_ratio_fail} FAIL cut); the depth-cost check
  (median {med.loc['bad_seeing', 'depth_cost_mag']:+.2f} mag) fired in
  support. Reasons: {_fmt_reasons(reasons['bad_seeing'])}.
* **bright sky** — caught by the sky-level outlier $z$-score (`sky_z`, median
  {med.loc['bright_sky', 'sky_z']:.0f} $>$ {thr.sky_z_review}; the enormous
  $z$ reflects the near-zero clean-frame scatter of the synthetic sky) plus
  the depth-cost check (median
  {med.loc['bright_sky', 'depth_cost_mag']:+.2f} mag). Its
  `sky_noise_ratio` stays at {med.loc['bright_sky', 'sky_noise_ratio']:.2f}
  because a genuinely brighter sky raises the CCD-equation expectation in
  step — the absolute-consistency and night-relative sky checks are
  complementary by design.
* **noisy readout (header lie)** — {noisy_text}.
* **cloud** — {cloud_text}

The {len(out)} production detections took {det_s:.0f} s of wall time
(single worker). Full per-frame metrics and decisions:
`validation/paper/data/frame_qc/frame_qc_night.csv`.
"""
    CAPDIR.mkdir(parents=True, exist_ok=True)
    cap_path = CAPDIR / "fig6_qc_validation.md"
    cap_path.write_text(caption, encoding="utf-8")
    return cap_path


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    t_start = time.perf_counter()
    reuse = "--reuse" in sys.argv
    DATADIR.mkdir(parents=True, exist_ok=True)
    metrics_path = DATADIR / "frame_metrics.csv"

    if reuse and metrics_path.exists():
        print(f"[fig6] --reuse: loading cached frame metrics from {metrics_path}")
        metrics = pd.read_csv(metrics_path)
    else:
        metrics = build_night_metrics()
        metrics.to_csv(metrics_path, index=False)
        print(f"[fig6] wrote frame metrics: {metrics_path}")

    out, thr = run_qc(metrics)
    out.to_csv(DATADIR / "frame_qc_night.csv", index=False)

    counts = (
        out.pivot_table(index="truth", columns="qc_status", aggfunc="size", fill_value=0)
        .reindex(index=TRUTH_ORDER, columns=DECISION_ORDER, fill_value=0)
        .astype(int)
    )
    sub = _subcheck_table(out, thr)

    fig = make_figure(out, counts, thr)
    paths = save_fig(fig, "fig6_qc_validation", OUTDIR)
    plt.close(fig)

    elapsed = time.perf_counter() - t_start
    cap_path = write_caption(out, counts, thr, sub)

    # ── verification report ──────────────────────────────────────────────────
    print("\n=== fig6 frame-QC validation ===")
    print("truth x decision matrix:")
    print(counts.to_string())
    print("\nper-class qc_reasons:")
    for truth in TRUTH_ORDER:
        print(f"  {truth:14s}: {dict(_reason_counter(out[out['truth'] == truth])) or '—'}")
    print("\nreview-level sub-checks fired (frames per class):")
    print(
        sub.groupby("truth").sum().reindex(TRUTH_ORDER).astype(int).to_string()
    )
    print("\nmedian diagnostics per class:")
    cols = ["n_sources", "fwhm_med", "sky_med", "sky_sigma", "fwhm_z", "sky_z",
            "sky_sigma_z", "nsrc_z", "fwhm_model_ratio", "sky_noise_ratio",
            "depth_cost_mag"]
    print(out.groupby("truth")[cols].median().reindex(TRUTH_ORDER).round(3).to_string())
    print(f"\ntotal runtime: {elapsed:.1f} s")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    print(f"wrote caption: {cap_path}  exists={cap_path.exists()}")

    n_clean_pass = int(counts.loc['clean', PASS])
    n_defect_flagged = int(
        counts.loc[[t for t in TRUTH_ORDER if t != 'clean'], [REVIEW, FAIL]].sum().sum()
    )
    print(f"clean->PASS: {n_clean_pass}/24 | defect->flagged: {n_defect_flagged}/20")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
