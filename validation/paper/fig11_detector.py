"""Figure 11 — Detector characterisation from the data (gain, read noise, dark).

The FITS header of these Moravian C3-61000 frames records the nominal MaxIm gain
(EGAIN = 0.0495 e-/ADU, the max-gain register value), but the actual delivered
gain differs. We measure it from first principles:

  (a) Photon-transfer curve: variance of same-level flat-pair differences vs
      signal, over 12 clean pairs. The slope gives gain = 0.681 +/- 0.014 e-/ADU
      (2x2 stored pixel); the header value (0.0495) would be a 14x steeper line
      (46 sigma away) and is decisively ruled out. Read noise from a bias-pair
      difference is 3.45 ADU = 2.35 e- (stored) / 1.18 e- (native), matching the
      IMX455 spec (Alarcon et al. 2023).
  (b) Dark current: source-free background vs exposure across a 10-480 s ladder,
      linear (R^2 = 0.998), slope = 0.008 e-/s at +5 C.

This is why gain must be measured (PTC), not read from the header, and anchors
the photometric error model to the detector's real physics.

Run: .venv-deploy\\Scripts\\python validation\\paper\\fig11_detector.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def main() -> int:
    ptc = json.loads((DATA / "detector_ptc.json").read_text())
    det = json.loads((DATA / "detector_characterization.json").read_text())

    S = np.array(ptc["signal_adu"]); V = np.array(ptc["var_half_adu2"])
    gain = ptc["gain_e_per_adu"]; gerr = ptc["gain_err"]
    rn_adu = ptc["read_noise_adu"]; egain = ptc["header_egain"]
    slope = 1.0 / gain

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.2))

    # (a) PTC
    flat = S > 0
    axa.scatter(S[flat], V[flat], s=16, color=C["data"], zorder=3, label="flat pairs")
    axa.scatter([0], [V[~flat][0] if (~flat).any() else rn_adu**2], s=22,
                marker="s", color=C["reference"], zorder=4, label="bias (read noise)")
    xline = np.linspace(0, S.max() * 1.05, 100)
    axa.plot(xline, slope * xline + (V[~flat][0] if (~flat).any() else rn_adu**2),
             color=C["model"], lw=1.5, zorder=2,
             label=f"fit: gain = {gain:.3f} e$^-$/ADU")
    # what the header gain would imply (14x steeper)
    axa.plot(xline, (1.0 / egain) * xline, color=C["bad"], lw=1.2, ls="--", zorder=1,
             label=f"header EGAIN {egain:.3f} (46$\\sigma$ off)")
    axa.set_xlabel("signal  $S$  (ADU)")
    axa.set_ylabel(r"$\frac{1}{2}\,\mathrm{Var}(\mathrm{flat}_1-\mathrm{flat}_2)$  (ADU$^2$)")
    axa.set_ylim(0, V[flat].max() * 1.25)
    axa.set_xlim(0, S.max() * 1.05)
    axa.legend(loc="upper left", fontsize=6.6)
    axa.set_title("(a) Photon-transfer curve → gain", loc="left")
    axa.text(0.97, 0.05,
             f"gain = {gain:.3f} $\\pm$ {gerr:.3f} e$^-$/ADU\n"
             f"read noise = {rn_adu*gain:.2f} e$^-$ (stored)",
             transform=axa.transAxes, va="bottom", ha="right", fontsize=6.8,
             bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                   "alpha": 0.85, "edgecolor": PALETTE["grey"]})

    # (b) dark current ladder
    t = np.array(det["dark"]["exptime"]); lvl = np.array(det["dark"]["level_dn"])
    dn_s = det["dark"]["dn_per_s"]; e_s = det["dark"]["e_per_s"]; r2 = det["dark"]["linearity_r2"]
    axb.scatter(t, lvl, s=18, color=C["data"], zorder=3)
    tl = np.linspace(0, t.max() * 1.05, 50)
    axb.plot(tl, dn_s * tl + det["dark"]["intercept_dn"], color=C["model"], lw=1.4,
             zorder=2, label=f"slope = {dn_s:.4f} DN/s")
    axb.set_xlabel("exposure  (s)")
    axb.set_ylabel("dark background above bias  (DN)")
    axb.legend(loc="upper left", fontsize=7.0)
    axb.set_title("(b) Dark current ladder", loc="left")
    axb.text(0.97, 0.05,
             f"dark = {e_s:.4f} e$^-$/s @ +5$^\\circ$C\n$R^2$ = {r2:.4f}",
             transform=axb.transAxes, va="bottom", ha="right", fontsize=6.8,
             bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                   "alpha": 0.85, "edgecolor": PALETTE["grey"]})

    fig.suptitle(
        "Moravian C3-61000 (Sony IMX455, 2×2): gain/read-noise/dark measured from data, "
        "not the header (which records the nominal max-gain value).",
        fontsize=7.2, y=1.02, color="#333333")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    paths = save_fig(fig, "fig11_detector", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig11_detector.md").write_text(
        f"""# Figure 11 — Detector characterisation from the data

**Figure 11.** Gain, read noise, and dark current of the Moravian C3-61000
(Sony IMX455, 2×2 binned) measured directly from APEX's own calibration frames.
**(a)** Photon-transfer curve: the variance of same-level flat-pair differences
(which cancel PRNU and vignette) versus signal, over {ptc['n_pairs']} clean pairs
spanning {S[S>0].min():.0f}–{S.max():.0f} ADU. The slope gives
gain = {gain:.3f} ± {gerr:.3f} e⁻/ADU (stored pixel); the read noise from a
bias-pair difference is {rn_adu*gain:.2f} e⁻ (stored) / {rn_adu*gain/2:.2f} e⁻
(native), consistent with the IMX455 laboratory value (Alarcón et al. 2023). The
FITS-header EGAIN ({egain:.4f} e⁻/ADU, the nominal MaxIm max-gain value) implies
a 14× steeper line and is ruled out at 46σ — so the gain must be measured, not
read from the header. **(b)** Dark current from the source-free background versus
exposure across a 10–480 s ladder: linear (R² = {r2:.4f}), slope
{e_s:.4f} e⁻/s at +5 °C. These measured values anchor APEX's photometric error
model to the detector's real physics.
""", encoding="utf-8")

    print("=== fig11 detector characterisation ===")
    print(f"gain = {gain:.4f} +/- {gerr:.4f} e-/ADU  (config 0.689 within {abs(0.689-gain)/gerr:.1f} sigma)")
    print(f"header EGAIN {egain} ruled out at {abs(egain-gain)/gerr:.0f} sigma")
    print(f"read noise {rn_adu*gain:.2f} e- stored / {rn_adu*gain/2:.2f} e- native")
    print(f"dark {e_s:.4f} e-/s, R2 {r2:.4f}")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
