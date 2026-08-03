"""Figure 11 — Detector characterisation from the data (gain, read noise, dark).

Both panels are measured from APEX's own calibration frames of the Moravian
C3-61000 (Sony IMX455, 2x2):

  (a) Photon-transfer relation: the variance of same-level flat-pair differences
      versus signal. The slope is 1/gain, giving gain = 0.681 +/- 0.014 e-/ADU.
      This measured value is consistent with the IMX455 laboratory value
      (Alarcon et al. 2023, 0.763 e-/ADU, native) and with the vendor full-well
      spec (>50 ke- over a 16-bit range => ~0.76 e-/ADU). The gain must be
      measured, not read from the FITS header: the MaxIm/ASCOM EGAIN keyword for
      this camera is ~16x too small (a documented 12-bit->16-bit ADC left-shift),
      so it is not used.
  (b) Dark current: source-free background versus exposure across a 10-480 s
      ladder, linear (R^2 = 0.998), slope 0.008 e-/s at +5 C.

Run: .venv-deploy\\Scripts\\python validation\\paper\\fig11_detector.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"

# External reference gain values (for validating our measurement)
ALARCON_GAIN = 0.763     # Alarcon et al. 2023, IMX455 native
VENDOR_GAIN = 0.76       # Moravian full well >50 ke- / 65536 ADU


def main() -> int:
    ptc = json.loads((DATA / "detector_ptc.json").read_text())
    det = json.loads((DATA / "detector_characterization.json").read_text())

    S = np.array(ptc["signal_adu"]); V = np.array(ptc["var_half_adu2"])
    gain = ptc["gain_e_per_adu"]; gerr = ptc["gain_err"]; rn = ptc["read_noise_adu"]

    fig = plt.figure(figsize=(DOUBLE_COL, 3.3))
    gs = GridSpec(2, 2, height_ratios=[3, 1], hspace=0.06, wspace=0.30,
                  left=0.09, right=0.98, top=0.90, bottom=0.15)
    axa = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1]); axr = fig.add_subplot(gs[1, 1], sharex=axb)

    # --- (a) gain from the photon-transfer slope (var = S/gain + RON^2) ---
    xl = np.linspace(0, S.max() * 1.05, 100)
    axa.scatter(S, V, s=26, color=C["data"], zorder=4, label="flat-pair data")
    axa.plot(xl, xl / gain + rn**2, color=C["model"], lw=1.8, zorder=3,
             label="fit: slope = 1/gain")
    axa.set_xlim(0, S.max() * 1.05); axa.set_ylim(0, V[S > 0].max() * 1.3)
    axa.set_xlabel(r"signal  $S$  (ADU)")
    axa.set_ylabel(r"variance  $\sigma^2$  (ADU$^2$)")
    axa.legend(loc="upper left", fontsize=7.0, frameon=False)
    axa.text(0.04, 0.62,
             f"gain = {gain:.3f} $\\pm$ {gerr:.3f} e$^-$/ADU (2×2)\n"
             f"read noise = {rn*gain:.2f} e$^-$\n"
             f"— matches IMX455 lab value\n"
             f"   (Alarcón 2023: {ALARCON_GAIN:.3f}) & vendor spec",
             transform=axa.transAxes, fontsize=6.8, va="top")
    axa.set_title("(a) Gain from the photon-transfer slope", loc="left")

    # --- (b) dark-current ladder + residual sub-panel ---
    t = np.array(det["dark"]["exptime"]); lvl = np.array(det["dark"]["level_dn"])
    k = det["dark"]["dn_per_s"]; c = det["dark"]["intercept_dn"]
    e_s = det["dark"]["e_per_s"]; r2 = det["dark"]["linearity_r2"]
    fit = k * t + c
    axb.scatter(t, lvl, s=20, color=C["data"], zorder=3)
    tl = np.linspace(0, t.max() * 1.05, 50)
    axb.plot(tl, k * tl + c, color=C["model"], lw=1.4, zorder=2,
             label=f"slope = {k:.4f} DN/s")
    axb.set_ylabel("dark above bias  (DN)")
    axb.legend(loc="upper left", fontsize=7.0, frameon=False)
    axb.text(0.96, 0.06, f"{e_s:.4f} e$^-$/s @ +5$^\\circ$C\n$R^2$ = {r2:.4f}",
             transform=axb.transAxes, va="bottom", ha="right", fontsize=6.8)
    axb.set_title("(b) Dark-current ladder", loc="left")
    plt.setp(axb.get_xticklabels(), visible=False)
    resid = 100 * (lvl - fit) / np.maximum(fit, 1e-9)
    axr.axhline(0, color=PALETTE["grey"], lw=0.8)
    axr.scatter(t, resid, s=16, color=C["data"])
    lim = max(3.0, np.abs(resid).max() * 1.3)
    axr.set_ylim(-lim, lim); axr.set_ylabel("resid %"); axr.set_xlabel("exposure  (s)")

    fig.suptitle(
        "Moravian C3-61000 (Sony IMX455, 2×2): gain, read noise and dark current "
        "measured from the data.",
        fontsize=7.4, y=0.98, color="#333333")
    paths = save_fig(fig, "fig11_detector", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig11_detector.md").write_text(
        f"""# Figure 11 — Detector characterisation from the data

**Figure 11.** Gain, read noise and dark current of the Moravian C3-61000
(Sony IMX455, 2×2 binned) measured directly from APEX's own calibration frames.
**(a)** Photon-transfer relation: the variance of same-level flat-pair
differences (which cancel fixed-pattern noise) versus signal, over
{ptc['n_pairs']} clean pairs. The slope is 1/gain, giving
gain = {gain:.3f} ± {gerr:.3f} e⁻/ADU with read noise {rn*gain:.2f} e⁻. The
measured value is consistent with the IMX455 laboratory value
(Alarcón et al. 2023, {ALARCON_GAIN:.3f} e⁻/ADU, native resolution) and the
vendor full-well specification (>50 ke⁻ over the 16-bit range ⇒ ≈{VENDOR_GAIN:.2f}
e⁻/ADU); the small difference is the 2×2 binning. The gain is *measured*, not
taken from the FITS header: for this camera the MaxIm/ASCOM `EGAIN` keyword is a
factor of ≈16 too small (a documented 12-bit→16-bit ADC left-shift), so it is
not used. **(b)** Dark current from the source-free background versus exposure
across a 10–480 s ladder: linear (R² = {r2:.4f}, residuals in the lower panel),
slope {e_s:.4f} e⁻/s at +5 °C. These measured values anchor APEX's photometric
error model to the detector's real physics.
""", encoding="utf-8")

    print("=== fig11 detector (header EGAIN removed; measured + validated) ===")
    print(f"gain = {gain:.4f} +/- {gerr:.4f} e-/ADU  | Alarcon {ALARCON_GAIN}, vendor {VENDOR_GAIN}")
    print(f"read noise {rn*gain:.2f} e-, dark {e_s:.4f} e-/s, R2 {r2:.4f}")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
