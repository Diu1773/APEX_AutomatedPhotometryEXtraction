"""Figure 13 — Cross-instrument + cross-pipeline validation on public LCO data.

APEX reduces raw frames from two different LCO cameras and is compared, pixel
for pixel, against the archive's BANZAI-processed product (an independent,
published pipeline) — a foreign-camera + foreign-pipeline check of the reduction:

  * QHY600 (CMOS, single-amplifier, 0.4 m): the whole frame agrees to a uniform
    ~0.06 e- offset (difference is featureless) — APEX reproduces BANZAI cleanly.
  * Sinistro (Fairchild CCD, 4 amplifiers, 1 m): the sky/source structure agrees
    to ~0.3%, but the difference shows a 4-quadrant pattern — the per-amplifier
    assembly (gain/overscan/cross-talk) that BANZAI does with dedicated Sinistro
    handling and a generic reduction does not. The calibration arithmetic
    generalises; multi-amp detector assembly is instrument-specific.

Data prepared by ``_make_lco_figdata.py`` (raw frames from archive.lco.global).
Run: .venv-deploy\\Scripts\\python validation\\paper\\fig13_cross_instrument.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO / "validation" / "paper"))
from apex_paper_style import apply_paper_style, save_fig, PALETTE, DOUBLE_COL

apply_paper_style()
DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def main() -> int:
    st = json.loads((DATA / "lco_crossinstrument.json").read_text())
    rows = [("qhy", "QHY600 · 0.4 m\n(CMOS, 1 amp)"),
            ("sinistro", "Sinistro · 1 m\n(CCD, 4 amps)")]
    z = ZScaleInterval()
    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE_COL, DOUBLE_COL * 0.72))
    col_titles = ["APEX-reduced", "BANZAI e91", "APEX − BANZAI"]

    for r, (key, label) in enumerate(rows):
        apex = np.load(DATA / f"lco_{key}_apex.npy")
        e91 = np.load(DATA / f"lco_{key}_e91.npy")
        diff = np.load(DATA / f"lco_{key}_diff.npy")
        for c, img in enumerate([apex, e91]):
            lo, hi = z.get_limits(img[np.isfinite(img)])
            axes[r, c].imshow(img, vmin=lo, vmax=hi, origin="lower", cmap="gray")
        dl = float(np.nanpercentile(np.abs(diff[np.isfinite(diff)]), 99))
        im = axes[r, 2].imshow(diff, vmin=-dl, vmax=dl, origin="lower", cmap="RdBu_r")
        fig.colorbar(im, ax=axes[r, 2], fraction=0.046, pad=0.03).set_label("e⁻", fontsize=6)
        for c in range(3):
            axes[r, c].set_xticks([]); axes[r, c].set_yticks([])
            if r == 0:
                axes[r, c].set_title(col_titles[c], fontsize=8.5)
        axes[r, 0].set_ylabel(label, fontsize=7.5)
        s = st[key]
        axes[r, 2].text(0.5, -0.09,
                        f"Δmed {s['delta_median']:+.3f} · σ {s['robust_sigma']:.3f} e⁻",
                        transform=axes[r, 2].transAxes, ha="center", va="top", fontsize=6.6)

    fig.suptitle("Cross-instrument + cross-pipeline: APEX vs LCO BANZAI on two public cameras",
                 fontsize=8.5, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    paths = save_fig(fig, "fig13_cross_instrument", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    q, s = st["qhy"], st["sinistro"]
    (CAPDIR / "fig13_cross_instrument.md").write_text(
        f"""# Figure 13 — Cross-instrument + cross-pipeline validation (LCO)

**Figure 13.** APEX reduces public raw frames from two different Las Cumbres
Observatory cameras and is compared pixel-for-pixel against the archive's
BANZAI-processed product — an independent, published pipeline — testing the
calibration on foreign detectors. **Top:** a QHY600 CMOS camera (single
amplifier, 0.4 m; Proxima Cen field). The whole frame agrees to a uniform
{q['delta_median']:+.3f} e⁻ offset (robust σ {q['robust_sigma']:.3f} e⁻); the
difference image is featureless. **Bottom:** a Sinistro CCD (four amplifiers,
1 m; NGC 5985 field). The sky and sources agree to ≈0.3 %, but the difference
shows a four-quadrant pattern (Δmedian {s['delta_median']:+.2f} e⁻,
σ {s['robust_sigma']:.2f} e⁻): the per-amplifier assembly — gain, overscan and
cross-talk — that BANZAI performs with dedicated Sinistro handling and a generic
reduction does not. The bias/dark/flat calibration arithmetic generalises across
cameras; multi-amplifier detector assembly is instrument-specific (APEX targets
single-CCD detectors). ZScale stretch; raw data from archive.lco.global.
""", encoding="utf-8")

    print("=== fig13 cross-instrument ===")
    for k, v in st.items():
        print(f"  {k}: Δmed {v['delta_median']:+.3f}  σ {v['robust_sigma']:.3f} e-")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
