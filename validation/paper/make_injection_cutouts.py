# -*- coding: utf-8 -*-
"""Inject a magnitude ladder of artificial stars into the REAL M13 frame using the
exact pipeline injection path (inject_flux_catalog + empirical PSF), then save
postage-stamp cutouts — the AutoPhOT Fig.12 / DES Balrog illustration."""
import sys
from pathlib import Path
REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
import numpy as np, pandas as pd
from astropy.io import fits
from apex.benchmark.psf_artificial_stars import inject_flux_catalog

FRAME = r"E:\APEX_validation\reprocess\M13\calibrated\20260515\pp_messier13-0001-V.fit"
PSF = REPO / "validation/paper/data_realframe_M13V/artificial_star/benchmark_run/empirical_psf.fits"
OUT = REPO / "validation/paper/data_realframe_M13V/injection_cutouts.npz"
GAIN = 0.689
ZP = 25.0
MAGS = [13.0, 14.0, 14.5, 15.0, 15.5, 16.0]
HALF = 20  # stamp half-size → 41x41

data = fits.getdata(FRAME).astype(float)
kernel = fits.getdata(PSF).astype(float)

# pick a calm region away from the crowded core: use a low-background patch.
# scan a coarse grid for the calmest 400x400 block (low median + low MAD).
best = None
H, W = data.shape
for yy in range(200, H - 600, 400):
    for xx in range(200, W - 600, 400):
        blk = data[yy:yy + 400, xx:xx + 400]
        med = np.median(blk); mad = np.median(np.abs(blk - med))
        score = med + 5 * mad
        if best is None or score < best[0]:
            best = (score, xx, yy, med, mad)
_, bx, by, bmed, bmad = best
print(f"calm region @ x={bx} y={by} median={bmed:.0f} mad={bmad:.1f}")

# lay the ladder horizontally, well separated, in the calm block
xs = [bx + 60 + i * 62 for i in range(len(MAGS))]
ys = [by + 200 for _ in MAGS]
cat = pd.DataFrame({
    "x_true": xs, "y_true": ys,
    "true_flux_e": [10 ** ((ZP - m) / 2.5) for m in MAGS],
    "mag": MAGS,
})
rng = np.random.default_rng(7)
injected, _, _, out = inject_flux_catalog(
    data, kernel, cat, gain_e_per_adu=GAIN, rng=rng, return_layers=False)

# shared stretch from the calm block background (so faint really looks faint)
lo = bmed - 1 * bmad
hi = bmed + 15 * bmad   # tuned for contrast: bright core saturates, faint = smudge
stamps = []
for x, y, m in zip(xs, ys, MAGS):
    s = injected[y - HALF:y + HALF + 1, x - HALF:x + HALF + 1].astype(float)
    stamps.append(s)
stamps = np.array(stamps)
np.savez(OUT, stamps=stamps, mags=np.array(MAGS), lo=lo, hi=hi,
         m50=14.9)
print("saved", OUT, "stamps", stamps.shape, "stretch", round(lo), round(hi))
# quick peak S/N sanity
for m, s in zip(MAGS, stamps):
    pk = s.max() - bmed
    print(f"  m={m}: peak_over_bkg={pk:.0f} ADU  (bkg_mad={bmad:.1f})")
