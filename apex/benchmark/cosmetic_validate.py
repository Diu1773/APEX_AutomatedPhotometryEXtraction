"""Synthetic validation of cosmetic correction (cosmic-ray + hot pixel).

Injects a KNOWN population of cosmic rays and hot pixels into a synthetic star
field, runs :func:`apex.analysis.cosmetic.clean_frame` (L.A.Cosmic via
astroscrappy, van Dokkum 2001), and quantifies:

* cosmic-ray removal completeness (fraction of injected CR pixels flagged),
* star protection — false-positive rate on star cores (must be ~0),
* photometric preservation — aperture flux of undisturbed stars before vs after
  cleaning (must be unchanged; cleaning must not eat real starlight),
* hot-pixel removal.

Because astroscrappy IS the reference L.A.Cosmic implementation, this benchmark
demonstrates that APEX's cosmetic step removes the injected artefacts without
harming photometry — the citable claim for the paper.

Usage::

    python -m apex.benchmark.cosmetic_validate
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict

import numpy as np

from apex.analysis.cosmetic import clean_frame, HAS_ASTROSCRAPPY


SIZE = 512
GAIN = 1.5
READNOISE = 6.0
SKY = 300.0
N_STARS = 80
N_CR = 60
N_HOT = 120
SEED = 20260706


def _gaussian(img, x, y, flux, sigma):
    h, w = img.shape
    x0, x1 = max(0, int(x - 6 * sigma)), min(w, int(x + 6 * sigma) + 1)
    y0, y1 = max(0, int(y - 6 * sigma)), min(h, int(y + 6 * sigma) + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    img[y0:y1, x0:x1] += flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma ** 2))


def _aper_flux(img, xy, r=6.0):
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    out = []
    for (x, y) in xy:
        m = (xx - x) ** 2 + (yy - y) ** 2 <= r * r
        vals = img[m]
        out.append(float(np.nansum(vals[np.isfinite(vals)])))
    return np.array(out)


def run(seed: int = SEED) -> Dict:
    if not HAS_ASTROSCRAPPY:
        return {"error": "astroscrappy not installed", "PASS": False}
    rng = np.random.default_rng(seed)
    fwhm = 3.2
    sigma = fwhm / 2.3548

    # clean star field (truth)
    truth = np.full((SIZE, SIZE), SKY, dtype=np.float64)
    margin = 20
    star_xy = rng.uniform(margin, SIZE - margin, size=(N_STARS, 2))
    # genuinely-detectable stars (S/N >~ 30) so preservation/false-positive
    # metrics reflect real sources, not near-noise blips.
    fluxes = 10 ** rng.uniform(3.3, 4.5, size=N_STARS)   # ~2000 .. 32000
    for (x, y), f in zip(star_xy, fluxes):
        _gaussian(truth, x, y, f, sigma)
    # realise noise once; this noisy-but-clean frame is the preservation truth
    clean_noisy = rng.poisson(np.clip(truth * GAIN, 0, None)) / GAIN
    clean_noisy += rng.normal(0, READNOISE / GAIN, size=truth.shape)

    # inject cosmic rays in EMPTY sky (away from stars) so preservation is testable
    def _far_from_stars(x, y, d=8):
        return np.all((star_xy[:, 0] - x) ** 2 + (star_xy[:, 1] - y) ** 2 > d * d)
    dirty = clean_noisy.copy()
    cr_mask_true = np.zeros_like(dirty, dtype=bool)
    placed = 0
    while placed < N_CR:
        x = int(rng.integers(margin, SIZE - margin)); y = int(rng.integers(margin, SIZE - margin))
        if not _far_from_stars(x, y):
            continue
        length = int(rng.integers(1, 4))           # 1-3 px streak
        for k in range(length):
            xx = min(SIZE - 1, x + k)
            dirty[y, xx] += rng.uniform(3000, 40000)
            cr_mask_true[y, xx] = True
        placed += 1

    # inject hot pixels
    hot_mask_true = np.zeros_like(dirty, dtype=bool)
    for _ in range(N_HOT):
        x = int(rng.integers(0, SIZE)); y = int(rng.integers(0, SIZE))
        dirty[y, x] += rng.uniform(5000, 30000)
        hot_mask_true[y, x] = True

    # run cosmetic (no dark-derived hot mask here; astroscrappy must catch them)
    cleaned, corrected, ncorr = clean_frame(
        dirty, gain=GAIN, readnoise=READNOISE, satlevel=65535.0,
        sigclip=4.5, objlim=5.0, hot_mask=None)

    # metrics
    cr_completeness = float(np.count_nonzero(corrected & cr_mask_true) / max(1, cr_mask_true.sum()))
    hot_completeness = float(np.count_nonzero(corrected & hot_mask_true) / max(1, hot_mask_true.sum()))

    # star-core false positives (pixels within r=2 of any star)
    yy, xx = np.mgrid[0:SIZE, 0:SIZE]
    core = np.zeros((SIZE, SIZE), dtype=bool)
    for (x, y) in star_xy:
        core |= (xx - x) ** 2 + (yy - y) ** 2 <= 4.0
    star_fp = float(np.count_nonzero(corrected & core) / max(1, core.sum()))

    # photometric preservation: aperture flux of stars, cleaned vs clean truth
    f_truth = _aper_flux(clean_noisy, star_xy)
    f_clean = _aper_flux(cleaned, star_xy)
    dflux = (f_clean - f_truth) / np.clip(f_truth - SKY * np.pi * 6 ** 2, 1.0, None)
    dmag = 2.5 / np.log(10) * np.abs(dflux)
    med_dmag = float(np.median(dmag))
    max_dmag = float(np.percentile(dmag, 95))

    report = {
        "n_stars": N_STARS, "n_cr": int(cr_mask_true.sum()), "n_hot": N_HOT,
        "cr_completeness": cr_completeness,
        "hot_completeness": hot_completeness,
        "star_core_false_positive": star_fp,
        "flux_preservation_med_mmag": med_dmag * 1000.0,
        "flux_preservation_p95_mmag": max_dmag * 1000.0,
        "pixels_corrected": ncorr,
    }
    report["PASS"] = bool(
        cr_completeness > 0.90 and hot_completeness > 0.90
        and star_fp < 1e-3 and med_dmag * 1000.0 < 5.0
    )
    return report


def _print(r: Dict):
    print("APEX cosmetic (cosmic-ray + hot pixel) validation - L.A.Cosmic/astroscrappy")
    print("-" * 60)
    if r.get("error"):
        print("  ERROR:", r["error"]); return
    print(f"  injected: {r['n_stars']} stars, {r['n_cr']} CR px, {r['n_hot']} hot px")
    print(f"  CR removal completeness   : {r['cr_completeness']*100:.1f} %")
    print(f"  hot-pixel completeness    : {r['hot_completeness']*100:.1f} %")
    print(f"  star-core false positives : {r['star_core_false_positive']*100:.4f} %  (target ~0)")
    print(f"  flux preservation (median): {r['flux_preservation_med_mmag']:.3f} mmag")
    print(f"  flux preservation (p95)   : {r['flux_preservation_p95_mmag']:.3f} mmag")
    print("-" * 60)
    print(f"  RESULT: {'PASS' if r['PASS'] else 'FAIL'}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Cosmetic (CR + hot pixel) validation")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    r = run()
    _print(r)
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "cosmetic_validation.json").write_text(json.dumps(r, indent=2))
    return 0 if r.get("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
