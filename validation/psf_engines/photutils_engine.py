"""Measure the same implanted stars with photutils' own PSF photometry.

The DAOPHOT comparison has a circularity problem that is worth naming. Finding
that ALLSTAR weights pixels differently, copying its parameter, and then
observing that APEX now agrees with ALLSTAR would prove nothing except that two
programs were made to behave alike. Two things keep this benchmark out of that
trap, and only the first is automatic:

* the metric is |measured - injected truth|, not |APEX - ALLSTAR|, so an
  improvement is an improvement in absolute accuracy;
* a *third* engine, developed independently of DAOPHOT, breaks the two-body
  comparison — if APEX, photutils and ALLSTAR disagree in a structured way,
  the structure is about the algorithms and not about one lineage.

photutils is the natural third engine here, and the choice is not neutral: APEX
already depends on it for detection and aperture photometry, so "why write a
PSF fitter at all instead of calling `photutils.psf`?" is a fair question that
this script answers with numbers rather than argument.

Fairness is enforced by giving photutils exactly what APEX had:

* **the same PSF** — APEX's own saved ePSF array, wrapped as an `ImagePSF` with
  the oversampling recorded in its header, so neither engine benefits from a
  better model;
* **the same positions** — the step-7 forced catalog that seeds APEX's step 8,
  passed as `init_params` with no finder, so this compares fitting and not
  detection;
* **the same fit footprint** — `fit_shape` set to APEX's window;
* **the same noise model** — the error image is built from the background RMS
  and gain APEX used.

What differs is only what is being tested: how each engine solves for flux in
the presence of neighbours.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))


def load_psf(path: Path):
    """APEX's saved ePSF as a photutils model, normalised to unit flux."""
    from photutils.psf import ImagePSF

    with fits.open(path) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        oversampling = int(hdul[0].header.get("OVERSAMPL", 1))
    # ImagePSF integrates data/oversampling**2; normalise so flux is in counts.
    total = data.sum() / oversampling ** 2
    if total > 0:
        data = data / total
    return ImagePSF(data, oversampling=oversampling)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--apex-run", default="apexfix",
                    help="run whose saved ePSF and step-7 seeds are reused")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--fit-shape", type=int, default=19)
    ap.add_argument("--gain", type=float, default=0.68)
    ap.add_argument("--background-rms", type=float, default=29.10280458)
    ap.add_argument("--grouper-fwhm", type=float, default=0.0,
                    help="group radius in FWHM; 0 disables simultaneous fitting")
    ap.add_argument("--fwhm-px", type=float, default=7.052391)
    args = ap.parse_args()

    from photutils.psf import PSFPhotometry, SourceGrouper

    work = Path(args.work)
    for trial in range(1, args.trials + 1):
        # Two layouts exist: the engine-variant runs are `<run>_trial<N>`, while
        # the injector writes `trial_000N`. Accept either so a benchmark built
        # on a new target does not need its own runner.
        candidates = [work / f"{args.apex_run}_trial{trial}" / "result" / "cmd_psf",
                      work / f"trial_{trial:04d}" / "result" / "cmd_psf"]
        psf_files: list[Path] = []
        for cmd_dir in candidates:
            psf_files = sorted(cmd_dir.glob("epsf_model_*.fits"))
            if psf_files:
                break
        frame_path = work / f"trial_{trial:04d}" / "data" / args.frame
        seeds = (work / f"trial_{trial:04d}" / "result" / "step7_forced_phot"
                 / f"photometry_{args.frame}.tsv")
        if not psf_files or not frame_path.exists() or not seeds.exists():
            print(f"[건너뜀] trial {trial}: 입력 없음")
            continue

        image = fits.getdata(frame_path).astype(float)
        # APEX fits a background-subtracted image; match that.
        from astropy.stats import sigma_clipped_stats
        _, sky, _ = sigma_clipped_stats(image, sigma=3.0, maxiters=5)
        data = image - sky
        error = np.sqrt(args.background_rms ** 2
                        + np.clip(data, 0.0, None) / max(args.gain, 1e-6))

        s7 = pd.read_csv(seeds, sep="\t")
        x = pd.to_numeric(s7["x_fit"], errors="coerce").to_numpy()
        y = pd.to_numeric(s7["y_fit"], errors="coerce").to_numpy()
        flux = pd.to_numeric(s7.get("flux_net_adu", s7.get("flux")),
                             errors="coerce").to_numpy()
        ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(flux) & (flux > 0)
        init = Table({"x": x[ok], "y": y[ok], "flux": flux[ok]})

        grouper = (SourceGrouper(min_separation=args.grouper_fwhm * args.fwhm_px)
                   if args.grouper_fwhm > 0 else None)
        phot = PSFPhotometry(load_psf(psf_files[0]), fit_shape=args.fit_shape,
                             finder=None, grouper=grouper, progress_bar=False)
        started = time.perf_counter()
        result = phot(data, error=error, init_params=init)
        elapsed = time.perf_counter() - started

        out = pd.DataFrame({
            "x_fit": np.asarray(result["x_fit"], dtype=float),
            "y_fit": np.asarray(result["y_fit"], dtype=float),
            "flux_fit": np.asarray(result["flux_fit"], dtype=float),
            "flux_err": np.asarray(result["flux_err"], dtype=float),
            "flags": np.asarray(result["flags"], dtype=int),
        })
        out["mag_psf"] = np.where(out["flux_fit"] > 0,
                                  -2.5 * np.log10(out["flux_fit"].clip(lower=1e-12)),
                                  np.nan)
        dest = work / f"photutils_trial{trial}.csv"
        out.to_csv(dest, index=False)
        print(f"  trial {trial}: {len(out)}개 · 유한 등급 "
              f"{int(np.isfinite(out['mag_psf']).sum())}개 · {elapsed:.1f}s "
              f"· grouper={'off' if grouper is None else f'{args.grouper_fwhm:g}FWHM'}"
              f" -> {dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
