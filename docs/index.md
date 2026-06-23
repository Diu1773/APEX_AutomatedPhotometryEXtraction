# APEX — Automated Photometry EXtraction

APEX is a Python application for astronomical aperture and PSF photometry. It
runs both as a **PyQt5 desktop app** and as a **headless command-line pipeline**,
covering the full path from raw FITS frames to calibrated science products.

## Two analysis modes

- **CMD mode** — cluster photometry: source detection → WCS → forced/PSF
  photometry → zeropoint calibration → color-magnitude diagram → PARSEC
  isochrone fitting.
- **LC mode** — light-curve analysis: differential photometry, detrending,
  multi-night merge, and Lomb–Scargle / PDM / BLS period analysis.

Both modes share Steps 1–7 (file selection, crop, sky/QC, source detection, WCS
solving, master catalog build, forced aperture photometry).

## Validation

APEX's photometry is validated *directly* against independent software on real
frames — SExtractor (`sep`) and IRAF/DAOPHOT agree with APEX to **~3 mmag** (robust
MAD), and Gaia DR3 synthetic photometry confirms **no colour-dependent bias**. See
[photometry cross-validation](validation_crosscheck.md) (the primary measurement
validation, `apex validate --suite crosscheck`), [real-data cluster reproduction](validation_realdata.md)
(distance to ~4 % across 0.85–7.8 kpc), and the [isochrone-fitting methodology](isochrone_fitting_methodology.md)
(§10 — the gri/BVR degeneracy, the u-band cure, and what is/isn't a code limit).

## Install

```bash
pip install -e .            # headless core (servers / CI / batch)
pip install -e ".[gui]"     # + desktop GUI
pip install -e ".[dev]"     # + tests, build, docs tooling
```

## Run

```bash
apex doctor                 # check your environment
apex run --mode cmd --steps 1-7   # headless pipeline
apex gui                    # desktop launcher
apex export --format aavso --input lc.csv --output report.txt
```

See the [command-line interface](cli.md) and [configuration](configuration.md)
guides to get started, the [multi-night merger](multi-night-merger.md) guide for
combining separate single-night runs into one light curve, and the design notes
under *Design & internals* for how the pipeline works.

## Citing APEX

If you use APEX in published work, please cite it — see `CITATION.cff` in the
repository. APEX is released under the MIT license.
