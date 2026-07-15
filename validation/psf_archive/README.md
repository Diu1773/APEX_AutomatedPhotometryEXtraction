# APEX PSF validation archive (2026-07-15)

This directory is the compact, commit-safe backup of the CPU-only Step 8 PSF
development and validation work. It keeps code-generated figures, aggregate
tables, compact source catalogues, ePSF models, and run metadata. The original
intermediate FITS/NPY products under `validation/real_gui_run` total about
18.6 GB and are intentionally not duplicated here.

## Contents

- `aperture_psf/`: M13 and NGC 6811 same-frame aperture/PSF agreement,
  reported-error diagnostics, matched catalogues, and binned statistics.
- `figures/real_gui_run/`: final and aggregate M3/M13/M5 figures, including
  fit-window A/B, ePSF contamination, artificial-star, repeatability, seeing,
  residual, core-policy, and GUI diagnostic figures. Per-trial duplicates and
  the invalid M60 cluster test are excluded.
- `figures/paper/`: PSF iteration, residual sequence, and NGC 6811 legacy
  repeatability figures in PNG/PDF form.
- `summaries/real_gui_run/`: aggregate artificial-star recovery tables,
  A/B summaries, diagnostic JSON, and seeing/error-bin tables.
- `runs/`: compact current-engine M13, M5, and NGC 6811 Step 4/Step 8
  catalogues, ePSF models, reference-star tables, residual metadata, and
  one-CPU run metadata. Full-size image/model/residual arrays are excluded.

## Current-engine policy

- One CPU worker; no GPU dependency.
- Per-frame ePSF with contamination-aware, spatially balanced references.
- Automatic fit window targeting 90% ePSF encircled energy.
- Local residual-source refit followed by a short full-catalogue final pass.
- `qfit/noise <= 3`, flux-scale correction off, grouper off, hard core cut off.
- Fit all usable sources and flag unresolved blends rather than masking the
  complete cluster core.

## Aperture versus PSF error result

The comparison removes only a constant offset measured from bright isolated
stars. Aperture photometry is not treated as truth. The reported combined
error is `sqrt(sigma_psf^2 + sigma_aperture^2)`; this is an approximation
because both measurements use the same pixels.

| Cluster | Type | N (SNR >= 5) | robust scatter | median combined error | scatter/error | 1-sigma coverage |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| M13 | globular | 672 | 0.0458 mag | 0.0435 mag | 1.05 | 57.0% |
| NGC 6811 | open | 815 | 0.0346 mag | 0.0280 mag | 1.24 | 47.9% |

For M13, the PSF-aperture scatter rises from 0.0241 mag at neighbor distances
of at least 6 FWHM to 0.1421 mag inside 1.5 FWHM. This confirms that the
asymmetric tail is primarily crowding/background sensitivity, not a global
ePSF flux offset. For NGC 6811, only four matched stars are inside 1.5 FWHM,
so no strong crowding conclusion is supported.

Both reported per-source error curves closely follow `1.0857 / SNR`. They do
not include a separate high-SNR systematic floor. In NGC 6811 the SNR >= 100
PSF-aperture scatter is about 3.2 times the quadrature error, showing that a
frame/method systematic term is still needed for precision applications. This
cross-method excess does not identify PSF photometry alone as the cause.

## Current one-CPU frame checks

- M13: 789/805 clean fits (98.0%), median qfit/noise 0.955, reduced chi-square
  1.095, hard core cut disabled.
- NGC 6811: 838/868 clean fits (96.5%), qfit/noise 0.981, reduced chi-square
  1.092, QC PASS, 54.5 s for the Step 8 frame worker.
- M5 good/mid/poor seeing: clean fractions 97.6%, 90.9%, and 90.2%; the PSF
  model remains stable while the adaptive footprint and CPU time grow with
  seeing.

The cross-cluster artificial-star tests found an information limit near
1.5 FWHM: recovery was effectively zero inside that separation and generally
92-100% outside it. High-precision calibration samples should use at least
6 FWHM neighbor separation.

## Verification

- Full repository suite: 529 passed.
- Final PSF/detection/QC subset after review fixes: 38 passed.
- `python -m compileall apex main.py`: passed.
- `git diff --check`: passed before archive creation.

Regenerate the new comparison with:

```powershell
python validation/compare_aperture_psf_errors.py
```

The default command records the external source paths in
`aperture_psf/aperture_psf_error_manifest.json`; the matched compact tables in
this archive allow the plotted values to be audited without the 18.6 GB raw
validation tree.
