# APEX paper validation figures

This directory contains the generators and final assets for the 17 figures in
the Korean manuscript. Figures use the production APEX calculation paths or
retained validation products. Synthetic experiments use fixed seeds. The
shared style is monochrome; categories are distinguished with line style,
marker shape, and hatching. Each final figure prints its data provenance.

| Fig. | Generator | Validation or application | Main input |
| ---: | --- | --- | --- |
| 1 | `fig_architecture.py` | Workflow and software layers | repository structure |
| 2 | `fig_calibration_step0.py` | Step-0 calibration data and effect | NGC 6811 calibration/science frames |
| 3 | `fig11_detector.py` | Gain, read noise, and dark current | retained detector-characterization data |
| 4 | `fig12_preproc_crosscheck.py` | APEX–ccdproc pixel arithmetic | retained NGC 6811 calibration products |
| 5 | `fig13_cross_instrument.py` | Cross-instrument calibration | LCO QHY600 and Sinistro products |
| 6 | `fig6_qc_validation.py` | Frame-QC decisions | fixed-seed synthetic night |
| 7 | `fig_completeness_realvssynth.py` | Detection completeness | artificial stars in seven real frames |
| 8 | `fig_detection_threshold.py` | Threshold and false-detection contamination | five real frames and Gaia checks |
| 9 | `fig_wcs_engines.py` | WCS engine comparison | eight M13 frames |
| 10 | `fig2_error_model.py` | Photometric uncertainty model | fixed-seed Monte Carlo data |
| 11 | `fig3_parameter_sweep.py` | Aperture, sky, and seeing sensitivity | fixed-seed parameter sweeps |
| 12 | `fig_photometry_crosschecks.py` | SEP and IRAF/DAOPHOT cross-checks | synthetic frame and NGC 6811 V frame |
| 13 | `fig_psf_validation.py` | PSF–aperture internal agreement | 67 frames from three cameras |
| 14 | `fig9_crowded_field.py` | Crowded-field internal agreement | retained M5 and M13 products |
| 15 | `fig_external_validation.py` | PS1 residuals and CMD comparison | NGC 6811 products and cached PS1 match |
| 16 | `fig_timeseries_validation.py` | PDM and SYSREM injection tests | fixed-seed synthetic time series |
| 17 | `fig_lc_yzboo.py` | YZ Boo end-to-end application | retained YZ Boo products |

Shared infrastructure:

- `apex_paper_style.py`: monochrome publication style and PNG/PDF saving.
- `_make_canonical_data.py`: canonical artificial-star data generator.
- `run_all.py`: runs the final generators in manuscript order and renders the
  Korean HTML preview.
- `render_preview.py`: maps manuscript figure numbers to final files and builds
  `MANUSCRIPT_ko_preview.html` and `MANUSCRIPT_ko_artifact.html`.

## How to run

Use the deploy environment and UTF-8 mode on Windows:

```powershell
.venv-deploy/Scripts/python.exe -X utf8 validation/paper/run_all.py
.venv-deploy/Scripts/python.exe -X utf8 validation/paper/run_all.py --only 4 12 15
.venv-deploy/Scripts/python.exe -X utf8 validation/paper/run_all.py --fast
```

`--fast` skips Figure 11's parameter sweep and uses the existing final asset
when rendering. The runner also passes UTF-8 mode to every child generator.

## Data and reproducibility notes

- Heavy `data*/` directories are regenerable or retained validation data and
  are gitignored. Scripts, captions, and final PNG/PDF assets are committed.
- Figures 6, 10, 11, and 16 use only synthetic data with fixed seeds. Figures
  2–5, 7–9, 12–15, and 17 require retained products or the external observation
  volume; a source-only checkout cannot rebuild them from raw data.
- Figure 12's IRAF comparison uses the retained fixed-coordinate benchmark.
  Recreating the IRAF measurements requires PyRAF in WSL.
- Figure 14 reads the retained Step-7/8 products. The raw FITS files are no
  longer present in the current archive, so the generator uses the known
  4800×3200 detector shape when needed. For M13 it falls back to
  `cmd_psf_backup_gui_20260729` if the primary PSF directory is empty.
- Figure 15 uses the cached PS1 cross-match at
  `data/ps1_match_ngc6811.csv`; no live catalogue query is needed.
- `validation/paper` is a junction. Generators use `Path.absolute()` for
  repository discovery; changing this to `resolve()` can redirect the inferred
  repository root to the junction target.
