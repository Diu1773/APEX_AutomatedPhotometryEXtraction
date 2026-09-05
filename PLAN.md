# APEX Current Pipeline Contract

This file records the live step numbering and output layout. Treat
`apex/gui/main_window.py` and `apex/utils/step_paths*.py` as the source of
truth; this table is a convenience copy and must be updated with them.

Last verified against the code: 2026-09-05.

## Off-chain Step 0

Detector calibration (bias / dark / flat) is an **optional off-chain step**. It
is a `ToolWindowBase`, not a numbered step window, and `main_window.py` renders
it above the chain. It turns raw frames into science frames before Step 1.

| Name | Module | Output |
|---|---|---|
| Detector Calibration | `step0_detector_calibration.py` | user-chosen output dir |

## Shared Steps (both modes)

| UI Step | Name | Module | Output |
|---:|---|---|---|
| 1 | File Selection | `step1_file_selection*` | `step1_file_selection/` |
| 2 | Image Crop | `step2_crop_selector.py` | `step2_crop/` |
| 3 | Sky Preview & QC | `step3_sky_preview.py` | `step3_sky_preview/` |
| 4 | Source Detection | `step4_source_detection.py` | `step4_detection/` |
| 5 | WCS Plate Solving | `step5_wcs_plate_solving.py` | `step5_wcs/` |
| 6 | Master Catalog Build | `step6_ref_build.py` | `step6_refbuild/` |
| 7 | Forced Aperture Phot | `step7_forced_aperture_phot.py` | `step7_forced_phot/` |
| 8 | PSF Photometry | `cmd/step8_psf_photometry.py` | `cmd_psf/` |

Step 8 is shared: required in CMD mode, optional and window-only in LC mode.
The two modes branch after it.

## CMD Steps

| UI Step | Name | Module | Output |
|---:|---|---|---|
| 9 | Master ID Editor | `cmd/step9_master_id_editor.py` | `cmd_selection/` |
| 10 | Zeropoint Calibration | `cmd/step10_zeropoint_calibration.py` | `cmd_zeropoint/` |
| 11 | CMD Plot | `cmd/step11_cmd_plot.py` | `cmd_plot/` |
| 12 | Isochrone Model | `cmd/step12_isochrone_model.py` | `cmd_isochrone/` |

## LC Steps

**The four LC file names are one lower than their UI step numbers.** An
optional PSF window was inserted at Step 8 on 2026-07-15 and the files were not
renamed. Read `step_index=` inside the file and add one; never infer the step
number from an LC file name.

| UI Step | Name | Module (`step_index`) | Output |
|---:|---|---|---|
| 9 | Target/Comparison Selection | `lc/step8_target_selection.py` (8) | `lc_selection/` |
| 10 | Light Curve Builder | `lc/step9_lightcurve_builder.py` (9) | `lc_lightcurve/` |
| 11 | Detrend & Night Merge | `lc/step10_detrend_merge.py` (10) | `lc_detrend/` |
| 12 | Period Analysis | `lc/step11_period_analysis.py` (11) | `lc_period/` |

Every other step window follows `stepN_*.py` ↔ `step_index = N − 1`, so
`step7_forced_aperture_phot.py` holds `step_index = 6`. The helper functions in
`step_paths_lc.py` (`step8_selection_dir()` and friends) also keep their
pre-2026-07-15 names as stable APIs; the module docstring carries the current
table.

## Step 7 Forced Photometry Outputs

`step7_forced_phot/` owns the science photometry table:

- `photometry_{fname}.tsv`
- `photometry_index.csv`
- `apcorr_summary.csv`
- `centering_stats.csv`
- `filter_frames.json`
- `master_sources.csv`
- `frame_stats.csv`

Do not reintroduce `step7_refbuild/`, `step6_wcs/`, `step_forced_phot/`,
`step5_aperture/`, or `step8_idmatch/` as current output paths. Legacy readers
may keep fallback support for old result folders, but new writes should use the
current names above.
