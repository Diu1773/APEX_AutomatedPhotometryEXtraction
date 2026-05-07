# APEX Current Pipeline Contract

This file records the live step numbering and output layout after the forced
photometry refactor. Treat `main_window.py` and `apex/utils/step_paths*.py` as
the source of truth.

## Shared Steps

| UI Step | Name | Module | Output |
|---:|---|---|---|
| 1 | File Selection | `step1_file_selection*` | `step1_file_selection/` |
| 2 | Image Crop | `step2_crop_selector.py` | `step2_crop/` |
| 3 | Sky Preview & QC | `step3_sky_preview.py` | `step3_sky_preview/` |
| 4 | Source Detection | `step4_source_detection.py` | `step4_detection/` |
| 5 | WCS Plate Solving | `step5_wcs_plate_solving.py` | `step5_wcs/` |
| 6 | Master Catalog Build | `step6_ref_build.py` | `step6_refbuild/` |
| 7 | Forced Aperture Phot | `step7_forced_aperture_phot.py` | `step7_forced_phot/` |

## CMD Steps

| UI Step | Name | Module | Output |
|---:|---|---|---|
| 8 | PSF Photometry | `cmd/step8_psf_photometry.py` | `cmd_psf/` |
| 9 | Master ID Editor | `cmd/step9_master_id_editor.py` | `cmd_selection/` |
| 10 | Zeropoint Calibration | `cmd/step10_zeropoint_calibration.py` | `cmd_zeropoint/` |
| 11 | CMD Plot | `cmd/step11_cmd_plot.py` | `cmd_plot/` |
| 12 | Isochrone Model | `cmd/step12_isochrone_model.py` | `cmd_isochrone/` |

## LC Steps

| UI Step | Name | Module | Output |
|---:|---|---|---|
| 8 | Target/Comparison Selection | `lc/step8_target_selection.py` | `lc_selection/` |
| 9 | Light Curve Builder | `lc/step9_lightcurve_builder.py` | `lc_lightcurve/` |
| 10 | Detrend & Night Merge | `lc/step10_detrend_merge.py` | `lc_detrend/` |
| 11 | Period Analysis | `lc/step11_period_analysis.py` | `lc_period/` |

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
