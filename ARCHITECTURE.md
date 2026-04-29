# APEX Architecture

## Package Overview

```
apex/
  core/                    Shared core: ProjectState, FileManager, InstrumentConfig
  config/
    parameters_cmd.py      Parameters class for CMD mode (TOML-backed)
    parameters_lc.py       Parameters class for LC mode (TOML-backed)
    schema.py              TOML schema validators
  utils/
    step_paths.py          Shared step dir helpers (steps 1-8)
    step_paths_cmd.py      CMD-specific step dir helpers (steps 9-13)
    step_paths_lc.py       LC-specific step dir helpers (steps 9-12)
    photometry_utils.py    Aperture photometry utilities
    astro_utils.py         Airmass, BJD, filter normalization
    io_utils.py            CSV/ECSV int64 safe readers
    cache_utils.py         Detection/WCS cache management
    param_loader.py        TOML parameter file resolution
    photometry_loader.py   Frame photometry CSV loader (LC)
    run_workspace.py       Multi-night workspace helpers (LC)
    ... (logging, constants, qc, header_cache, common_helpers)
  analysis/
    light_curve/           LC science: lightcurve, detrend, period, eclipse, asteroid
    merge/                 Multi-night workspace scan/build/id-match (LC)
    cmd/                   CMD science: isochrone_fitter, isochrone_fitter_v2
  gui/
    main_window.py         UNIFIED main window (mode="cmd" or "lc")
    widgets/
      image_viewer.py      Zoomable FITS image viewer
    workflow/
      step_window_base.py  Base class for all step windows
      step2_crop_selector.py   }
      step3_sky_preview.py     }  Shared steps 2-8
      step4_source_detection.py}  (same window for both modes)
      step5_aperture_worker.py }
      step5_aperture_photometry.py}
      step6_wcs_plate_solving.py }
      step7_ref_build.py        }
      step8_star_id_matching.py }
      cmd/
        step1_file_selection.py   CMD file selection
        step6_psf_photometry.py   PSF photometry (CMD step 5, index 5)
        step10_master_id_editor.py
        step11_zeropoint_calibration.py
        step12_cmd_plot.py
        step13_isochrone_model.py
      lc/
        step1_file_selection.py   LC file selection (multi-night aware)
        step9_target_selection.py
        step10_lightcurve_builder.py
        step11_detrend_merge.py
        step12_period_analysis.py
    tools/
      extinction_fit.py     Bouguer extinction / zeropoint fitting
      iraf_photometry.py    IRAF/DAOPHOT integration
      iraf_comparison.py    IRAF comparison photometry
      qa_report.py          QA / publication validation report
      airmass_debug.py      Airmass header diagnostics
      aperture_overlay.py   Aperture overlay visualizer
      variable_star.py      Variable star classification (LC)
      multi_night_merger.py Multi-night LC merge tool (LC)
      transit_tool.py       Exoplanet transit fitting (LC)
      eb_tool.py            Eclipsing binary fitting (LC)
      gaia_3d_viewer.py     Gaia 3D cluster viewer (CMD)
      cmd_iso_tool.py       CMD + isochrone from results (CMD)
      cluster_structure/    Cluster structure analysis (CMD)
```

## Mode Concept

Both modes launch from `MainWindowWorkflow(mode=...)`:

- **Shared steps (1-8)** use the same window classes regardless of mode.
- **CMD-only steps**: PSF photometry (step 5), master ID editor, zeropoint, CMD plot, isochrone.
- **LC-only steps**: target/comparison selection, light curve builder, detrend/merge, period analysis.

Step index → file dispatch is in `_open_step_window(step_index)` of `main_window.py`.

## Adding a New Shared Step

1. Create `apex/gui/workflow/step_N_xxx.py` inheriting `StepWindowBase`.
2. Add a `stepN_xxx_dir()` helper to `apex/utils/step_paths.py`.
3. Wire up in `main_window._open_step_window()` for both modes.

## Adding a New Mode-Specific Step

1. Create `apex/gui/workflow/cmd/stepN_xxx.py` or `apex/gui/workflow/lc/stepN_xxx.py`.
2. Add path helpers to `step_paths_cmd.py` or `step_paths_lc.py`.
3. Add the step name to the appropriate `step_names` list and wire dispatch in `main_window.py`.

## Data Flow

```
params.P.data_dir/          Raw FITS files
params.P.result_dir/
  step1_file_selection/     FITS scan manifest
  step2_crop/               Crop region + cropped images
  step3_sky_preview/        Sky QC metadata
  step4_detection/          Source catalogs + frame QC
  step5_aperture/           Aperture photometry CSVs
  step6_wcs/                WCS-solved FITS headers
  step7_refbuild/           Master star catalog, Gaia IDs
  step8_idmatch/            Per-frame star ID matches
  [cmd_*/ or lc_*/]         Mode-specific outputs
  cache/                    Intermediate caches (header scan, detect, WCS)
```

## Step Directory Conventions

Each step writes to a named subdirectory of `result_dir`, defined in `step_paths*.py`.
- Shared names (step1-8): defined in `step_paths.py`.
- CMD names (cmd_psf, cmd_selection, cmd_zeropoint, cmd_plot, cmd_isochrone): in `step_paths_cmd.py`.
- LC names (lc_selection, lc_lightcurve, lc_detrend, lc_period): in `step_paths_lc.py`.
