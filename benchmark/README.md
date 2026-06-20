# APEX Automated Artificial-Star Benchmark

This benchmark inserts low-density artificial stars into a copy of an existing
FITS observation and reruns APEX production code.

```powershell
python benchmark/run_benchmark.py `
  --config benchmark/configs/baseline.toml `
  --input "E:\path\to\science.fit"
```

The runner performs the following operations automatically:

1. Run the production Step 4 detector on the unmodified image.
2. Build an empirical PSF by median-stacking isolated real stars.
3. Create deterministic subpixel injections over multiple low-density trials.
4. Add source Poisson noise without adding sky or read noise a second time.
5. Run the same production Step 4 detector on every injected image.
6. Measure original and injected images with Step 7's shared
   `phot_vectorized` aperture implementation.
7. Produce truth matches, completeness, position error, differential flux
   bias, plots, and a version/hash manifest.

`baseline_confounded` stars are excluded from the primary completeness number
because a source already existed within the matching radius before injection.
They remain in `stars.csv` for separate crowded-field analysis.

The primary photometric bias uses the flux measured directly from the injected
observation. `differential_mag_error`, which subtracts the original-frame
measurement, is retained only as a diagnostic because that subtraction also
cancels real background and confusion residuals. Non-isolated stars are not
used to build the empirical PSF unless the configuration explicitly permits
the fallback.

Artificial-star flux follows the same convention as APEX Step 7 and Step 10:

```text
mag_inst = -2.5 log10(total source electrons)
mag_cal  = mag_inst + frame zeropoint
```

The runner does not multiply source flux by exposure time. By default it uses
the matching row from `cmd_zeropoint/frame_zeropoint.csv` and assumes a zero
color term for injected stars. A missing calibrated zeropoint is an error
unless the engineering-only initial-ZP fallback is explicitly enabled.

For a precision run, set `magnitude_grid` and
`stars_per_magnitude_per_trial`. Each trial then receives the same number of
stars at every magnitude. A binomial logistic curve estimates the 50%
completeness magnitude, and a trial-level bootstrap provides its 95%
confidence interval. Set `save_injected_fits = false` to remove large
intermediate FITS files after their result tables have been saved.

This test estimates completeness and measurement bias for the selected image
and configuration. Absolute calibration, color terms, and cross-instrument
accuracy require standard-star and independent-pipeline validation.

## CMD observing-condition batch

For a completed CMD project, the batch runner joins the Step 7 photometry index
to the Step 10 frame zeropoints and selects one unique frame nearest the best,
median, and worst FWHM quantiles in each filter.

Inspect the selection without running injections:

```powershell
python benchmark/run_cmd_batch.py `
  --project-root "E:\observed_Analysis\NGC457\pp" `
  --output "benchmark/runs/ngc457_cmd_selection" `
  --select-only
```

Remove `--select-only` to run the full batch. Each frame receives a coarse pass
over a common injected-electron range, followed by a precision pass centered
on that frame's fitted 50% completeness magnitude. Existing non-empty frame
outputs are rejected rather than silently reused.

## CMD validation package

After a CMD batch has finished, build a paper-facing validation package:

```powershell
python benchmark/run_cmd_validation.py `
  --batch-root "benchmark\runs\ngc457_cmd_batch_v2" `
  --project-root "E:\observed_Analysis\NGC457\pp" `
  --output "benchmark\runs\ngc457_cmd_validation_v1"
```

The report reuses existing outputs and does not inject new stars. It writes
tables and plots for completeness, photometric bias/scatter, placement
crowding, false positives, repeated-frame photometric repeatability, and
Step 10 zeropoint residuals.

Combine completed validation reports:

```powershell
python benchmark/run_cmd_validation_combined.py `
  --input "NGC457:seeing variation=benchmark\runs\ngc457_cmd_validation_v2" `
  --input "M5:crowded-field stress=benchmark\runs\m5_cmd_validation_v1" `
  --output "benchmark\runs\cmd_validation_combined_v1"
```

The combined report is the paper-facing aggregation layer. It does not rerun
photometry or injection; it only compares completed validation packages.

## IRAF/PyRAF cross-check

For an independent reference sanity check, compare a completed Step 7 frame
against IRAF/DAOPHOT through PyRAF:

```powershell
python benchmark/run_iraf_crosscheck.py `
  --input "E:\observed_Analysis\NGC457\pp\pp_-0016-gfilter_20240907.fit" `
  --step7 "E:\observed_Analysis\NGC457\pp\result\step7_forced_phot\photometry_pp_-0016-gfilter_20240907.fit.tsv" `
  --output "benchmark\runs\ngc457_iraf_crosscheck_g0016_v1" `
  --mode both
```

`phot_fixed_coords` runs IRAF `phot` at the same Step 7 positions and reports
median-removed magnitude scatter. This is the cleanest implementation
cross-check because APEX aperture correction and electron-flux conventions can
introduce a constant offset relative to IRAF.

The fixed-coordinate and matched-detector outputs also include
`iraf_mag_cal_apcorr_zp`, an IRAF calibrated-equivalent magnitude. It converts
IRAF `phot` output from `zmag - 2.5 log10(flux_ADU)` onto the APEX Step 7/10
scale by applying the matching gain, aperture correction, and frame zero point:

```text
iraf_mag_inst_apcorr_equiv =
    iraf_mag - zmag - 2.5 log10(gain_e_per_adu) - 2.5 log10(apcorr)

iraf_mag_cal_apcorr_zp = iraf_mag_inst_apcorr_equiv + zp_frame
```

This is a fair pipeline-scale comparison, not a claim that IRAF lacks aperture
correction or zero-point tools.

`daofind_phot` starts at high `findpars.threshold` values and moves downward.
If a threshold exceeds the configured source-count guard, lower thresholds are
skipped to prevent runaway detections on noisy or crowded frames.

Run the IRAF detector comparison over every calibrated CMD frame:

```powershell
python benchmark/run_iraf_crosscheck_batch.py `
  --project-root "E:\observed_Analysis\M67\pp" `
  --output "benchmark\runs\m67_iraf_daofind_all_v1" `
  --threshold-grid 12,9,7,5
```

The batch runner writes one per-frame IRAF cross-check package under
`frames/<filter>/<frame>/`, then aggregates threshold-dependent source counts,
APEX match fractions, ZP-aligned residual scatter, recommended thresholds, and
summary plots at the batch root. Repeatability tables compare APEX calibrated
magnitudes, raw IRAF magnitudes, IRAF apcorr+ZP-equivalent magnitudes, and IRAF
median frame-aligned magnitudes separately so raw photometry and full-pipeline
effects are not mixed.
