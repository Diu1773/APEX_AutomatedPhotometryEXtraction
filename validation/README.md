# APEX Validation

This directory is for science validation that is larger than unit tests.

Use two layers:

1. Synthetic truth checks: deterministic, small, and fast enough to run after
   Step12 fitter changes.
2. Real-data checks: documented runs on known clusters, compared against
   literature ranges and inspected plots.

## Step12 Synthetic Truth Check

Run:

```bash
python3 validation/cmd_step12_synthetic.py
```

Optional JSON report:

```bash
python3 validation/cmd_step12_synthetic.py --output validation/reports/cmd_step12_synthetic.json
```

What it validates:

- `IsochroneFitterV2` can recover a known synthetic CMD input.
- The recovered `log_age`, `[M/H]`, distance modulus, and `E(color)` stay within
  explicit tolerances.
- Grid scan metadata and best-fit result are serializable for release reports.

This is not a real astrophysical calibration. It is a regression guard for the
fitter mechanics and Step12 science plumbing.

## Real-Data Validation Template

For a known cluster such as M13 or M38:

1. Run the normal CMD workflow through Step10 zeropoint calibration.
2. Open Step12 and select the intended color/magnitude axes.
3. Run color-color `E(B-V)` fit if at least 3 calibrated bands exist.
4. Run Auto Fit or Grid Scan.
5. Export Step12 results.
6. Save the final CMD screenshot with the isochrone overlay.
7. Compare against a documented literature range for:
   - `log_age` or age in Gyr
   - `[M/H]`
   - distance modulus
   - `E(B-V)` or the selected `E(color)`
8. Record pass/fail plus any manual adjustment needed.

Recommended acceptance style:

- Synthetic truth: strict numeric tolerances, automated.
- Known-cluster real data: range-based acceptance plus plot inspection.
- GUI behavior: manual click-through note, separate from compile/test status.

## Step12 Real-Data Runner

Run against an APEX result directory:

```bash
python validation/cmd_step12_realdata.py E:/observed_Analysis/M5/light/result \
  --color g-r --mag g \
  --age-bounds 10.00 10.12 \
  --mh-bounds -1.70 -1.00 \
  --ecolor-bounds 0.00 0.10 \
  --output validation/reports/M5_cmd_step12_realdata.json \
  --output-plot validation/reports/M5_cmd_step12_realdata.png
```

The runner reads `cmd_zeropoint/median_by_ID_filter_wide_cmd.csv`, loads a
compact PARSEC subset, runs the Step12 grid scan, and writes:

- JSON report with selected bands, fit parameters, and broad sanity checks.
- CMD overlay plot for manual inspection.
- Reusable isochrone subset cache under `validation/reports/iso_cache/`.

Use the JSON pass/fail as a first gate only. The plot still has to be inspected,
especially for sparse turnoff coverage, horizontal-branch dominated CMDs, or
strong field contamination.

## Multi-Exposure Photometry Validation

Real-data check of the multi-exposure CMD support (count-rate magnitudes,
union master catalog, dynamic-range bridging). There is no standalone runner —
the diagnostics are written by Steps 6–10 themselves; the validation reads them
back.

- **Report:** [`reports/NGC6811_multiexposure_validation.md`](reports/NGC6811_multiexposure_validation.md)
  (dataset, evidence tables, verdict) + `reports/NGC6811_multiexposure_nonlinearity.png`.
- **What it confirms:** `zp_frame` stays constant across a 16× exposure range
  (count-rate normalization), bridge-star Δmag is flat (linear detector,
  exposures agree < 7 mmag), the union master spans `n_det_frames` 1→N, and
  short/long exposures cover the bright/faint ends respectively.
- **Pipeline artifacts to inspect:** `cmd_zeropoint/frame_zeropoint.csv`,
  `cmd_zeropoint/nonlinearity_summary.csv` (+ `nonlinearity_check.png`),
  `step6_refbuild/ref_catalog.tsv` (`n_det_frames`).
- **Unit coverage:** `tests/test_step6_union_master.py`,
  `tests/test_step10_nonlinearity.py`.
