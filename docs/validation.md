# Validation

APEX ships a reproducible validation harness. Run it with:

```bash
apex validate --suite all --output validation/report
```

This generates `validation/report/validation_summary.md` (human-readable report) and `validation/report/validation_manifest.json` (machine-readable metrics with the git commit and package versions).

## What it measures

- **Artificial-star completeness & photometry** — injects an empirical PSF into a reference frame (synthetic when none is supplied), reruns the production detector, and reports the 50% completeness magnitude, photometric bias/scatter, and position RMSE.
- **IRAF cross-check** (optional) — compares APEX aperture magnitudes against IRAF `phot` on the same fixed coordinates.
- **Known-target reproduction** — a zero-data synthetic isochrone-recovery regression guard (plus optional real-cluster checks).

## Latest headline numbers

_Run `apex validate` to populate the latest numbers._
