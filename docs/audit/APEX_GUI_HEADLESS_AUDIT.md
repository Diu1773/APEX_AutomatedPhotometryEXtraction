# APEX GUI/headless architecture audit

## Observed architecture

| Area | Actual path | Status |
|---|---|---|
| GUI shell | `apex/gui/main_window.py`, PyQt5 workflow windows | Desktop orchestration, parameter editing, previews, review decisions and worker presentation. |
| Qt-free shared core | `apex/analysis/{calibration,detection,crop,refbuild,forced_photometry,wcs_solve}.py`; `apex/pipeline/*` | Steps 0–7 compute paths are largely callable without importing `apex.gui` or PyQt5. Detection, WCS and forced photometry explicitly document this boundary. |
| Qt adapter | `step4_source_detection.DetectionWorker`, `step5_wcs_plate_solving.*Worker`, `step6_ref_build.RefBuildWorker`, `step7_forced_aperture_phot.ForcedPhotWorker` | Thin QThread wrappers re-emit progress/error callbacks around shared analysis functions. This is a real shared path for these stages. |
| Headless runner | `apex/cli.py:_cmd_run` → `PipelineRunner` | Runs registered shared Steps 1–7; Step 0 is explicitly separate. It does not automatically run CMD 8–12 or LC 8–11. |
| Step 8 PSF | `apex/gui/workflow/cmd/step8_psf_photometry.py:Step6PSFWorker` | Main ePSF/PSF fitting still lives in a GUI workflow module. `scripts/run_step8_headless.py` imports PyQt5/QCoreApplication and runs the worker synchronously, so it is off-screen Qt rather than a Qt-free headless service. |
| Step 10 zeropoint | `apex/gui/workflow/cmd/step10_zeropoint_calibration.py:ZeropointCalibrationWorker`; `scripts/run_step10_headless.py` | Same off-screen Qt pattern. Supporting catalogue/fit services are reusable, but orchestration is not fully GUI-independent. |
| Step 12 CMD | `apex/analysis/cmd/isochrone_fit_service.py`; `scripts/run_step12_headless.py` | Genuine Qt-free service/runner; GUI `McmcFitWorker` delegates to it. This is the strongest parity example. |
| LC services | `apex/analysis/light_curve/*` plus GUI workers in `step8_target_selection.py`, `step9_lightcurve_builder.py`, `step10_detrend_merge.py`, `step11_period_analysis.py` | Scientific primitives are reusable, but full LC workflow/review state remains GUI-centric; no single `apex run --mode lc` equivalent covers all LC steps. |

## Claim audit

| Manuscript/document wording | Source trace | Assessment |
|---|---|---|
| “GUI and headless execute the same shared computation” for Steps 0–7 | `analysis/*` docstrings + GUI worker delegation + `pipeline/registry.py` | Defensible for the shared functions, subject to identical parameter/config and input path. |
| “The whole pipeline runs headlessly without Qt” | `scripts/run_step8_headless.py`, `run_step10_headless.py`, `run_step7_headless.py` | Not defensible as written. Step 8/10/legacy Step 7 scripts import PyQt5/QCoreApplication. Use “scriptable/off-screen worker execution” until refactored. |
| “Step 8 is independently validated by the same code as a user GUI run” | `step8_psf_photometry.py` and `run_step8_headless.py` | Partly defensible: same worker path, but it is not a separate Qt-free implementation and therefore is not independent of the GUI module. |
| “Pipeline is reproducible” | `PipelineRunner` manifest, cache signatures, seeded benchmark/MCMC paths | Only conditional. Manifest lacks full config/environment; network/solver availability and external databases change outputs. |
| “Parallel processing” | Thread pools in detection/WCS/forced photometry/period bootstrap/PSF | Defensible as frame/bootstrap-level threading. Do not imply process parallelism or linear scaling. |

## Refactor needed for a strong architecture claim

1. Move Step 8 fit orchestration (not just policy/diagnostics) into `apex/analysis/psf_service.py` with explicit callbacks and a serialisable request/result contract.
2. Move Step 10 worker orchestration into a Qt-free calibration service; have GUI and CLI adapters call it.
3. Add a first-class LC pipeline registry and headless runner, or narrow the paper claim to shared analysis services rather than end-to-end LC execution.
4. Make each run manifest include resolved parameters, package versions, input signatures, solver selection/fallback, seed, worker counts, and skipped external suites.
5. Add parity tests that run GUI adapters and headless services on the same fixture and compare files/metrics, not merely imports.

## GUI-specific scientific risk

GUI fields can encode user judgement (crop, engine, QC acceptance, reference-star selection, priors). That judgement is scientifically valuable only if the selected values and review outcome are persisted in the project state/manifest. A screenshot or a GUI default is not provenance. The paper should describe “recorded human-in-the-loop choices” only where the corresponding state file and run manifest are available.
