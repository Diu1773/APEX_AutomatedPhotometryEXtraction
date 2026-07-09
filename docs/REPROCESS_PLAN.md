# Full-APEX headless reprocessing — build & batch plan

Goal: reprocess every cluster / variable-star target **from RAW, entirely with
APEX** (Step 0 detector calibration → photometry → CMD/LC science), headless, so
the paper's figures rest on APEX raw→science (no AIPPI dependency). Isochrone
fitting (CMD step 12) is left for the user to do interactively.

## HARD SAFETY RAILS (never violate)
1. **Non-destructive.** Never modify or delete `E:\observe_raw_Analysis`,
   `E:\observe_DSY`, or `E:\observed_Analysis`. All output goes to
   `E:\APEX_validation\reprocess\<target>\`. (`observed_Analysis` is backed up to
   `D:\APEX_backup\observed_Analysis`.)
2. **Disk: free space by deleting the backed-up old output (one clean delete),
   then keep frames.** Sequence: (a) `observed_Analysis` (248 GB) finishes copying
   to `D:\APEX_backup\observed_Analysis`; (b) **STRICT verify** — file count AND
   total bytes on D: match E: exactly (E: = 29 659 files / 248.39 GB) + spot-check
   a few files open; (c) only then delete `E:\observed_Analysis`, freeing ~248 GB
   (E: → ~362 GB free). After that the reprocess can keep all calibrated frames —
   no per-target streaming-delete needed. **The 248 GB delete is gated on (b); if
   counts/bytes don't match, DO NOT delete — report.** Still never let E: free
   space drop < 20 GB. Deleting is only ever `observed_Analysis` (the old,
   backed-up, regenerable output) — never raw.
3. **Verify before batch.** Build each headless runner, then verify it on ONE
   target end-to-end and confirm results match the existing analysis (they should,
   since APEX Step 0 ≡ AIPPI to sub-DN — a mismatch means a bug). Only batch after
   the single-target chain is green.
4. **Bias/dark** come from `E:\bias` (+ each target folder); scan picks them up.

## FOUNDATION (already working — do not rebuild)
- Step 0 (detector calibration): headless via `_ScanWorker` + `_CalibrationWorker._run()`
  (see `scratchpad/reprocess_cluster.py`). NGC6811 done (8 min / 2.94 GB / 38 lights).
- Steps 1-7 (scan→crop→sky→detect→wcs→refbuild→forcedphot): `apex run --mode cmd
  --steps 1-7 --config <t>.toml --force` (PipelineRunner). Verified resolvable.
- CMD 8 (PSF phot): `scripts/run_step8_headless.py`. CMD 10 (zeropoint):
  `scripts/run_step10_headless.py`. Both drive the GUI QThread workers headless.

## BUILD (missing pieces — the actual work)
Each GUI step's real work lives in a QThread worker; drive it headless with a
`QCoreApplication` exactly as run_step8/10_headless.py and Step 0 do.
1. **CMD 9 headless runner** — Master ID assignment. Use the *auto* path
   (`_auto_add_detections_to_master` / `save_master_ids` in
   `apex/gui/workflow/cmd/step9_master_id_editor.py`); no interactive editing.
2. **CMD 11 headless runner** — CMD plot/product from
   `apex/gui/workflow/cmd/step11_cmd_plot.py` (worker → CMD table + figure).
3. **LC 8-11 headless runners** — target selection, lightcurve builder, detrend/
   merge, period analysis (`apex/gui/workflow/lc/step8_target_selection.py`,
   `step9_lightcurve_builder.py`, `step10_detrend_merge.py`, `step11_period_analysis.py`).
4. **Per-target Step 0 → photometry reorganisation.** Step 0 writes
   `calibrated/<night>/pp_*.fit` mixing targets (e.g. 20260611 has NGC6811 + NGC3231).
   Photometry is per-target: split calibrated frames into per-object folders
   (by filename stem or OBJECT header) before pointing `data_dir` at them.
5. **Per-target config generation.** Copy `parameters.toml`, set `[io].data_dir` →
   the target's APEX Step-0 output, `[io].result_dir` → `reprocess/<target>/result`,
   `[target].ra_deg/dec_deg` from the object (existing: parameters_M13/M5.toml).
6. **Batch orchestrator** `scripts/reprocess_batch.py` — for each target:
   Step 0 (headless) → split per object → gen config → `apex run 1-7` → CMD 8/9/10/11
   (or LC 8-11) → keep results **and** calibrated frames (space freed by the
   observed_Analysis delete, rail #2) → log. Resumable (skip targets already done).
   Disk guard: stop if E: free < 20 GB.

## TARGETS
- **CMD clusters (RAW in observe_raw_Analysis):** NGC6811(20260611), M13(20260515),
  M67(20260208), M3. (M5 raw is NOT present → Fig 9 M5 cannot be redone; note it.)
- **LC variable stars:** AE UMa, YZ Boo (multi-night).
- **observe_DSY:** mostly nebulae/galaxies (M27, M51, M81, M57, …) — NOT star
  clusters → Step 0 + photometry catalog only, NO CMD/isochrone. Only run CMD on
  the actual clusters among them (e.g. NGC457, M16 region). Preprocess the rest to
  catalogs if time/disk allow; else skip (paper doesn't need them).
- Complete calib set required (light + dark + flat; bias from E:\bias). Skip
  incomplete targets with a logged reason.

## VERIFY / COMPARE
After a target's chain, compare the new APEX-raw→science products to the existing
`observed_Analysis` result (photometry MAD, CMD ridge, LC period). Expect a match
(validates the reprocess). Log any deviation for review — do NOT silently accept.

## BATCH ORDER (fast → slow, verify early)
NGC6811 (proof, config exists) → M13 → M67 → M3 → AE UMa (LC) → YZ Boo (LC) → DSY.

## PROGRESS LOG
Append per-target status to `E:\APEX_validation\reprocess\PROGRESS.md`
(target, steps done, disk after, result vs existing, issues).
