# Cache Manager Foundation

This document defines the first shared cache layer for APEX. Existing step
payloads such as `detect_*.json`, `photometry_*.csv`, and WCS outputs are not
renamed in this phase. The manager adds a standard manifest next to future or
migrated payloads so cache validity can be checked consistently.

## Manifest Fields

Each manifest is JSON and contains:

- `manifest_version`: version of the shared manifest envelope.
- `cache_schema_version`: step-specific payload schema version.
- `step_id`: stable step identifier such as `step4_detection`.
- `created_at`: UTC ISO timestamp.
- `parameter_hash`: optional parameter hash used for invalidation.
- `input_files`: source file signatures with path, size, mtime, and optional
  sha1.
- `payload_paths`: named output payload files.
- `dependency_versions`: optional package/tool version metadata.
- `extra`: step-specific metadata that should not be parsed by the generic
  manager.

## Validation Reasons

The first shared manager reports explicit invalidation reasons:

- `missing_manifest`
- `manifest_version_mismatch`
- `step_id_mismatch`
- `cache_schema_mismatch`
- `parameter_hash_mismatch`
- `input_count_mismatch`
- `input_signature_mismatch`
- `missing_payload:<name>`
- `payload_not_found:<name>`
- `ok`

## Migration Strategy

1. Keep all existing payload filenames readable.
2. Add manifests for one step at a time, starting with Step 4 Detection.
3. Keep legacy cache readers active while writing the new manifest.
4. Use explicit invalidation reasons in logs instead of silent cache misses.
5. Only after Step 4 is stable, migrate the current shared chain in order:
   Step 5 WCS, Step 6 Reference Build, Step 7 Forced Aperture Photometry, and
   downstream steps.

## Current Cache Controls

The shared cache layer is still transitional. These are the user-facing cache
controls that should be treated as active behavior today:

| Step | UI control | Runtime parameter | Meaning | Current payloads |
| --- | --- | --- | --- | --- |
| Step 4 Source Detection | `Use detection cache` | `resume_mode`, inverse of `force_redetect` | Reuse compatible `detect_*.json` payloads and skip completed frames. Disabled means detect every selected frame again. | `cache/detect_*.json`, `cache/detect_*.csv`, mirrored Step 4 outputs |
| Step 5 WCS / astrometry.net | `Use Cached Outputs` | `astnet_local_use_cache` | Reuse compatible local solve-field sidecars instead of running solve-field again. | local astrometry.net sidecars and Step 5 WCS outputs |
| Step 5 Gaia fallback | `Gaia cache miss` | `gaia_allow_no_cache` | Not a cache reuse toggle. This only decides whether online Gaia query is allowed when the local Gaia cache is missing. | Gaia query cache |
| Step 6 Reference Build | `Use existing output if complete` | inverse of `force_master_build` | Reuse the existing reference-build summary when complete unless rebuilding is forced. | `step6_refbuild/ref_build_meta.json`, master/reference catalogs |
| Step 7 Forced Aperture Photometry | `Use existing output if complete` | inverse of `force_rephot` | Reuse the complete Step 7 output set when every selected frame has an OK index row and TSV payload. Incomplete output triggers a full recompute. | `step7_forced_phot/photometry_index.csv`, `photometry_*.tsv`, `frame_stats.csv` |
| CMD Step 8 PSF Photometry | `Use existing output if complete` | `project_state.psf_photometry.use_existing_psf_output` | Reuse Step 8 only when the saved signature matches selected frames, Step 4/7 input mtimes, crop mode, and PSF parameters, and all expected per-frame outputs exist. Disabled means clear stale Step 8 products and recompute. | `cmd_psf/psf_output_signature.json`, `photometry_index.csv`, `photometry_*.tsv`, `epsf_model_*.fits`, residual metadata/FITS/NPY products |

The old `parallel.step9_use_cache` parameter was removed from the active CMD
parameter map because Step 9 does not implement step-level output reuse. If it
appears in an older TOML file, it should be treated as an ignored legacy key
until a real Step 9 reuse path is designed.

## Shared Cache/Parameter UI Helpers

Workflow parameter buttons and the migrated cache/reuse controls should use
`apex/gui/workflow/ui_helpers.py` rather than local styles. The active helpers
are:

- `create_parameter_button()` for workflow `... Parameters` buttons.
- `configure_parameter_dialog()` for shared parameter-dialog title, size, and
  visual styling.
- `create_collapsible_section()` for scrollable parameter-dialog sections.
- `create_detection_cache_checkbox()` for Step 4 detection-cache reuse.
- `create_output_reuse_checkbox()` for complete step-output reuse controls.
- `create_cache_checkbox()` for other cache-like options that need a shared
  style but do not fit the two standard labels.

Current migrated section dialogs include Step 4 detection parameters, Step 5
WCS/Astrometry.net parameters, and CMD Step 8 PSF parameters. Smaller dialogs
can stay flat, but their button/dialog styling should still use the shared
helpers.

## Other Cache-Like Paths

These paths use the word cache or reuse previous artifacts, but they are not the
same as a step-level recompute toggle:

| Area | Type | Current behavior | UI toggle? | Policy |
| --- | --- | --- | --- | --- |
| Step 1 File Selection | project setup | Creates `result/cache` and records file-selection state. | No | Keep as project/session state, not step cache. |
| Step 2 Crop | persistent artifact reuse | Reuses an existing cropped FITS only when source mtime, crop rectangle, and output dimensions still match. | No separate toggle | Treat as output reuse. If crop rectangle changes after Step 4, rerun Step 4 and later steps. |
| Step 3 Sky Preview | in-memory preview cache | LRU FITS data/header cache for fast frame switching. | No | Never persist or expose as recompute cache. Clear when frame list/project changes. |
| CMD Step 8 PSF Photometry in-run ePSF cache | in-run model cache | May reuse an in-memory shared ePSF during one run when `psf_shared_filter_epsf` is enabled. Step-level output reuse is handled separately by the Step 8 signature check above. | No separate toggle | Keep this as an internal speed path; do not treat it as persisted cache validity. |
| CMD Step 9 Master ID Editor | in-memory UI cache | LRU FITS and photometry-table caches for interactive overlay work. Output validity is `master_star_ids.csv`. | No | Keep as UI performance cache, not a recompute cache. Add a new output-reuse toggle only after implementing signature checks around `master_star_ids.csv` and ROI outputs. |
| CMD Step 10 Zeropoint | in-memory plot cache | Stores scatter-plot pick data for current view. | No | Keep internal; no user cache toggle needed. |
| CMD Step 11 CMD Plot | output restore | Restores existing plot/output state when files exist. | No | No cache toggle unless a slow build step is added. |
| CMD Step 12 Isochrone | persistent external-data cache + in-memory data cache | Merges an isochrone folder into `.apex_cache/combined_isochrones.dat` with a file signature, then caches loaded CMD/isochrone arrays in memory. | No | This is an external data parsing cache; keep separate from workflow step cache. Add clear/rebuild only in the isochrone UI if users need it. |
| LC Step 8 Target Selection | in-memory UI/network cache | FITS/display/source caches and SIMBAD type cache for interactive target selection. | No | Keep as UI performance cache. Network cache should have its own refresh action if needed. |
| LC Step 9 Lightcurve Builder | in-memory series cache + persisted comp QC summary | Caches photometry/header/diff-series data in memory. `comp_qc_summary.csv` can be loaded when its JSON signature matches target, comp IDs, date filter, exclusions, thresholds, and photometry-index mtime. | No global toggle | Keep the QC summary cache automatic because it is signature guarded. Add an explicit `Recompute QC` action if users report stale QC. |
| LC Step 10 Detrend / Step 11 Period | output restore | Reads saved products and GUI state. | No | Treat as step outputs, not cache toggles, unless slow recomputation is introduced. |

## Cache Policy

- Default to cache reuse only when the payload is complete and signature- or
  schema-compatible.
- Rebuild/off toggles must recompute the complete step output and avoid mixing
  stale rows with new rows.
- Partial reuse is allowed only after the step has cache-aware index merging and
  stale-output cleanup. Step 7 therefore reuses all existing output or recomputes
  all selected frames; it does not currently merge partial forced-photometry
  cache.
- Options named `allow_no_cache` mean "continue by querying/computing without a
  local cache" and should not be labeled like cache reuse toggles.
- New or migrated caches should use `StepCacheManager` for manifest validation.
  Existing step payloads can remain in place while manifests are added beside
  them.

## First Integration Target

Step 4 currently stores source path, size, mtime, cache schema, and detection
engine in `detect_*.json`. That should become:

- payload: existing `detect_*.json`, `detect_*.csv`, optional
  `detect_peak_*.csv`
- manifest: `_manifests/<frame>.manifest.json`
- `step_id`: `step4_detection`
- `cache_schema_version`: current Step 4 schema
- `parameter_hash`: existing parameter hash
- `extra.detect_engine`: selected detection engine

The old fields can remain in the payload during transition for backward
compatibility.
