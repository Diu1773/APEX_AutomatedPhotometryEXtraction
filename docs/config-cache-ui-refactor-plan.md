# APEX Config, Cache, Legacy, and Shared UI Refactor Plan

## Goal

Stabilize APEX configuration and workflow infrastructure before larger workflow
changes. The current codebase has long-lived legacy parameter names, duplicated
CMD/LC loaders, step-specific cache layouts, and repeated log/parameter UI
patterns. This refactor should make those foundations explicit, versioned, and
safe to migrate incrementally.

## Non-Goals

- Do not rewrite scientific algorithms as part of the foundation work.
- Do not remove legacy readers until equivalent migration coverage exists.
- Do not convert every step UI in one patch.
- Do not break existing `params.P.<flat_name>` access during the first phase.

## Phase 0: Inventory

Create a complete inventory before changing behavior.

- List every key in `parameters.toml` and `parameters.example.toml`.
- List every `params.P.<name>` attribute used under `apex/`.
- List every key written by `save_toml()`.
- List every cache file pattern, output directory, and JSON/CSV schema by step.
- Mark whether each value belongs in configuration, cache metadata, or
  `project_state`.
- Identify legacy aliases such as `pix` vs `px`, old flat keys, and fallback
  output paths.

Deliverables:

- `docs/parameter-inventory.md`
- `docs/cache-inventory.md`
- A small script or test that reports unmapped `params.P` attributes.

## Phase 1: Canonical Parameter Schema

Use the nested TOML schema as the canonical user-facing format.

- Add a top-level `schema_version`.
- Move shared TOML path mappings into one module, for example
  `apex/config/parameter_map.py`.
- Keep `params.P.<flat_name>` as a compatibility runtime facade for now.
- Read legacy keys as aliases, but save only canonical TOML keys.
- Normalize naming rules:
  - pixel units use `_px`
  - angular units use `_arcsec`
  - ADU units use `_adu`
  - electrons use `_e`
  - seconds use `_s`
  - booleans use consistent prefixes such as `use_`, `force_`, `require_`, or
    `enabled`
- Replace duplicated CMD/LC mapping drift with shared definitions plus
  mode-specific extensions.
- Decide what happens to catch-all sections such as `[parameters]`.

Deliverables:

- Shared parameter map.
- Canonical `parameters.example.toml`.
- Loader tests for current TOML and selected legacy aliases.
- Save/load round-trip test that preserves canonical keys.

## Phase 2: Unified Cache Structure

Standardize cache metadata after parameter hashing and schema versioning are
available.

- Add a cache manager module, for example `apex/core/cache_manager.py` or
  `apex/utils/cache_manager.py`.
- Define standard cache metadata:
  - `cache_schema_version`
  - `step_id`
  - `input_files`
  - input file mtime/hash policy
  - parameter hash
  - dependency versions where relevant
  - `created_at`
  - payload paths
- Keep legacy cache readers active during migration.
- Use explicit invalidation reasons instead of silently accepting stale files.
- Start with Step 4 detection because it already has cache schema handling.
- Then migrate Step 5/Forced Aperture, WCS, MasterBuild, and downstream steps.

Deliverables:

- Cache manifest format.
- Step 4 cache writer/reader using the shared manager.
- Legacy cache compatibility tests.

## Phase 3: Legacy Cleanup

Remove or quarantine old behavior only after the new loaders are proven.

- Convert old TOML aliases to read-only compatibility paths.
- Add warnings for deprecated keys and cache layouts.
- Remove duplicated parameter maps from CMD/LC loaders.
- Move project/session state out of parameters where appropriate.
- Retire obsolete fallback output paths only when downstream steps no longer
  depend on them.

Deliverables:

- Deprecation report in logs.
- Tests for the supported legacy migration cases.
- Reduced duplicated loader code.

## Phase 4: Shared Log and Parameter UI

Build common GUI utilities after canonical parameter metadata exists.

Suggested modules:

- `apex/gui/common/collapsible.py`
- `apex/gui/common/parameter_dialog.py`
- `apex/gui/common/log_panel.py`
- `apex/gui/common/progress_panel.py`

Parameter UI behavior:

- Build sections from parameter metadata.
- Use scrollable, collapsible sections by default.
- Support common field types: int, float, bool, choice, string, path, list.
- Allow custom section hooks for dynamic cases such as per-filter sigma.
- Save through the shared parameter facade.

Log/progress UI behavior:

- Standard timestamped log append.
- Optional worker progress rows.
- Clear/copy/save actions.
- Auto-scroll with line limit.
- Step/tool adapters rather than rewriting every worker at once.

Migration order:

1. Step 4 Detection
2. Step 5 Aperture and Forced Aperture
3. Step 6 WCS
4. Step 7 MasterBuild
5. PSF Photometry
6. Zeropoint/CMD/Isochrone
7. Tools

## Validation Strategy

At minimum for every phase:

- `python -m compileall -q apex main.py`
- TOML load/save smoke test for CMD and LC parameter classes.
- Existing current `parameters.toml` can be loaded.
- A canonical example TOML can be loaded.
- Legacy alias fixtures can be loaded.

For cache changes:

- Run a small Step 4 detection smoke test.
- Verify old cache files are either accepted with a legacy path or rejected with
  a clear reason.

For GUI changes:

- Instantiate target dialogs where possible.
- Manually open the migrated step and verify save/cancel behavior.

## First Patch Scope

Keep the first patch intentionally small:

1. Add parameter and cache inventory docs/scripts.
2. Add shared parameter map skeleton.
3. Add `schema_version` support without changing runtime behavior.
4. Add tests that prove current TOML still loads.

Do not migrate all steps or remove legacy keys in the first patch.
