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
5. Only after Step 4 is stable, migrate Step 5/Forced Aperture, WCS, and
   downstream steps.

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
