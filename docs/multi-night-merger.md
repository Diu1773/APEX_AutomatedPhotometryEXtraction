# Multi-Night Light Curve Merger

The **Multi-Night Light Curve Merger** combines several independently-processed
single-night workspaces (`RESULT_*`) into one merged workspace (`MERGED_*`) whose
light curves span all of the input nights. It is a standalone tool in **LC mode**,
reachable from **Tools → Multi-Night Light Curve Merger** (`Ctrl+M`).

> **Not the same as Step 10 "Detrend & Night Merge".** Step 10 merges the nights
> that already live *inside one* workspace during detrending (`lc_detrend/`). This
> tool merges *separate* `RESULT_*` workspaces that were each reduced on their own
> — e.g. you processed 2025-04-29 and 2025-04-30 as two runs and now want a single
> combined light curve. Use Step 10 for one multi-night run; use this tool to join
> runs.

## Prerequisites

Each input folder must be a workspace that has been processed at least through the
**Light Curve Builder (LC Step 9)**, i.e. it must contain:

| Requirement | Canonical location | Legacy location (auto-detected) |
|---|---|---|
| Forced photometry + index | `step7_forced_phot/photometry_index.csv` | `step5_photometry/` |
| Target/comparison selection | `lc_selection/master_catalog_*.tsv` + `selection_*.json` | `step9_selection/` |
| Differential light curves | `lc_lightcurve/lightcurve_*.csv` | `step10_lightcurve/` |

The Step 1 scan reports each folder as **merge-ready** only when all three are
present. All inputs must share the same target/label, and they are merged
per-filter (filters that appear in only some nights are unioned).

> **Legacy workspaces are supported.** Folders produced *before* the step
> renumbering — which use `step5_photometry` / `step9_selection` /
> `step10_lightcurve` and store per-frame identity as `det_uid` plus a separate
> `cache/idmatch/idmatch_<frame>.csv` map instead of an embedded `source_id` — are
> detected and read automatically. You do **not** need to re-run the pipeline to
> merge old runs.

## Workflow

The tool is a six-page wizard. The first two pages are the merge itself; the
remaining four embed the live LC step windows so you can continue analysis on the
merged workspace without leaving the tool.

1. **Step 1 · 폴더 선택 (Select folders).** Add the `RESULT_*` (or earlier
   `MERGED_*`) folders to combine, then **폴더 스캔** to validate them. Confirm the
   output path — the default is `MERGED_<target>_<start>_<end>/` next to the
   inputs.
2. **Step 2 · ID 매칭 (ID match).** Pick a **position match radius** (arcsec) and
   run **ID 매칭 실행**. Sources are reconciled across nights with this priority:
   **Gaia `source_id` → existing canonical `source_id` → positional match within
   the radius.** Unmatched stars become new canonical entries. The merge then
   materializes the combined workspace; this runs on a **background thread**, so
   the window stays responsive (progress is logged in the panel). A status line
   reports `exact / positional / new` counts per night and filter.
3. **Steps 3–6 · 선택 / Light Curve / Detrend / Period.** These are the embedded LC
   Step 8–11 windows operating on the merged workspace (target/comparison
   selection, light-curve builder, detrend & ensemble correction, period
   analysis). The base night's selection is carried over as the default.

## Output

A new `MERGED_<target>_<start>_<end>/` workspace in the **current** step layout:

```text
MERGED_.../
  step1_file_selection/   merged headers, night assignments, file path map
  step7_forced_phot/      per-frame merged photometry (source_id + unified ID),
                          photometry_index.csv, master_sources.csv, filter_frames.json
  lc_selection/           merged master_catalog_<filter>.tsv + selection_<filter>.json
  lc_lightcurve/  lc_detrend/  lc_period/   (populated as you run Steps 9–11)
  merge_manifest.json     inputs, filters, folder tags
  merge_id_map.csv        per-source local→merged ID/source_id mapping
```

Night numbering is reassigned so every input night is distinct in the merged
workspace. The merged workspace is a normal LC workspace — you can also open it
directly in LC mode afterwards.

## Tips & limitations

- **Same target only.** The tool refuses to merge runs whose labels differ.
- **Frames without a recoverable identity are skipped.** For legacy inputs, a
  frame is only merged if its `cache/idmatch/idmatch_<frame>.csv` exists; frames
  that failed WCS/idmatch in the original run are dropped (reported in the log).
- **Match radius.** Start at the default (2.0″). Too large risks blending close
  neighbours into one canonical star; too small splits a star that drifted
  between nights into duplicates.
