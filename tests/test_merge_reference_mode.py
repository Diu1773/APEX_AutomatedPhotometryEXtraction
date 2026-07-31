"""Reference-mode merged workspaces (P7).

A merged workspace normally copies every per-frame photometry TSV; reference
mode records where each frame came from instead. The two modes must be
indistinguishable to every consumer that goes through the loader.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from apex.analysis.merge import reference_store
from apex.analysis.merge.id_match import reconcile_workspace_catalogs
from apex.analysis.merge.workspace_build import materialize_merged_workspace
from apex.utils.photometry_loader import load_frame_photometry
from apex.utils.step_paths_lc import step7_forced_phot_dir, step8_selection_dir


def _make_input_workspace(root: Path, name: str, source_ids, mags, night_id=1):
    """A minimal but real RESULT_* layout the merger accepts."""
    rd = root / name
    s1 = rd / "step1_file_selection"
    s7 = rd / "step7_forced_phot"
    s8 = rd / "step8_selection"
    for d in (s1, s7, s8):
        d.mkdir(parents=True, exist_ok=True)

    fname = f"{name}_frame.fit"
    pd.DataFrame([
        {"det_uid": i, "source_id": sid, "ID": i + 1, "FILTER": "V",
         "x": 100.0 + i, "y": 200.0 + i, "mag": mag, "mag_err": 0.01,
         "file": fname}
        for i, (sid, mag) in enumerate(zip(source_ids, mags))
    ]).to_csv(s7 / f"photometry_{fname}.tsv", sep="\t", index=False)
    pd.DataFrame([{"file": fname, "filter": "V", "night_id": night_id}]).to_csv(
        s7 / "photometry_index.csv", index=False)

    pd.DataFrame([
        {"ID": i + 1, "source_id": sid, "ra_deg": 250.0 + i * 0.01, "dec_deg": 36.0}
        for i, sid in enumerate(source_ids)
    ]).to_csv(s8 / "master_catalog_V.tsv", sep="\t", index=False)
    (s8 / "selection_V.json").write_text(json.dumps({
        "filter": "V", "target_source_id": int(source_ids[0]),
        "comparison_source_ids": [int(s) for s in source_ids[1:]],
        "check_source_id": None,
    }), encoding="utf-8")

    (s1 / "night_assignments.json").write_text(
        json.dumps({"night_assignments": {fname: night_id}}), encoding="utf-8")
    (rd / "run_manifest.json").write_text(json.dumps({
        "run_type": "result", "label": "TGT",
        "date_start": f"2026010{night_id}", "date_end": f"2026010{night_id}",
    }), encoding="utf-8")
    return rd, fname


def _build(tmp_path, mode):
    a, fa = _make_input_workspace(tmp_path, "RESULT_a", [1001, 1002, 1003],
                                  [15.0, 16.0, 17.0], night_id=1)
    b, fb = _make_input_workspace(tmp_path, "RESULT_b", [1001, 1002, 1004],
                                  [15.1, 16.1, 18.0], night_id=2)
    folders = [a, b]
    tags = {str(a): "F01_a", str(b): "F02_b"}
    cats = {
        str(f): {"V": pd.read_csv(f / "step8_selection" / "master_catalog_V.tsv", sep="\t")}
        for f in folders
    }
    rec = reconcile_workspace_catalogs(folders, cats, tags, 2.0)
    out = tmp_path / f"MERGED_{mode}"
    info = materialize_merged_workspace(
        out_dir=out, folders=folders, folder_tags=tags,
        local_id_maps=rec["local_id_maps"],
        merged_catalogs=rec["canonical_by_filter"],
        selection_target_by_filter={"V": 1001},
        selection_comp_by_filter={"V": {1002}},
        selection_check_by_filter={"V": None},
        match_records=rec["match_records"],
        storage_mode=mode,
    )
    return out, info, [f"F01_a__{fa}", f"F02_b__{fb}"]


def test_reference_mode_writes_no_per_frame_tsv(tmp_path):
    full_dir, _, _ = _build(tmp_path, "full")
    ref_dir, _, _ = _build(tmp_path, "reference")

    full_tsvs = list(step7_forced_phot_dir(full_dir).glob("photometry_*.tsv"))
    ref_tsvs = list(step7_forced_phot_dir(ref_dir).glob("photometry_*.tsv"))
    assert len(full_tsvs) == 2
    assert ref_tsvs == []                                  # nothing duplicated
    assert (ref_dir / reference_store.REFERENCE_INDEX).exists()
    assert not (full_dir / reference_store.REFERENCE_INDEX).exists()

    assert reference_store.is_reference_workspace(ref_dir)
    assert not reference_store.is_reference_workspace(full_dir)


def test_both_modes_load_identical_photometry(tmp_path):
    """The point of the mode: consumers must not be able to tell."""
    full_dir, _, frames = _build(tmp_path, "full")
    ref_dir, _, _ = _build(tmp_path, "reference")

    for fname in frames:
        a = load_frame_photometry(full_dir, fname, "V")
        b = load_frame_photometry(ref_dir, fname, "V")
        assert a is not None and b is not None, fname
        cols = ["ID", "source_id", "mag", "mag_err", "file",
                "source_folder", "original_file"]
        left = a[cols].sort_values("source_id").reset_index(drop=True)
        right = b[cols].sort_values("source_id").reset_index(drop=True)
        pd.testing.assert_frame_equal(left, right, check_dtype=False)


def test_reference_mode_keeps_the_merged_identity(tmp_path):
    ref_dir, _, frames = _build(tmp_path, "reference")
    df = load_frame_photometry(ref_dir, frames[1], "V")     # second night
    assert set(df["file"]) == {frames[1]}
    assert set(df["source_folder"]) == {"F02_b"}
    # merged source_ids, and every row carries a merged display ID
    assert df["source_id"].notna().all()
    assert df["ID"].notna().all()


def test_reference_mode_still_writes_a_complete_frame_index(tmp_path):
    """Load-bearing: tools that fall back to globbing photometry_*.tsv (the QA
    report does, in three places) only do so when photometry_index.csv is
    missing or unreadable. Reference mode must therefore keep it complete."""
    full_dir, _, frames = _build(tmp_path, "full")
    ref_dir, _, _ = _build(tmp_path, "reference")

    full_idx = pd.read_csv(step7_forced_phot_dir(full_dir) / "photometry_index.csv")
    ref_idx = pd.read_csv(step7_forced_phot_dir(ref_dir) / "photometry_index.csv")
    assert sorted(ref_idx["file"]) == sorted(full_idx["file"]) == sorted(frames)
    assert set(ref_idx["filter"]) == {"V"}
    assert sorted(ref_idx["night_id"]) == sorted(full_idx["night_id"])

    # The catalogs the merged workspace serves are written either way.
    assert (step8_selection_dir(ref_dir) / "master_catalog_V.tsv").exists()
    assert (step8_selection_dir(ref_dir) / "selection_V.json").exists()


def test_manifests_record_the_mode(tmp_path):
    ref_dir, info, _ = _build(tmp_path, "reference")
    assert info["storage_mode"] == "reference"
    merge_manifest = json.loads((ref_dir / "merge_manifest.json").read_text(encoding="utf-8"))
    run_manifest = json.loads((ref_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert merge_manifest["storage_mode"] == "reference"
    assert run_manifest["storage_mode"] == "reference"

    full_dir, info, _ = _build(tmp_path, "full")
    assert info["storage_mode"] == "full"
    assert json.loads(
        (full_dir / "run_manifest.json").read_text(encoding="utf-8"))["storage_mode"] == "full"


def test_unknown_storage_mode_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="storage_mode"):
        _build(tmp_path, "sometimes")


def test_reference_workspace_breaks_loudly_if_the_input_moves(tmp_path):
    """The documented trade-off: no silent wrong answer, just no data."""
    ref_dir, _, frames = _build(tmp_path, "reference")
    assert load_frame_photometry(ref_dir, frames[0], "V") is not None
    (tmp_path / "RESULT_a").rename(tmp_path / "RESULT_a_moved")
    reference_store._INDEX_CACHE.clear()
    assert load_frame_photometry(ref_dir, frames[0], "V") is None
