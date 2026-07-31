"""Regression test for the multi-night merge pipeline.

Locks the behaviour of ``reconcile_workspace_catalogs`` (id_match) and
``materialize_merged_workspace`` (workspace_build) before/after the H3/H4
performance refactors. Covers all three match outcomes — exact source_id,
positional, and brand-new source — across two folders and one filter.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from apex.analysis.merge.id_match import reconcile_workspace_catalogs
from apex.analysis.merge.workspace_build import materialize_merged_workspace
from apex.analysis.merge.workspace_scan import (
    folder_tag,
    load_master_catalogs_by_filter,
)
from apex.utils.io_utils import read_csv_int64_source_id
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_lc import step8_selection_dir


def _write_master_catalog(folder: Path, flt: str, rows: list[dict]) -> None:
    sdir = step8_selection_dir(folder)
    sdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(sdir / f"master_catalog_{flt}.tsv", sep="\t", index=False)


def _write_frame(folder: Path, fname: str, flt: str, night_id: int, rows: list[dict]) -> None:
    fdir = step7_forced_phot_dir(folder)
    fdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(fdir / f"photometry_{fname}.tsv", sep="\t", index=False)
    idx_path = fdir / "photometry_index.csv"
    new_row = pd.DataFrame([{"file": fname, "filter": flt, "night_id": night_id}])
    if idx_path.exists():
        prev = pd.read_csv(idx_path)
        new_row = pd.concat([prev, new_row], ignore_index=True)
    new_row.to_csv(idx_path, index=False)


@pytest.fixture
def merge_workspace(tmp_path):
    """Two folders, one filter, exercising exact/positional/new matches."""
    folder_a = tmp_path / "RESULT_A"
    folder_b = tmp_path / "RESULT_B"

    # Base folder A: three Gaia-matched stars.
    _write_master_catalog(folder_a, "V", [
        {"ID": 1, "source_id": 1001, "ra_deg": 10.0000, "dec_deg": 20.0000, "x_ref": 100.0, "y_ref": 100.0},
        {"ID": 2, "source_id": 1002, "ra_deg": 10.0010, "dec_deg": 20.0010, "x_ref": 110.0, "y_ref": 110.0},
        {"ID": 3, "source_id": 1003, "ra_deg": 10.0020, "dec_deg": 20.0020, "x_ref": 120.0, "y_ref": 120.0},
    ])
    _write_frame(folder_a, "frameA1.fits", "V", 1, [
        {"ID": 1, "source_id": 1001, "mag": 15.0},
        {"ID": 2, "source_id": 1002, "mag": 16.0},
        {"ID": 3, "source_id": 1003, "mag": 17.0},
    ])

    # Folder B:
    #  ID 1 -> exact source_id match (1001)
    #  ID 2 -> no gaia (source_id -5), positional match to A's star 1002
    #  ID 3 -> new star with its own positive source_id 2001 (far away)
    _write_master_catalog(folder_b, "V", [
        {"ID": 1, "source_id": 1001, "ra_deg": 10.00000, "dec_deg": 20.00000, "x_ref": 101.0, "y_ref": 101.0},
        {"ID": 2, "source_id": -5, "ra_deg": 10.00101, "dec_deg": 20.00101, "x_ref": 111.0, "y_ref": 111.0},
        {"ID": 3, "source_id": 2001, "ra_deg": 10.50000, "dec_deg": 20.50000, "x_ref": 500.0, "y_ref": 500.0},
    ])
    _write_frame(folder_b, "frameB1.fits", "V", 1, [
        {"ID": 1, "source_id": 1001, "mag": 15.1},
        {"ID": 2, "source_id": -5, "mag": 16.1},
        {"ID": 3, "source_id": 2001, "mag": 17.1},
    ])

    folders = [folder_a, folder_b]
    folder_tags = {str(f): folder_tag(i, f) for i, f in enumerate(folders)}
    catalogs_by_folder = {str(f): load_master_catalogs_by_filter(f) for f in folders}
    return tmp_path, folders, folder_tags, catalogs_by_folder


def _reconcile(folders, catalogs_by_folder, folder_tags):
    return reconcile_workspace_catalogs(
        folders, catalogs_by_folder, folder_tags, pos_tol_arcsec=2.0,
    )


def test_reconcile_matches_and_ids(merge_workspace):
    _, folders, folder_tags, catalogs_by_folder = merge_workspace
    res = _reconcile(folders, catalogs_by_folder, folder_tags)

    canon = res["canonical_by_filter"]["V"].sort_values("ID").reset_index(drop=True)
    # 3 base stars + 1 new from folder B
    assert canon["ID"].tolist() == [1, 2, 3, 4]
    assert canon["source_id"].astype("int64").tolist() == [1001, 1002, 1003, 2001]

    # Folder B summary: 1 exact, 1 positional, 1 new
    b_summary = [r for r in res["match_summary_rows"]
                 if r["folder"] == folders[1].name and r["filter"] == "V"][0]
    assert (b_summary["exact"], b_summary["pos"], b_summary["new"]) == (1, 1, 1)

    # Folder B local_id -> merged mapping
    bmap = res["local_id_maps"][str(folders[1])]["V"]
    assert bmap[1] == {"merged_id": 1, "merged_source_id": 1001}     # exact
    assert bmap[2] == {"merged_id": 2, "merged_source_id": 1002}     # positional -> A's sid
    assert bmap[3] == {"merged_id": 4, "merged_source_id": 2001}     # new


def test_materialize_writes_merged_photometry(merge_workspace):
    tmp_path, folders, folder_tags, catalogs_by_folder = merge_workspace
    res = _reconcile(folders, catalogs_by_folder, folder_tags)

    out_dir = tmp_path / "MERGED"
    info = materialize_merged_workspace(
        out_dir=out_dir,
        folders=folders,
        folder_tags=folder_tags,
        local_id_maps=res["local_id_maps"],
        merged_catalogs=res["canonical_by_filter"],
        selection_target_by_filter={"V": 1001},
        selection_comp_by_filter={"V": {1002}},
        selection_check_by_filter={"V": None},
        match_records=res["match_records"],
    )

    # Merged index has both frames.
    idx = pd.read_csv(step7_forced_phot_dir(out_dir) / "photometry_index.csv")
    merged_files = sorted(idx["file"].tolist())
    assert len(merged_files) == 2
    assert any("frameA1.fits" in f for f in merged_files)
    assert any("frameB1.fits" in f for f in merged_files)

    # Folder B's merged photometry: local source_ids remapped to merged identity.
    b_tag = folder_tags[str(folders[1])]
    b_phot = read_csv_int64_source_id(
        step7_forced_phot_dir(out_dir) / f"photometry_{b_tag}__frameB1.fits.tsv", sep="\t"
    ).sort_values("ID").reset_index(drop=True)
    assert b_phot["ID"].astype("int64").tolist() == [1, 2, 4]
    assert b_phot["source_id"].astype("int64").tolist() == [1001, 1002, 2001]

    # Two distinct merged nights (one per folder).
    assert len(set(info["night_assignments"].values())) == 2


def test_reconcile_positional_collision_and_tolerance():
    """Two rows nearest the same base star -> first matches, second is new;
    a row beyond tolerance is new. Locks the used_sid / tolerance edges that
    the vectorized positional matcher must preserve."""
    folder_a = Path("RESULT_A")
    folder_b = Path("RESULT_B")
    folders = [folder_a, folder_b]
    folder_tags = {str(f): folder_tag(i, f) for i, f in enumerate(folders)}

    df_a = pd.DataFrame([
        {"ID": 1, "source_id": 1001, "ra_deg": 10.0, "dec_deg": 20.0},
    ])
    # Row order matters: ID 10 is processed first and claims star 1001;
    # ID 11 is also nearest 1001 but it's taken -> new; ID 12 is far -> new.
    df_b = pd.DataFrame([
        {"ID": 10, "source_id": -1, "ra_deg": 10.00001, "dec_deg": 20.00001},
        {"ID": 11, "source_id": -2, "ra_deg": 10.00002, "dec_deg": 20.00002},
        {"ID": 12, "source_id": 3003, "ra_deg": 10.10000, "dec_deg": 20.10000},
    ])
    catalogs_by_folder = {
        str(folder_a): {"V": df_a},
        str(folder_b): {"V": df_b},
    }

    res = reconcile_workspace_catalogs(
        folders, catalogs_by_folder, folder_tags, pos_tol_arcsec=2.0,
    )

    b_summary = [r for r in res["match_summary_rows"]
                 if r["folder"] == folder_b.name and r["filter"] == "V"][0]
    assert (b_summary["exact"], b_summary["pos"], b_summary["new"]) == (0, 1, 2)

    bmap = res["local_id_maps"][str(folder_b)]["V"]
    assert bmap[10]["merged_source_id"] == 1001          # positional winner
    assert bmap[11]["merged_id"] != bmap[10]["merged_id"]  # collision -> distinct new id
    assert bmap[12]["merged_source_id"] == 3003          # far star keeps its own sid
    # canonical now has the base star plus the two new ones
    assert sorted(res["canonical_by_filter"]["V"]["ID"].astype(int).tolist()) == [1, 2, 3]


def test_compute_selection_defaults_filters_to_available():
    pytest.importorskip("PyQt5")  # multi_night_merger imports Qt at module load
    from apex.gui.tools.multi_night_merger import compute_selection_defaults

    merged = {"V": pd.DataFrame({"source_id": [1001, 1002, 1003]})}
    base = {"V": {"target_source_id": 1001,
                  "comparison_source_ids": [1002, 9999],   # 9999 absent -> dropped
                  "check_source_id": 1003}}
    t, c, k = compute_selection_defaults(merged, base)
    assert t["V"] == 1001
    assert c["V"] == {1002}
    assert k["V"] == 1003

    base_missing = {"V": {"target_source_id": 8888,
                          "comparison_source_ids": [],
                          "check_source_id": None}}
    t2, c2, k2 = compute_selection_defaults(merged, base_missing)
    assert t2["V"] is None and c2["V"] == set() and k2["V"] is None


def test_legacy_workspace_read_compat(tmp_path):
    """Pre-renumbering RESULT_* layout (step5_photometry / step9_selection /
    step10_lightcurve, with per-frame source_id recoverable only via
    cache/idmatch/) must still be recognized and loadable by the merger."""
    from apex.utils.photometry_loader import load_frame_photometry
    from apex.analysis.merge.workspace_scan import scan_merge_input_workspace

    rd = tmp_path / "RESULT_legacy"
    (rd / "step5_photometry").mkdir(parents=True)
    (rd / "step9_selection").mkdir(parents=True)
    (rd / "step10_lightcurve").mkdir(parents=True)
    (rd / "cache" / "idmatch").mkdir(parents=True)

    fname = "frameA.fit"
    # legacy forced-phot TSV: det_uid only, no source_id / ID
    pd.DataFrame([
        {"det_uid": 0, "FILTER": "V", "mag": 15.0},
        {"det_uid": 1, "FILTER": "V", "mag": 16.0},
    ]).to_csv(rd / "step5_photometry" / f"{fname}_photometry.tsv", sep="\t", index=False)
    pd.DataFrame([{"file": fname, "filter": "V", "night_id": 1}]).to_csv(
        rd / "step5_photometry" / "photometry_index.csv", index=False)
    # legacy idmatch: det_idx -> source_id
    pd.DataFrame([
        {"det_idx": 0, "source_id": 1001, "file": fname, "filter": "V"},
        {"det_idx": 1, "source_id": 1002, "file": fname, "filter": "V"},
    ]).to_csv(rd / "cache" / "idmatch" / f"idmatch_{fname}.csv", index=False)
    # selection master: source_id -> stable ID
    pd.DataFrame([
        {"ID": 5, "source_id": 1001},
        {"ID": 6, "source_id": 1002},
    ]).to_csv(rd / "step9_selection" / "master_catalog_V.tsv", sep="\t", index=False)
    (rd / "step9_selection" / "selection_V.json").write_text(
        '{"filter":"V","target_source_id":1001,"comparison_source_ids":[1002],"check_source_id":null}',
        encoding="utf-8")
    pd.DataFrame([{"ID": 5, "mjd": 1.0, "mag": 15.0}]).to_csv(
        rd / "step10_lightcurve" / "lightcurve_ID5_raw.csv", index=False)

    info = scan_merge_input_workspace(rd)
    assert info["merge_ready"]
    assert info["has_step7"] and info["has_step8"] and info["has_step9"]

    df = load_frame_photometry(rd, fname, "V")
    assert df is not None
    assert sorted(df["source_id"].astype("int64").tolist()) == [1001, 1002]  # via cache/idmatch
    assert sorted(df["ID"].astype("int64").tolist()) == [5, 6]               # via selection master

    # Step 9 is not a merge input: the merged workspace rebuilds the light
    # curves, so requiring it per night was busywork.
    for lc in (rd / "step10_lightcurve").glob("lightcurve_*.csv"):
        lc.unlink()
    info = scan_merge_input_workspace(rd)
    assert info["merge_ready"]
    assert not info["has_step9"]

    # Step 7 and Step 8 remain hard requirements.
    (rd / "step9_selection" / "master_catalog_V.tsv").unlink()
    info = scan_merge_input_workspace(rd)
    assert not info["merge_ready"]
    assert info["has_step7"] and not info["has_step8"]


def test_merge_worker_runs_full_chain_synchronously(merge_workspace):
    """Drive _MergeWorker.run() directly (no event loop): exercises the exact
    reconcile -> selection-defaults -> materialize chain the GUI now offloads."""
    pytest.importorskip("PyQt5.QtWidgets")
    from PyQt5.QtWidgets import QApplication
    from apex.gui.tools.multi_night_merger import _MergeWorker

    _app = QApplication.instance() or QApplication([])

    tmp_path, folders, folder_tags, catalogs_by_folder = merge_workspace
    out_dir = tmp_path / "MERGED_WORKER"

    captured: dict = {}
    failures: list = []
    worker = _MergeWorker(
        folders=list(folders),
        catalogs_by_folder=catalogs_by_folder,
        folder_tags=dict(folder_tags),
        match_radius=2.0,
        base_selection_by_filter={
            "V": {"target_source_id": 1001,
                  "comparison_source_ids": [1002],
                  "check_source_id": None},
        },
        out_dir=out_dir,
        cached_reconcile=None,
    )
    worker.finished_ok.connect(captured.update)
    worker.failed.connect(failures.append)
    worker.run()  # synchronous; signals deliver inline on this thread

    assert not failures, failures
    assert captured, "worker did not emit finished_ok"
    assert captured["sel_t"]["V"] == 1001
    assert captured["sel_c"]["V"] == {1002}
    assert captured["out_dir"] == out_dir
    assert (step7_forced_phot_dir(out_dir) / "photometry_index.csv").exists()
    # one merged night per input folder
    assert len(set(captured["build_info"]["night_assignments"].values())) == 2
