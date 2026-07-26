from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd

from apex.analysis.light_curve.check_star_io import load_check_star_csv
from apex.analysis.light_curve.lightcurve_output_service import (
    save_combined_raw_outputs,
    save_dataset_raw_outputs,
)
from apex.gui.workflow.lc.step9_lightcurve_builder import (
    LightCurveBuilderWindow,
    _active_comparison_ids_for_filter,
    _load_selection_ids_by_filter,
)
from apex.gui.workflow.lc.step10_detrend_merge import (
    DetrendNightMergeWindow,
    _load_step9_comparison_ids_by_filter,
)
from apex.utils.step_paths_lc import (
    step8_selection_dir,
    step9_lc_dir,
    step10_current_lc_path,
)


def _write_selection(
    result_dir,
    *,
    target_id: int = 1,
    comparison_ids: list[int] | None = None,
    comparison_source_ids: list[int] | None = None,
    check_id: int | None = 7,
    check_source_id: int | None = 7007,
) -> None:
    selection_dir = step8_selection_dir(result_dir)
    selection_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "target_id": target_id,
        "target_source_id": 1001,
        "comparison_ids": comparison_ids or [],
        "comparison_source_ids": comparison_source_ids or [],
        "check_id": check_id,
        "check_source_id": check_source_id,
    }
    (selection_dir / "selection_g.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )


def _raw_curve(jd: list[float], dataset: str = "night") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "JD": jd,
            "filter": ["g"] * len(jd),
            "diff_mag_raw": [0.01 * i for i in range(len(jd))],
            "dataset": [dataset] * len(jd),
        }
    )


def test_selection_recovery_preserves_id_source_pair_order(tmp_path):
    _write_selection(
        tmp_path,
        comparison_ids=[10, 20],
        comparison_source_ids=[2002, 2001],
    )
    pd.DataFrame(
        {
            "ID": [10, 20],
            "source_id": [2001, 2002],
        }
    ).to_csv(
        step8_selection_dir(tmp_path) / "master_catalog_g.tsv",
        sep="\t",
        index=False,
    )

    selection = _load_selection_ids_by_filter(tmp_path)["g"]

    assert selection["comparison_ids"] == [20, 10]
    assert selection["comparison_source_ids"] == [2002, 2001]
    assert _active_comparison_ids_for_filter(selection, [10, 20, 30]) == [20, 10]
    assert _load_step9_comparison_ids_by_filter(tmp_path) == {"g": {10, 20}}


def test_dataset_save_removes_stale_check_outputs(tmp_path):
    output_dir = step9_lc_dir(tmp_path)
    output_dir.mkdir(parents=True)
    stale = output_dir / "lightcurve_check_g_ID99_raw.csv"
    _raw_curve([1.0]).to_csv(stale, index=False)

    save_dataset_raw_outputs(
        result_dir=tmp_path,
        target_id=1,
        raw_df=_raw_curve([1.0]),
        check_df=pd.DataFrame(),
    )

    assert not stale.exists()
    assert (output_dir / "lightcurve_ID1_raw.csv").exists()


def test_builder_reads_legacy_forced_photometry_index(tmp_path):
    legacy_dir = tmp_path / "step5_photometry"
    legacy_dir.mkdir()
    pd.DataFrame(
        {"file": ["legacy.fits"], "filter": ["g"]}
    ).to_csv(legacy_dir / "photometry_index.csv", index=False)
    window = LightCurveBuilderWindow.__new__(LightCurveBuilderWindow)
    window.params = SimpleNamespace(P=SimpleNamespace(result_dir=tmp_path))
    window.project_state = None
    window._photometry_source_cache = {}
    window._force_aperture_for_datasets = False

    loaded = window._load_active_photometry_index(tmp_path)

    assert loaded["file"].tolist() == ["legacy.fits"]


def test_combined_check_curve_keeps_all_datasets(tmp_path):
    _write_selection(tmp_path)
    check_a = _raw_curve([1.0], "night-a").assign(
        check_id=7, check_source_id=7007
    )
    check_b = _raw_curve([2.0], "night-b").assign(
        check_id=8, check_source_id=8008
    )

    save_combined_raw_outputs(
        base_result_dir=tmp_path,
        target_id=1,
        combined_raw=[_raw_curve([1.0]), _raw_curve([2.0])],
        single_dataset_mode=False,
        comp_candidate_ids=[10, 20, 30],
        active_comp_ids=[10, 20, 30],
        combined_check=[check_a, check_b],
        check_ids_by_filter={"g": 7},
    )

    check_id, loaded = load_check_star_csv(tmp_path, filt="g")

    assert check_id is None
    assert loaded["JD"].tolist() == [1.0, 2.0]
    assert loaded["dataset"].tolist() == ["night-a", "night-b"]


def test_builder_validation_requires_a_usable_check_curve(tmp_path):
    _write_selection(tmp_path, comparison_ids=[10, 20, 30])
    output_dir = step9_lc_dir(tmp_path)
    output_dir.mkdir(parents=True)
    _raw_curve([1.0]).to_csv(output_dir / "lightcurve_ID1_raw.csv", index=False)
    window = SimpleNamespace(
        params=SimpleNamespace(P=SimpleNamespace(result_dir=tmp_path)),
        datasets=[("base", tmp_path)],
    )

    assert LightCurveBuilderWindow.validate_step(window) is False

    _raw_curve([1.0]).assign(check_id=7, check_source_id=7007).to_csv(
        output_dir / "lightcurve_check_combined_raw.csv", index=False
    )
    assert LightCurveBuilderWindow.validate_step(window) is True


def test_detrend_validation_ignores_previous_step_raw_output(tmp_path):
    output_dir = step9_lc_dir(tmp_path)
    output_dir.mkdir(parents=True)
    _raw_curve([1.0]).to_csv(output_dir / "lightcurve_ID1_raw.csv", index=False)
    window = SimpleNamespace(
        params=SimpleNamespace(P=SimpleNamespace(result_dir=tmp_path)),
    )

    assert DetrendNightMergeWindow.validate_step(window) is False

    current_path = step10_current_lc_path(tmp_path, 1)
    current_path.parent.mkdir(parents=True, exist_ok=True)
    _raw_curve([1.0]).assign(diff_mag_corr=[0.0]).to_csv(current_path, index=False)
    assert DetrendNightMergeWindow.validate_step(window) is True
