"""Tests for apex.analysis.light_curve.global_ensemble — inhomogeneous ensemble solver."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.light_curve.global_ensemble import (
    _robust_sigma,
    comparison_network_diagnostics,
    generate_synthetic_data,
    run_synthetic_test,
    select_comparisons_from_qc,
    solve_global_ensemble,
)


# ---------------------------------------------------------------------------
# _robust_sigma
# ---------------------------------------------------------------------------

def test_robust_sigma_gaussian():
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, 1000)
    sig = _robust_sigma(x)
    assert abs(sig - 1.0) < 0.1


def test_robust_sigma_with_outliers():
    rng = np.random.default_rng(42)
    # Base distribution near 0 with small spread, plus two huge outliers.
    # MAD stays non-zero because base values have small random scatter.
    x = np.concatenate([rng.normal(0, 0.05, 200), [100.0, -100.0]])
    sig = _robust_sigma(x)
    # MAD-based sigma ≈ 0.05 * 1.4826, well below 1.0
    assert sig < 1.0, f"_robust_sigma={sig:.4f} expected < 1.0 for outlier-dominated data"


def test_robust_sigma_constant():
    sig = _robust_sigma(np.ones(50))
    assert sig == 0.0 or sig < 1e-10


# ---------------------------------------------------------------------------
# solve_global_ensemble — validation on synthetic data
# ---------------------------------------------------------------------------

def _make_simple_df(n_frames=30, n_comps=10, noise=0.01, seed=7):
    rng = np.random.default_rng(seed)
    target_id = 1
    comp_ids = list(range(2, 2 + n_comps))
    jd = np.linspace(2450000.0, 2450000.0 + 2.0, n_frames)
    Z_true = rng.normal(0, 0.05, n_frames)
    M_comp = rng.normal(12.0, 0.3, n_comps)
    rows = []
    for t, (tid, jd_t) in enumerate(zip([f"f{t:03d}" for t in range(n_frames)], jd)):
        for i, sid in enumerate(comp_ids):
            rows.append({"time_id": tid, "jd": jd_t, "filter": "V", "star_id": sid,
                         "mag_inst": M_comp[i] + Z_true[t] + rng.normal(0, noise),
                         "err": noise})
        rows.append({"time_id": tid, "jd": jd_t, "filter": "V", "star_id": target_id,
                     "mag_inst": 13.0 + Z_true[t] + rng.normal(0, noise), "err": noise})
    return pd.DataFrame(rows), target_id, comp_ids, Z_true


def test_solve_returns_required_keys():
    df, tid, cids, _ = _make_simple_df()
    result = solve_global_ensemble(df, tid, cids, min_comps=3, n_iter=1)
    assert "zp_df" in result
    assert "mean_df" in result
    assert "lc_df" in result
    assert "diagnostics" in result


def test_solve_zp_correlation_with_injected():
    """Recovered frame ZPs should correlate strongly with the injected ones."""
    df, tid, cids, Z_true = _make_simple_df(n_frames=60, n_comps=15, noise=0.005)
    result = solve_global_ensemble(df, tid, cids, min_comps=5, n_iter=2)
    zp_df = result["zp_df"]
    z_rec = pd.to_numeric(zp_df["Z"], errors="coerce").dropna().to_numpy()
    ok = np.isfinite(z_rec)
    assert ok.sum() > 10
    corr = np.corrcoef(z_rec[ok], Z_true[ok])[0, 1]
    assert corr > 0.90, f"ZP correlation too low: {corr:.3f}"


def test_solve_lc_lower_scatter_than_raw():
    df, tid, cids, _ = _make_simple_df(n_frames=50, n_comps=10, noise=0.01)
    result = solve_global_ensemble(df, tid, cids, min_comps=5, n_iter=2)
    lc = result["lc_df"]
    assert not lc.empty

    raw  = pd.to_numeric(lc.get("diff_mag_raw",  pd.Series(dtype=float)), errors="coerce").dropna()
    corr = pd.to_numeric(lc.get("diff_mag_corr", pd.Series(dtype=float)), errors="coerce").dropna()
    if len(raw) > 5 and len(corr) > 5:
        assert np.std(corr) <= np.std(raw) + 1e-9, (
            f"Corrected scatter {np.std(corr):.5f} > raw {np.std(raw):.5f}"
        )


def test_solve_mean_df_has_most_comp_ids():
    """mean_df must contain most comp IDs; rms_clip_pct may remove a small fraction."""
    df, tid, cids, _ = _make_simple_df()
    result = solve_global_ensemble(df, tid, cids, min_comps=3, n_iter=1, rms_clip_pct=0.0)
    mean_df = result["mean_df"]
    recovered_ids = set(mean_df["star_id"].tolist())
    for cid in cids:
        assert cid in recovered_ids, f"comp {cid} missing from mean_df (rms_clip_pct=0)"


def test_solve_zp_df_has_all_frame_ids():
    df, tid, cids, _ = _make_simple_df(n_frames=20, n_comps=8)
    result = solve_global_ensemble(df, tid, cids, min_comps=3, n_iter=1)
    zp_df = result["zp_df"]
    all_tids = set(df["time_id"].unique())
    assert set(zp_df["time_id"].tolist()) == all_tids


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def test_empty_df_raises():
    with pytest.raises((ValueError, Exception)):
        solve_global_ensemble(pd.DataFrame(), target_id=1, comp_ids=[2, 3])


def test_missing_columns_raises():
    df = pd.DataFrame({"time_id": ["a"], "star_id": [1], "mag_inst": [12.0]})
    with pytest.raises((ValueError, Exception)):
        solve_global_ensemble(df, target_id=1, comp_ids=[2])


def test_empty_comp_ids_raises():
    df, tid, cids, _ = _make_simple_df()
    with pytest.raises((ValueError, Exception)):
        solve_global_ensemble(df, target_id=tid, comp_ids=[])


def test_too_few_comps_after_clip_raises():
    df, tid, cids, _ = _make_simple_df(n_comps=3)
    with pytest.raises(Exception):
        solve_global_ensemble(df, target_id=tid, comp_ids=cids, min_comps=10)


# ---------------------------------------------------------------------------
# Per-filter mode
# ---------------------------------------------------------------------------

def test_per_filter_produces_one_block_per_filter():
    df, tid, cids, _ = _make_simple_df()
    df2 = df.copy()
    df2["filter"] = "B"
    combined = pd.concat([df, df2], ignore_index=True)
    result = solve_global_ensemble(combined, tid, cids, min_comps=3, per_filter=True)
    zp = result["zp_df"]
    assert set(zp["filter"].unique()) == {"V", "B"}


# ---------------------------------------------------------------------------
# generate_synthetic_data + run_synthetic_test
# ---------------------------------------------------------------------------

def test_generate_synthetic_data_shape():
    df, tid, cids, truth = generate_synthetic_data(n_frames=50, n_comps=10)
    n_rows = 50 * (10 + 1)
    assert len(df) == n_rows
    assert set(df.columns) >= {"time_id", "jd", "filter", "star_id", "mag_inst", "err"}


def test_run_synthetic_test_high_correlation():
    result = run_synthetic_test()
    corr = result["corr_Z"]
    assert np.isfinite(corr), "correlation is non-finite"
    assert corr > 0.85, f"ZP recovery correlation too low: {corr:.3f}"


def _make_inhomogeneous_network(disconnected: bool = False, seed: int = 22):
    rng = np.random.default_rng(seed)
    target_id = 1
    rows = []
    for frame in range(12):
        first_block = frame < 6
        if first_block:
            stars = [(2, 12.0), (3, 12.5), (4, 13.0)]
        else:
            bridge_id = 7 if disconnected else 4
            stars = [(bridge_id, 13.0), (5, 13.5), (6, 14.0)]
        zeropoint = 0.08 * np.sin(frame / 2.0) + (0.35 if not first_block else 0.0)
        for star_id, mean_mag in stars:
            rows.append(
                {
                    "time_id": f"f{frame:02d}",
                    "jd": 2450000.0 + frame,
                    "filter": "V",
                    "star_id": star_id,
                    "mag_inst": mean_mag + zeropoint + rng.normal(0, 0.002),
                    "err": 0.005,
                }
            )
        rows.append(
            {
                "time_id": f"f{frame:02d}",
                "jd": 2450000.0 + frame,
                "filter": "V",
                "star_id": target_id,
                "mag_inst": 11.2 + zeropoint + rng.normal(0, 0.002),
                "err": 0.005,
            }
        )
    comp_ids = [2, 3, 4, 5, 6] + ([7] if disconnected else [])
    return pd.DataFrame(rows), target_id, comp_ids


def test_disconnected_comparison_network_is_rejected():
    df, target_id, comp_ids = _make_inhomogeneous_network(disconnected=True)
    comp_df = df[df["star_id"].isin(comp_ids)]
    diagnostic = comparison_network_diagnostics(comp_df)
    assert diagnostic["connected"] is False
    assert diagnostic["n_components"] == 2
    with pytest.raises(ValueError, match="disconnected"):
        solve_global_ensemble(
            df,
            target_id,
            comp_ids,
            min_comps=3,
            n_iter=1,
            rms_clip_pct=0,
        )


def test_shared_comparison_connects_nights_and_preserves_target_baseline():
    df, target_id, comp_ids = _make_inhomogeneous_network(disconnected=False)
    result = solve_global_ensemble(
        df,
        target_id,
        comp_ids,
        min_comps=3,
        n_iter=1,
        rms_clip_pct=0,
    )
    lc = result["lc_df"].sort_values("jd")
    corrected = lc["diff_mag_corr"].to_numpy(float)
    assert np.nanstd(corrected) < 0.01
    assert abs(np.nanmean(corrected[:6]) - np.nanmean(corrected[6:])) < 0.01
    network = result["diagnostics"]["filters"]["V"]["network"]
    assert network["connected"] is True


def test_missing_comparison_errors_are_imputed_and_reported():
    df, target_id, comp_ids, _ = _make_simple_df(n_frames=12, n_comps=5)
    df.loc[df["star_id"].isin(comp_ids), "err"] = np.nan
    result = solve_global_ensemble(
        df, target_id, comp_ids, min_comps=3, n_iter=1, rms_clip_pct=0
    )
    assert result["zp_df"]["Z"].notna().all()
    diagnostic = result["diagnostics"]["filters"]["V"]
    assert diagnostic["n_errors_imputed"] == 12 * 5


def test_rms_clip_percentage_is_a_cap_not_forced_removal():
    df, target_id, comp_ids, _ = _make_simple_df(n_frames=30, n_comps=10)
    result = solve_global_ensemble(
        df,
        target_id,
        comp_ids,
        min_comps=3,
        n_iter=3,
        sigma=100.0,
        rms_clip_pct=50.0,
    )
    assert result["diagnostics"]["filters"]["V"]["removed_comps"] == []


def test_high_rms_comparison_is_removed_and_reason_is_reported():
    df, target_id, comp_ids, _ = _make_simple_df(
        n_frames=50, n_comps=10, noise=0.003
    )
    bad_id = comp_ids[-1]
    mask = df["star_id"] == bad_id
    df.loc[mask, "mag_inst"] += 0.08 * np.sin(np.linspace(0, 8 * np.pi, mask.sum()))
    result = solve_global_ensemble(
        df,
        target_id,
        comp_ids,
        min_comps=5,
        n_iter=3,
        sigma=3.0,
        rms_clip_pct=20.0,
    )
    diagnostic = result["diagnostics"]["filters"]["V"]
    assert bad_id in diagnostic["removed_comps"]
    details = {item["star_id"]: item for item in diagnostic["removed_comp_details"]}
    assert "relative_rms" in details[bad_id]["reasons"]


def test_meanz0_gauge_propagates_nonzero_reference_frame_error():
    df, target_id, comp_ids, _ = _make_simple_df(n_frames=20, n_comps=8)
    result = solve_global_ensemble(
        df,
        target_id,
        comp_ids,
        min_comps=3,
        n_iter=1,
        rms_clip_pct=0,
        gauge="meanZ0",
    )
    assert np.all(result["zp_df"]["Z_err"].to_numpy(float) > 0)


def test_qc_selection_uses_absolute_thresholds_when_enough_pass():
    qc = pd.DataFrame(
        {
            "comp_id": [2, 3, 4, 5],
            "n": [100, 100, 100, 100],
            "rms": [0.008, 0.009, 0.011, 0.08],
            "sigma_nights": [0.004, 0.005, 0.006, 0.2],
            "outlier_frac": [0.01, 0.02, 0.01, 0.3],
            "coverage_fraction": [1.0, 1.0, 1.0, 1.0],
        }
    )
    result = select_comparisons_from_qc(qc, 0.02, 0.1, 10)
    assert result["method"] == "absolute"
    assert result["selected_ids"] == [2, 3, 4]


def test_qc_selection_relative_fallback_rejects_unstable_low_coverage_comps():
    qc = pd.DataFrame(
        {
            "comp_id": [1, 6, 7, 8, 9, 12, 14, 15, 16],
            "n": [518, 527, 527, 527, 527, 527, 363, 363, 527],
            "rms": [0.123, 0.559, 0.223, 0.246, 0.217, 0.243, 0.204, 0.191, 0.227],
            "sigma_nights": [0.178, 0.761, 0.062, 0.182, 0.058, 0.061, 0.050, 0.012, 0.063],
            "outlier_frac": [0.328, 0.347, 0.021, 0.017, 0.023, 0.017, 0.063, 0.091, 0.019],
            "coverage_fraction": [0.983, 1.0, 1.0, 1.0, 1.0, 1.0, 0.689, 0.689, 1.0],
        }
    )
    result = select_comparisons_from_qc(qc, 0.02, 0.1, 10)
    assert result["method"] == "relative_fallback"
    assert len(result["selected_ids"]) == 3
    assert 1 not in result["selected_ids"]
    assert 6 not in result["selected_ids"]
    assert 14 not in result["selected_ids"]
    assert 15 not in result["selected_ids"]
