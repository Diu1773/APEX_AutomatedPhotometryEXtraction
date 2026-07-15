import numpy as np
import pandas as pd

from apex.benchmark.psf_repeatability import _pairwise_products, analyze_psf_repeatability, fit_additive_model, match_psf_to_step7


def _series(n_stars=1000, zps=(0.0, .31, -.18), sigma=.01, sigmas=None, seed=4):
    rng = np.random.default_rng(seed); rows = []
    mags = np.linspace(14, 19, n_stars)
    for j, zp in enumerate(zps):
        for star, mag in enumerate(mags, 1):
            noise = sigmas[j] if sigmas is not None else sigma
            rows.append({"master_id": star, "frame": f"f{j}", "mag_psf": mag + zp + rng.normal(0, noise), "mag_psf_err": noise,
                         "snr_psf": 10.0 + 30.0 * (star % 4), "flags_psf": 0})
    return pd.DataFrame(rows)


def test_exact_cv_n3_has_expected_noise_scale_and_unit_pull():
    sigma = .01; result = fit_additive_model(_series(sigma=sigma, seed=12))
    cv = result["exact_leave_one_out_cv_residual"]; pull = result["exact_leave_one_out_cv_pull"]
    expected = np.sqrt(1.0 + 1.0 / 2.0) * sigma
    measured = 1.4826 * np.median(np.abs(cv - np.median(cv)))
    pull_scatter = 1.4826 * np.median(np.abs(pull - np.median(pull)))
    assert np.isclose(measured, expected, rtol=.12)
    assert np.isclose(pull_scatter, 1.0, rtol=.12)
    assert not np.allclose(result["raw_in_sample_residual"], cv)


def test_frame_zeropoints_are_recovered():
    result = fit_additive_model(_series(n_stars=300, seed=8))
    offsets = result["frame_zp"]
    assert np.allclose([offsets["f0"], offsets["f1"], offsets["f2"]], [0, .31, -.18], atol=.035)


def test_pair_scatter_and_nnls_recover_heteroscedastic_frame_noise():
    sigmas = (.01, .03, .06); data = _series(n_stars=2500, sigmas=sigmas, seed=21)
    fit = fit_additive_model(data); data["corrected_mag"] = fit["corrected_mag"]; data["frame_fwhm_px"] = 5.; data["frame_fwhm_arcsec"] = 2.
    data["star_min_snr"] = data.groupby("master_id")["snr_psf"].transform("min")
    pair_sets = []; noise_sets = []
    for cut in (0, 20, 50, 100):
        pairs, noise = _pairwise_products(data, cut); pair_sets.append(pairs); noise_sets.append(noise)
    pairs = pair_sets[0]; noise = noise_sets[0]
    expected_pairs = {"f0__f1": np.hypot(.01, .03), "f0__f2": np.hypot(.01, .06), "f1__f2": np.hypot(.03, .06)}
    for label, expected in expected_pairs.items():
        measured = pairs.loc[pairs.pair_label == label, "robust_mad"].iloc[0]
        assert np.isclose(measured, expected, rtol=.10)
    robust = noise[noise.noise_metric == "robust"].set_index("frame").frame_scatter
    assert np.allclose([robust["f0"], robust["f1"], robust["f2"]], sigmas, atol=.012)
    assert set(noise.noise_metric) == {"robust", "rmse"}
    expected_counts = [2500, 1875, 1250, 625]
    assert [int(pairs.n_stars.iloc[0]) for pairs in pair_sets] == expected_counts
    assert all(set(p.n_stars) == {n} for p, n in zip(pair_sets, expected_counts))
    assert all(set(n.n_stars) == {count} for n, count in zip(noise_sets, expected_counts))


def test_robust_fit_survives_variable_outliers():
    data = _series(n_stars=300); data.loc[(data.master_id == 3) & (data.frame == "f1"), "mag_psf"] += 2
    data.loc[(data.master_id == 4) & (data.frame == "f2"), "mag_psf"] += 1
    result = fit_additive_model(data)
    assert abs(result["frame_zp"]["f1"] - .31) < .06 and abs(result["frame_zp"]["f2"] + .18) < .06


def test_min_frames_excludes_single_observation_stars(tmp_path):
    result = tmp_path / "result"; psf = result / "cmd_psf"; forced = result / "step7_forced_phot"; psf.mkdir(parents=True); forced.mkdir()
    for frame in range(3):
        rows = []
        for star in range(1, 5):
            rows.append({"x_fit": float(star), "y_fit": 0., "master_id": star, "mag_psf": 15 + frame * .1 + star / 10, "mag_psf_err": .01, "flags_psf": 0})
        if frame == 2: rows.append({"x_fit": 20., "y_fit": 0., "master_id": 99, "mag_psf": 15, "mag_psf_err": .01, "flags_psf": 0})
        psf_df = pd.DataFrame(rows).drop(columns="master_id"); forced_df = pd.DataFrame(rows)[["master_id", "x_fit", "y_fit"]]
        name = f"photometry_f{frame}.tsv"; psf_df.to_csv(psf / name, sep="\t", index=False); forced_df.to_csv(forced / name, sep="\t", index=False)
    out = analyze_psf_repeatability(result, tmp_path / "out", min_frames=3)
    by_star = pd.read_csv(tmp_path / "out" / "psf_repeatability_by_star.csv")
    assert set(by_star.master_id) == {1, 2, 3, 4} and out["n_stars"] == 4


def test_invalid_and_flagged_rows_are_excluded():
    data = _series(n_stars=50); data.loc[0, "flags_psf"] = 8; data.loc[1, "mag_psf"] = np.nan
    good = data[(data.flags_psf == 0) & np.isfinite(data.mag_psf)]
    result = fit_additive_model(good)
    assert len(result["exact_leave_one_out_cv_residual"]) == len(good)


def test_step7_position_match_is_one_to_one_and_reports_distance():
    psf = pd.DataFrame({"x_fit": [10.1, 20.2, 20.3], "y_fit": [10, 20, 20], "mag_psf": [1, 2, 3]})
    forced = pd.DataFrame({"master_id": [7, 8], "x_fit": [10, 20], "y_fit": [10, 20]})
    out, meta = match_psf_to_step7(psf, forced, 1.0)
    assert out.master_id.notna().sum() == 2 and out.master_id.dropna().nunique() == 2
    assert meta["match_fraction"] == 2 / 3 and meta["match_distance_median_px"] < .3
