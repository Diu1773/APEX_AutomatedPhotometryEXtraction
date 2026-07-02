"""Generate the canonical paper-quality artificial-star dataset (run once).

Writes stars.csv / magnitude_points.csv / magnitude_bins.csv / completeness_fit.json
under validation/paper/data/artificial_star/benchmark_run/ for the figure scripts
to consume. Kept short (under the repo) to stay clear of the Windows MAX_PATH limit.
"""
import sys, time, json, shutil
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

from apex.benchmark.validate import run_artificial_star_suite

OUT = REPO / "validation" / "paper" / "data"
shutil.rmtree(OUT, ignore_errors=True)

t0 = time.perf_counter()
res = run_artificial_star_suite(
    OUT,
    reference_frame=None,     # self-contained synthetic frame
    config=None,
    seed=20260702,
    trials=12,
    stars_per_trial=70,
    bootstrap_samples=500,
)
dt = time.perf_counter() - t0

summary = {
    "elapsed_s": round(dt, 1),
    "status": res.get("status"),
    "n_injected": res.get("n_injected"),
    "m50": res.get("completeness_m50"),
    "m90": res.get("completeness_m90"),
    "m10": res.get("completeness_m10"),
    "m50_ci95": res.get("completeness_m50_ci95"),
    "width_mag": res.get("completeness_width_mag"),
    "rms_bright": res.get("photometric_rms_bright_mag"),
    "rms_faint": res.get("photometric_rms_faint_mag"),
    "bias_median": res.get("photometric_bias_median_mag"),
    "pos_rmse_px": res.get("position_rmse_px"),
    "false_positive_rate": res.get("false_positive_rate"),
    "non_degenerate": res.get("completeness_non_degenerate"),
    "run_dir": res.get("benchmark_run_dir"),
}
(OUT / "canonical_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
