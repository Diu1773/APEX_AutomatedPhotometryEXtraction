"""Does PSF photometry undo the crowding penalty that aperture photometry pays?

The zero-point fit in M13 scatters 41 mmag in B while M67 scatters 13 mmag in g,
and the two data sets have almost the same seeing (2.72" vs 2.67"). What differs
is how close the stars are: median FWHM / nearest-neighbour distance is 0.307 in
M13 against 0.119 in M67. Swapping the reference system for Gaia DR3 synthetic
photometry did not shrink the scatter, so the reference is not the cause.

Aperture photometry has no way to separate two overlapping stars — every
neighbour inside the aperture is counted as target flux, and every neighbour
inside the sky annulus inflates the background. PSF photometry fits a model to
each source simultaneously, so the blend is decomposed rather than summed. If
crowding is really what drives the scatter, PSF should weaken the dependence on
neighbour distance; if something else is responsible, the dependence survives.

The measurement is a single quantity computed identically for both photometry
sources: the Spearman correlation between |zero-point residual| and distance to
the nearest catalogued neighbour, plus the ratio of the crowded quartile's
median |residual| to the isolated quartile's.

Residuals come from step10's own fit path (inverse-variance weights,
sigma-clipped, quadratic colour term), not a re-implementation, so the numbers
are comparable to the published zero points. Neighbour distances use the full
master catalogue, not just the calibrators: a calibrator's nearest neighbour is
usually a star too faint to calibrate on, and restricting the search to
calibrators would report every star as isolated.

The crowding statistics keep the sigma-clipped stars. The published scatter is
an inlier scatter by construction, but the stars the clip throws away are the
ones a crowding hypothesis predicts, so dropping them would remove the effect
being measured. `n_inliers` (fit) and `n_crowding` (statistics) therefore
differ, and both are reported.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

REPO = Path(__file__).absolute().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (  # noqa: E402
    robust_weighted_polyfit,
)
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

SNR_CUT, CLIP, ITERS, MIN_N = 20.0, 3.0, 5, 10

BANDS = {
    "B": ("B", "V"),
    "V": ("B", "V"),
    "R": ("V", "R"),
    "g": ("g", "r"),
    "r": ("g", "r"),
    "i": ("r", "i"),
}


def neighbour_distance(master: pd.DataFrame, x: np.ndarray,
                       y: np.ndarray) -> np.ndarray:
    """Distance in pixels from each calibrator to its nearest catalogued star.

    `master_sources.csv` repeats a star once per filter, so positions are
    de-duplicated first — otherwise every star's nearest neighbour is itself at
    zero separation. k=2 then skips the self-match.
    """
    pos = (master[["x_ref", "y_ref"]]
           .apply(pd.to_numeric, errors="coerce")
           .dropna()
           .drop_duplicates()
           .to_numpy(float))
    dist, _ = cKDTree(pos).query(np.column_stack([x, y]), k=2)
    return dist[:, 1]


def residuals(cal: pd.DataFrame, band: str,
              color: tuple[str, str]) -> tuple[np.ndarray, np.ndarray, dict]:
    """Run step10's zero-point fit and return the surviving stars' residuals."""
    ref = pd.to_numeric(cal[f"ref_{band}"], errors="coerce").to_numpy(float)
    inst = pd.to_numeric(cal[f"mag_inst_{band}"], errors="coerce").to_numpy(float)
    err = pd.to_numeric(cal[f"mag_inst_err_{band}"], errors="coerce").to_numpy(float)
    snr = pd.to_numeric(cal.get(f"snr_{band}"), errors="coerce").to_numpy(float)
    ca = pd.to_numeric(cal[f"ref_{color[0]}"], errors="coerce").to_numpy(float)
    cb = pd.to_numeric(cal[f"ref_{color[1]}"], errors="coerce").to_numpy(float)

    delta, x = ref - inst, ca - cb
    w = np.where(np.isfinite(err) & (err > 0), 1.0 / err ** 2, np.nan)
    qual, _ = gaia_quality_report(cal, cstar_nsigma=None)
    m = (np.isfinite(delta) & np.isfinite(x) & np.isfinite(w)
         & np.isfinite(snr) & (snr >= SNR_CUT) & qual)
    if int(m.sum()) < MIN_N:
        return np.array([]), np.array([]), {}

    coeffs, n_inliers, scatter = robust_weighted_polyfit(
        x[m], delta[m], w=w[m], degree=2,
        clip_sigma=CLIP, iters=ITERS, min_n=MIN_N)
    if coeffs is None:
        return np.array([]), np.array([]), {}

    resid = delta - np.polyval(coeffs, x)
    info = {"n_fit": int(m.sum()), "n_inliers": int(n_inliers),
            "zp": float(coeffs[2]), "scatter_mmag": float(scatter) * 1000}
    return resid, m, info


def analyse(cal_path: Path, master: pd.DataFrame, label: str) -> dict:
    cal = pd.read_csv(cal_path)
    bands = [b for b in BANDS if f"mag_inst_{b}" in cal.columns]
    nn = neighbour_distance(
        master,
        pd.to_numeric(cal["x_pix"], errors="coerce").to_numpy(float),
        pd.to_numeric(cal["y_pix"], errors="coerce").to_numpy(float))

    out: dict = {"label": label, "path": str(cal_path), "bands": {}}
    for band in bands:
        source = cal.get(f"photometry_source_{band}")
        out["bands"].setdefault(band, {})["photometry_source"] = (
            str(source.dropna().iloc[0]) if source is not None
            and source.notna().any() else "unknown")

        resid, mask, info = residuals(cal, band, BANDS[band])
        if resid.size == 0:
            out["bands"][band]["error"] = "fit failed"
            continue

        keep = mask & np.isfinite(resid) & np.isfinite(nn)
        a, d = np.abs(resid[keep]), nn[keep]
        if a.size < 40:
            out["bands"][band].update(info, error="too few stars for quartiles")
            continue

        rho, p = spearmanr(d, a)
        q1, q3 = np.percentile(d, [25, 75])
        crowded = float(np.median(a[d <= q1])) * 1000
        isolated = float(np.median(a[d >= q3])) * 1000
        out["bands"][band].update(
            info,
            n_crowding=int(a.size),
            spearman_r=float(rho), spearman_p=float(p),
            q1_px=float(q1), q3_px=float(q3),
            crowded_mmag=crowded, isolated_mmag=isolated,
            ratio=crowded / isolated if isolated > 0 else float("nan"),
        )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result-dir", required=True)
    ap.add_argument("--aperture-zp", default="cmd_zeropoint_APERTURE_BACKUP")
    ap.add_argument("--psf-zp", default="cmd_zeropoint")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    root = Path(args.result_dir)
    master = pd.read_csv(root / "step7_forced_phot" / "master_sources.csv")

    runs = []
    for label, sub in (("aperture", args.aperture_zp), ("psf", args.psf_zp)):
        path = root / sub / "gaia_sdss_calibrator_by_ID.csv"
        if not path.exists():
            print(f"[skip] {label}: {path} 없음")
            continue
        runs.append(analyse(path, master, label))

    head = (f"{'band':>5}{'source':>10}{'N':>6}{'zp':>10}{'scatter':>10}"
            f"{'rho':>8}{'crowded':>10}{'isolated':>10}{'ratio':>8}")
    print(head)
    print("-" * len(head))
    for run in runs:
        for band, b in run["bands"].items():
            if "ratio" not in b:
                print(f"{band:>5}{b.get('photometry_source','?'):>10}"
                      f"   {b.get('error','')}")
                continue
            print(f"{band:>5}{b['photometry_source']:>10}{b['n_inliers']:>6}"
                  f"{b['zp']:>10.4f}{b['scatter_mmag']:>9.1f}m"
                  f"{b['spearman_r']:>8.3f}{b['crowded_mmag']:>9.1f}m"
                  f"{b['isolated_mmag']:>9.1f}m{b['ratio']:>8.2f}")

    out = Path(args.output) if args.output else (
        REPO / "validation" / f"crowding_aperture_vs_psf_{root.parent.name}.json")
    out.write_text(json.dumps(
        {"result_dir": str(root), "snr_cut": SNR_CUT, "runs": runs},
        indent=1), encoding="utf-8")
    print(f"\nsaved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
