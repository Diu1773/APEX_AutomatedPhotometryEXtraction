"""
Shared aperture photometry worker — per-frame growth curve analysis and
aperture correction (apcorr) determination.

ApertureWorker outputs:
    aperture_by_frame.csv   – optimal r_ap/r_in/r_out per frame
    apcorr_summary.csv      – per-frame apcorr as flux ratio (ref_flux / opt_flux)
    growth_curve.csv        – full growth curve data

The ``apcorr`` column is a flux ratio.  Apply it as::

    flux_corrected = flux_measured * apcorr

or equivalently in magnitudes::

    mag_corrected = mag_measured - 2.5 * log10(apcorr)
"""
from __future__ import annotations

import json
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.spatial import cKDTree

from PyQt5.QtCore import QThread, pyqtSignal

from apex.utils.step_paths import (
    step2_cropped_dir,
    step4_dir,
    step5_aperture_dir,
    crop_is_active,
)
from apex.utils.photometry_utils import phot_one_star as _phot_one_star
from apex.utils.qc_utils import resolve_frame_quality_path


# ── Scalar helpers ─────────────────────────────────────────────────────────────

def _to_float(val, default):
    try:
        if val is None:
            return float(default)
        out = float(val)
        return out if np.isfinite(out) else float(default)
    except Exception:
        return float(default)


def _to_int(val, default):
    try:
        if val is None:
            return int(default)
        return int(float(val))
    except Exception:
        return int(default)


def _build_scale_grid(start, stop, step):
    start, stop, step = float(start), float(stop), float(step)
    if step <= 0:
        step = 0.25
    if stop < start:
        start, stop = stop, start
    n = int(np.floor((stop - start) / step + 1e-9)) + 1
    n = int(np.clip(n, 1, 200))
    vals = [start + k * step for k in range(n)]
    if not vals or vals[-1] < stop - 1e-9:
        vals.append(stop)
    return [float(round(v, 6)) for v in vals]


def _load_fwhm_from_meta(
    fname: str, cache_dir: Path, result_dir: Path, params_fwhm_guess: float = 6.0
) -> float:
    candidates = [
        cache_dir / f"detect_{fname}.json",
        step4_dir(result_dir) / f"detect_{fname}.json",
    ]
    existing = sorted(
        [c for c in candidates if c.exists()],
        key=lambda q: q.stat().st_mtime_ns,
        reverse=True,
    )
    for p in existing:
        try:
            meta = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        for k in ("fwhm_med_rad_px", "fwhm_med_px", "fwhm_px", "fwhm_med"):
            v = meta.get(k)
            if v is not None:
                try:
                    v = float(v)
                    if np.isfinite(v) and v > 0:
                        return v
                except Exception:
                    continue
    return float(params_fwhm_guess)


def _load_detect_positions(fname: str, cache_dir: Path, result_dir: Path):
    """Load x, y (and det_uid if present) from detect_{fname}.csv.

    Returns a DataFrame with columns [det_uid, x, y], or None if not found.
    """
    candidates = [
        cache_dir / f"detect_{fname}.csv",
        step4_dir(result_dir) / f"detect_{fname}.csv",
        result_dir / f"detect_{fname}.csv",
    ]
    for p in candidates:
        if not p.exists() or p.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(p)
            x_col = "x" if "x" in df.columns else ("xcenter" if "xcenter" in df.columns else None)
            y_col = "y" if "y" in df.columns else ("ycenter" if "ycenter" in df.columns else None)
            if x_col is None or y_col is None:
                continue
            out = pd.DataFrame({
                "x": pd.to_numeric(df[x_col], errors="coerce"),
                "y": pd.to_numeric(df[y_col], errors="coerce"),
            })
            if "det_uid" in df.columns:
                out["det_uid"] = (
                    pd.to_numeric(df["det_uid"], errors="coerce")
                    .fillna(pd.Series(range(len(df))))
                    .astype(int)
                )
            else:
                out["det_uid"] = range(len(df))
            out = out.dropna(subset=["x", "y"])
            return out[["det_uid", "x", "y"]]
        except Exception:
            continue
    return None


# ── ApertureWorker ─────────────────────────────────────────────────────────────

class ApertureWorker(QThread):
    """Per-frame aperture sizing and aperture correction via growth curve analysis.

    The ``apcorr`` column written to *apcorr_summary.csv* is a **flux ratio**
    (ref_aperture_flux / optimal_aperture_flux).  Downstream consumers apply it as::

        flux_corrected = flux_measured * apcorr
    """

    progress = pyqtSignal(int, int, str)
    file_done = pyqtSignal(str, dict)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)
    log = pyqtSignal(str)

    def __init__(
        self,
        file_list,
        params,
        data_dir,
        result_dir,
        cache_dir,
        use_cropped=False,
        output_dir=None,
    ):
        super().__init__()
        self.file_list = list(file_list)
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.use_cropped = use_cropped
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self._stop_requested = False

    def stop(self):
        self._stop_requested = True

    def _resolve_fits_path(self, fname: str):
        if self.use_cropped and crop_is_active(self.result_dir):
            cdir = step2_cropped_dir(self.result_dir)
            p = cdir / fname
            if p.exists():
                return p
        try:
            p = Path(self.params.get_file_path(fname))
        except Exception:
            p = self.data_dir / fname
        return p if p.exists() else None

    def run(self):  # noqa: C901
        try:
            P = self.params.P
            output_dir = (
                self.output_dir
                if self.output_dir is not None
                else step5_aperture_dir(self.result_dir)
            )
            output_dir.mkdir(parents=True, exist_ok=True)

            ap_scale          = _to_float(getattr(P, "phot_aperture_scale",    1.0),   1.0)
            ann_in_scale      = _to_float(getattr(P, "fitsky_annulus_scale",   4.0),   4.0)
            ann_out_scale     = _to_float(getattr(P, "fitsky_dannulus_scale",  2.0),   2.0)
            fwhm_px_min       = _to_float(getattr(P, "fwhm_px_min",            3.5),   3.5)
            fwhm_px_max       = _to_float(getattr(P, "fwhm_px_max",            8.0),   8.0)
            min_r_ap_px       = _to_float(getattr(P, "min_r_ap_px",            4.0),   4.0)
            ann_gap           = _to_float(getattr(P, "annulus_min_gap_px",     6.0),   6.0)
            apcorr_apply      = bool(getattr(P, "apcorr_apply",                True))
            apcorr_use_min_n  = _to_int(getattr(P,  "apcorr_use_min_n",        20),    20)
            apcorr_scatter_max = _to_float(getattr(P, "apcorr_scatter_max",    0.05),  0.05)
            ann_sigma         = _to_float(getattr(P, "annulus_sigma_clip",     3.0),   3.0)
            ann_maxiter       = _to_int(getattr(P,  "fitsky_max_iter",          5),    5)
            gc_scale_min      = max(0.3, _to_float(getattr(P, "apcorr_scale_min", 0.5), 0.5))
            gc_scale_max      = max(gc_scale_min + 0.5, _to_float(getattr(P, "apcorr_scale_max", 5.0), 5.0))
            gc_scale_step     = max(0.1, _to_float(getattr(P, "apcorr_scale_step", 0.25), 0.25))
            gc_large_ref_scale = _to_float(getattr(P, "apcorr_large_ref_scale", 5.0), 5.0)
            gc_isolation_factor = max(1.0, _to_float(getattr(P, "apcorr_isolation_factor", 2.0), 2.0))
            apcorr_max_sources = max(30, _to_int(getattr(P, "apcorr_max_sources", 250), 250))
            gain              = _to_float(getattr(P, "gain_e_per_adu",         1.0),   1.0)
            rn_e              = _to_float(getattr(P, "rdnoise_e",              7.5),   7.5)
            ps                = _to_float(getattr(P, "pixel_scale_arcsec",     np.nan), np.nan)
            fwhm_guess        = _to_float(getattr(P, "fwhm_pix_guess",         6.0),   6.0)

            gc_scales = _build_scale_grid(gc_scale_min, gc_scale_max, gc_scale_step)

            # Optional QC filtering
            phot_use_qc_pass_only = bool(getattr(P, "phot_use_qc_pass_only", False))
            files = list(self.file_list)
            if phot_use_qc_pass_only:
                qpath = resolve_frame_quality_path(self.result_dir)
                if qpath is not None and qpath.exists():
                    try:
                        dfq = pd.read_csv(qpath)
                        good = set(dfq.loc[dfq["passed"] == True, "file"].astype(str).tolist())
                        files = [f for f in files if f in good]
                    except Exception:
                        pass

            rows_ap, rows_apcorr, rows_gc = [], [], []
            total = len(files)

            for i, fname in enumerate(files, 1):
                if self._stop_requested:
                    break

                fwhm_med = float(
                    _load_fwhm_from_meta(fname, self.cache_dir, self.result_dir, fwhm_guess)
                )
                fwhm_used = float(np.clip(fwhm_med, fwhm_px_min, fwhm_px_max))

                radii_px = np.array([max(s * fwhm_used, min_r_ap_px) for s in gc_scales])
                _, unique_idx = np.unique(np.round(radii_px, 4), return_index=True)
                unique_idx = np.sort(unique_idx)
                radii_px = radii_px[unique_idx]
                scales_used = np.array(gc_scales)[unique_idx]

                r_large_ref = max(
                    gc_large_ref_scale * fwhm_used,
                    radii_px[-1] if len(radii_px) else min_r_ap_px,
                )
                r_in = max(ann_in_scale * fwhm_used, r_large_ref + ann_gap)
                r_out = r_in + ann_out_scale * fwhm_used
                r_ap_default = max(ap_scale * fwhm_used, min_r_ap_px)

                apc_row = dict(
                    file=fname, fwhm_med=fwhm_med, fwhm_used=fwhm_used,
                    optimal_scale=ap_scale, r_optimal=r_ap_default,
                    r_large_ref=r_large_ref, n_used=0, apcorr=np.nan,
                    mag_err_optimal=np.nan, snr_optimal=np.nan, apply=False,
                    apply_reason="apcorr_not_run" if apcorr_apply else "apcorr_disabled",
                )

                if apcorr_apply:
                    det_df = _load_detect_positions(fname, self.cache_dir, self.result_dir)
                    img_path = self._resolve_fits_path(fname)
                    if det_df is not None and len(det_df) and img_path is not None:
                        try:
                            img = fits.getdata(img_path).astype(float)
                            h, w = img.shape
                            xy_all = det_df[["x", "y"]].to_numpy(float)
                            valid = np.isfinite(xy_all[:, 0]) & np.isfinite(xy_all[:, 1])
                            xy_all = xy_all[valid]

                            if len(xy_all):
                                vals = img[
                                    xy_all[:, 1].astype(int).clip(0, h - 1),
                                    xy_all[:, 0].astype(int).clip(0, w - 1),
                                ]
                                xy_all = xy_all[np.argsort(vals)[::-1]][:apcorr_max_sources]

                            if len(xy_all):
                                edge_pad = int(np.ceil(max(r_large_ref, r_out) + 2.0))
                                edge_ok = (
                                    (xy_all[:, 0] >= edge_pad)
                                    & (xy_all[:, 0] <= w - 1 - edge_pad)
                                    & (xy_all[:, 1] >= edge_pad)
                                    & (xy_all[:, 1] <= h - 1 - edge_pad)
                                )
                                xy_all = xy_all[edge_ok]

                            if len(xy_all) >= 3:
                                tree = cKDTree(xy_all)
                                dists, _ = tree.query(xy_all, k=min(2, len(xy_all)), workers=1)
                                nn_dists = dists[:, 1] if dists.ndim > 1 else dists
                                xy_iso = xy_all[nn_dists > gc_isolation_factor * r_large_ref]
                                if len(xy_iso) < 3:
                                    xy_iso = xy_all

                                # Growth curve
                                gc_mags, gc_errs, gc_snrs, gc_n = [], [], [], []
                                for r_ap in radii_px:
                                    mags, errs, snrs = [], [], []
                                    for xc, yc in xy_iso:
                                        cut_r = int(np.ceil(r_out + 2))
                                        xi, yi = int(round(xc)), int(round(yc))
                                        x0 = max(0, xi - cut_r); x1 = min(w, xi + cut_r + 1)
                                        y0 = max(0, yi - cut_r); y1 = min(h, yi + cut_r + 1)
                                        try:
                                            fe, se, snr, *_ = _phot_one_star(
                                                img[y0:y1, x0:x1], xc - x0, yc - y0,
                                                r_ap, r_in, r_out,
                                                sigma_clip_val=ann_sigma, maxiters=ann_maxiter,
                                                gain=gain, rn_param_e=rn_e,
                                            )
                                            if np.isfinite(fe) and fe > 0 and np.isfinite(se) and se > 0:
                                                mags.append(-2.5 * np.log10(fe))
                                                errs.append(2.5 / np.log(10) * se / fe)
                                                snrs.append(snr)
                                        except Exception:
                                            pass
                                    gc_mags.append(float(np.nanmedian(mags)) if mags else np.nan)
                                    gc_errs.append(float(np.nanmedian(errs)) if errs else np.nan)
                                    gc_snrs.append(float(np.nanmedian(snrs)) if snrs else np.nan)
                                    gc_n.append(len(mags))

                                # Find optimal aperture (minimum mag_err)
                                gc_errs_arr = np.array(gc_errs, dtype=float)
                                valid_err = np.isfinite(gc_errs_arr)
                                if np.any(valid_err):
                                    valid_pos = np.where(valid_err)[0]
                                    opt_pos = valid_pos[int(np.argmin(gc_errs_arr[valid_err]))]
                                    r_optimal    = float(radii_px[opt_pos])
                                    sc_optimal   = float(scales_used[opt_pos])
                                    err_optimal  = float(gc_errs_arr[opt_pos])
                                    snr_optimal  = float(gc_snrs[opt_pos]) if np.isfinite(gc_snrs[opt_pos]) else np.nan

                                    # Measure optimal fluxes for apcorr ratio
                                    opt_fluxes = []
                                    for xc, yc in xy_iso:
                                        cut_r = int(np.ceil(r_out + 2))
                                        xi, yi = int(round(xc)), int(round(yc))
                                        x0 = max(0, xi - cut_r); x1 = min(w, xi + cut_r + 1)
                                        y0 = max(0, yi - cut_r); y1 = min(h, yi + cut_r + 1)
                                        try:
                                            fe, *_ = _phot_one_star(
                                                img[y0:y1, x0:x1], xc - x0, yc - y0,
                                                r_optimal, r_in, r_out,
                                                sigma_clip_val=ann_sigma, maxiters=ann_maxiter,
                                                gain=gain, rn_param_e=rn_e,
                                            )
                                            opt_fluxes.append(fe if (np.isfinite(fe) and fe > 0) else np.nan)
                                        except Exception:
                                            opt_fluxes.append(np.nan)
                                    opt_arr = np.array(opt_fluxes, dtype=float)

                                    # Sweep candidate large-ref scales; pick minimum rel_scatter.
                                    cand_scales_raw = list(np.arange(1.5, max(gc_large_ref_scale, 2.0) + 0.01, 0.5))
                                    if gc_large_ref_scale not in cand_scales_raw:
                                        cand_scales_raw.append(gc_large_ref_scale)
                                    cand_scales = sorted(
                                        s for s in cand_scales_raw
                                        if s * fwhm_used > r_optimal * 1.1 + 1.0
                                    )
                                    if not cand_scales:
                                        cand_scales = [gc_large_ref_scale]

                                    apcorr = np.nan
                                    rel_sc = np.nan
                                    best_rel_sc = np.inf
                                    best_r_ref = gc_large_ref_scale * fwhm_used
                                    n_used = 0

                                    for cand_scale in cand_scales:
                                        r_ref = max(cand_scale * fwhm_used, r_optimal + 2.0)
                                        ref_fluxes = []
                                        for xc, yc in xy_iso:
                                            cut_r = int(np.ceil(r_out + 2))
                                            xi, yi = int(round(xc)), int(round(yc))
                                            x0 = max(0, xi - cut_r); x1 = min(w, xi + cut_r + 1)
                                            y0 = max(0, yi - cut_r); y1 = min(h, yi + cut_r + 1)
                                            try:
                                                fe, *_ = _phot_one_star(
                                                    img[y0:y1, x0:x1], xc - x0, yc - y0,
                                                    r_ref, r_in, r_out,
                                                    sigma_clip_val=ann_sigma, maxiters=ann_maxiter,
                                                    gain=gain, rn_param_e=rn_e,
                                                )
                                                ref_fluxes.append(
                                                    fe if (np.isfinite(fe) and fe > 0) else np.nan
                                                )
                                            except Exception:
                                                ref_fluxes.append(np.nan)

                                        ref_arr = np.array(ref_fluxes[: len(opt_arr)], dtype=float)
                                        valid_c = (
                                            np.isfinite(ref_arr) & np.isfinite(opt_arr)
                                            & (ref_arr > 0) & (opt_arr > 0)
                                        )
                                        if np.sum(valid_c) < 3:
                                            continue
                                        ratios = ref_arr[valid_c] / opt_arr[valid_c]
                                        ratios = ratios[(ratios > 0.5) & (ratios < 10.0)]
                                        if len(ratios) < 3:
                                            continue
                                        r_med = float(np.nanmedian(ratios))
                                        r_sig = 1.4826 * float(np.nanmedian(np.abs(ratios - r_med)))
                                        if np.isfinite(r_sig) and r_sig > 0:
                                            keep = np.abs(ratios - r_med) <= 3.0 * r_sig
                                            if np.sum(keep) >= 3:
                                                ratios = ratios[keep]
                                        apcorr_c = float(np.nanmedian(ratios))
                                        mad_c = float(np.nanmedian(np.abs(ratios - apcorr_c)))
                                        rel_sc_c = 1.4826 * mad_c / apcorr_c if apcorr_c > 0 else np.nan

                                        if np.isfinite(rel_sc_c) and rel_sc_c < best_rel_sc:
                                            best_rel_sc = rel_sc_c
                                            apcorr = apcorr_c
                                            rel_sc = rel_sc_c
                                            best_r_ref = r_ref
                                            n_used = len(ratios)

                                    apply_reasons = []
                                    if n_used < apcorr_use_min_n:
                                        apply_reasons.append(f"n_used<{apcorr_use_min_n}")
                                    if not np.isfinite(apcorr):
                                        apply_reasons.append("apcorr_nan")
                                    elif not (0.8 <= apcorr <= 5.0):
                                        apply_reasons.append(f"apcorr_out_of_range({apcorr:.3f})")
                                    if not np.isfinite(rel_sc):
                                        apply_reasons.append("rel_scatter_nan")
                                    elif rel_sc > apcorr_scatter_max:
                                        apply_reasons.append(
                                            f"rel_scatter>{apcorr_scatter_max:.3f}({rel_sc:.3f})"
                                        )
                                    apply_flag = len(apply_reasons) == 0
                                    apply_reason = "ok" if apply_flag else "|".join(apply_reasons)
                                    apc_row.update(
                                        optimal_scale=sc_optimal, r_optimal=r_optimal,
                                        r_large_ref=best_r_ref, n_used=n_used, apcorr=apcorr,
                                        rel_scatter=rel_sc, mag_err_optimal=err_optimal,
                                        snr_optimal=snr_optimal, apply=apply_flag,
                                        apply_reason=apply_reason,
                                    )
                                    self.log.emit(
                                        f"[APCORR] {fname} apply={apply_flag} reason={apply_reason} | "
                                        f"n={n_used} apcorr={apcorr:.4f}(flux ratio) "
                                        f"rel_scatter={rel_sc:.4f} "
                                        f"r_ref={best_r_ref:.1f}px"
                                        f"(×{best_r_ref/fwhm_used:.1f}FWHM)"
                                    )
                                    for k, (r_px, scale, m_mag, m_err, m_snr, n_s) in enumerate(
                                        zip(radii_px, scales_used, gc_mags, gc_errs, gc_snrs, gc_n)
                                    ):
                                        rows_gc.append(dict(
                                            file=fname, r_px=float(r_px), scale=float(scale),
                                            median_mag=m_mag, median_mag_err=m_err,
                                            median_snr=m_snr, n_stars=n_s,
                                            selected=(k == opt_pos),
                                        ))
                        except Exception as exc:
                            self.log.emit(f"[Apcorr] {fname}: {exc}")

                r_ap_out = float(apc_row.get("r_optimal", r_ap_default))
                r_in_out = max(ann_in_scale * fwhm_used, r_ap_out + ann_gap)
                r_out_out = r_in_out + ann_out_scale * fwhm_used
                row = dict(
                    file=fname,
                    fwhm_med=fwhm_med,
                    fwhm_used=fwhm_used,
                    r_ap=r_ap_out,
                    r_in=r_in_out,
                    r_out=r_out_out,
                )
                if np.isfinite(ps) and ps > 0:
                    row.update(
                        fwhm_med_arcsec=fwhm_med * ps,
                        fwhm_used_arcsec=fwhm_used * ps,
                        r_ap_arcsec=r_ap_out * ps,
                    )
                rows_ap.append(row)
                rows_apcorr.append(apc_row)
                try:
                    (self.cache_dir / f"apcorr_{fname}.json").write_text(
                        json.dumps(apc_row, indent=2), encoding="utf-8"
                    )
                except Exception:
                    pass
                self.file_done.emit(fname, row)
                self.progress.emit(i, total, fname)

            pd.DataFrame(rows_ap).to_csv(output_dir / "aperture_by_frame.csv", index=False)
            pd.DataFrame(rows_apcorr).to_csv(output_dir / "apcorr_summary.csv", index=False)
            if rows_gc:
                pd.DataFrame(rows_gc).to_csv(output_dir / "growth_curve.csv", index=False)

            self.finished.emit({"total": len(rows_ap), "frames": len(rows_ap)})
        except Exception as e:
            self.error.emit("WORKER", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")
            self.finished.emit({})
