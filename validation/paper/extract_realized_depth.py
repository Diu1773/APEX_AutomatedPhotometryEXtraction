"""Per-frame REALIZED detection depth from real stars (the ATLAS/Kessler-style half
of the QC-depth figure; the PREDICTED half comes from the predict-m50 QC gate).

For every step7 frame of every reprocessed target: take the master-catalog stars on
the frame, assign each star its median count-rate instrumental magnitude over ALL
frames of the same filter (exposure-invariant, and independent of this frame's own
measurement — kills the Eddington bias that distorts single-frame binning), convert
to the injection scale (total electrons: m_inj = mag_rate - 2.5 log10 t), bin the
frame's detected_flag against it, and read off the 50% crossing.

Output: validation/paper/data_qc_depth/realized_m50.csv
    .venv-deploy\\Scripts\\python.exe validation\\paper\\extract_realized_depth.py
"""
from __future__ import annotations
import sys
from pathlib import Path
REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

import numpy as np
import pandas as pd

REPRO = Path(r"E:\APEX_validation\reprocess")
TARGETS = ["M13", "M67", "NGC6811"]
OUT = REPO / "validation" / "paper" / "data_qc_depth"
OUT.mkdir(exist_ok=True)


def read_off_50(x, frac):
    o = np.argsort(x); x, frac = x[o], frac[o]
    for i in range(len(x) - 1):
        if frac[i] >= 0.5 >= frac[i + 1]:
            d = frac[i] - frac[i + 1]
            f = (frac[i] - 0.5) / d if d else 0.0
            return float(x[i] + f * (x[i + 1] - x[i]))
    return float("nan")


def main() -> int:
    rows = []
    for tgt in TARGETS:
        base = REPRO / tgt / "result" / "step7_forced_phot"
        fs = pd.read_csv(base / "frame_stats.csv")
        fs = fs[fs["status"] == "ok"]
        # load all frames per filter once, to build per-star median rate-mags
        per_filter: dict[str, list[tuple[str, float, float, pd.DataFrame]]] = {}
        for _, fr in fs.iterrows():
            tsv = base / f"photometry_{fr['file']}.tsv"
            if not tsv.exists():
                continue
            d = pd.read_csv(tsv, sep="\t",
                            usecols=["master_id", "mag_inst", "detected_flag",
                                     "off_frame_flag", "exptime"])
            d = d[(d.off_frame_flag == False) & np.isfinite(d.mag_inst)]
            per_filter.setdefault(fr["filter"], []).append(
                (fr["file"], float(d.exptime.iloc[0]), float(fr["fwhm_px"]), d))
        for filt, frames in per_filter.items():
            med = pd.concat(
                [f[3].set_index("master_id")["mag_inst"] for f in frames], axis=1
            ).median(axis=1)   # per-star median count-rate mag (exposure-invariant)
            # master faint limit in this filter (count-rate mags). A frame whose
            # 50% crossing approaches this limit is testing stars that only exist
            # in the master BECAUSE frames like it detected them → circular, flag.
            mlim_rate = float(np.nanpercentile(med.to_numpy(float), 90))
            for fname, expt, fwhm, d in frames:
                m_true = med.reindex(d.master_id).to_numpy(float)
                ok = np.isfinite(m_true)
                m_inj = m_true[ok] - 2.5 * np.log10(expt)   # injection scale (total e-)
                det = d.detected_flag.to_numpy(bool)[ok]
                bw = 0.3
                e = np.arange(np.floor(m_inj.min() / bw) * bw, m_inj.max() + bw, bw)
                xs, cs = [], []
                for lo, hi in zip(e[:-1], e[1:]):
                    k = (m_inj >= lo) & (m_inj < hi)
                    if k.sum() >= 12:
                        xs.append(0.5 * (lo + hi)); cs.append(det[k].mean())
                m50 = read_off_50(np.array(xs), np.array(cs)) if len(xs) > 2 else float("nan")
                mlim_inj = mlim_rate - 2.5 * np.log10(expt)
                margin = mlim_inj - m50 if np.isfinite(m50) else float("nan")
                valid = bool(np.isfinite(m50) and margin > 0.7)
                rows.append(dict(target=tgt, file=fname, filter=filt, exptime=expt,
                                 fwhm_px=fwhm, n_stars=int(ok.sum()),
                                 n_detected=int(det.sum()), realized_m50=m50,
                                 master_limit_inj=round(mlim_inj, 2),
                                 margin_mag=round(margin, 2) if np.isfinite(margin) else np.nan,
                                 depth_valid=valid))
                flag = "" if valid else "  ⚠ near master limit (circular)"
                print(f"[{tgt:8s}] {fname:34s} {filt:2s} t={expt:5.0f}s "
                      f"n={ok.sum():4d} det={det.sum():4d} m50={m50:6.2f}{flag}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "realized_m50.csv", index=False)
    good = out[np.isfinite(out.realized_m50)]
    print(f"\nframes: {len(out)}  with readable m50: {len(good)}  "
          f"range {good.realized_m50.min():.2f}–{good.realized_m50.max():.2f}")
    print("wrote", OUT / "realized_m50.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
