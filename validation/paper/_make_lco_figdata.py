"""Prepare Fig 13 data: reduce the two public LCO datasets with APEX-style
calibration and save downsampled (APEX, BANZAI, difference) arrays + agreement
stats to data/ (git-ignored, regenerable). Needs the raw LCO frames under
E:\\APEX_validation\\external\\ (downloaded from archive.lco.global).

Run: .venv-deploy\\Scripts\\python validation\\paper\\_make_lco_figdata.py
"""
from __future__ import annotations
import json, re
from pathlib import Path
import numpy as np
from astropy.io import fits

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
EXT = Path(r"E:\APEX_validation\external")
DATA = REPO / "validation" / "paper" / "data"
DATA.mkdir(parents=True, exist_ok=True)


def _sci(fz):
    with fits.open(fz) as h:
        for hd in h:
            if hd.name == "SCI" and hd.data is not None:
                return np.asarray(hd.data, float), hd.header, h[0].header
        return np.asarray(h[1].data, float), h[1].header, h[0].header


def _block(a, n):
    """Block-mean downsample to <= n on the long axis (for display only)."""
    f = max(1, a.shape[0] // n)
    h, w = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    return a[:h, :w].reshape(h // f, f, w // f, f).mean(axis=(1, 3))


def _rstat(d):
    d = d[np.isfinite(d)]; m = float(np.median(d)); s = float(1.4826 * np.median(np.abs(d - m)))
    return m, s


def sinistro():
    L = EXT / "LCO_sinistro"
    def sec(s): return tuple(int(x) for x in re.findall(r"-?\d+", s))
    with fits.open(L / "raw_e00.fits.fz") as h:
        mos = np.zeros((4096, 4096))
        for i in range(1, 5):
            hd = h[i].header; d = h[i].data.astype(float)
            x1, x2, y1, y2 = sec(hd["DATASEC"]); bx1, bx2, by1, by2 = sec(hd["BIASSEC"])
            over = np.median(d[by1-1:by2, bx1-1:bx2])
            sci = (d[y1-1:y2, x1-1:x2] - over) * float(hd["GAIN"])
            dx1, dx2, dyy1, dyy2 = sec(hd["DETSEC"])
            if dx1 > dx2: sci = np.fliplr(sci); dx1, dx2 = dx2, dx1
            if dyy1 > dyy2: sci = np.flipud(sci); dyy1, dyy2 = dyy2, dyy1
            mos[dyy1-1:dyy2, dx1-1:dx2] = sci
        lexp = float(h[0].header["EXPTIME"])
    mb, *_ = _sci(L/"master_bias.fits.fz"); md, mdh, _ = _sci(L/"master_dark.fits.fz")
    mf, *_ = _sci(L/"master_flat.fits.fz"); e91, *_ = _sci(L/"banzai_e91.fits.fz")
    apex = (mos - mb - md*(lexp/float(mdh["EXPTIME"]))) / np.where(np.abs(mf) < 1e-6, np.nan, mf)
    return apex, e91


def qhy():
    Q = EXT / "LCO_qhy600_0m4"
    raw, rh, _ = _sci(Q/"raw_e00.fits.fz"); mb, *_ = _sci(Q/"master_bias.fits.fz")
    md, mdh, _ = _sci(Q/"master_dark.fits.fz"); mf, *_ = _sci(Q/"master_flat.fits.fz")
    e91, e91h, _ = _sci(Q/"banzai_e91.fits.fz")
    bl = float(e91h["BIASLVL"]); ratio = float(e91h["EXPTIME"])/float(mdh["EXPTIME"])
    apex = (raw*float(rh["GAIN"]) - bl - mb - md*ratio) / np.where(np.abs(mf) < 1e-6, np.nan, mf)
    return apex, e91


stats = {}
for name, fn, cam in [("qhy", qhy, "QHY600 (CMOS, single-amp)"),
                      ("sinistro", sinistro, "Sinistro (CCD, 4-amp)")]:
    apex, e91 = fn()
    d = apex - e91
    m, s = _rstat(d)
    np.save(DATA / f"lco_{name}_apex.npy", _block(apex, 700).astype(np.float32))
    np.save(DATA / f"lco_{name}_e91.npy", _block(e91, 700).astype(np.float32))
    np.save(DATA / f"lco_{name}_diff.npy", _block(d, 700).astype(np.float32))
    stats[name] = {"camera": cam, "delta_median": m, "robust_sigma": s,
                   "sky_median": float(np.nanmedian(e91))}
    print(f"{name}: Δmedian={m:+.4f} robustσ={s:.4f} e-  (sky {np.nanmedian(e91):.2f})")
(DATA / "lco_crossinstrument.json").write_text(json.dumps(stats, indent=2))
print("wrote lco_*.npy + lco_crossinstrument.json to", DATA)
