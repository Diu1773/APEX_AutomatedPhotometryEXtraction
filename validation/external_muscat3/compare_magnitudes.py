"""다른 기기로 잰 같은 성단의 등급이 얼마나 다른가 — 정량 비교.

**이것이 목표 둘의 완료 조건이다** (STATUS.md). 「다른 기기에서 돈다」가 아니라
**「다른 기기에서 잰 등급이 얼마나 다른가」**를 수치로 낸다.

## 무엇과 무엇을 견주나

    A. MuSCAT3 (LCO 2 m, ogg)  vs  Moravian C3-61000 (사용자 망원경)
       같은 성단 M67, 같은 필터 g·r·i, 둘 다 APEX 로 처리.
       **기기 이식성을 직접 재는 축이다.**

    B. 각각  vs  Pan-STARRS1 (PS1 DR2, VizieR II/349/ps1)
       바깥 기준. 어느 쪽이 어긋났는지 가른다 — A 만 보면 둘 중 누가
       틀렸는지 알 수 없다.

## 별을 어떻게 잇나

**`gaia_source_id` 로 잇는다.** 워크스페이스가 다르면 `source_id` 는 서로 다른
뜻이므로 쓰면 안 된다(`Main/FAILURES.md` F-285). Gaia 번호가 없는 별은 하늘
좌표로 1 초각 안에서 잇는다.

## 무엇을 내나

밴드마다 이것들을 낸다. 절대값과 흩어짐을 **함께** 낸다 — 「상관 제거」는 전부
균일하게 나쁘게 만들어도 달성되므로 절대값을 늘 같이 봐야 한다.

    N          이어진 별 수
    중앙 차이  Δ = 이쪽 − 저쪽 의 중앙값 (계통 어긋남)
    산포       Δ 의 MAD (별마다의 흩어짐)
    밝기 의존  등급 구간별 Δ — 어두운 쪽에서 휘는지
    색 의존    (g−i) 에 대한 Δ 의 기울기 — 색항이 남았는지

실행:
    python -X utf8 validation/external_muscat3/compare_magnitudes.py \
        --muscat3 <job>/results --moravian E:/APEX_validation/reprocess/M67/result
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).absolute().parents[2]))

from apex.utils.io_utils import (  # noqa: E402
    normalize_id_columns,
    read_csv_int64_source_id,
)

BANDS = ("g", "r", "i")


def mad(v) -> float:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def load_workspace(result_dir: Path, label: str) -> pd.DataFrame:
    """등급 + 하늘 좌표 + Gaia 번호를 한 표로."""
    zp = result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide.csv"
    if not zp.exists():
        raise SystemExit(f"{label}: 등급 표가 없다 — {zp}")
    mags = read_csv_int64_source_id(zp)

    cat = None
    for name in ("ref_catalog.tsv", "master_catalog.tsv"):
        p = result_dir / "step6_refbuild" / name
        if p.exists():
            cat = read_csv_int64_source_id(p, sep="\t")
            break
    if cat is None:
        raise SystemExit(f"{label}: 마스터 목록이 없다 — {result_dir/'step6_refbuild'}")

    keep = [c for c in ("ID", "source_id", "gaia_source_id", "ra_deg", "dec_deg")
            if c in cat.columns]
    out = mags.merge(cat[keep], on="ID", how="left")
    normalize_id_columns(out)
    have = [b for b in BANDS if f"mag_cal_{b}" in out.columns]
    print(f"[{label}] {len(out)} 별 · 밴드 {have} · "
          f"Gaia 번호 있는 별 {int(out.get('gaia_source_id', pd.Series(dtype='Int64')).notna().sum())}")
    return out


def join(a: pd.DataFrame, b: pd.DataFrame, tol_arcsec: float = 1.0) -> pd.DataFrame:
    """Gaia 번호 먼저, 없으면 위치로."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ga = a["gaia_source_id"] if "gaia_source_id" in a.columns else pd.Series(dtype="Int64")
    gb = b["gaia_source_id"] if "gaia_source_id" in b.columns else pd.Series(dtype="Int64")
    pairs: list[tuple[int, int]] = []
    used_b: set[int] = set()
    if len(ga) and len(gb):
        idx_b = {int(v): i for i, v in enumerate(gb) if pd.notna(v) and int(v) > 0}
        for i, v in enumerate(ga):
            if pd.notna(v) and int(v) in idx_b:
                j = idx_b[int(v)]
                pairs.append((i, j))
                used_b.add(j)
    n_by_id = len(pairs)

    left = [i for i in range(len(a)) if i not in {p[0] for p in pairs}]
    right = [j for j in range(len(b)) if j not in used_b]
    if left and right:
        ca = SkyCoord(a["ra_deg"].to_numpy(float)[left],
                      a["dec_deg"].to_numpy(float)[left], unit="deg")
        cb = SkyCoord(b["ra_deg"].to_numpy(float)[right],
                      b["dec_deg"].to_numpy(float)[right], unit="deg")
        k, sep, _ = ca.match_to_catalog_sky(cb)
        for n, (kk, ss) in enumerate(zip(k, sep.arcsec)):
            if ss <= tol_arcsec and right[int(kk)] not in used_b:
                pairs.append((left[n], right[int(kk)]))
                used_b.add(right[int(kk)])
    print(f"이어진 별 {len(pairs)} (Gaia 번호 {n_by_id} · 위치 {len(pairs)-n_by_id})")

    ia = [p[0] for p in pairs]
    ib = [p[1] for p in pairs]
    out = a.iloc[ia].reset_index(drop=True).add_suffix("_A")
    out = pd.concat([out, b.iloc[ib].reset_index(drop=True).add_suffix("_B")], axis=1)
    return out


def compare(joined: pd.DataFrame, name_a: str, name_b: str) -> list[dict]:
    rows = []
    print()
    print(f"=== {name_a}  −  {name_b} ===")
    print(f"{'밴드':<5}{'N':>6}{'중앙 차이':>11}{'산포(MAD)':>11}"
          f"{'밝은쪽':>10}{'어두운쪽':>10}{'색 기울기':>11}")
    print("-" * 66)
    gi = None
    if all(f"mag_cal_{b}_A" in joined.columns for b in ("g", "i")):
        gi = joined["mag_cal_g_A"] - joined["mag_cal_i_A"]
    for band in BANDS:
        ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
        if ca not in joined.columns or cb not in joined.columns:
            continue
        d = joined[ca] - joined[cb]
        m = np.isfinite(d) & np.isfinite(joined[ca])
        if m.sum() < 5:
            continue
        dv, mv = d[m].to_numpy(float), joined[ca][m].to_numpy(float)
        cut = np.nanmedian(mv)
        bright, faint = dv[mv <= cut], dv[mv > cut]
        slope = float("nan")
        if gi is not None:
            gm = m & np.isfinite(gi)
            if gm.sum() > 10:
                slope = float(np.polyfit(gi[gm].to_numpy(float),
                                         d[gm].to_numpy(float), 1)[0])
        row = dict(band=band, n=int(m.sum()),
                   median=float(np.median(dv)), mad=mad(dv),
                   median_bright=float(np.median(bright)) if bright.size else float("nan"),
                   median_faint=float(np.median(faint)) if faint.size else float("nan"),
                   mag_split=float(cut), colour_slope=slope)
        rows.append(row)
        print(f"{band:<5}{row['n']:>6}{row['median']:>+11.4f}{row['mad']:>11.4f}"
              f"{row['median_bright']:>+10.4f}{row['median_faint']:>+10.4f}"
              f"{slope:>+11.4f}")
    return rows


def query_ps1(ra, dec, radius_deg: float) -> pd.DataFrame | None:
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        from astroquery.vizier import Vizier
    except Exception as exc:  # noqa: BLE001
        print(f"PS1 조회 불가 (astroquery 없음): {exc}")
        return None
    v = Vizier(columns=["RAJ2000", "DEJ2000", "gmag", "rmag", "imag",
                        "e_gmag", "e_rmag", "e_imag"],
               column_filters={"gmag": "<21"}, row_limit=-1)
    try:
        res = v.query_region(SkyCoord(ra, dec, unit="deg"),
                             radius=radius_deg * u.deg, catalog="II/349/ps1")
    except Exception as exc:  # noqa: BLE001
        print(f"PS1 조회 실패: {exc}")
        return None
    if not res:
        print("PS1 결과 없음")
        return None
    t = res[0].to_pandas().dropna(subset=["gmag", "rmag", "imag"])
    print(f"PS1 {len(t)} 줄")
    return t.rename(columns={"RAJ2000": "ra_deg", "DEJ2000": "dec_deg",
                             "gmag": "mag_cal_g", "rmag": "mag_cal_r",
                             "imag": "mag_cal_i"})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="다른 기기로 잰 같은 성단의 등급 비교")
    ap.add_argument("--muscat3", required=True, help="MuSCAT3 결과 폴더")
    ap.add_argument("--moravian", default="", help="Moravian 결과 폴더 (선택)")
    ap.add_argument("--out", default="", help="요약 JSON 경로")
    ap.add_argument("--skip-ps1", action="store_true")
    a = ap.parse_args(argv)

    summary: dict = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "comparisons": []}
    m3 = load_workspace(Path(a.muscat3), "MuSCAT3")

    if a.moravian:
        mv = load_workspace(Path(a.moravian), "Moravian")
        j = join(m3, mv)
        summary["comparisons"].append(
            {"pair": "MuSCAT3 - Moravian", "n_joined": len(j),
             "bands": compare(j, "MuSCAT3", "Moravian C3-61000")})

    if not a.skip_ps1:
        ra0, dec0 = float(np.nanmedian(m3["ra_deg"])), float(np.nanmedian(m3["dec_deg"]))
        span = float(np.nanmax(np.hypot(
            (m3["ra_deg"] - ra0) * np.cos(np.deg2rad(dec0)), m3["dec_deg"] - dec0)))
        ps1 = query_ps1(ra0, dec0, span + 0.02)
        if ps1 is not None:
            ps1["gaia_source_id"] = pd.NA
            summary["comparisons"].append(
                {"pair": "MuSCAT3 - PS1", "n_joined": None,
                 "bands": compare(join(m3, ps1), "MuSCAT3", "PS1 DR2")})
            if a.moravian:
                summary["comparisons"].append(
                    {"pair": "Moravian - PS1", "n_joined": None,
                     "bands": compare(join(mv, ps1), "Moravian", "PS1 DR2")})

    out = Path(a.out) if a.out else Path(a.muscat3) / "magnitude_comparison.json"
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
