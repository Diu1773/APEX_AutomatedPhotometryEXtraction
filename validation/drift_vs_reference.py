"""기준을 갈아 끼워도 밝기 치우침이 남는가 — Gaia 를 빼고 Pan-STARRS1 로 잰다.

**여태 기준을 한 번도 안 바꿔 봤다.** kb26 의 밝기 치우침이 기기 탓이라고 말해
온 근거는 「같은 Gaia 기준에 견준 Moravian·MuSCAT3 는 평평하다」였는데, 그것만으로는
기준이 무죄라고 못 한다. 세 기기가 **같은 변환식을 같은 방식으로** 쓰는 것이
아니기 때문이다.

    kb26      r 을 `r−i` 색으로 맞춘다   (gp 를 안 찍어서 g−r 을 못 쓴다)
    Moravian  r 을 `g−r` 색으로 맞춘다
    MuSCAT3   r 을 `g−r` 색으로 맞춘다

같은 Jordi+2010 이라도 **색축이 다르면 다른 변환**이다. 그러니 이 변인을 잡으려면
기준 자체를 갈아야 한다.

## 무엇으로 가나

**Pan-STARRS1 (PS1 DR2)** 이다. Gaia 와 완전히 독립인 측광 목록이고, r·i 는
변환 없이 거의 그대로 견줄 수 있어 **변환식이라는 변인이 통째로 사라진다.**
M67 은 적위 +11.8 이라 PS1 이 덮는다.

## 두 번에 나눠 본다

1 부는 기기 등급에서 PS1 등급을 뺀 값을 PS1 등급으로 나눠 본다. 다만 이것만으로는
아직 못 가른다. PS1 자체가 밝은 별에서 포화해 실제보다 어둡게 적히기 때문에,
**기준의 결함과 기기의 결함이 한 표에 섞여 들어오기 때문이다.**

**2 부가 그것을 가른다.** 같은 별에서 차이를 둘로 쪼갠다.

    기기 − PS1  =  (기기 − Gaia 기준)  +  (Gaia 기준 − PS1)
                     기기·변환 쪽            두 기준 사이

오른쪽 둘째 항은 기기가 무엇이든 **같은 별에서 같은 값**이어야 한다. 그러므로
이 항이 기기마다 다르면 변환 탓이고, 네 기기에서 같은 모양이면 기준 탓이며,
첫째 항만 kb26 에서 크면 기기 탓이다.

실행:
    python -X utf8 validation/drift_vs_reference.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "external_muscat3"))

from compare_magnitudes import load_workspace  # noqa: E402
from apex.utils.io_utils import read_csv_int64_source_id  # noqa: E402

WORKSPACES = {
    "kb26": (REPO / "validation/external_kb26/results", ("r", "i")),
    "kb27": (REPO / "validation/external_kb27/results", ("r", "i")),
    "Moravian": (Path("E:/APEX_validation/reprocess/M67/result"), ("r", "i")),
    "MuSCAT3": (REPO / "validation/external_muscat3/results", ("r", "i")),
}
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0]

#: PS1 이 포화하기 시작하는 등급. 이보다 밝은 별은 PS1 쪽이 실제보다 어둡게
#: 적혀 있어 **어느 기기로 재든 같은 방향으로 어긋난다.** 그래서 「PS1 이 성한
#: 구간」의 폭을 따로 낸다 — 표의 맨 왼쪽 세 칸이 그 영향권이다.
PS1_SAFE_MAG = 14.0
MATCH_ARCSEC = 1.0
MIN_CELL = 8
CACHE = REPO / "validation/ps1_m67.csv"


def ps1_for_m67(ra: float, dec: float, radius_deg: float = 0.45) -> pd.DataFrame:
    """PS1 목록. 한 번 받아 두고 그 다음부터는 파일에서 읽는다."""
    if CACHE.exists():
        out = pd.read_csv(CACHE)
        print(f"PS1 {len(out)} 줄 (저장해 둔 것)")
        return out
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    v = Vizier(columns=["RAJ2000", "DEJ2000", "gmag", "rmag", "imag",
                        "e_gmag", "e_rmag", "e_imag", "Nd"],
               column_filters={"rmag": "<20"}, row_limit=-1)
    res = v.query_region(SkyCoord(ra, dec, unit="deg"),
                         radius=radius_deg * u.deg, catalog="II/349/ps1")
    if not res:
        raise SystemExit("PS1 결과가 없다")
    t = res[0].to_pandas().dropna(subset=["rmag", "imag"])
    t = t.rename(columns={"RAJ2000": "ra_deg", "DEJ2000": "dec_deg",
                          "gmag": "ps1_g", "rmag": "ps1_r", "imag": "ps1_i"})
    t.to_csv(CACHE, index=False)
    print(f"PS1 {len(t)} 줄 (새로 받아 {CACHE.name} 에 저장)")
    return t


def bin_medians(diff: np.ndarray, mag: np.ndarray) -> tuple[list, list]:
    """등급 구간마다 중앙값과 별 수. 별이 적은 칸은 비운다."""
    meds, ns = [], []
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        k = np.isfinite(diff) & np.isfinite(mag) & (mag >= lo) & (mag < hi)
        ns.append(int(k.sum()))
        meds.append(float(np.median(diff[k])) if k.sum() >= MIN_CELL else np.nan)
    return meds, ns


def span_of(meds: list, lo_edge: float = -np.inf) -> float:
    """구간 중앙값의 최대에서 최소를 뺀 값. `lo_edge` 보다 밝은 칸은 뺀다."""
    got = [m for m, lo in zip(meds, MAG_EDGES[:-1])
           if np.isfinite(m) and lo >= lo_edge]
    return float(max(got) - min(got)) if len(got) >= 2 else float("nan")


def cells(meds: list) -> str:
    return " ".join((f"{m:+.3f}".rjust(8) if np.isfinite(m) else "·".rjust(8))
                    for m in meds)


def header(first: str) -> str:
    hdr = "{:<11}{:<5}{:>7}{:>8}{:>8}".format(first, "밴드", "별수", "폭", "14등↓")
    hdr += "   " + " ".join(f"{lo:g}~{hi:g}".rjust(8)
                            for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]))
    return hdr


def line(first: str, band: str, n: int, full: float, safe: float,
         meds: list) -> str:
    return "{:<11}{:<5}{:>7}{:>8.3f}{:>8.3f}   {}".format(
        first, band, n, full, safe, cells(meds))


def gaia_reference(result_dir: Path) -> pd.DataFrame | None:
    """영점 보정이 실제로 쓴 Gaia 변환 기준 등급. 없으면 None."""
    p = Path(result_dir) / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv"
    if not p.exists():
        return None
    tab = read_csv_int64_source_id(p)
    keep = ["ID"] + [c for c in tab.columns if c.startswith("ref_")]
    return tab[keep] if len(keep) > 1 else None


def main() -> int:
    from astropy.coordinates import SkyCoord

    ws = {}
    for name, (rd, bands) in WORKSPACES.items():
        try:
            tab = load_workspace(Path(rd), name)
        except SystemExit:
            print(f"[{name}] 등급 표가 없어 건너뛴다")
            continue
        ref = gaia_reference(Path(rd))
        if ref is None:
            print(f"[{name}] Gaia 기준 등급 표가 없어 2 부에서 빠진다")
        else:
            tab = tab.merge(ref, on="ID", how="left")
        ws[name] = (tab, bands)

    any_tab = next(iter(ws.values()))[0]
    ps1 = ps1_for_m67(float(np.nanmedian(any_tab["ra_deg"])),
                      float(np.nanmedian(any_tab["dec_deg"])))
    cps = SkyCoord(ps1["ra_deg"].to_numpy(float),
                   ps1["dec_deg"].to_numpy(float), unit="deg")

    # 기기마다 PS1 을 한 번만 붙여 두고 1·2 부가 같은 짝을 쓴다.
    joined: dict[str, tuple[pd.DataFrame, tuple, np.ndarray]] = {}
    for name, (tab, bands) in ws.items():
        c = SkyCoord(pd.to_numeric(tab["ra_deg"], errors="coerce").to_numpy(float),
                     pd.to_numeric(tab["dec_deg"], errors="coerce").to_numpy(float),
                     unit="deg")
        k, sep, _ = c.match_to_catalog_sky(cps)
        ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
        joined[name] = (tab, bands, np.where(ok, np.asarray(k), -1))

    rows: list[dict] = []

    print()
    print("=== 1 부. 기준을 PS1 로 갈면 밝기 치우침이 남는가 ===")
    print("기기 등급에서 PS1 등급을 뺀 값을 PS1 등급으로 나눈 구간 중앙값.")
    print(f"「폭」은 전 구간이고, 「14등↓」은 PS1 이 포화하지 않는 "
          f"{PS1_SAFE_MAG:g} 등급보다 어두운 칸만 쓴 폭이다.")
    print()
    hdr = header("기기")
    print(hdr)
    print("-" * len(hdr))

    for name, (tab, bands, kk) in joined.items():
        for band in bands:
            col = f"mag_cal_{band}"
            if col not in tab.columns:
                continue
            mine = pd.to_numeric(tab[col], errors="coerce").to_numpy(float)
            ref = np.full(len(tab), np.nan)
            hit = kk >= 0
            ref[hit] = ps1[f"ps1_{band}"].to_numpy(float)[kk[hit]]
            diff = mine - ref
            meds, ns = bin_medians(diff, ref)
            full, safe = span_of(meds), span_of(meds, PS1_SAFE_MAG)
            if not np.isfinite(full):
                continue
            n_used = int((np.isfinite(diff) & np.isfinite(ref)).sum())
            print(line(name, band, n_used, full, safe, meds))
            rows.append(dict(part=1, instrument=name, band=band, n=n_used,
                             span_vs_ps1=full, span_vs_ps1_unsaturated=safe,
                             bin_medians=meds, bin_counts=ns))

    print()
    print("=== 2 부. 그 치우침이 기기 쪽인가 기준 쪽인가 ===")
    print("같은 별에서  (기기 − PS1) = (기기 − Gaia기준) + (Gaia기준 − PS1)  으로 쪼갠다.")
    print("둘째 항은 기기와 무관하므로, 네 기기에서 같은 모양이면 기준 탓이다.")
    print()
    hdr2 = header("무엇")
    print(hdr2)
    print("-" * len(hdr2))

    for name, (tab, bands, kk) in joined.items():
        shown = False
        for band in bands:
            col, rcol = f"mag_cal_{band}", f"ref_{band}"
            if col not in tab.columns or rcol not in tab.columns:
                continue
            mine = pd.to_numeric(tab[col], errors="coerce").to_numpy(float)
            gref = pd.to_numeric(tab[rcol], errors="coerce").to_numpy(float)
            pref = np.full(len(tab), np.nan)
            hit = kk >= 0
            pref[hit] = ps1[f"ps1_{band}"].to_numpy(float)[kk[hit]]
            # **같은 별만 쓴다.** 한쪽에만 있는 별이 섞이면 두 항의 합이
            # 왼쪽 변과 안 맞아서 쪼갠 뜻이 없어진다.
            both = np.isfinite(mine) & np.isfinite(gref) & np.isfinite(pref)
            if int(both.sum()) < MIN_CELL * 3:
                continue
            if not shown:
                print(f"  [{name}]")
                shown = True
            axis = np.where(both, pref, np.nan)
            for what, d in (("기기−Gaia", mine - gref), ("Gaia−PS1", gref - pref)):
                meds, ns = bin_medians(np.where(both, d, np.nan), axis)
                full, safe = span_of(meds), span_of(meds, PS1_SAFE_MAG)
                print(line("  " + what, band, int(both.sum()), full, safe, meds))
                rows.append(dict(part=2, instrument=name, band=band, term=what,
                                 n=int(both.sum()), span=full,
                                 span_unsaturated=safe,
                                 bin_medians=meds, bin_counts=ns))

    out = REPO / "validation/drift_vs_reference.json"
    out.write_text(json.dumps(rows, ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
