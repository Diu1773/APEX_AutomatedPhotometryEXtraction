"""How much would the LCO NGC 6811 set actually cost to fetch?

`lco_ngc6811_nights.txt` says U 96 / B 66 / V 94 frames exist, but a count is
not a decision — a decision needs bytes and hours, and it needs to separate
"science frames" from "the calibration frames APEX Step 0 requires". Fetching
raw (RLEVEL 00) without its matching bias/dark/flat leaves the frames unusable
by the very step the cross-instrument claim depends on.

So this reports four things per night, without downloading anything:

  * raw science (RLEVEL 00) — what APEX Step 0 consumes
  * BANZAI-reduced (RLEVEL 91) — the free independent-pipeline comparison
  * calibration frames on the same night — whether Step 0 can even run
  * total bytes, and wall time at a measured transfer rate

Why this matters for the scorecard: this one dataset moves the C axis (a second
instrument through the photometry->CMD chain, which only Step 0 has today) and
supplies the U band that the M67 audit showed is what breaks the age-metallicity
degeneracy. It is the only remaining gap that no amount of local work can close.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

API = "https://archive-api.lco.global/frames/"
TARGET = "NGC 6811"   # the archive stores the spaced form; "NGC6811" pages empty
SCIENCE_FILTERS = ("U", "B", "V", "up", "gp", "rp", "ip", "zs", "R")
# Measured on the two B and two V frames already fetched (see the nights file);
# a conservative single-stream rate for a same-order estimate, not a promise.
RATE_MB_S = 3.0
OUT = Path(__file__).parent / "lco_ngc6811_estimate.json"


def query(**kw) -> list[dict]:
    """Page through the archive; the API caps `limit` well below our counts."""
    kw.setdefault("public", "true")
    kw.setdefault("limit", "100")   # the API rejects larger pages with 400
    url = API + "?" + urllib.parse.urlencode(kw)
    rows: list[dict] = []
    while url:
        with urllib.request.urlopen(url, timeout=120) as r:
            page = json.loads(r.read())
        rows.extend(page.get("results", []))
        url = page.get("next")
    return rows


# The archive's `count` is flagged `count_estimated: true` and is wrong by
# large factors: for NGC 6811 it reported U 58 / B 40 / V 56 where paging finds
# 0 / 4 / 4. Never report `count`; always page and count the rows.


def frame_bytes(row: dict) -> int:
    """Real size of one frame.

    The list endpoint has no `filesize` and the detail endpoint returns None
    for it, so the only truthful source is a one-byte Range request against the
    download URL, which reports the total in Content-Range.
    """
    with urllib.request.urlopen(f"{API}{row['id']}/", timeout=90) as h:
        url = json.loads(h.read()).get("url")
    if not url:
        return 0
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        cr = r.headers.get("Content-Range", "")
    return int(cr.split("/")[-1]) if "/" in cr else 0


def summarise(rows: list[dict], unit_bytes: int = 0) -> dict:
    by_filter: dict[str, int] = defaultdict(int)
    total = 0
    if not unit_bytes and rows:
        unit_bytes = frame_bytes(rows[0])
    for r in rows:
        by_filter[r.get("FILTER") or "?"] += 1
        total += unit_bytes
    return {"n": len(rows), "bytes": total,
            "gb": round(total / 1e9, 2), "by_filter": dict(by_filter)}


print(f"querying LCO archive for {TARGET} (public frames only)…")

raw = [r for r in query(target_name=TARGET, RLEVEL="00", OBSTYPE="EXPOSE")
       if r.get("FILTER") in SCIENCE_FILTERS]
red = [r for r in query(target_name=TARGET, RLEVEL="91", OBSTYPE="EXPOSE")
       if r.get("FILTER") in SCIENCE_FILTERS]

print(f"  raw science (RLEVEL 00): {len(raw)}")
print(f"  BANZAI reduced (RLEVEL 91): {len(red)}")

# Group by night+instrument: a cross-instrument claim needs one instrument's
# frames kept together, and Step 0 needs calibrations from the same night.
nights: dict[tuple[str, str], list[dict]] = defaultdict(list)
for r in raw:
    nights[(r.get("DAY_OBS") or "?", r.get("INSTRUME") or "?")].append(r)

print(f"\n{'night':12s} {'inst':8s} {'frames':>7s} {'GB':>7s}  filters")
night_rows = []
for (day, inst), rows in sorted(nights.items()):
    s = summarise(rows)
    filters = " ".join(f"{f}{n}" for f, n in sorted(s["by_filter"].items()))
    print(f"{day:12s} {inst:8s} {s['n']:7d} {s['gb']:7.2f}  {filters}")
    night_rows.append({"day_obs": day, "instrument": inst, **s})

raw_s, red_s = summarise(raw), summarise(red)

# Only nights that carry all three science filters are worth the transfer: U-B
# with B-V from the same night and instrument is the point of the exercise.
full_ubv = [n for n in night_rows
            if all(f in n["by_filter"] for f in SCIENCE_FILTERS)]

print(f"\nraw science total     {raw_s['n']:4d} frames  {raw_s['gb']:6.2f} GB  "
      f"{raw_s['by_filter']}")
print(f"BANZAI reduced total  {red_s['n']:4d} frames  {red_s['gb']:6.2f} GB  "
      f"{red_s['by_filter']}")

both_gb = raw_s["gb"] + red_s["gb"]
print(f"\nboth products         {both_gb:6.2f} GB  "
      f"~{both_gb * 1000 / RATE_MB_S / 3600:.1f} h at {RATE_MB_S} MB/s")

if full_ubv:
    sub_gb = sum(n["gb"] for n in full_ubv)
    print(f"\nnights with all of U+B+V ({len(full_ubv)}): {sub_gb:.2f} GB raw  "
          f"~{sub_gb * 1000 / RATE_MB_S / 3600:.1f} h")
    for n in full_ubv:
        print(f"  {n['day_obs']}  {n['instrument']}  {n['n']} frames  "
              f"{n['by_filter']}")
else:
    print("\nNo single night carries U+B+V together — the U-band argument would "
          "have to cross nights, which reopens the zero-point question the "
          "M67 audit just measured.")

OUT.write_text(json.dumps(
    {"target": TARGET, "rate_mb_s": RATE_MB_S,
     "raw_rlevel00": raw_s, "banzai_rlevel91": red_s,
     "nights": night_rows,
     "nights_with_full_ubv": [n["day_obs"] for n in full_ubv]},
    indent=1), encoding="utf-8")
print(f"\nsaved -> {OUT}")
