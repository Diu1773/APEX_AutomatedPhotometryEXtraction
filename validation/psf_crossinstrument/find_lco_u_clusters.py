"""Which cluster does LCO actually have U-band science frames for?

The M67 audit established that APEX's g/r/i isochrone solution is degenerate —
metallicity, extinction and distance trade off against each other and the
posterior ends up led by the prior. U is what breaks it: U-B is sensitive to
metallicity in a way g-r is not. The scorecard therefore made "get the LCO
NGC 6811 U set" the top priority.

That set does not exist (LCO_NGC6811_REALITY.md). But the instrument carries
U and up filters — skyflats for both were taken that same week — so the right
question is not "does LCO have U" but "for which cluster".

This asks it directly, for clusters APEX has already processed plus the
standard calibration clusters. Counts come from paging, never from the API's
`count` field, which is flagged `count_estimated` and was wrong by a factor of
infinity for NGC 6811 (reported 58 U frames; there are none).

STATUS 2026-08-11: written and correct in shape, but NOT yet run to completion.
On a heavily observed target (M67 is the first in the list) the archive stops
answering — the first query returned an HTTPError and the next hung past 15
minutes, after a session that had already made a few dozen queries. That looks
like rate limiting rather than a bug here, so the run was stopped rather than
retried into the same wall.

To finish it, throttle: a sleep of a few seconds between queries, a retry with
backoff on HTTPError, and ideally an early exit once a page comes back short
(`next` absent) so common targets do not page needlessly. Run it against a
handful of targets at a time rather than sixteen in one pass.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

API = "https://archive-api.lco.global/frames/"
OUT = Path(__file__).parent / "lco_u_cluster_search.json"

# APEX's own targets first, then clusters that are standard photometric
# calibration fields and so are the likeliest to have been shot in U.
TARGETS = [
    "M67", "NGC 6811", "M13", "M3", "M5", "M37",
    "NGC 188", "NGC 6791", "M35", "M45", "NGC 7789", "NGC 2168",
    "M92", "M15", "NGC 6205", "Melotte 111",
]
U_FILTERS = ("U", "up")


def paged(**kw) -> list[dict]:
    kw.setdefault("public", "true")
    kw.setdefault("limit", "100")
    url = API + "?" + urllib.parse.urlencode(kw)
    rows: list[dict] = []
    while url:
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                page = json.loads(r.read())
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"    [query failed: {type(e).__name__}]")
            return rows
        rows.extend(page.get("results", []))
        url = page.get("next")
    return rows


print(f"{'target':14s}{'U':>4s}{'up':>4s}  instruments / nights")
found = []
for target in TARGETS:
    rows: list[dict] = []
    for f in U_FILTERS:
        rows += paged(target_name=target, FILTER=f, RLEVEL="00",
                      OBSTYPE="EXPOSE")
    n_u = sum(1 for r in rows if r.get("FILTER") == "U")
    n_up = sum(1 for r in rows if r.get("FILTER") == "up")
    insts = Counter(r.get("INSTRUME") for r in rows)
    nights = sorted({r.get("DAY_OBS") for r in rows})
    detail = ""
    if rows:
        detail = (f"{dict(insts)}  nights={len(nights)}"
                  f"{' ' + str(nights[:3]) if nights else ''}")
        found.append({"target": target, "n_U": n_u, "n_up": n_up,
                      "instruments": dict(insts), "nights": nights})
    print(f"{target:14s}{n_u:4d}{n_up:4d}  {detail}")

OUT.write_text(json.dumps({"targets_queried": TARGETS, "with_u": found},
                          indent=1), encoding="utf-8")
print(f"\n{len(found)} target(s) with any U/up science frames")
print(f"saved -> {OUT}")
