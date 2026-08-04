import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

API = "https://archive-api.lco.global/frames/"
OUT = Path(__file__).with_name("lco_ngc6811_nights.txt")


def q(**kw):
    kw.setdefault("public", "true")
    last = None
    for i in range(6):
        try:
            url = API + "?" + urllib.parse.urlencode(kw)
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.loads(r.read())
        except Exception as exc:
            last = exc
            time.sleep(10 * (i + 1))
    raise RuntimeError(f"query failed after retries: {last}")


def q_all(**kw):
    """Paginated fetch (API limit is 100/page); follows 'next' links."""
    kw["limit"] = 100
    res = q(**kw)
    frames = list(res["results"])
    next_url = res.get("next")
    while next_url:
        for i in range(6):
            try:
                with urllib.request.urlopen(next_url, timeout=180) as r:
                    page = json.loads(r.read())
                break
            except Exception:
                time.sleep(10 * (i + 1))
        else:
            break
        frames.extend(page["results"])
        next_url = page.get("next")
    return res["count"], frames


lines = []
days = defaultdict(lambda: defaultdict(int))
expo = defaultdict(set)
for f in ("U", "B", "V"):
    count, frames = q_all(target_name="NGC 6811", primary_optical_element=f,
                          reduction_level=0)
    lines.append(f"{f}: count={count} fetched={len(frames)}")
    for fr in frames:
        key = (fr["observation_day"], fr["instrument_id"])
        days[key][f] += 1
        expo[key].add((f, fr["exposure_time"]))

for k in sorted(days):
    lines.append(f"{k[0]} {k[1]} {dict(days[k])} {sorted(expo[k])}")

OUT.write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines))
