"""Download and normalize BaSTI Johnson-Cousins isochrones for APEX.

The BaSTI-IAC precomputed archives contain one file per age and metallicity.
APEX Step 12 expects a single numeric table whose first columns are compatible
with its existing grid loader. This script downloads a compact near-solar grid
and converts it without extracting the archives to disk.
"""

from __future__ import annotations

import argparse
import io
import math
import re
import tarfile
import urllib.request
from pathlib import Path

import numpy as np


BASTI_BASE_URL = (
    "http://basti-iac.oa-abruzzo.inaf.it/PREISOCS/P00O1D0E1Y247"
)
ARCHIVES = {
    -0.20: "isocz102y260p00o1d0e1.isc_john.tar.gz",
    -0.08: "isocz132y264p00o1d0e1.isc_john.tar.gz",
    0.06: "isocz172y269p00o1d0e1.isc_john.tar.gz",
}
EXPECTED_SIZES = {
    "isocz102y260p00o1d0e1.isc_john.tar.gz": 26_994_421,
    "isocz132y264p00o1d0e1.isc_john.tar.gz": 29_726_787,
    "isocz172y269p00o1d0e1.isc_john.tar.gz": 33_395_783,
}

_AGE_RE = re.compile(r"/(\d+)z")
_META_RE = re.compile(
    r"\[M/H\]\s*=\s*(?P<mh>[+-]?\d+(?:\.\d+)?)"
    r".*?Z\s*=\s*(?P<z>\d+(?:\.\d+)?)"
    r".*?Age \(Myr\)\s*=\s*(?P<age>\d+(?:\.\d+)?)"
)

HEADER = (
    "# Isochrone model: BaSTI-IAC updated models\n"
    "# Model grid: solar-scaled, overshooting=yes, diffusion=no, eta=0.3, "
    "primordial Y=0.247\n"
    "# Photometric system: Johnson-Cousins / Bessell (BaSTI-IAC)\n"
    "# Source: http://basti-iac.oa-abruzzo.inaf.it/isocs.html\n"
    "# APEX normalized columns follow. Placeholder columns preserve the common "
    "Step 12 grid layout.\n"
    "# Zini MH logAge Mini int_IMF Mass logL logTe logg label "
    "Umag BXmag Bmag Vmag Rmag Imag Jmag Hmag Kmag Lpmag Lmag Mmag\n"
)


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "APEX-BaSTI-downloader/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        with temporary.open("wb") as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
    temporary.replace(destination)


def ensure_archives(source_dir: Path) -> list[tuple[float, Path]]:
    """Download the three near-solar BaSTI archives when absent."""
    results = []
    for metallicity, filename in ARCHIVES.items():
        path = source_dir / filename
        expected_size = EXPECTED_SIZES[filename]
        if not path.exists() or path.stat().st_size != expected_size:
            print(f"Downloading {filename} ({expected_size / 1024**2:.1f} MiB)...")
            _download(f"{BASTI_BASE_URL}/{filename}", path)
        results.append((metallicity, path))
    return results


def _member_age_myr(name: str) -> int | None:
    match = _AGE_RE.search(name.replace("\\", "/"))
    return int(match.group(1)) if match else None


def _parse_member(data: bytes) -> tuple[dict[str, float], np.ndarray]:
    text = data.decode("ascii", errors="replace")
    metadata = None
    numeric_lines = []
    for line in text.splitlines():
        if metadata is None and "[M/H]" in line and "Age (Myr)" in line:
            match = _META_RE.search(line)
            if match:
                metadata = {
                    key: float(value) for key, value in match.groupdict().items()
                }
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            numeric_lines.append(stripped)

    if metadata is None:
        raise ValueError("BaSTI metadata line was not found")
    if not numeric_lines:
        raise ValueError("BaSTI member contains no numeric rows")

    raw = np.loadtxt(io.StringIO("\n".join(numeric_lines)), ndmin=2)
    if raw.shape[1] != 16:
        raise ValueError(f"Expected 16 BaSTI columns, found {raw.shape[1]}")
    return metadata, raw


def normalize_archives(
    archives: list[tuple[float, Path]],
    output_path: Path,
    age_min_myr: int,
    age_max_myr: int,
) -> tuple[int, int]:
    """Write selected BaSTI members in the common APEX numeric layout."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".part")
    row_count = 0
    grid_count = 0

    with temporary.open("w", encoding="ascii", newline="\n") as output:
        output.write(HEADER)
        for expected_mh, archive_path in archives:
            with tarfile.open(archive_path, mode="r:gz") as archive:
                members = []
                for member in archive.getmembers():
                    age_myr = _member_age_myr(member.name)
                    if (
                        member.isfile()
                        and age_myr is not None
                        and age_min_myr <= age_myr <= age_max_myr
                    ):
                        members.append((age_myr, member))

                for age_myr, member in sorted(members):
                    extracted = archive.extractfile(member)
                    if extracted is None:
                        continue
                    metadata, raw = _parse_member(extracted.read())
                    if not math.isclose(
                        metadata["mh"], expected_mh, abs_tol=0.011
                    ):
                        raise ValueError(
                            f"{member.name}: expected [M/H]={expected_mh}, "
                            f"found {metadata['mh']}"
                        )

                    n_rows = raw.shape[0]
                    normalized = np.zeros((n_rows, 22), dtype=float)
                    normalized[:, 0] = metadata["z"]
                    normalized[:, 1] = metadata["mh"]
                    normalized[:, 2] = math.log10(
                        metadata["age"] * 1_000_000.0
                    )
                    normalized[:, 3] = raw[:, 0]
                    normalized[:, 5] = raw[:, 1]
                    normalized[:, 6] = raw[:, 2]
                    normalized[:, 7] = raw[:, 3]
                    normalized[:, 9] = 1.0
                    normalized[:, 10:] = raw[:, 4:]
                    np.savetxt(output, normalized, fmt="%.8g")
                    row_count += n_rows
                    grid_count += 1
                    print(
                        f"Converted [M/H]={metadata['mh']:+.2f}, "
                        f"age={age_myr:4d} Myr ({n_rows} rows)"
                    )

    temporary.replace(output_path)
    return grid_count, row_count


def build_binary_cache(output_path: Path) -> Path:
    """Create the binary sidecar used by Step 12 for fast loading."""
    data = np.loadtxt(output_path, comments="#")
    cache_path = output_path.with_suffix(output_path.suffix + ".apex.npy")
    np.save(cache_path, data, allow_pickle=False)
    return cache_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("isochrone") / "BaSTI" / "johnson",
    )
    parser.add_argument("--age-min-myr", type=int, default=300)
    parser.add_argument("--age-max-myr", type=int, default=3200)
    args = parser.parse_args()

    archives = ensure_archives(args.output_dir / "source")
    output_path = args.output_dir / "basti_p00_o1_d0_e1_y247_johnson.dat"
    grid_count, row_count = normalize_archives(
        archives,
        output_path,
        args.age_min_myr,
        args.age_max_myr,
    )
    cache_path = build_binary_cache(output_path)
    print(
        f"Wrote {grid_count} age/metallicity grids and {row_count} rows:\n"
        f"  {output_path}\n"
        f"  {cache_path}"
    )


if __name__ == "__main__":
    main()
