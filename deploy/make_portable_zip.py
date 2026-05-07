"""Create the APEX portable release zip."""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


def make_zip(source: Path, output: Path) -> int:
    source = source.resolve()
    output = output.resolve()
    if not source.is_dir():
        print(f"[ERROR] Source directory not found: {source}")
        return 1

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()

    files = [path for path in source.rglob("*") if path.is_file()]
    with ZipFile(output, "w", compression=ZIP_DEFLATED, compresslevel=6) as zf:
        for index, path in enumerate(files, start=1):
            arcname = path.relative_to(source)
            zf.write(path, arcname.as_posix())
            if index % 500 == 0:
                print(f"[ZIP] {index}/{len(files)} files")
    print(f"[OK] Portable zip: {output}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    return make_zip(args.source, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
