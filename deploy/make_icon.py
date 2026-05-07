"""Build the Windows .ico file from the APEX SVG logo."""

from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

from PIL import Image
from PyQt5.QtCore import QSize
from PyQt5.QtGui import QGuiApplication
from PyQt5.QtGui import QImage, QPainter
from PyQt5.QtSvg import QSvgRenderer


def render_svg(svg_path: Path, png_path: Path, size: int) -> None:
    renderer = QSvgRenderer(str(svg_path))
    image = QImage(QSize(size, size), QImage.Format_ARGB32)
    image.fill(0)
    painter = QPainter(image)
    renderer.render(painter)
    painter.end()
    if not image.save(str(png_path), "PNG"):
        raise RuntimeError(f"Failed to render icon PNG: {png_path}")


def make_icon(svg_path: Path, output_path: Path) -> int:
    app = QGuiApplication.instance() or QGuiApplication([])
    svg_path = svg_path.resolve()
    output_path = output_path.resolve()
    if not svg_path.is_file():
        print(f"[ERROR] SVG not found: {svg_path}")
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        png_path = Path(tmp) / "apex_256.png"
        render_svg(svg_path, png_path, 256)
        image = Image.open(png_path).convert("RGBA")
        image.save(
            output_path,
            format="ICO",
            sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)],
        )
    print(f"[OK] Icon: {output_path}")
    app.quit()
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--svg", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    return make_icon(args.svg, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
