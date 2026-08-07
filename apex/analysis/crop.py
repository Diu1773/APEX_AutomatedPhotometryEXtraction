"""Qt-free image crop (headless port of ``CropSelectorWindow.apply_crop``).

This module is a VERBATIM RELOCATION of the crop compute body of
``apex.gui.workflow.step2_crop_selector.CropSelectorWindow.apply_crop`` (the
``crop_one`` closure and the ThreadPoolExecutor loop) plus the cache-validation
helper ``_is_crop_cache_valid``. Only the Qt coupling is substituted:

* ``self.crop_x0/y0/x1/y1`` -> the ``crop_x0/y0/x1/y1`` parameters
* ``self.params`` -> the ``params`` parameter, ``self.file_manager.filenames``
  -> the ``file_list`` parameter
* ``QProgressDialog`` / ``QApplication.processEvents`` / ``QMessageBox`` ->
  optional ``progress_cb`` callback and a ``should_stop`` predicate

The crop slice, the FITS header keywords, the cache freshness rules, the worker
count, and the parallelism are unchanged.

LAYER RULE: this module must NOT import apex.gui or PyQt5.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

from apex.utils.constants import get_parallel_workers
from apex.utils.step_paths import step2_cropped_dir


def _is_crop_cache_valid(original_path: Path, cropped_path: Path,
                         crop_x0: int, crop_y0: int, crop_x1: int, crop_y1: int) -> bool:
    """Return True if existing cropped file can be reused for current crop rectangle."""
    if not cropped_path.exists():
        return False
    try:
        src_stat = original_path.stat()
        dst_stat = cropped_path.stat()
        if int(dst_stat.st_mtime_ns) < int(src_stat.st_mtime_ns):
            return False

        header = fits.getheader(cropped_path, ext=0)
        if int(header.get("CROP_X0", -1)) != int(crop_x0):
            return False
        if int(header.get("CROP_Y0", -1)) != int(crop_y0):
            return False
        if int(header.get("CROP_X1", -1)) != int(crop_x1):
            return False
        if int(header.get("CROP_Y1", -1)) != int(crop_y1):
            return False

        expected_w = int(crop_x1 - crop_x0)
        expected_h = int(crop_y1 - crop_y0)
        if int(header.get("NAXIS1", -1)) != expected_w:
            return False
        if int(header.get("NAXIS2", -1)) != expected_h:
            return False
        return True
    except Exception:
        return False


def run_crop(file_list, params, crop_x0, crop_y0, crop_x1, crop_y1, *,
             progress_cb=None, should_stop=None, logger=None):
    """Headless image crop over ``file_list`` for a fixed crop rectangle.

    Verbatim relocation of ``CropSelectorWindow.apply_crop``'s worker body:
    same crop slice, same FITS header keywords, same cache-reuse rules, same
    worker count and ThreadPoolExecutor logic. Qt dialogs/signals are replaced
    by an optional ``progress_cb(done, total, filename)`` callback and a
    ``should_stop()`` predicate.

    Returns ``(n_total, errors)`` where ``errors`` is a list of
    ``(filename, message)`` tuples (empty on full success), mirroring the GUI's
    ``errors`` accumulator so the caller can decide on the skip marker.
    """
    crop_x0 = int(crop_x0)
    crop_y0 = int(crop_y0)
    crop_x1 = int(crop_x1)
    crop_y1 = int(crop_y1)

    def _should_stop():
        return should_stop() if should_stop else False

    # Create cropped directory only
    cropped_dir = step2_cropped_dir(params.P.result_dir)
    cropped_dir.mkdir(parents=True, exist_ok=True)

    files = list(file_list)

    def crop_one(filename: str) -> str:
        original_path = params.get_file_path(filename)
        cropped_path = cropped_dir / filename

        if _is_crop_cache_valid(original_path, cropped_path,
                                crop_x0, crop_y0, crop_x1, crop_y1):
            return filename  # reuse existing cropped file

        with fits.open(original_path, memmap=False) as hdul:
            data = None
            header = None
            for hdu in hdul:
                if hdu.data is not None:
                    data = hdu.data
                    header = hdu.header
                    break
            if data is None:
                raise RuntimeError(f"No image data found in FITS file: {original_path.name}")
            header = header.copy() if header is not None else hdul[0].header.copy()

        cropped_data = data[crop_y0:crop_y1, crop_x0:crop_x1]

        header['CROPPED'] = (True, 'Image has been cropped')
        header['CROP_X0'] = (crop_x0, 'Crop region start X')
        header['CROP_Y0'] = (crop_y0, 'Crop region start Y')
        header['CROP_X1'] = (crop_x1, 'Crop region end X')
        header['CROP_Y1'] = (crop_y1, 'Crop region end Y')

        fits.PrimaryHDU(data=cropped_data, header=header).writeto(
            cropped_path, overwrite=True
        )
        return filename

    max_workers = get_parallel_workers(params, stage="crop")
    max_workers = min(max_workers, len(files)) if files else 1

    errors = []
    if max_workers <= 1 or len(files) <= 1:
        for i, filename in enumerate(files, 1):
            if _should_stop():
                break
            (progress_cb(i, len(files), filename) if progress_cb else None)
            try:
                crop_one(filename)
            except Exception as e:
                errors.append((filename, str(e)))
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(crop_one, fn): fn for fn in files}
            for i, future in enumerate(as_completed(futures), 1):
                if _should_stop():
                    for f in futures:
                        f.cancel()
                    break
                filename = futures[future]
                try:
                    future.result()
                except Exception as e:
                    errors.append((filename, str(e)))
                (progress_cb(i, len(files), filename) if progress_cb else None)

    return len(files), errors
