"""Moved to `apex.analysis.header_target` (2026-08-16).

It never used Qt, but living here made the headless scan step import
`apex.gui` for a function that has nothing to do with a window. Kept as a
re-export so existing imports keep working.
"""

from apex.analysis.header_target import select_header_target  # noqa: F401
