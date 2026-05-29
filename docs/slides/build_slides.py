#!/usr/bin/env python3
"""Render the slide deck.

Executes ``docs/slides/tour.ipynb`` end-to-end and writes
``docs/slides/tour.slides.html`` (reveal.js).

Usage::

    python docs/slides/build_slides.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent  # docs/slides/
NB = HERE / "tour.ipynb"


def main() -> int:
    if not NB.exists():
        print(f"error: notebook not found at {NB}", file=sys.stderr)
        return 1
    if shutil.which("jupyter") is None:
        print("error: 'jupyter' not on PATH. Install the docs extras: pip install -e \".[docs]\"", file=sys.stderr)
        return 1
    print(f"Rendering slides from {NB} ...", flush=True)
    rc = subprocess.call(
        [
            "jupyter", "nbconvert",
            "--to", "slides",
            "--execute",
            "--ExecutePreprocessor.timeout=300",
            "--output", "tour",
            "--output-dir", str(HERE),
            str(NB),
        ],
    )
    if rc != 0:
        return rc
    out = HERE / "tour.slides.html"
    if not out.exists():
        print(f"error: nbconvert reported success but {out} is missing", file=sys.stderr)
        return 1
    size_kb = out.stat().st_size // 1024
    print(f"OK: wrote {out} ({size_kb} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
