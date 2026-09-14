#!/usr/bin/env python3
"""Render the Marimo tour to a static HTML page."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NB = HERE / "tour.py"
OUT = HERE / "tour.html"


def main() -> int:
    if not NB.exists():
        print(f"error: tour notebook not found at {NB}", file=sys.stderr)
        return 1
    if shutil.which("marimo") is None:
        print("error: 'marimo' not on PATH. Install the docs extras: pip install -e \".[docs]\"", file=sys.stderr)
        return 1
    print(f"Rendering tour from {NB} ...", flush=True)
    rc = subprocess.call(
        [
            "marimo", "export", "html",
            "--include-code",
            str(NB),
            "-o", str(OUT),
        ],
    )
    if rc != 0:
        return rc
    if not OUT.exists():
        print(f"error: marimo reported success but {OUT} is missing", file=sys.stderr)
        return 1
    size_kb = OUT.stat().st_size // 1024
    print(f"OK: wrote {OUT} ({size_kb} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
