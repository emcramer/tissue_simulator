"""Regression test: the MkDocs Material site builds cleanly under --strict.

This test guards the docs-site scaffolding (mkdocs.yml, docs/index.md,
docs/changelog.md, docs/reference/*.md, etc.) against regressions:

- broken cross-references inside pages the site owns
- pages dropped from the navigation
- mkdocstrings autodoc failing to import a tissue_simulator submodule
- the include-markdown plugin losing track of ../CHANGELOG.md

If any of the docs toolchain isn't installed (e.g. on a minimal CI image),
the test is skipped rather than failed - it's a docs-site regression check,
not a hard prerequisite for shipping tissue_simulator itself.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


def _repo_root() -> Path:
    """Walk up from this file until we find mkdocs.yml."""
    here = Path(__file__).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "mkdocs.yml").is_file():
            return candidate
    raise RuntimeError(
        "Could not find mkdocs.yml by walking up from "
        f"{here}; is the docs site checked out?"
    )


def test_mkdocs_strict_build(tmp_path: Path) -> None:
    # Skip cleanly if any of the docs deps isn't importable in this env.
    pytest.importorskip("mkdocs")
    pytest.importorskip("material")
    pytest.importorskip("mkdocstrings")
    pytest.importorskip("mkdocs_include_markdown_plugin")

    repo_root = _repo_root()
    site_dir = tmp_path / "site"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mkdocs",
            "build",
            "--strict",
            "--site-dir",
            str(site_dir),
        ],
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=240,
    )

    assert result.returncode == 0, (
        "mkdocs build --strict failed (exit "
        f"{result.returncode}).\n"
        f"--- stdout ---\n{result.stdout}\n"
        f"--- stderr ---\n{result.stderr}\n"
    )

    # Sanity-check that the rendered site contains the pages we care
    # about. Any one of these missing means the nav or the autodoc
    # silently dropped a page.
    sentinels = [
        site_dir / "index.html",
        site_dir / "quickstart" / "index.html",
        site_dir / "changelog" / "index.html",
        site_dir / "reference" / "tissue" / "index.html",
        site_dir / "reference" / "graph_coloring" / "index.html",
    ]
    missing = [str(p.relative_to(site_dir)) for p in sentinels if not p.is_file()]
    assert not missing, (
        "mkdocs build succeeded but expected pages are missing from the "
        f"rendered site: {missing}"
    )
