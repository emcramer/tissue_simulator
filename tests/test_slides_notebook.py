"""Regression test: the slide-deck notebook is well-formed.

Guards ``docs/slides/tour.ipynb`` against silent breakage:

- the JSON schema is valid (nbformat.validate)
- the notebook is non-empty
- every cell carries ``cell.metadata["slideshow"]["slide_type"]`` so
  nbconvert --to slides knows where each cell goes on the deck
- kernelspec and language_info are present (so nbconvert can pick a
  kernel without prompting)

Skipped if ``nbformat`` isn't importable (minimal CI image).
"""

from __future__ import annotations

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


def test_slides_notebook_is_valid() -> None:
    nbformat = pytest.importorskip("nbformat")

    nb_path = _repo_root() / "docs" / "slides" / "tour.ipynb"
    assert nb_path.is_file(), f"slide-deck notebook not found at {nb_path}"

    nb = nbformat.read(str(nb_path), as_version=4)

    # Schema check - raises ValidationError on bad JSON.
    nbformat.validate(nb)

    assert nb.cells, "tour.ipynb has no cells"

    # Every cell must declare its slide_type so nbconvert --to slides
    # can place it. Forgetting this on a newly added cell silently
    # drops it from the deck.
    missing = []
    for i, cell in enumerate(nb.cells):
        slide_type = cell.metadata.get("slideshow", {}).get("slide_type")
        if not slide_type:
            missing.append(i)
    assert not missing, (
        "cells missing metadata.slideshow.slide_type: "
        f"{missing} (every cell needs one of slide/subslide/fragment/notes/skip/-)"
    )

    # nbconvert needs a kernel to execute. If kernelspec is missing it
    # falls back to prompting, which breaks CI.
    assert nb.metadata.get("kernelspec"), "tour.ipynb missing metadata.kernelspec"
    assert nb.metadata.get("language_info"), "tour.ipynb missing metadata.language_info"
