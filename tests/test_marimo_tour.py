"""Regression test: docs/slides/tour.py is a valid Marimo notebook.

Guards the Marimo-rendered tour (the replacement for the old reveal.js
slide deck) against the kinds of breakage that would only surface in CI:

- a SyntaxError in tour.py
- an ImportError from a top-level cell
- Marimo's own MultipleDefinitionError, raised when a top-level name is
  bound in more than one cell (Marimo's reactive-graph invariant)
- accidental deletion of all the @app.cell decorators
- the title cell's headline going missing

If Marimo isn't installed in this environment (e.g. on a minimal CI image
that skips the docs extras), the test is skipped cleanly rather than
failed: it's a tour-source-validity check, not a hard prerequisite for
shipping tissue_simulator itself.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

marimo = pytest.importorskip("marimo")

# Marimo's MultipleDefinitionError lives at marimo._ast.errors as of v0.23.
# Be defensive across versions: if the import fails, fall back to a sentinel
# that the except clause will simply never match (so other errors still fire).
try:
    from marimo._ast.errors import MultipleDefinitionError as _MDE
except Exception:  # pragma: no cover - depends on marimo internals
    class _MDE(Exception):  # type: ignore[no-redef]
        """Fallback when marimo's MultipleDefinitionError isn't importable."""


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


def _tour_path() -> Path:
    path = _repo_root() / "docs" / "slides" / "tour.py"
    assert path.is_file(), f"docs/slides/tour.py is missing at {path}"
    return path


def test_tour_loads_as_marimo_app() -> None:
    """tour.py imports cleanly and exposes a marimo.App named `app`."""
    path = _tour_path()
    spec = importlib.util.spec_from_file_location("tour", path)
    assert spec is not None and spec.loader is not None, (
        f"could not build an import spec for {path}"
    )
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except SyntaxError as exc:
        pytest.fail(f"docs/slides/tour.py has a SyntaxError: {exc}")
    except _MDE as exc:
        pytest.fail(
            "docs/slides/tour.py violates Marimo's top-level-uniqueness "
            f"invariant (MultipleDefinitionError): {exc}"
        )
    except ImportError as exc:
        pytest.fail(f"docs/slides/tour.py raised ImportError on load: {exc}")

    assert hasattr(module, "app"), "tour.py must define a top-level `app`"
    assert isinstance(module.app, marimo.App), (
        f"tour.py's `app` must be a marimo.App; got {type(module.app)!r}"
    )


def test_tour_has_expected_cell_count() -> None:
    """tour.py has at least 20 cells (the Driver reported 23-24)."""
    path = _tour_path()

    # Prefer the Marimo introspection path when available; fall back to a
    # simple textual count of @app.cell decorators if the API shape shifts.
    cell_count: int
    try:
        spec = importlib.util.spec_from_file_location("tour_for_count", path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        cell_count = len(list(module.app._cell_manager.cell_ids()))
    except Exception:
        cell_count = path.read_text(encoding="utf-8").count("@app.cell")

    assert cell_count >= 20, (
        f"tour.py has only {cell_count} cells; expected the full tour "
        "(>= 20). Did a section get accidentally deleted?"
    )


def test_tour_title_cell_present() -> None:
    """The headline 'tissue_simulator: an end-to-end tour' is in the source."""
    src = _tour_path().read_text(encoding="utf-8")
    assert "tissue_simulator: an end-to-end tour" in src, (
        "tour.py is missing its title-cell headline; "
        "the landing-page title would render empty."
    )
