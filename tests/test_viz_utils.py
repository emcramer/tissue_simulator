"""Tests for ``tissue_simulator._viz_utils.make_color_map`` and downstream
deterministic color assignment in ``TissueSection.visualize``.

The CPython interpreter randomizes its per-process hash seed
(``PYTHONHASHSEED``) at startup, which leaks into the iteration order of
``set()`` / unsorted ``dict`` objects. Prior to v0.1.12 our ``visualize_*``
functions built their cell-type → color maps by iterating an unsorted set, so
two independent Python processes could assign different colors to the same
cell type.  These tests pin down the new deterministic behaviour.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap

import pytest

from tissue_simulator._viz_utils import make_color_map


# ---------------------------------------------------------------------------
# Unit tests on the helper itself (in-process — these run fast)
# ---------------------------------------------------------------------------


def test_make_color_map_basic():
    """Keys come back in alphabetical order."""
    result = make_color_map(["cancer", "immune", "stroma"])
    assert list(result.keys()) == ["cancer", "immune", "stroma"]


def test_make_color_map_deterministic_across_input_order():
    """Same set of types in different orders → identical mapping."""
    forward = make_color_map(["cancer", "immune", "stroma"])
    reverse = make_color_map(["stroma", "immune", "cancer"])
    assert forward.keys() == reverse.keys()
    for key in forward:
        assert tuple(forward[key]) == tuple(reverse[key])


def test_make_color_map_deduplicates():
    """Duplicate inputs collapse to one entry per unique type."""
    result = make_color_map(["a", "b", "a", "b", "a"])
    assert set(result.keys()) == {"a", "b"}
    assert len(result) == 2


def test_make_color_map_empty():
    """Empty input → empty dict, no exception."""
    result = make_color_map([])
    assert result == {}


def test_make_color_map_palette_passthrough():
    """A user-supplied palette is consumed in order against sorted types."""
    palette = [(1.0, 0.0, 0.0, 1.0), (0.0, 1.0, 0.0, 1.0)]
    result = make_color_map(["a", "b"], palette=palette)
    assert result == {
        "a": (1.0, 0.0, 0.0, 1.0),
        "b": (0.0, 1.0, 0.0, 1.0),
    }


def test_make_color_map_palette_too_short_raises():
    """ValueError when the palette has fewer entries than unique types."""
    with pytest.raises(ValueError):
        make_color_map(["a", "b", "c"], palette=[(1, 0, 0, 1)])


# ---------------------------------------------------------------------------
# Subprocess-based determinism tests
# ---------------------------------------------------------------------------
# CPython only honours ``PYTHONHASHSEED`` at interpreter startup. Setting
# ``os.environ["PYTHONHASHSEED"]`` after import has no effect, so these tests
# spawn fresh interpreters.


_HELPER_SCRIPT = textwrap.dedent(
    """
    import json
    from tissue_simulator._viz_utils import make_color_map

    result = make_color_map(["cancer", "immune", "stroma"])
    print(json.dumps(list(result.keys())))
    """
).strip()


def _run_with_hashseed(script: str, hashseed: str) -> str:
    """Run ``script`` in a fresh Python with the given PYTHONHASHSEED."""
    env = {**os.environ, "PYTHONHASHSEED": hashseed}
    proc = subprocess.run(
        [sys.executable, "-c", script],
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    return proc.stdout.strip()


def test_make_color_map_pythonhashseed_independent():
    """Two interpreters with different PYTHONHASHSEED → identical keys."""
    out_a = _run_with_hashseed(_HELPER_SCRIPT, "0")
    out_b = _run_with_hashseed(_HELPER_SCRIPT, "1")
    assert out_a == out_b
    # And the value really is the sorted ordering, not just consistently
    # arbitrary.
    assert json.loads(out_a) == ["cancer", "immune", "stroma"]


_TISSUE_SCRIPT = textwrap.dedent(
    """
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from tissue_simulator import TissueSection

    try:
        tissue = TissueSection(
            150, 150, 40,
            {"cancer": (5, 8), "immune": (3, 6), "stroma": (8, 12)},
            seed=42,
        )
        tissue.generate_cells(max_attempts=300)

        # ``visualize()`` calls ``plt.show()`` at the end. Under the Agg
        # backend that's a no-op, but we still want the figure object so we
        # can inspect the legend handles. ``plt.gcf()`` after the call gives
        # us the most recently created figure.
        tissue.visualize()
        fig = plt.gcf()
        ax = fig.axes[0]

        # ``TissueSection.visualize`` adds its legend with an explicit
        # ``handles=`` argument (the Line2D proxies built from the color
        # map), which does NOT populate ``get_legend_handles_labels`` since
        # those Line2D objects aren't attached as plot artists.  Read the
        # binding straight off the Legend object instead.
        legend = ax.get_legend()
        assert legend is not None, "TissueSection.visualize did not draw a legend"

        binding = []
        for text, handle in zip(legend.get_texts(), legend.legend_handles):
            getter = getattr(handle, "get_markerfacecolor", None)
            color = getter() if getter is not None else handle.get_color()
            binding.append((text.get_text(), str(tuple(color))))

        print(json.dumps(binding))
    finally:
        plt.close("all")
    """
).strip()


def test_tissue_visualize_color_map_deterministic():
    """TissueSection.visualize's legend → color binding is HASHSEED-stable."""
    pytest.importorskip("matplotlib")
    out_a = _run_with_hashseed(_TISSUE_SCRIPT, "0")
    out_b = _run_with_hashseed(_TISSUE_SCRIPT, "1")
    assert out_a == out_b, (
        f"Visualization color bindings differ between PYTHONHASHSEED runs:\n"
        f"  seed=0 → {out_a}\n  seed=1 → {out_b}"
    )
    binding = json.loads(out_a)
    # Sanity: the three cell types should all be represented (order is
    # alphabetical because make_color_map sorts).
    labels = [label for label, _ in binding]
    assert labels == sorted(labels)
    assert set(labels) == {"cancer", "immune", "stroma"}
