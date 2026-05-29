"""Internal matplotlib-visualization helpers.

Centralizes the cell-type → color mapping used by ``TissueSection.visualize``,
``TissueSlicer.visualize_slice_2d`` / ``visualize_slice_in_3d``,
``SpatialNetworkAnalyzer.visualize_network``, ``visualize_colored_graph``, and
``visualize_graph_comparison`` so the assignment is deterministic regardless
of ``PYTHONHASHSEED`` or input iteration order.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np


def make_color_map(
    cell_types: Iterable[str],
    palette: Optional[Sequence] = None,
) -> Dict[str, Tuple[float, float, float, float]]:
    """Return a deterministic ``{cell_type: rgba_color}`` mapping.

    Cell types are deduplicated and sorted alphabetically before assignment,
    so the same set of types always yields the same mapping regardless of
    the order they were iterated in (which CPython's set/dict iteration
    leaks via the per-interpreter hash randomization seed).

    Args:
        cell_types: Iterable of cell-type names. Duplicates are ignored.
        palette: Optional pre-built sequence of RGBA colors. When omitted,
            the function uses ``matplotlib.pyplot.cm.tab10`` sampled at
            evenly spaced points. When provided, the first ``len(types)``
            entries are zipped with the sorted types.

    Returns:
        Ordered mapping (sorted by key) from cell type to RGBA tuple.
    """
    types_sorted = sorted(set(cell_types))
    n = max(len(types_sorted), 1)
    if palette is None:
        import matplotlib.pyplot as plt  # imported lazily so the package import doesn't pull matplotlib

        palette = plt.cm.tab10(np.linspace(0, 1, n))
    else:
        palette = list(palette)
        if len(palette) < n:
            raise ValueError(
                f"palette has {len(palette)} entries but {n} are needed for cell types {types_sorted}"
            )
    return {ct: tuple(palette[i]) for i, ct in enumerate(types_sorted)}
