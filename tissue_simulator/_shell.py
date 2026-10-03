"""Short-range (first-shell) spacing measures on size-normalised pair distances.

``s = d / (r_i + r_j)``; touching circles sit at ``s = 1``. All functions are
2-D (x, y), deterministic (no RNG) and use a KD-tree pair search, so memory is
O(pairs within ``s_max * 2 * r_max``).
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree

_SAMPLE_PAIRS = 20000
_PEAK_MARGIN = 0.05


def _pairs(points, radii, s_max):
    """Unordered pairs within ``s_max (r_i + r_j)``: ``(i, j, s)``."""
    empty = np.zeros(0, dtype=int)
    if len(points) < 2:
        return empty, empty, np.zeros(0)
    pairs = cKDTree(points).query_pairs(s_max * 2.0 * radii.max(), output_type="ndarray")
    if pairs.size == 0:
        return empty, empty, np.zeros(0)
    i, j = pairs[:, 0], pairs[:, 1]
    s = np.linalg.norm(points[i] - points[j], axis=1) / np.maximum(radii[i] + radii[j], 1e-12)
    keep = s < s_max
    return i[keep], j[keep], s[keep]


def pair_profile(points, radii, width, height, s_max=2.2, bin_width=0.05):
    """Binned size-normalised pair statistics of one window.

    Returns a dict with ``edges`` (s bin edges), ``pairs_per_cell`` (per bin,
    2 x unordered pairs / n: mean neighbours per cell, NOT edge-corrected, for
    comparing against a replicate in the same window) and ``g`` (edge-corrected
    size-normalised pair correlation).

    ``g_k = 2 C_k / (rho^2 A_k)``: ``C_k`` sums ``1 / ((W - |dx|)(H - |dy|))``
    over unordered pairs in bin k (translation correction), ``rho = n / (W H)``
    and ``A_k = pi (s_hi^2 - s_lo^2) <(r_i + r_j)^2>``, the mean over 20000
    deterministic random cell pairs (fixed-seed local generator, not global
    RNG), i.e. the annulus area in d of an uncorrelated pattern with the same
    radii. An uncorrelated pattern gives g ~ 1.
    """
    points = np.asarray(points, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n = len(points)
    edges = np.arange(0.0, s_max + 0.5 * bin_width, bin_width)
    nb = len(edges) - 1
    if n < 2:
        return {"edges": edges, "pairs_per_cell": np.zeros(nb), "g": np.zeros(nb)}
    i, j, s = _pairs(points, radii, s_max)
    k = np.minimum(np.searchsorted(edges, s, side="right") - 1, nb - 1)
    counts = np.bincount(k, minlength=nb).astype(float)
    dx = np.abs(points[i, 0] - points[j, 0])
    dy = np.abs(points[i, 1] - points[j, 1])
    w = 1.0 / np.maximum((width - dx) * (height - dy), 1e-6 * width * height)
    corrected = np.bincount(k, weights=w, minlength=nb)
    rng = np.random.default_rng(0)
    a, b = rng.integers(0, n, _SAMPLE_PAIRS), rng.integers(0, n, _SAMPLE_PAIRS)
    ok = a != b
    mean_sq = float(np.mean((radii[a[ok]] + radii[b[ok]]) ** 2))
    rho = n / (width * height)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2) * mean_sq
    g = 2.0 * corrected / (rho ** 2 * area)
    return {"edges": edges, "pairs_per_cell": 2.0 * counts / n, "g": g}


def shell_edge(profile, smooth=0.1, s_min=1.0, s_cap=2.1):
    """s of the first minimum after the first peak of the smoothed ``g``.

    The peak is searched over the whole range (in dense tissue with overlapping
    segmentation radii it can sit below s = 1); only the minimum must be at
    ``s >= s_min``. Returns None when there is no peak-then-minimum in ``[s_min, s_cap]`` or the
    minimum is not at least 5% below the peak (Poisson / RSA-like patterns).
    """
    edges = np.asarray(profile["edges"], dtype=float)
    g = np.asarray(profile["g"], dtype=float)
    centres = 0.5 * (edges[1:] + edges[:-1])
    width = float(edges[1] - edges[0])
    gs = gaussian_filter1d(g, smooth / width, mode="nearest")
    idx = np.flatnonzero(centres <= s_cap)
    if idx.size < 3:
        return None
    lo, hi = idx[0], idx[-1]
    peak = None
    for t in range(lo + 1, hi):
        if gs[t] > gs[t - 1] and gs[t] >= gs[t + 1]:
            peak = t
            break
    if peak is None:
        return None
    for t in range(peak + 1, hi):
        if centres[t] >= s_min and gs[t] < gs[t - 1] and gs[t] <= gs[t + 1]:
            if gs[t] <= (1.0 - _PEAK_MARGIN) * gs[peak]:
                return float(centres[t])
            return None
    return None


def first_shell_summary(points, radii, width, height, factor=1.5):
    """Scalar first-shell descriptors of one window (no edge correction)."""
    points = np.asarray(points, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n = len(points)
    if n == 0:
        return {"mean_degree": 0.0, "overlap_pairs_per_cell": 0.0, "median_radius": 0.0,
                "area_fraction": 0.0, "n_cells": 0}
    _, _, s = _pairs(points, radii, max(factor, 1.0))
    return {"mean_degree": float(2.0 * np.count_nonzero(s <= factor) / n),
            "overlap_pairs_per_cell": float(np.count_nonzero(s < 1.0) / n),
            "median_radius": float(np.median(radii)),
            "area_fraction": float(np.pi * np.sum(radii ** 2) / (width * height)),
            "n_cells": int(n)}


def first_shell_ratios(replicate_summary, source_summary):
    """Replicate / source ratio of each summary value (NaN when source is 0)."""
    out = {}
    for key in ("mean_degree", "overlap_pairs_per_cell", "median_radius", "area_fraction"):
        src = float(source_summary[key])
        out[key] = float(replicate_summary[key]) / src if src != 0.0 else float("nan")
    return out
