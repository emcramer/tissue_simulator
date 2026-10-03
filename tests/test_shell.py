"""Tests for the first-shell spacing measures (tissue_simulator._shell)."""

import numpy as np
from scipy.spatial import Voronoi

from tissue_simulator import _shell


def _hex_lattice(a=10.0, nx=20, ny=20):
    pts = [(a * (i + 0.5 * (j % 2)), a * np.sqrt(3) / 2 * j) for j in range(ny) for i in range(nx)]
    pts = np.array(pts)
    return pts, nx * a, ny * a * np.sqrt(3) / 2


def _lloyd(pts, size, iters=15, seed=0):
    """Lloyd relaxation in a box using mirrored points to bound the cells."""
    for _ in range(iters):
        mirrors = [pts]
        for axis in (0, 1):
            for edge in (0.0, size):
                m = pts.copy()
                m[:, axis] = 2 * edge - m[:, axis]
                mirrors.append(m)
        allp = np.vstack(mirrors)
        vor = Voronoi(allp)
        new = pts.copy()
        for k in range(len(pts)):
            verts = vor.vertices[vor.regions[vor.point_region[k]]]
            new[k] = verts.mean(axis=0)  # vertex mean: cheap centroid proxy
        pts = np.clip(new, 0, size)
    return pts


def test_hex_lattice_profile_and_summary():
    pts, w, h = _hex_lattice(10.0, 40, 40)
    r = np.full(len(pts), 4.1)  # spacing 10 > 2r; s = 1.22 is off a bin edge
    prof = _shell.pair_profile(pts, r, w, h, bin_width=0.05)
    centres = 0.5 * (prof["edges"][1:] + prof["edges"][:-1])
    k = np.argmax(np.where(centres < 1.6, prof["pairs_per_cell"], 0))
    assert abs(centres[k] - 10.0 / 8.2) < 0.05
    n_int = np.count_nonzero(
        (pts[:, 0] > 20) & (pts[:, 0] < w - 20) & (pts[:, 1] > 20) & (pts[:, 1] < h - 20))
    assert n_int > 100
    assert 5.0 < prof["pairs_per_cell"][k] <= 6.0
    summ = _shell.first_shell_summary(pts, r, w, h, factor=1.5)
    assert 5.0 < summ["mean_degree"] <= 6.0
    assert summ["overlap_pairs_per_cell"] == 0.0
    assert summ["n_cells"] == len(pts)
    tight = _shell.first_shell_summary(pts, np.full(len(pts), 6.0), w, h)
    assert tight["overlap_pairs_per_cell"] > 2.0


def test_relaxed_dense_pattern_has_shell_edge():
    rng = np.random.default_rng(3)
    size, n = 400.0, 900
    pts = _lloyd(rng.random((n, 2)) * size, size)
    pts = np.clip(pts + rng.normal(0, 0.8, pts.shape), 0, size)
    spacing = np.sqrt(size * size / n)
    r = np.full(n, 0.55 * spacing)
    prof = _shell.pair_profile(pts, r, size, size)
    edge = _shell.shell_edge(prof)
    assert edge is not None and 1.2 <= edge <= 1.8


def test_poisson_has_no_shell_edge_and_g_near_one():
    rng = np.random.default_rng(1)
    size, n = 300.0, 1500
    pts = rng.random((n, 2)) * size
    r = rng.uniform(2.0, 3.0, n)
    prof = _shell.pair_profile(pts, r, size, size)
    assert _shell.shell_edge(prof) is None
    assert abs(np.mean(prof["g"][prof["edges"][:-1] >= 0.5]) - 1.0) < 0.15


def test_ratios_of_self_are_one_and_nan_safe():
    pts, w, h = _hex_lattice()
    summ = _shell.first_shell_summary(pts, np.full(len(pts), 4.0), w, h)
    ratios = _shell.first_shell_ratios(summ, summ)
    assert ratios["mean_degree"] == 1.0 and ratios["median_radius"] == 1.0
    assert ratios["area_fraction"] == 1.0
    assert np.isnan(ratios["overlap_pairs_per_cell"])  # source 0
