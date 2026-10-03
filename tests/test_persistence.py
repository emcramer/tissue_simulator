"""Tests for tissue_simulator._persistence (deterministic, small grids)."""
import numpy as np

from tissue_simulator._persistence import (
    barcode_summary, persistence_distance, superlevel_h0, superlevel_h1)

N = 61
YY, XX = np.mgrid[0:N, 0:N].astype(float)


def blobs(noise=0.0):
    f = np.zeros((N, N))
    for (r, c, h) in [(15, 15, 1.0), (15, 45, 0.7), (45, 30, 0.4)]:
        f += h * np.exp(-((YY - r) ** 2 + (XX - c) ** 2) / (2 * 4.0 ** 2))
    if noise:
        f = f + noise * np.random.default_rng(0).uniform(-1, 1, f.shape) * f.max()
    return f


def ring(cr, cc, R=12.0, w=3.0, shape_n=N):
    return np.exp(-(np.hypot(YY - cr, XX - cc) - R) ** 2 / (2 * w ** 2))


def test_flat_one_essential():
    bars = superlevel_h0(np.ones((10, 10)))
    assert len(bars) == 1 and bars[0]["persistence"] == 0.0
    assert bars[0]["area_at_death"] == 100
    assert superlevel_h1(np.ones((10, 10))) == []


def test_three_blobs():
    bars = superlevel_h0(blobs())
    top = bars[:3]
    assert [b["peak"] for b in top] == [(15, 15), (15, 45), (45, 30)]
    pers = [b["persistence"] for b in top]
    assert pers[0] > pers[1] > pers[2] > 0.3


def test_ring_h1():
    f = ring(30, 30)
    bars = superlevel_h1(f)
    assert len(bars) == 1
    b = bars[0]
    assert b["basin_min"] == (30, 30)
    ridge = min(f[30, 30 + 12], f[30 + 12, 30], f[30, 30 - 12], f[30 - 12, 30])
    assert abs(b["persistence"] - (ridge - f[30, 30])) < 0.05
    assert b["basin_area"] > 100


def test_two_rings():
    f = np.maximum(ring(15, 15, 8, 2), ring(45, 45, 8, 2))
    assert len([b for b in superlevel_h1(f) if b["persistence"] > 0.05]) == 2


def test_ring_touching_border():
    f = ring(0, 30, R=12)  # basin opens through the array edge
    assert superlevel_h1(f) == []


def test_mask():
    f = blobs()
    mask = np.ones_like(f, bool)
    mask[:30, 30:] = False  # removes the blob at (15, 45)
    peaks = [b["peak"] for b in superlevel_h0(f, mask)]
    assert (15, 45) not in peaks[:2] and peaks[0] == (15, 15)
    # ring whose basin reaches the mask edge is not counted
    g = ring(30, 30)
    m = np.ones_like(g, bool)
    m[:, 30:] = False
    assert superlevel_h1(g, m) == []
    assert len(superlevel_h1(g)) == 1


def test_distance():
    a, b = superlevel_h0(blobs()), superlevel_h0(blobs() * 0.5)
    assert persistence_distance(a, a) == 0.0
    d = persistence_distance(a, b)
    assert d > 0 and abs(d - persistence_distance(b, a)) < 1e-12
    s = barcode_summary(a, 2)
    assert s["count"] == len(a) and len(s["top"]) == 2


def test_noise_robustness():
    clean = superlevel_h0(blobs())
    noisy = superlevel_h0(blobs(noise=0.01))
    pm = clean[0]["persistence"]
    for c, n in zip(clean[:3], noisy[:3]):
        assert abs(n["persistence"] - c["persistence"]) / c["persistence"] < 0.05
    assert all(b["persistence"] < 0.02 * pm for b in noisy[3:])


def test_flat_peak_first_pixel():
    assert superlevel_h0(np.ones((5, 5)))[0]["peak"] == (0, 0)


def test_min_persistence_and_drop_essential():
    f = blobs(noise=0.01)
    kept = superlevel_h0(f, min_persistence=0.05)
    assert len(kept) == 3 and sum(b["essential"] for b in kept) == 1
    assert all(b["persistence"] >= 0.05 for b in superlevel_h1(ring(30, 30), min_persistence=0.05))
    a, b = superlevel_h0(blobs()), superlevel_h0(blobs() * 0.5)
    assert persistence_distance(a, b, drop_essential=True) < persistence_distance(a, b)
