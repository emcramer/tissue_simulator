import numpy as np
import pytest

from tissue_simulator import SpherePacker
from tissue_simulator._voids import infer_voids, place_voids, summarize_components

STEP = 4.0


def _points():
    t = SpherePacker((200, 200, 0), {"c": (3.5, 5.0)}, min_spacing=0.3, seed=0).pack(max_attempts=300)
    cells = t.cells if hasattr(t, "cells") else t
    return (np.array([c.center[0] for c in cells]), np.array([c.center[1] for c in cells]),
            np.array([c.radius for c in cells]))


def _run(keep):
    x, y, r = _points()
    k = keep(x, y)
    return infer_voids(x[k], y[k], r[k], None, STEP, 200, 200)


def _near(x, y, c, rad):
    return np.hypot(x - c[0], y - c[1]) <= rad


def test_uniform_has_no_holes():
    _, info = _run(lambda x, y: np.ones(x.size, bool))
    assert info["n_holes"] == 0


def test_single_hole():
    mask, info = _run(lambda x, y: ~_near(x, y, (100, 100), 25))
    assert info["n_holes"] == 1
    assert abs(info["equivalent_diameters"][0] - 50) <= 10
    cx, cy = info["centroids"][0]
    assert np.hypot(cx - 100, cy - 100) < 8
    assert info["elongations"][0] < 1.3
    assert not mask[25, 25]


def test_two_holes_separation():
    _, info = _run(lambda x, y: ~(_near(x, y, (60, 60), 20) | _near(x, y, (140, 140), 20)))
    assert info["n_holes"] == 2
    assert info["min_center_separation"] == pytest.approx(113, abs=15) or \
        info["min_center_separation"] >= 90
    json_ok = __import__("json").dumps(info)
    assert json_ok


def test_edge_touching_excluded_from_pool():
    _, info = _run(lambda x, y: ~(_near(x, y, (0, 100), 40) | _near(x, y, (140, 100), 30)))
    assert info["n_edge_touching"] >= 1
    assert info["n_holes"] == 1
    assert len(info["equivalent_diameters"]) == 1
    assert not info["edge_touching"][0]


def _info(n, d=50.0, sep=None):
    return {"n_holes": n, "equivalent_diameters": [d], "min_center_separation": sep, "tau": 5.0}


def test_place_respects_separation_and_window():
    shape = (50, 50)
    tm = np.ones(shape, bool)
    out, p = place_voids(np.random.default_rng(1), shape, STEP, _info(3, 30.0, 60.0), tm)
    assert p["n_placed"] == 3 and tm.all() and not out.all()
    c = np.array(p["centers"])
    for (x, y), r in zip(c, p["radii"]):
        assert r <= x <= 200 - r and r <= y <= 200 - r
    for i in range(3):
        for j in range(i):
            assert np.hypot(*(c[i] - c[j])) >= 60


def test_place_anchor():
    info = {"n_holes": 2, "equivalent_diameters": [20.0, 40.0], "min_center_separation": 10.0, "tau": 5}
    _, p = place_voids(np.random.default_rng(2), (50, 50), STEP, info, np.ones((50, 50), bool),
                       anchor=(100, 100))
    assert p["anchored"] and p["centers"][0] == [100.0, 100.0]
    assert p["radii"][0] == max(p["radii"])


def test_place_shortfall():
    _, p = place_voids(np.random.default_rng(0), (15, 15), STEP, _info(3, 50.0, 50.0),
                       np.ones((15, 15), bool))
    assert p["n_placed"] < p["n_requested"] == 3


def test_place_reproducible():
    a = place_voids(np.random.default_rng(5), (50, 50), STEP, _info(3, 30.0), np.ones((50, 50), bool))
    b = place_voids(np.random.default_rng(5), (50, 50), STEP, _info(3, 30.0), np.ones((50, 50), bool))
    assert a[1] == b[1] and (a[0] == b[0]).all()
