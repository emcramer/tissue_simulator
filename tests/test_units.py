"""Persistent-homology unit detection (blobs, follicles, glomeruli, arteries)."""
import json
import time

import numpy as np
import pytest

from tissue_simulator.density import DensityModel
from tissue_simulator.packing import SpherePacker

SIZE = 300.0
NULLS = 9


def _points(seed):
    cells = SpherePacker((SIZE, SIZE, 0.0), {'c': (3.5, 5.0)}, min_spacing=0.3,
                         seed=seed).pack(max_attempts=300)
    return np.array([c.center[:2] for c in cells]), np.array([c.radius for c in cells])


def _distance(xy, centers):
    return np.linalg.norm(xy[:, None] - np.asarray(centers)[None], axis=2).min(axis=1)


def _fit(x, y, r, t, **kw):
    t0 = time.perf_counter()
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, SIZE, SIZE), strategy="adaptive",
                             n_null=NULLS, seed=0, **kw)
    model.fit_seconds = time.perf_counter() - t0
    return model


def _stroma(xy, rng, p_b=0.0):
    return np.where(rng.random(len(xy)) < p_b, 'B', 'A')


BLOBS = [(70, 70), (230, 80), (80, 230), (220, 225)]
FOLLICLES = [(70, 70), (225, 85), (150, 225)]
GLOMERULI = [(75, 75), (225, 80), (150, 225)]


@pytest.fixture(scope="module")
def blobs():
    xy, r = _points(1)
    rng = np.random.default_rng(1)
    t = np.where(_distance(xy, BLOBS) < 38, 'B', _stroma(xy, rng, 0.03))
    return _fit(xy[:, 0], xy[:, 1], r, t)


@pytest.fixture(scope="module")
def follicles():
    xy, r = _points(2)
    rng = np.random.default_rng(2)
    d = _distance(xy, FOLLICLES)
    t = np.where(d < 16, 'C', np.where(d < 38, 'B', _stroma(xy, rng, 0.03)))
    return _fit(xy[:, 0], xy[:, 1], r, t)


@pytest.fixture(scope="module")
def glomeruli():
    xy, r = _points(3)
    rng = np.random.default_rng(3)
    d = _distance(xy, GLOMERULI)
    keep = (d >= 18) & ((d < 40) | (rng.random(len(r)) < 1.1))
    t = np.where(d < 40, 'B', _stroma(xy, rng, 0.03))
    return _fit(xy[keep, 0], xy[keep, 1], r[keep], t[keep])


@pytest.fixture(scope="module")
def artery():
    xy, r = _points(4)
    rng = np.random.default_rng(4)
    d = _distance(xy, [(150, 150)])
    keep = d >= 30
    t = np.where(d < 65, 'B', _stroma(xy, rng, 0.03))
    return _fit(xy[keep, 0], xy[keep, 1], r[keep], t[keep])


def _check_units(model, kind, centers, radius, rtol=0.25, tol=10.0, core=None):
    u = model.units
    assert u, "no units detected"
    assert u["n_units"] == len(centers) and u["kinds"] == {kind: len(centers)}
    found = np.array([x["center"] for x in u["units"]])
    for c in centers:
        assert np.linalg.norm(found - np.array(c), axis=1).min() < tol
    for x in u["units"]:
        assert abs(x["outer_radius"] - radius) <= rtol * radius
        if core:
            assert x["core_type"] == core
    assert model.estimation["units"]["n_units"] == len(centers)


def test_blobs(blobs):
    _check_units(blobs, "blob", BLOBS, 38, tol=10)


def test_follicles(follicles):
    _check_units(follicles, "ring_core", FOLLICLES, 38, core='C')


def test_glomeruli(glomeruli):
    _check_units(glomeruli, "ring_lumen", GLOMERULI, 40)


def test_artery(artery):
    u = artery.units
    assert u and u["kinds"] == {"ring_lumen": 1}
    x = u["units"][0]
    assert np.linalg.norm(np.array(x["center"]) - 150) < 10
    assert abs(x["outer_radius"] - 65) <= 0.25 * 65 and abs(x["inner_radius"] - 30) <= 8


def test_uniform_control():
    xy, r = _points(5)
    t = np.random.default_rng(5).choice(['A', 'B', 'C'], len(r))
    model = _fit(xy[:, 0], xy[:, 1], r, t)
    assert model.units == {}


def test_gradient_control():
    xy, r = _points(6)
    rng = np.random.default_rng(6)
    t = np.where(rng.random(len(r)) < 0.1 + 0.8 * xy[:, 0] / SIZE, 'A', 'B')
    model = _fit(xy[:, 0], xy[:, 1], r, t)
    assert model.units == {}


def test_serialization_round_trip(blobs):
    data = json.loads(json.dumps(blobs.to_dict()))
    again = DensityModel.from_dict(data)
    assert again.units == blobs.units
    data.pop("units")
    assert DensityModel.from_dict(data).units == {}


def test_legacy_and_units_none_have_no_units(blobs):
    xy, r = _points(7)
    t = np.where(_distance(xy, BLOBS) < 38, 'B', 'A')
    legacy = DensityModel.fit(xy[:, 0], xy[:, 1], r, t, bounds=(0, 0, SIZE, SIZE),
                              n_null=0, seed=0)
    assert legacy.units == {} and "units" not in legacy.estimation
    off = DensityModel.fit(xy[:, 0], xy[:, 1], r, t, bounds=(0, 0, SIZE, SIZE),
                           strategy="adaptive", n_null=0, seed=0, units="none")
    assert off.units == {}
    with pytest.raises(ValueError):
        DensityModel.fit(xy[:, 0], xy[:, 1], r, t, units="bogus")


def test_fallback_without_null_draws():
    xy, r = _points(1)
    t = np.where(_distance(xy, BLOBS) < 38, 'B', 'A')
    m = DensityModel.fit(xy[:, 0], xy[:, 1], r, t, bounds=(0, 0, SIZE, SIZE),
                         strategy="adaptive", n_null=0, seed=0)
    assert m.units.get("null_draws", 0) == 0


NESTS = np.array([[70.0, 70.0], [220.0, 80.0], [80.0, 225.0], [215.0, 220.0]])


@pytest.fixture(scope="module")
def nests():
    """Full-density tumor nests in stroma thinned to 25 % (as in test_density)."""
    rng = np.random.default_rng(0)
    xy, r = _points(0)
    inside = (np.linalg.norm(xy[:, None] - NESTS[None], axis=2) < 45).any(axis=1)
    keep = inside | (rng.random(len(r)) < 0.25)
    t = np.where(inside, np.where(rng.random(len(r)) < 0.9, 'Tumor', 'CD8'),
                 np.where(rng.random(len(r)) < 0.8, 'Stroma', 'CD8'))
    return _fit(xy[keep, 0], xy[keep, 1], r[keep], t[keep])


def test_nest_region_has_exactly_four_units(nests):
    assert nests.units["n_units"] == 4
    assert nests.units["kinds"] == {"blob": 4}
    assert "barcodes" not in nests.units
    assert "barcode_top" in nests.estimation["units"]


def _potts_labels(xy, seed, j=0.5, sweeps=8, radius=14.0):
    """Gibbs sweeps of a 2-state Potts-like model on a 14 um neighbourhood."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    nbrs = cKDTree(xy).query_ball_point(xy, radius)
    lab = rng.integers(0, 2, len(xy))
    for _ in range(sweeps):
        for i in rng.permutation(len(xy)):
            nb = [k for k in nbrs[i] if k != i]
            same1 = lab[nb].sum()
            e = j * np.array([len(nb) - same1, same1])
            p = np.exp(e - e.max())
            lab[i] = rng.random() < p[1] / p.sum()
    return np.where(lab == 1, 'A', 'B')


@pytest.mark.parametrize("seed", [11, 12, 13])
def test_potts_clustering_control(seed):
    xy, r = _points(seed)
    model = _fit(xy[:, 0], xy[:, 1], r, _potts_labels(xy, seed))
    n = model.units.get("n_units", 0)
    print("potts", seed, n)
    assert n <= 1


def test_layered_trend_gives_no_units():
    """Residual structure only: an accepted planar trend must not yield units."""
    from tests.test_density import _layered_region
    x, y, r, t = _layered_region(2)
    model = DensityModel.fit(x, y, r, t, bounds=(0, 0, 240.0, 120.0), strategy="adaptive",
                             seed=0)
    assert model.organization["model"].startswith("planar")
    assert model.units == {}
