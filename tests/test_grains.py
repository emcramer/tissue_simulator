import numpy as np
import pytest

from tissue_simulator import SpherePacker
from tissue_simulator import _grains as G

CENTERS = [(70, 70), (220, 90), (150, 220)]


def _points():
    cs = SpherePacker((300, 300, 0), {"c": (3.5, 5.0)}, min_spacing=0.3, seed=0).pack()
    return np.array([[c.center[0], c.center[1], c.radius] for c in cs])


def _units(kind, inner, outer, n=3):
    return {"units": [{"kind": kind, "center": list(c), "outer_radius": outer,
                       "inner_radius": inner} for c in CENTERS[:n]]}


def _types(xy, core_r, ring_r, lumen=False):
    d = np.min([np.hypot(xy[:, 0] - c[0], xy[:, 1] - c[1]) for c in CENTERS], axis=0)
    t = np.zeros(len(xy), int)  # A stroma
    t[d < ring_r] = 1  # B ring
    if not lumen:
        t[d < core_r] = 2  # C core
    return t


def test_follicle_profiles():
    p = _points()
    ti = _types(p[:, :2], 16, 38)
    prof = G.fit_unit_profiles(p[:, 0], p[:, 1], ti, 3, p[:, 2],
                               _units("ring_core", 16, 38), 300, 300, 8)
    pr = prof["ring_core"]
    comp = np.array(pr["composition"])
    s = np.array(pr["s_knots"])
    assert np.allclose(comp.sum(1), 1)
    assert np.all(np.isfinite(pr["density_rel"]))
    assert np.all(comp[(s > 0.2) & (s < 0.42)].argmax(1) == 2)
    assert np.all(comp[(s > 0.42) & (s < 1)].argmax(1) == 1)
    assert pr["inner_ratio"] == pytest.approx(16 / 38)


def test_glomerulus_lumen_zero():
    p = _points()
    d = np.min([np.hypot(p[:, 0] - c[0], p[:, 1] - c[1]) for c in CENTERS], axis=0)
    keep = d >= 18
    ti = np.where(d[keep] < 40, 1, 0)
    prof = G.fit_unit_profiles(p[keep, 0], p[keep, 1], ti, 2, p[keep, 2],
                               _units("ring_lumen", 18, 40), 300, 300, 8)["ring_lumen"]
    s = np.array(prof["s_knots"])
    dens = np.array(prof["density_rel"])
    assert np.all(dens[s < 18 / 40] == 0)
    assert dens[(s > 0.5) & (s < 1)].min() > 0


SHAPE = (60, 60)


def _place(seed=1, **kw):
    args = dict(units_summary=_units("blob", 0, 20), min_separation=60.0)
    args.update(kw)
    return G.place_grains(np.random.default_rng(seed), SHAPE, 5.0, **args)


def test_place_separation_and_repro():
    a, b = _place(), _place()
    assert a == b
    assert a["n_placed"] == 3 == a["n_requested"]
    c = np.array(a["centers"])
    for i in range(3):
        for j in range(i):
            assert np.hypot(*(c[i] - c[j])) >= 60.0
    assert (c - 20 >= 0).all() and (c + 20 <= 300).all()


def test_place_shortfall_and_anchor():
    s = G.place_grains(np.random.default_rng(0), (20, 20), 5.0,
                       _units("blob", 0, 20), 90.0)
    assert s["n_placed"] < s["n_requested"]
    a = _place(anchor=(100.0, 120.0))
    assert a["anchored"] and a["centers"][0] == [100.0, 120.0]


def _profiles():
    k = [0.1, 0.5, 1.0, 1.4]
    return {"ring_lumen": {"s_knots": k, "density_rel": [0.0, 2.0, 2.0, 1.0],
                           "composition": [[0.0, 1.0]] * 4, "inner_ratio": 0.4},
            "blob": {"s_knots": k, "density_rel": [3.0] * 4,
                     "composition": [[1.0, 0.0]] * 4, "inner_ratio": 0.0}}


def test_rasterize():
    placed = {"centers": [[100.0, 100.0]], "outer_radii": [50.0], "inner_radii": [20.0],
              "kinds": ["ring_lumen"]}
    f, comp, g, v = G.rasterize_grains(placed, _profiles(), SHAPE, 5.0, [0.5, 0.5])
    py, px = np.mgrid[0:60, 0:60]
    s = np.hypot((px + .5) * 5 - 100, (py + .5) * 5 - 100) / 50
    assert np.array_equal(v, s < 0.4)
    assert (f[v] == 0).all()
    assert (f[~g] == 1).all() and not g[s > 1.5].any()
    assert np.allclose(comp.sum(0), 1)
    assert np.allclose(comp[:, ~g], 0.5)


def test_rasterize_nearest_grain():
    placed = {"centers": [[100.0, 100.0], [140.0, 100.0]], "outer_radii": [30.0, 30.0],
              "inner_radii": [0.0, 0.0], "kinds": ["blob", "ring_lumen"]}
    placed["inner_radii"][1] = 12.0
    f, comp, g, v = G.rasterize_grains(placed, _profiles(), SHAPE, 5.0, [0.5, 0.5])
    # pixel (row 20, col 25): center (127.5,102.5) is nearer grain 2 -> ring_lumen comp type 1
    assert comp[1, 20, 25] == pytest.approx(1.0)
    # pixel (row 20, col 18): (92.5,102.5) nearer grain 1 -> blob comp type 0
    assert comp[0, 20, 18] == pytest.approx(1.0)


def test_edge_unit_density_and_empty():
    p = _points()
    ti = np.zeros(len(p), int)
    def dens(c):
        u = {"units": [{"kind": "blob", "center": list(c), "outer_radius": 38.0,
                        "inner_radius": 0.0}]}
        pr = G.fit_unit_profiles(p[:, 0], p[:, 1], ti, 1, p[:, 2], u, 300, 300)["blob"]
        k = np.array(pr["s_knots"])
        return np.mean(np.array(pr["density_rel"])[k < 1])
    assert dens((20, 150)) == pytest.approx(dens((150, 150)), rel=0.15)
    assert G.fit_unit_profiles(p[:, 0], p[:, 1], ti, 1, p[:, 2], {}, 300, 300) == {}


def test_floor_and_shortfall_key():
    u = {"units": [{"kind": "blob", "center": [10.0, 10.0], "outer_radius": 5.0,
                    "inner_radius": 0.0}]}
    pr = G.fit_unit_profiles([200.0], [200.0], [0], 1, [3.0], u, 300, 300)["blob"]
    assert min(pr["density_rel"]) >= G.DEFAULT_N_KNOTS * 0 + G.DENSITY_FLOOR
    s = G.place_grains(np.random.default_rng(0), (20, 20), 5.0, _units("blob", 0, 20), 90.0)
    assert s["shortfall"] == s["n_requested"] - s["n_placed"] > 0
    a = _place(anchor=(-50.0, 1e4))
    assert a["centers"][0] == [0.0, 300.0]
