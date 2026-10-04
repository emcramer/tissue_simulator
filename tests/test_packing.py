"""Tests for the spatial hash grid and the packers in ``tissue_simulator.packing``."""

import numpy as np
import pytest

from tissue_simulator import Cell, _shell
from tissue_simulator.density import Layout, RadiusMarks
from tissue_simulator.packing import (
    InhomogeneousPacker, SpatialHashGrid, SpherePacker, _stochastic_round,
    separation_diagnostics,
)


BOUNDS = (120.0, 150.0, 20.0)

# Recorded from the exhaustive O(n^2) packer before the hash grid was added.
GOLDEN = {
    (0, True): dict(config={'a': (3.0, 6.0)}, n=352, checksum=53268.561613,
                    first=[40.468007065, 4.916822872, 0.330552711, 4.910885062],
                    last=[149.668430413, 48.938280667, 9.053423892, 3.039236415],
                    nboundary=225),
    (7, False): dict(config={'a': (2.0, 4.0), 'b': (5.0, 8.0)}, n=469,
                     checksum=68906.057986,
                     first=[112.111902598, 31.252078308, 9.077424252, 7.691641403],
                     last=[56.47422928, 111.051664595, 3.284418253, 2.072910671],
                     nboundary=0, types='bbaabbbabbab'),
}


def _brute_force_pack(packer, max_attempts):
    """The original exhaustive-collision RSA loop, kept as a reference."""
    cells, failed = [], 0
    while failed < max_attempts:
        cell_type, radius = packer._select_cell_type_and_radius()
        position = packer._generate_random_position(radius)
        candidate = Cell(center=position, radius=radius, cell_type=cell_type)
        if (not packer._is_valid_placement(candidate)
                or packer._check_collision(candidate, cells)):
            failed += 1
            continue
        candidate.is_boundary = not candidate.is_within_bounds(packer.bounds)
        cells.append(candidate)
        failed = 0
    return cells


def _as_array(cells):
    return np.array([[*c.center, c.radius] for c in cells])


def test_grid_neighbors_superset_of_true_neighbors():
    rng = np.random.default_rng(1)
    points = rng.uniform(0, 50, size=(300, 3))
    grid = SpatialHashGrid(4.0)
    for i, p in enumerate(points):
        grid.insert(i, p)
    for q in rng.uniform(0, 50, size=(40, 3)):
        near = set(np.flatnonzero(np.linalg.norm(points - q, axis=1) < 4.0))
        assert near <= set(grid.neighbors(q))


def test_grid_move_and_remove():
    grid = SpatialHashGrid(2.0)
    grid.insert(0, (1.0, 1.0, 1.0))
    grid.move(0, (1.0, 1.0, 1.0), (9.0, 9.0, 1.0))
    assert list(grid.neighbors((1.0, 1.0, 1.0))) == []
    assert list(grid.neighbors((9.5, 9.5, 1.0))) == [0]
    grid.remove(0, (9.0, 9.0, 1.0))
    assert list(grid.neighbors((9.5, 9.5, 1.0))) == []
    with pytest.raises(ValueError):
        SpatialHashGrid(0.0)


@pytest.mark.parametrize("seed,config,allow,spacing", [
    (3, {'a': (3.0, 6.0)}, True, 0.5),
    (11, {'a': (2.0, 4.0), 'b': (5.0, 8.0)}, False, 0.5),
    (5, {'x': (4.0, 4.5)}, True, 0.0),
    (8, {'x': (4.0, 5.0), 'y': (3.0, 4.0)}, True, -1.0),
])
def test_grid_pack_matches_exhaustive_reference(seed, config, allow, spacing):
    packed = SpherePacker(BOUNDS, config, min_spacing=spacing,
                          allow_boundary_cells=allow, seed=seed).pack(max_attempts=150)
    reference = _brute_force_pack(
        SpherePacker(BOUNDS, config, min_spacing=spacing,
                     allow_boundary_cells=allow, seed=seed), 150)
    assert len(packed) == len(reference)
    np.testing.assert_array_equal(_as_array(packed), _as_array(reference))
    assert [c.cell_type for c in packed] == [c.cell_type for c in reference]
    assert [c.is_boundary for c in packed] == [c.is_boundary for c in reference]


@pytest.mark.parametrize("key", list(GOLDEN))
def test_pack_matches_recorded_golden_output(key):
    seed, allow = key
    expected = GOLDEN[key]
    cells = SpherePacker(BOUNDS, expected['config'], min_spacing=0.5,
                         allow_boundary_cells=allow, seed=seed).pack(max_attempts=200)
    arr = _as_array(cells)
    assert len(cells) == expected['n']
    assert arr.sum() == pytest.approx(expected['checksum'], abs=1e-5)
    np.testing.assert_allclose(arr[0], expected['first'], atol=1e-8)
    np.testing.assert_allclose(arr[-1], expected['last'], atol=1e-8)
    assert sum(c.is_boundary for c in cells) == expected['nboundary']
    if 'types' in expected:
        assert ''.join(c.cell_type for c in cells[:12]) == expected['types']


def test_pack_with_progress_matches_pack():
    calls = []
    kwargs = dict(bounds=BOUNDS, cell_radii_config={'a': (3.0, 6.0)}, seed=2)
    plain = SpherePacker(**kwargs).pack(max_attempts=100)
    progress = SpherePacker(**kwargs).pack_with_progress(
        max_attempts=100, callback=lambda n, t: calls.append((n, t)))
    np.testing.assert_array_equal(_as_array(plain), _as_array(progress))
    assert calls and all(n % 10 == 0 for n, _ in calls)


def _step_layout(dense=0.0149, sparse=0.002, kappa=0.8, size=200.0, step=5.0):
    """Left half near hard-core saturation, right half sparse."""
    n = int(size / step)
    intensity = np.full((n, n), sparse)
    intensity[:, : n // 2] = dense
    marks = RadiusMarks(np.zeros(0), [np.linspace(3.5, 4.5, 11)])
    return Layout(width=size, height=size, grid_step=step, intensity=intensity,
                  composition=np.ones((1, n, n)), compartment=np.zeros((n, n), dtype=int),
                  cell_types=('A',), marks=marks, kappa=kappa,
                  n_target=int(round(intensity.sum() * step ** 2)),
                  target_overlap_fraction=0.02, mode='copy', bandwidth=20.0)


def _pack_layout(layout, seed=0, **kwargs):
    packer = InhomogeneousPacker((layout.height, layout.width, 0.0), layout, seed=seed, **kwargs)
    return packer.pack(), packer.report


def test_stochastic_round_is_exact_and_unbiased():
    counts = _stochastic_round(np.ones(100), 50, np.random.default_rng(0))
    assert counts.sum() == 50 and set(np.unique(counts)) <= {0, 1}
    assert counts[:50].sum() < 50


def test_inhomogeneous_packer_places_target_count_and_follows_layout():
    layout = _step_layout()
    cells, report = _pack_layout(layout)
    assert len(cells) == report.n_placed == layout.n_target
    assert report.n_rsa + report.n_inserted == report.n_placed
    assert report.bin_correlation > 0.9
    x = np.array([c.center[0] for c in cells])
    expected_left = layout.intensity[:, :20].sum() / layout.intensity.sum()
    assert np.mean(x < 100.0) == pytest.approx(expected_left, abs=0.03)
    assert report.max_displacement <= 0.5 * layout.marks.median_radius + 1e-9
    # The dense half is near hard-core saturation and relaxation is
    # displacement-capped, so some overlap remains by design.
    assert report.overlap_fraction < 0.2


def test_inhomogeneous_packer_is_seeded():
    layout = _step_layout()
    a, _ = _pack_layout(layout, seed=4)
    b, _ = _pack_layout(layout, seed=4)
    c, _ = _pack_layout(layout, seed=5)
    np.testing.assert_array_equal(_as_array(a), _as_array(b))
    assert not np.array_equal(_as_array(a), _as_array(c))


def test_inhomogeneous_packer_avoids_zero_intensity():
    layout = _step_layout(sparse=0.0)
    cells, _ = _pack_layout(layout)
    assert max(c.center[0] for c in cells) <= 100.0 + 0.5 * layout.marks.median_radius


def test_inhomogeneous_packer_keeps_cells_inside_when_requested():
    cells, _ = _pack_layout(_step_layout(dense=0.008), allow_boundary_cells=False)
    for c in cells:
        assert c.radius <= c.center[0] <= 200.0 - c.radius
        assert c.radius <= c.center[1] <= 200.0 - c.radius


def test_inhomogeneous_packer_rejects_mismatched_bounds():
    with pytest.raises(ValueError):
        InhomogeneousPacker((100.0, 200.0, 0.0), _step_layout())


def test_tissue_generate_cells_routes_layout_to_inhomogeneous_packer():
    from tissue_simulator import TissueSection

    layout = _step_layout(dense=0.008)
    tissue = TissueSection(200.0, 200.0, 0.0, {'A': (3.5, 4.5)}, seed=3)
    assert tissue.generate_cells(layout=layout) == layout.n_target
    assert tissue.packing_report.n_placed == layout.n_target


def test_report_legacy_quota_and_diagnostics():
    layout = _step_layout()
    _, report = _pack_layout(layout)
    assert report.bin_size == max(layout.bandwidth / 2, 10.0)
    assert report.quota_floor == report.bin_size and report.quota_floor_source == 'legacy'
    assert report.bin_shortfall.shape == report.bin_targets.shape
    assert set(report.clearance_quantiles) == {'p5', 'p50', 'p95'}
    assert 0.0 <= report.dense_bin_fraction_short <= 1.0
    d = report.to_dict()
    assert isinstance(d['bin_shortfall'], list) and isinstance(d['clearance_quantiles'], dict)


def test_layout_quota_scale_sets_bin_size():
    from types import SimpleNamespace
    layout = _step_layout()
    ns = SimpleNamespace(**{f: getattr(layout, f) for f in layout.__dataclass_fields__})
    ns.quota_scale = 15.0
    packer = InhomogeneousPacker((200.0, 200.0, 0.0), ns, seed=0)
    packer.pack()
    assert packer.bin_size == 15.0
    assert packer.report.quota_floor_source == 'layout'
    # explicit bin_size wins over the layout scale
    assert InhomogeneousPacker((200.0, 200.0, 0.0), ns, seed=0, bin_size=20.0).bin_size == 20.0


def test_overdense_layout_reports_shortfall_and_overlap():
    layout = _step_layout(dense=0.0149 * 4, sparse=0.002)
    _, report = _pack_layout(layout)
    assert report.dense_bin_fraction_short > 0
    assert report.clearance_quantiles['p50'] <= 0


def test_separation_diagnostics_two_cells():
    cells = [Cell(center=(0.0, 0.0, 0.0), radius=2.0, cell_type='a'),
             Cell(center=(9.0, 0.0, 0.0), radius=1.0, cell_type='a')]
    out = separation_diagnostics(cells)
    for k in ('p5', 'p50', 'p95'):
        assert out['clearance_quantiles'][k] == pytest.approx(6.0)
        assert out['normalized_distance_quantiles'][k] == pytest.approx(3.0)
    assert out['n_nearest_neighbour'] == 2
    assert separation_diagnostics(cells[:1])['n_nearest_neighbour'] == 0


def _wide_layout(strategy, dense_radii=None, sparse_radii=None, dense=0.012, sparse=0.003):
    """Crowded left half, sparse right half, two radius-mark bins."""
    layout = _step_layout(dense=dense, sparse=sparse)
    sparse_radii = np.linspace(2.0, 8.0, 41) if sparse_radii is None else sparse_radii
    dense_radii = np.linspace(2.0, 8.0, 41) if dense_radii is None else dense_radii
    layout.marks = RadiusMarks(np.array([0.5 * (dense + sparse)]), [sparse_radii, dense_radii])
    layout.strategy = strategy
    return layout


def test_radius_deck_matches_source_marginal_where_legacy_is_biased():
    pool = np.linspace(2.0, 8.0, 41)
    adaptive, rep_a = _pack_layout(_wide_layout("adaptive"))
    legacy, rep_l = _pack_layout(_wide_layout("legacy"))
    ra = np.array([c.radius for c in adaptive])
    rl = np.array([c.radius for c in legacy])
    grid = np.unique(pool)
    ks = np.max(np.abs(np.searchsorted(np.sort(ra), grid, side='right') / ra.size
                       - np.searchsorted(pool, grid, side='right') / pool.size))
    assert ks < 0.06
    assert np.median(ra) == pytest.approx(np.median(pool), rel=0.01)
    assert abs(np.median(ra) - np.median(pool)) < abs(np.median(rl) - np.median(pool))
    assert np.median(rl) < np.median(pool)
    assert rep_a.radius_assignment == "deck" and rep_l.radius_assignment == "per_candidate"
    assert rep_a.to_dict()["radius_assignment"] == "deck"


def test_radius_deck_keeps_density_conditioning():
    layout = _wide_layout("adaptive", dense_radii=np.linspace(2.0, 5.0, 31),
                          sparse_radii=np.linspace(5.0, 8.0, 31))
    cells, _ = _pack_layout(layout)
    left = [c.radius for c in cells if c.center[0] < 100.0]
    right = [c.radius for c in cells if c.center[0] >= 100.0]
    assert np.mean(left) < np.mean(right)


def test_radius_deck_is_seeded():
    layout = _wide_layout("adaptive")
    a, _ = _pack_layout(layout, seed=4)
    b, _ = _pack_layout(layout, seed=4)
    c, _ = _pack_layout(layout, seed=5)
    np.testing.assert_array_equal(_as_array(a), _as_array(b))
    assert not np.array_equal(_as_array(a), _as_array(c))


# -- first-shell refinement ---------------------------------------------------

_SIZE = 220.0
_RADIUS = 4.35


def _confluent_source(seed=0):
    """Jittered stretched-triangular lattice: neighbours sit near s = 1.15."""
    rng = np.random.default_rng(seed)
    pts = np.array([(10.0 * (i + 0.5 * (j % 2)) + 2.5, 8.8 * j + 4.4)
                    for j in range(25) for i in range(22)])
    pts = np.clip(pts + rng.normal(0.0, 0.8, pts.shape), 0.0, _SIZE)
    return pts, np.full(len(pts), _RADIUS)


def _refine_layout(strategy="adaptive", with_shell=True):
    pts, radii = _confluent_source()
    n = len(pts)
    grid = int(_SIZE / 5.0)
    marks = RadiusMarks(np.zeros(0), [np.linspace(0.95 * _RADIUS, 1.05 * _RADIUS, 11)])
    layout = Layout(width=_SIZE, height=_SIZE, grid_step=5.0,
                    intensity=np.full((grid, grid), n / _SIZE ** 2),
                    composition=np.ones((1, grid, grid)),
                    compartment=np.zeros((grid, grid), dtype=int), cell_types=('A',),
                    marks=marks, kappa=0.8, n_target=n, target_overlap_fraction=0.02,
                    mode='copy', bandwidth=20.0)
    layout.strategy = strategy
    if with_shell:
        prof = _shell.pair_profile(pts, radii, _SIZE, _SIZE)
        layout.shell = {"edges": prof["edges"].tolist(),
                        "pairs_per_cell": prof["pairs_per_cell"].tolist(),
                        "g": prof["g"].tolist(), "edge": _shell.shell_edge(prof),
                        "summary": _shell.first_shell_summary(pts, radii, _SIZE, _SIZE),
                        "s_floor": 0.8, "width": _SIZE, "height": _SIZE}
    return layout, pts, radii


def _xyzr(cells):
    return np.array([[*c.center, c.radius] for c in cells])


def _degree(cells):
    a = _xyzr(cells)
    return _shell.first_shell_summary(a[:, :2], a[:, 3], _SIZE, _SIZE)["mean_degree"]


def _pairs_below(cells, kappa):
    a = _xyzr(cells)
    d = np.linalg.norm(a[:, None, :3] - a[None, :, :3], axis=2)
    s = d / (a[:, None, 3] + a[None, :, 3])
    return int(np.count_nonzero(np.triu(s < kappa, 1)))


def test_refinement_moves_replicate_toward_source_shell():
    layout, pts, radii = _refine_layout()
    source = _shell.first_shell_summary(pts, radii, _SIZE, _SIZE)["mean_degree"]
    plain, rep_off = _pack_layout(layout, seed=3, refine_shell=False)
    refined, rep_on = _pack_layout(layout, seed=3)
    assert rep_off.refinement is None and rep_on.refinement is not None
    info = rep_on.refinement
    assert info["energy_after"] < info["energy_before"]
    assert 1 <= info["sweeps"] <= 20 and info["n_moved"] > 0
    assert 0.0 < info["accept_rate"] <= 1.0 and info["seconds"] >= 0.0
    assert info["max_displacement"] <= layout.marks.median_radius + 1e-9
    assert abs(_degree(refined) - source) < abs(_degree(plain) - source)
    assert rep_on.to_dict()["refinement"]["s_target"] == pytest.approx(2.0)
    assert info["taper_end"] == pytest.approx(layout.shell["edges"][-1] if layout.shell["edges"][-1] < 3.5 else 3.5)
    np.testing.assert_array_equal(rep_on.bin_achieved, rep_off.bin_achieved)
    assert rep_on.n_placed == rep_off.n_placed == layout.n_target


def test_refinement_respects_cap_floor_bins_and_z():
    layout, _, _ = _refine_layout()
    cap = 0.5 * layout.marks.median_radius
    plain, _ = _pack_layout(layout, seed=4, refine_shell=False)
    refined, rep = _pack_layout(layout, seed=4, refine_cap=cap)
    a, b = _xyzr(plain), _xyzr(refined)
    assert np.max(np.linalg.norm(b[:, :2] - a[:, :2], axis=1)) <= cap + 1e-9
    np.testing.assert_array_equal(a[:, 2:], b[:, 2:])
    assert _pairs_below(refined, layout.kappa) <= _pairs_below(plain, layout.kappa)
    assert np.all((b[:, :2] >= 0.0) & (b[:, :2] <= _SIZE))
    assert rep.refinement["max_displacement"] <= cap + 1e-9
    # Each cell stays in its own 10 um quota bin.
    assert np.array_equal(np.floor(a[:, :2] / 10.0), np.floor(b[:, :2] / 10.0))


def test_refinement_keeps_cells_out_of_voids():
    layout, _, _ = _refine_layout()
    layout.intensity[:, 20:24] = 0.0
    layout.n_target = int(round(layout.intensity.sum() * 25.0))
    refined, rep = _pack_layout(layout, seed=5)
    a = _xyzr(refined)
    assert rep.refinement is not None
    assert not np.any((a[:, 0] >= 100.0) & (a[:, 0] < 120.0))


def test_refinement_is_deterministic_and_off_switch_is_exact():
    layout, _, _ = _refine_layout()
    first, _ = _pack_layout(layout, seed=6)
    second, _ = _pack_layout(layout, seed=6)
    np.testing.assert_array_equal(_xyzr(first), _xyzr(second))
    off, rep = _pack_layout(layout, seed=6, refine_shell=False)
    stripped, _ = _pack_layout(_refine_layout(with_shell=False)[0], seed=6)
    np.testing.assert_array_equal(_xyzr(off), _xyzr(stripped))
    assert rep.refinement is None
    assert not np.array_equal(_xyzr(off), _xyzr(first))


def test_refinement_is_skipped_for_legacy_thick_small_or_per_candidate():
    layout, _, _ = _refine_layout()
    assert _pack_layout(_refine_layout(strategy="legacy")[0], seed=7)[1].refinement is None
    assert _pack_layout(layout, seed=7, radius_assignment="per_candidate")[1].refinement is None
    packer = InhomogeneousPacker((_SIZE, _SIZE, 3.0 * _RADIUS), layout, seed=7)
    packer.pack()
    assert packer.report.refinement is None
    small, _, _ = _refine_layout()
    small.n_target = 40
    assert _pack_layout(small, seed=7)[1].refinement is None
    forced = _pack_layout(layout, seed=7, radius_assignment="per_candidate", refine_shell=True)[1]
    assert forced.refinement is not None


def test_report_describes_refined_positions():
    layout, _, _ = _refine_layout()
    cells, rep = _pack_layout(layout, seed=8)
    a = _xyzr(cells)
    from scipy.spatial import cKDTree
    d, idx = cKDTree(a[:, :3]).query(a[:, :3], k=2)
    norm = d[:, 1] / (a[:, 3] + a[idx[:, 1], 3])
    assert rep.normalized_distance_quantiles["p50"] == pytest.approx(np.percentile(norm, 50))
    full = np.linalg.norm(a[:, None, :3] - a[None, :, :3], axis=2) / (a[:, None, 3] + a[None, :, 3])
    np.fill_diagonal(full, np.inf)
    assert rep.overlap_fraction == pytest.approx(np.mean(full.min(axis=1) < layout.kappa))


# -- first-shell report -------------------------------------------------------

_SHELL_KEYS = {"mean_degree", "overlap_pairs_per_cell", "median_radius", "area_fraction", "n_cells"}


@pytest.mark.parametrize("strategy", ["adaptive", "legacy"])
def test_first_shell_report_for_thin_slab(strategy):
    layout, _, _ = _refine_layout(strategy=strategy)
    cells, rep = _pack_layout(layout, seed=9)
    block = rep.first_shell
    assert block["factor"] == 1.5 and block["same_window"] is True
    assert set(block["replicate"]) == set(block["source"]) == _SHELL_KEYS
    assert set(block["ratios"]) == _SHELL_KEYS - {"n_cells"}
    a = _xyzr(cells)
    expected = _shell.first_shell_summary(a[:, :2], a[:, 3], _SIZE, _SIZE)
    assert block["replicate"] == expected
    assert block["ratios"]["mean_degree"] == pytest.approx(
        expected["mean_degree"] / layout.shell["summary"]["mean_degree"])
    import json
    assert json.loads(json.dumps(rep.to_dict()))["first_shell"]["same_window"] is True


def test_first_shell_report_none_for_thick_slab_or_no_summary():
    layout, _, _ = _refine_layout()
    packer = InhomogeneousPacker((_SIZE, _SIZE, 3.0 * _RADIUS), layout, seed=7)
    packer.pack()
    assert packer.report.first_shell is None and packer.report.to_dict()["first_shell"] is None
    assert _pack_layout(_refine_layout(with_shell=False)[0], seed=7)[1].first_shell is None


def test_first_shell_same_window_flag():
    layout, _, _ = _refine_layout()
    layout.shell = dict(layout.shell, width=2 * _SIZE)
    assert _pack_layout(layout, seed=9)[1].first_shell["same_window"] is False
    layout.shell = {k: v for k, v in layout.shell.items() if k not in ("width", "height")}
    assert _pack_layout(layout, seed=9)[1].first_shell["same_window"] is None
