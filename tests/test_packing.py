"""Tests for the spatial hash grid and the packers in ``tissue_simulator.packing``."""

import numpy as np
import pytest

from tissue_simulator import Cell
from tissue_simulator.density import Layout, RadiusMarks
from tissue_simulator.packing import (
    InhomogeneousPacker, SpatialHashGrid, SpherePacker, _stochastic_round,
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
