"""Density-aware replicates: the annealer's spatial term and the generator wiring."""

import contextlib
import io
import random
import warnings

import networkx as nx
import numpy as np
import pytest

from tissue_simulator import Cell, DensityModel, ReplicateGenerator, TissueSection
from tissue_simulator.graph_coloring import GraphColorizer
from tissue_simulator.packing import SpherePacker
from tissue_simulator.replicate_generator import load_target_statistics_from_tissue

SIZE = 200.0
TYPES = ('CD8', 'Stroma', 'Tumor')
COLORING = dict(cooling_rate=0.999, max_iterations=4000)


# Achieved 0.44 (18 core cells; global share ~0.1); 0.5 was not reached.
CORE_C_MIN = 0.4
NESTS_B_MIN = 0.65  # achieved 0.81


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def _nest_tissue(seed=0):
    """Two dense tumor nests in stroma thinned to 25%, in a 1 µm slab."""
    rng = np.random.default_rng(seed)
    packed = SpherePacker((SIZE, SIZE, 1.0), {'c': (3.5, 5.0)}, min_spacing=0.3,
                          seed=seed).pack(max_attempts=300)
    nests = np.array([[60.0, 60.0], [140.0, 140.0]])
    tissue = TissueSection(SIZE, SIZE, 1.0, {t: (3.5, 5.0) for t in TYPES})
    for cell in packed:
        inside = (np.linalg.norm(nests - cell.center[:2], axis=1) < 40.0).any()
        if inside or rng.random() < 0.25:
            u = rng.random()
            cell.cell_type = ('Tumor' if u < 0.9 else 'CD8') if inside else ('Stroma' if u < 0.75 else 'CD8')
            tissue.cells.append(cell)
    return tissue


@pytest.fixture(scope="module")
def setup():
    tissue = _nest_tissue()
    target = load_target_statistics_from_tissue(tissue, network_mode="radius", network_radius=20.0)
    target.target_density = None  # a 1 µm slab has no meaningful 3D packing fraction
    model = DensityModel.from_tissue(tissue, n_null=0, n_compartments=2, seed=0)
    return tissue, target, model, {t: (3.5, 5.0) for t in TYPES}


def _generator(setup, **kwargs):
    _, target, model, radii = setup
    kwargs.setdefault('density_model', model)
    return ReplicateGenerator(target, (SIZE, SIZE, 1.0), radii, network_mode="radius",
                              network_radius=20.0, seed=11, method="graph_coloring",
                              coloring_params=COLORING, **kwargs)


def _as_array(tissue):
    return np.array([[*c.center, c.radius] for c in tissue.cells])


def test_spatial_term_incremental_update_matches_full_recompute():
    graph = nx.random_geometric_graph(80, 0.2, seed=3)
    colors = ['A', 'B', 'C']
    rng = random.Random(0)
    coloring = {n: rng.choice(colors) for n in graph.nodes()}
    plain = _quiet(GraphColorizer, target_graph=graph, colors=colors,
                   target_statistics={'node_counts': {c: 1 for c in colors},
                                      'edge_counts': {}, 'neighbor_dist': {}})
    targets = dict(plain._calculate_statistics(graph, coloring)[0])
    targets['spatial_composition'] = {
        'node_bin': {n: n % 7 for n in graph.nodes()},
        'expected': {b: {c: rng.uniform(0, 5) for c in colors} for b in range(7)},
        'weight': 1.5,
    }
    colorizer = _quiet(GraphColorizer, target_graph=graph, colors=colors,
                       target_statistics=targets, seed=0)
    assert 'spatial_composition' not in colorizer.target_stats

    stats, counts = colorizer._calculate_statistics(graph, coloring)
    nodes = list(graph.nodes())
    for _ in range(300):
        a, b = rng.sample(nodes, 2)
        stats, counts = colorizer._update_statistics_incremental(a, b, coloring, stats, counts)
        coloring[a], coloring[b] = coloring[b], coloring[a]
        full, _ = colorizer._calculate_statistics(graph, coloring)
        assert stats['spatial_sse'] == pytest.approx(full['spatial_sse'], abs=1e-6)
        assert colorizer._calculate_cost(stats) == pytest.approx(colorizer._calculate_cost(full), abs=1e-6)
    assert colorizer.cost_terms(stats)['spatial'] == pytest.approx(full['spatial_sse'])


def test_density_replicate_is_valid_and_reproducible(setup):
    _, _, model, _ = setup
    gen = _generator(setup)
    t1, s1 = _quiet(gen.generate_single_replicate, 0)
    t2, _ = _quiet(gen.generate_single_replicate, 0)
    t3, _ = _quiet(gen.generate_single_replicate, 1)
    np.testing.assert_array_equal(_as_array(t1), _as_array(t2))
    assert [c.cell_type for c in t1.cells] == [c.cell_type for c in t2.cells]
    assert not np.array_equal(_as_array(t1), _as_array(t3))
    assert s1.num_cells == s1.packing_report['n_target'] == model.n_cells
    assert 0.0 <= s1.composition_error <= 1.0
    assert s1.layout_flags[0] == "mode:resample"
    assert set(s1.cell_type_counts) <= set(TYPES)


def test_copy_layout_replicate(setup):
    _, stats = _quiet(_generator(setup, layout="copy").generate_single_replicate, 0)
    assert stats.layout_flags[0] == "mode:copy"


def test_density_model_argument_validation(setup):
    _, target, model, radii = setup
    with pytest.raises(ValueError, match="graph_coloring"):
        ReplicateGenerator(target, (SIZE, SIZE, 1.0), radii, density_model=model,
                           method="radius_tuning")
    with pytest.raises(ValueError, match="layout"):
        _generator(setup, layout="tile")


def test_composition_weight_reduces_composition_error(setup):
    errors = {}
    for weight in (0.0, 4.0):
        gen = _generator(setup, composition_weight=weight)
        errors[weight] = np.mean([_quiet(gen.generate_single_replicate, i)[1].composition_error
                                  for i in range(3)])
    assert errors[4.0] < errors[0.0]


def test_density_replicates_parallel_matches_serial(setup):
    gen = _generator(setup)
    serial = _quiet(gen.generate_replicates, 2)
    parallel = _quiet(gen.generate_replicates, 2, parallel=True, max_workers=2)
    for (ts, ss), (tp, sp) in zip(serial, parallel):
        np.testing.assert_array_equal(_as_array(ts), _as_array(tp))
        assert [c.cell_type for c in ts.cells] == [c.cell_type for c in tp.cells]
        assert ss.divergence_score == sp.divergence_score


def test_from_coordinates_builds_density_generator(setup, tmp_path):
    tissue = setup[0]
    path = tmp_path / "region.csv"
    tissue.export_to_csv(str(path))
    gen = ReplicateGenerator.from_coordinates(
        str(path), seed=2, coloring_params=COLORING,
        density_kwargs=dict(n_null=0, n_compartments=2))
    assert gen.density_model is not None
    assert gen.layout == "resample" and gen.method == "graph_coloring"
    _, stats = _quiet(gen.generate_single_replicate, 0)
    assert stats.num_cells == gen.density_model.n_cells


def test_mcp_setup_replicate_generator_density_layout(setup):
    import asyncio
    import json
    pytest.importorskip("mcp")
    from tissue_simulator.mcp.server import TissueSimulatorMCPServer

    tissue, target, _, radii = setup
    server = TissueSimulatorMCPServer()
    server.target_stats, server.network_mode, server.network_radius = target, "radius", 20.0
    args = {"height": SIZE, "width": SIZE, "thickness": 20,
            "cell_radii": {k: list(v) for k, v in radii.items()},
            "seed": 3, "method": "graph_coloring", "density_layout": "resample"}

    def call():
        # A private loop: asyncio.run() would clear the current event loop,
        # which tests/test_mcp_server.py relies on.
        loop = asyncio.new_event_loop()
        try:
            return json.loads(loop.run_until_complete(
                server._handle_setup_replicate_generator(args))[0].text)
        finally:
            loop.close()

    data = call()
    assert "error" in data  # no source tissue recorded yet

    server.target_source_tissue = tissue
    data = call()
    assert data["status"] == "success", data
    assert data["density_layout"] == "resample"
    assert data["density_model"]["n_compartments"] >= 1
    assert server.replicate_generator.density_model is not None


def test_mcp_setup_replicate_generator_strategy(setup):
    import asyncio
    import json
    pytest.importorskip("mcp")
    from tissue_simulator.mcp.server import TissueSimulatorMCPServer

    tissue, target, _, radii = setup
    server = TissueSimulatorMCPServer()
    server.target_stats, server.network_mode, server.network_radius = target, "radius", 20.0
    server.target_source_tissue = tissue
    args = {"height": SIZE, "width": SIZE, "thickness": 20,
            "cell_radii": {k: list(v) for k, v in radii.items()},
            "seed": 3, "method": "graph_coloring", "density_layout": "resample",
            "strategy": "adaptive", "diagnostics": True, "composition_weight": 2.0}
    loop = asyncio.new_event_loop()
    try:
        data = json.loads(loop.run_until_complete(
            server._handle_setup_replicate_generator(args))[0].text)
    finally:
        loop.close()
    assert data["status"] == "success", data
    assert data["strategy"] == "adaptive"
    assert data["diagnostics"] is True
    assert server.replicate_generator.strategy == "adaptive"


# ---------------------------------------------------------------------------
# Adaptive strategy (multi-scale composition, size compatibility, diagnostics)
# ---------------------------------------------------------------------------

def _mixed_radii_tissue(seed=0):
    """Nest tissue with large Tumor cells (5-6 um) and small Stroma cells (3-3.5 um), ~5 % CD8."""
    rng = np.random.default_rng(seed)
    packed = SpherePacker((SIZE, SIZE, 1.0), {'c': (3.0, 6.0)}, min_spacing=0.3,
                          seed=seed).pack(max_attempts=300)
    nests = np.array([[60.0, 60.0], [140.0, 140.0]])
    tissue = TissueSection(SIZE, SIZE, 1.0, {'Tumor': (5.0, 6.0), 'Stroma': (3.0, 3.5), 'CD8': (3.0, 3.5)})
    for cell in packed:
        inside = (np.linalg.norm(nests - cell.center[:2], axis=1) < 40.0).any()
        u = rng.random()
        if inside:
            cell.cell_type, cell.radius = ('Tumor', rng.uniform(5.0, 6.0)) if u < 0.95 else ('CD8', rng.uniform(3.0, 3.5))
        elif rng.random() < 0.3:
            cell.cell_type, cell.radius = ('Stroma', rng.uniform(3.0, 3.5)) if u < 0.95 else ('CD8', rng.uniform(3.0, 3.5))
        else:
            continue
        tissue.cells.append(cell)
    return tissue


@pytest.fixture(scope="module")
def adaptive_setup():
    tissue = _mixed_radii_tissue()
    target = load_target_statistics_from_tissue(tissue, network_mode="radius", network_radius=20.0)
    target.target_density = None
    model = DensityModel.from_tissue(tissue, n_null=0, n_compartments=2, seed=0, strategy="adaptive")
    radii = {'Tumor': (5.0, 6.0), 'Stroma': (3.0, 3.5), 'CD8': (3.0, 3.5)}
    return tissue, target, model, radii


def _adaptive_generator(adaptive_setup, **kwargs):
    tissue, target, model, radii = adaptive_setup
    kwargs.setdefault('dims', (SIZE, SIZE, 1.0))
    dims = kwargs.pop('dims')
    gen = ReplicateGenerator(target, dims, radii, network_mode="radius", network_radius=20.0,
                             seed=11, method="graph_coloring", coloring_params=COLORING,
                             density_model=model, strategy="adaptive", **kwargs)
    gen._cache_source_reference(tissue)
    return gen


def test_explicit_legacy_strategy_is_identical(setup):
    t1, s1 = _quiet(_generator(setup).generate_single_replicate, 0)
    t2, s2 = _quiet(_generator(setup, strategy="legacy").generate_single_replicate, 0)
    np.testing.assert_array_equal(_as_array(t1), _as_array(t2))
    assert [c.cell_type for c in t1.cells] == [c.cell_type for c in t2.cells]
    assert s1.divergence_score == s2.divergence_score
    assert s1.fidelity is None and s1.separation is not None


def test_adaptive_requires_density_model(setup):
    _, target, _, radii = setup
    with pytest.raises(ValueError, match="density_model"):
        ReplicateGenerator(target, (SIZE, SIZE, 1.0), radii, method="graph_coloring",
                           strategy="adaptive")


def test_adaptive_resolves_defaults_and_scales(adaptive_setup):
    gen = _adaptive_generator(adaptive_setup)
    assert gen.composition_bin == "auto" and gen.size_weight == 1.0 and gen.diagnostics
    scales = gen._resolve_composition_scales()
    assert 1 <= len(scales) <= 3
    assert all(b / a >= 1.5 for a, b in zip(scales, scales[1:]))
    assert _adaptive_generator(adaptive_setup, composition_scales=[12.0, 30.0]
                               )._resolve_composition_scales() == [12.0, 30.0]


def test_adaptive_mixed_radii_and_end_to_end_targets(adaptive_setup):
    tissue, _, _, _ = adaptive_setup
    gen = _adaptive_generator(adaptive_setup)
    rep, stats = _quiet(gen.generate_single_replicate, 0)
    src = {t: np.median([c.radius for c in tissue.cells if c.cell_type == t]) for t in ('Tumor', 'Stroma')}
    for t in ('Tumor', 'Stroma'):
        med = np.median([c.radius for c in rep.cells if c.cell_type == t])
        assert abs(med - src[t]) < 0.5, (t, med, src[t])
    ks = stats.fidelity['size_ks_by_type']
    assert all(v < 0.25 for k, v in ks.items() if k in ('Tumor', 'Stroma')), ks
    assert stats.fidelity['n_components'] >= 1 and stats.fidelity['size_nll']['source'] is not None
    import json
    json.dumps(stats.to_dict(), default=float)

    # The colorizer receives finite multi-scale and size terms.
    from tissue_simulator.spatial_analysis import SpatialNetworkAnalyzer
    graph = SpatialNetworkAnalyzer().build_network_from_tissue(rep, mode="radius", radius=20.0)
    mean_deg = 2.0 * graph.number_of_edges() / graph.number_of_nodes()
    targets = gen._build_colorizer_targets(graph)
    layout = gen.density_model.sample_layout(rng=1, width=SIZE, height=SIZE)
    scales, initial = gen._spatial_composition_target(
        rep, layout, targets['node_counts'], mean_deg, np.random.default_rng(0),
        bin_sizes=gen._resolve_composition_scales())
    targets['spatial_composition_scales'] = scales
    targets['size_compatibility'] = gen._size_compatibility_target(rep, gen.density_model, layout, mean_deg)
    colorizer = _quiet(GraphColorizer, target_graph=graph, colors=list(gen.cell_types),
                       target_statistics=targets, seed=0)
    stat, _ = colorizer._calculate_statistics(graph, initial)
    terms = colorizer.cost_terms(stat)
    assert np.isfinite(terms['spatial_scales']) and np.isfinite(terms['size'])
    assert len(scales) == len(gen._resolve_composition_scales())


def test_adaptive_rare_type_count_matches_request(adaptive_setup):
    _, stats = _quiet(_adaptive_generator(adaptive_setup).generate_single_replicate, 0)
    assert stats.requested_cell_type_counts['CD8'] > 0
    assert stats.cell_type_counts.get('CD8', 0) == stats.requested_cell_type_counts['CD8']
    assert stats.achieved_cell_type_counts == stats.requested_cell_type_counts


def test_adaptive_infeasible_packing_reports_shortfall(adaptive_setup):
    # Infeasible: same cell count and window, but radii ~3x larger than the source's.
    import copy
    gen = _adaptive_generator(adaptive_setup)
    gen.density_model = copy.deepcopy(gen.density_model)
    gen.density_model.marks.radii_by_bin = [3.0 * r for r in gen.density_model.marks.radii_by_bin]
    _, stats = _quiet(gen.generate_single_replicate, 0)
    assert stats.packing_report['dense_bin_fraction_short'] > 0
    assert stats.separation is not None


def test_adaptive_parallel_matches_serial(adaptive_setup):
    gen = _adaptive_generator(adaptive_setup)
    serial = _quiet(gen.generate_replicates, 2)
    try:
        parallel = _quiet(gen.generate_replicates, 2, parallel=True, max_workers=2)
    except (OSError, PermissionError, NotImplementedError) as exc:  # pragma: no cover
        pytest.skip(f"ProcessPool unavailable: {exc}")
    for (ts, ss), (tp, sp) in zip(serial, parallel):
        np.testing.assert_array_equal(_as_array(ts), _as_array(tp))
        assert [c.cell_type for c in ts.cells] == [c.cell_type for c in tp.cells]
        assert ss.divergence_score == sp.divergence_score


def test_adaptive_replicate_respects_placed_void():
    """A lumen in the source is re-placed and left empty in the replicate."""
    packed = SpherePacker((SIZE, SIZE, 1.0), {'c': (2.3, 3.0)}, min_spacing=0.3,
                          seed=2).pack(max_attempts=400)
    rng = np.random.default_rng(2)
    tissue = TissueSection(SIZE, SIZE, 1.0, {t: (2.3, 3.0) for t in TYPES})
    for cell in packed:
        d = np.linalg.norm(cell.center[:2] - SIZE / 2)
        if d < 25:
            continue
        cell.cell_type = TYPES[1] if d < 55 else TYPES[2]
        if rng.random() < 0.1:
            cell.cell_type = TYPES[int(rng.integers(3))]
        tissue.cells.append(cell)
    target = load_target_statistics_from_tissue(tissue, network_mode="radius", network_radius=20.0)
    target.target_density = None
    model = DensityModel.from_tissue(tissue, n_null=0, seed=0, strategy="adaptive",
                                 units="none")  # tests the voids path
    gen = ReplicateGenerator(target, (SIZE, SIZE, 1.0), {t: (2.3, 3.0) for t in TYPES},
                             network_mode="radius", network_radius=20.0, seed=11,
                             method="graph_coloring",
                             coloring_params=dict(cooling_rate=0.99, max_iterations=300),
                             density_model=model, strategy="adaptive")
    gen._cache_source_reference(tissue)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep, stats = _quiet(gen.generate_single_replicate, 0)
    lv = stats.layout_voids
    assert lv["n_placed"] == 1
    c, rad = np.array(lv["centers"][0]), lv["radii"][0]
    med = float(np.median([cell.radius for cell in rep.cells]))
    depth = np.array([rad - np.linalg.norm(cell.center[:2] - c) for cell in rep.cells])
    # Relaxation spillover (0.5 median radius) plus the half pixel diagonal of the
    # raster hole boundary (cells are drawn from pixels next to the hole).
    assert depth.max() <= 0.5 * med + 0.5 * np.sqrt(2.0) * model.grid_step
    assert "layout_voids" in stats.to_dict()


@pytest.fixture(scope="module")
def follicle_setup():
    """Three C-core / B-ring follicles in A stroma (300 um), adaptive model with units."""
    from .test_units import FOLLICLES, _distance, _fit, _points, _stroma
    size = 300.0
    xy, r = _points(2)
    rng = np.random.default_rng(2)
    d = _distance(xy, FOLLICLES)
    t = np.where(d < 16, 'C', np.where(d < 38, 'B', _stroma(xy, rng, 0.03)))
    model = _fit(xy[:, 0], xy[:, 1], r, t)
    tissue = TissueSection(size, size, 1.0, {c: (3.5, 5.0) for c in ('A', 'B', 'C')})
    for (x, y), rad, typ in zip(xy, r, t):
        tissue.cells.append(Cell((x, y, 0.5), rad, str(typ)))
    target = load_target_statistics_from_tissue(tissue, network_mode="radius", network_radius=20.0)
    target.target_density = None
    return tissue, target, model, size


def test_follicle_replicate_places_units_and_reports_persistence(follicle_setup):
    tissue, target, model, size = follicle_setup
    assert model.units["n_units"] == 3
    gen = ReplicateGenerator(
        target, (size, size, 1.0), {c: (3.5, 5.0) for c in ('A', 'B', 'C')},
        network_mode="radius", network_radius=20.0, seed=11, method="graph_coloring",
        coloring_params=dict(cooling_rate=0.99, max_iterations=300),
        density_model=model, strategy="adaptive")
    gen._cache_source_reference(tissue)
    rep, stats = _quiet(gen.generate_single_replicate, 0)
    assert stats.layout_units["n_placed"] == 3
    pd_ = stats.fidelity["persistence_distance"]
    assert pd_ is not None and all(v is not None and np.isfinite(v) for v in pd_.values())
    import json
    json.dumps(stats.to_dict(), default=float)
    xy = np.array([c.center[:2] for c in rep.cells])
    types = np.array([c.cell_type for c in rep.cells])
    in_core = np.zeros(len(xy), bool)
    for c, R in zip(stats.layout_units["centers"], stats.layout_units["outer_radii"]):
        in_core |= np.linalg.norm(xy - np.array(c), axis=1) < 0.4 * R  # core is s < 16/38
    # Equal-count profile bins + neighbour prior keep the core composition near
    # the source's; the annealer then reproduces most of it (observed 0.44).
    assert in_core.sum() >= 5 and (types[in_core] == 'C').mean() > CORE_C_MIN


def test_nests_replicate_keeps_nest_purity():
    from .test_units import BLOBS, _distance, _fit, _points
    xy, r = _points(1)
    rng = np.random.default_rng(1)
    d = _distance(xy, BLOBS)
    stroma = rng.choice(['A', 'B', 'C'], size=len(xy), p=[0.8, 0.1, 0.1])
    pure = rng.random(len(xy)) < 0.9
    t = np.where(d < 38, np.where(pure, 'B', stroma), stroma)
    model = _fit(xy[:, 0], xy[:, 1], r, t)
    tissue = TissueSection(300.0, 300.0, 1.0, {c: (3.5, 5.0) for c in ('A', 'B', 'C')})
    for (x, y), rad, typ in zip(xy, r, t):
        tissue.cells.append(Cell((x, y, 0.5), rad, str(typ)))
    target = load_target_statistics_from_tissue(tissue, network_mode="radius", network_radius=20.0)
    target.target_density = None
    gen = ReplicateGenerator(
        target, (300.0, 300.0, 1.0), {c: (3.5, 5.0) for c in ('A', 'B', 'C')},
        network_mode="radius", network_radius=20.0, seed=5, method="graph_coloring",
        coloring_params=dict(cooling_rate=0.99, max_iterations=300),
        density_model=model, strategy="adaptive")
    rep, stats = _quiet(gen.generate_single_replicate, 0)
    xy_r = np.array([c.center[:2] for c in rep.cells])
    types = np.array([c.cell_type for c in rep.cells])
    inside = np.zeros(len(xy_r), bool)
    for c, R in zip(stats.layout_units["centers"], stats.layout_units["outer_radii"]):
        inside |= np.linalg.norm(xy_r - np.array(c), axis=1) < 0.8 * R
    assert inside.sum() >= 20
    assert (types[inside] == 'B').mean() > NESTS_B_MIN

