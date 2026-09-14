"""Density-aware replicates: the annealer's spatial term and the generator wiring."""

import contextlib
import io
import random

import networkx as nx
import numpy as np
import pytest

from tissue_simulator import DensityModel, ReplicateGenerator, TissueSection
from tissue_simulator.graph_coloring import GraphColorizer
from tissue_simulator.packing import SpherePacker
from tissue_simulator.replicate_generator import load_target_statistics_from_tissue

SIZE = 200.0
TYPES = ('CD8', 'Stroma', 'Tumor')
COLORING = dict(cooling_rate=0.999, max_iterations=4000)


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
