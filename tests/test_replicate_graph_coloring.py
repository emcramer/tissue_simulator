"""
Tests for the graph-coloring replicate mode, adaptive stopping, multi-restart,
parallelism, the consistency report, and the differential-evolution radius
optimizer added to ReplicateGenerator.
"""

from collections import Counter

import numpy as np
import pytest
import networkx as nx

from tissue_simulator import (
    TissueSection,
    SpatialNetworkAnalyzer,
    GraphColorizer,
    ReplicateGenerator,
    TargetStatistics,
    load_target_statistics_from_tissue,
)
from tissue_simulator.graph_coloring import color_graph_to_targets
from tissue_simulator.convergence import find_convergence_time


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

RADII = {'cancer': (8, 12), 'immune': (5, 8), 'fibroblast': (6, 10)}


def _reference_target(seed=1):
    ref = TissueSection(180, 180, 40, cell_radii=RADII, seed=seed)
    ref.generate_cells(max_attempts=500, seed=seed)
    return load_target_statistics_from_tissue(
        ref, network_mode="radius", network_radius=30.0
    )


def _make_generator(target, method="graph_coloring", **kwargs):
    params = dict(
        target_stats=target,
        tissue_dimensions=(180, 180, 40),
        base_cell_radii=RADII,
        network_mode="radius",
        network_radius=30.0,
        seed=7,
        method=method,
    )
    params.update(kwargs)
    return ReplicateGenerator(**params)


# ---------------------------------------------------------------------------
# Bridge: InteractionStatistics -> GraphColorizer target_statistics
# ---------------------------------------------------------------------------

def test_round_proportions_to_counts_sums_to_total():
    props = {'a': 0.5, 'b': 0.3, 'c': 0.2}
    for total in (0, 1, 7, 100, 999):
        counts = ReplicateGenerator._round_proportions_to_counts(props, total)
        assert sum(counts.values()) == total
        assert all(v >= 0 for v in counts.values())


def test_build_colorizer_targets_roundtrip_zero_cost():
    """The source coloring must score ~0 cost against its own derived targets."""
    g = nx.random_geometric_graph(120, 0.18, seed=1)
    pos = nx.get_node_attributes(g, 'pos')
    colors = ['cancer', 'immune', 'stroma']
    truth = {n: colors[min(2, int(pos[n][0] * 3))] for n in g.nodes()}
    for n in g.nodes():
        g.nodes[n]['cell_type'] = truth[n]
    for u, v in g.edges():
        g.edges[u, v]['distance'] = 1.0

    analyzer = SpatialNetworkAnalyzer()
    analyzer.graph = g
    istats = analyzer.compute_interaction_statistics()

    cnt = Counter(truth.values())
    n_total = sum(cnt.values())
    props = {c: cnt[c] / n_total for c in colors}
    target = TargetStatistics(interaction_stats=istats, cell_type_proportions=props)

    gen = _make_generator(target, method="graph_coloring",
                          base_cell_radii={c: (5, 8) for c in colors})
    targets = gen._build_colorizer_targets(g)

    # Node counts round-trip exactly and sum to N.
    assert sum(targets['node_counts'].values()) == g.number_of_nodes()

    gc = GraphColorizer(target_graph=g, colors=colors, target_statistics=targets)
    src_stats, _ = gc._calculate_statistics(g, truth)
    assert gc._calculate_cost(src_stats) == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# graph_coloring replicate mode
# ---------------------------------------------------------------------------

def test_graph_coloring_mode_produces_valid_replicates():
    target = _reference_target()
    gen = _make_generator(target, coloring_params={'max_iterations': 2000})
    reps = gen.generate_replicates(num_replicates=2, max_attempts=500)

    valid = set(gen.cell_types)
    for tissue, stats in reps:
        assert len(tissue.cells) > 0
        assert all(c.cell_type in valid for c in tissue.cells)
        assert not np.isnan(stats.divergence_score)
        # Proportions are locked: realized type counts == rounded target counts.
        expected = ReplicateGenerator._round_proportions_to_counts(
            target.cell_type_proportions, len(tissue.cells)
        )
        assert dict(Counter(c.cell_type for c in tissue.cells)) == expected


def test_graph_coloring_mode_reproducible():
    target = _reference_target()

    def run():
        gen = _make_generator(target, coloring_params={'max_iterations': 2000})
        reps = gen.generate_replicates(num_replicates=3, max_attempts=500)
        return [s.divergence_score for _, s in reps]

    assert run() == run()


def test_graph_coloring_lower_divergence_than_radius_tuning():
    """The graph-coloring path should match interaction targets at least as
    well as radius-tuning (typically much better)."""
    target = _reference_target()
    cp = {'max_iterations': 6000, 'initial_temp': 50.0, 'cooling_rate': 0.999}

    gc = _make_generator(target, method="graph_coloring", coloring_params=cp)
    rt = _make_generator(target, method="radius_tuning")

    gc_div = np.array([s.divergence_score for _, s in
                       gc.generate_replicates(num_replicates=3, max_attempts=500)])
    rt_div = np.array([s.divergence_score for _, s in
                       rt.generate_replicates(num_replicates=3, max_attempts=500,
                                              max_iterations=5)])

    assert np.nanmean(gc_div) < np.nanmean(rt_div)


def test_multi_restart_runs_and_stays_low():
    """n_restarts keeps the best (lowest SA-cost) of k runs per replicate.

    The selection is by SA cost, which correlates with -- but is not identical
    to -- the reported normalized-interaction divergence, so divergence is not
    strictly monotone in n_restarts. We assert the result stays valid and that
    multi-restart does not meaningfully worsen the mean divergence.
    """
    target = _reference_target()
    cp = {'max_iterations': 2000}

    one = _make_generator(target, n_restarts=1, coloring_params=cp)
    three = _make_generator(target, n_restarts=3, coloring_params=cp)

    d1 = np.nanmean([s.divergence_score for _, s in
                     one.generate_replicates(num_replicates=2, max_attempts=500)])
    d3 = np.nanmean([s.divergence_score for _, s in
                     three.generate_replicates(num_replicates=2, max_attempts=500)])
    assert np.isfinite(d3)
    assert d3 <= d1 + 0.02


def test_invalid_method_and_optimizer_raise():
    target = _reference_target()
    with pytest.raises(ValueError):
        _make_generator(target, method="not_a_method")
    with pytest.raises(ValueError):
        _make_generator(target, method="radius_tuning", radius_optimizer="bogus")


# ---------------------------------------------------------------------------
# consistency_report
# ---------------------------------------------------------------------------

def test_consistency_report_compares_methods():
    target = _reference_target()
    gc = _make_generator(target, coloring_params={'max_iterations': 2000})
    rt = _make_generator(target, method="radius_tuning")

    gc_reps = gc.generate_replicates(num_replicates=3, max_attempts=500)
    rt_reps = rt.generate_replicates(num_replicates=3, max_attempts=500, max_iterations=5)

    report = gc.consistency_report({"graph_coloring": gc_reps, "radius_tuning": rt_reps})
    assert set(report["per_method"]) == {"graph_coloring", "radius_tuning"}
    for stats in report["per_method"].values():
        assert "mean" in stats and "std" in stats and "cv" in stats
    assert len(report["pairwise"]) == 1

    # Single-method form also works.
    single = gc.consistency_report(gc_reps)
    assert "replicates" in single["per_method"]


# ---------------------------------------------------------------------------
# Parallel execution (determinism)
# ---------------------------------------------------------------------------

def test_parallel_matches_serial():
    target = _reference_target()

    def run(parallel):
        gen = _make_generator(target, coloring_params={'max_iterations': 1500})
        reps = gen.generate_replicates(num_replicates=3, max_attempts=500,
                                       parallel=parallel)
        return [round(s.divergence_score, 6) for _, s in reps]

    serial = run(False)
    try:
        parallel = run(True)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"ProcessPoolExecutor unavailable in this environment: {exc}")
    assert serial == parallel


# ---------------------------------------------------------------------------
# Adaptive stopping + cost history on colorize
# ---------------------------------------------------------------------------

def _small_colorizer(seed=5):
    g = nx.gnm_random_graph(60, 200, seed=3)
    ts = {
        'node_counts': {'a': 20, 'b': 20, 'c': 20},
        'edge_counts': {'a-a': 33, 'a-b': 33, 'a-c': 33,
                        'b-b': 33, 'b-c': 33, 'c-c': 33},
        'neighbor_dist': {'a': {'a': 2, 'b': 1, 'c': 1},
                          'b': {'a': 1, 'b': 2, 'c': 1},
                          'c': {'a': 1, 'b': 1, 'c': 2}},
    }
    return GraphColorizer(target_graph=g, colors=['a', 'b', 'c'],
                          target_statistics=ts, seed=seed)


def test_colorize_return_history_monotone():
    coloring, history = _small_colorizer().colorize(
        max_iterations=3000, verbose=False, return_history=True
    )
    assert isinstance(coloring, dict)
    assert len(history) > 0
    # best_cost history is monotone non-increasing.
    assert all(history[i + 1] <= history[i] for i in range(len(history) - 1))
    # convergence diagnostics accept the trajectory.
    t = find_convergence_time(history, window=50, cv_threshold=0.05,
                              require_stationary=False)
    assert t is None or t >= 0


def test_colorize_patience_stops_early():
    """With slow cooling + a high iteration cap, patience must stop before the cap."""
    _, history = _small_colorizer().colorize(
        initial_temp=10.0, final_temp=1e-6, cooling_rate=0.99999,
        max_iterations=200000, verbose=False, patience=1500, return_history=True
    )
    assert len(history) < 200000


def test_colorize_default_return_type_unchanged():
    """Without return_history, colorize still returns a plain dict (byte-compat)."""
    result = _small_colorizer().colorize(max_iterations=500, verbose=False)
    assert isinstance(result, dict)


def test_color_graph_to_targets_return_cost():
    g = nx.gnm_random_graph(40, 120, seed=1)
    ts = {
        'node_counts': {'a': 20, 'b': 20},
        'edge_counts': {'a-a': 30, 'a-b': 30, 'b-b': 30},
        'neighbor_dist': {'a': {'a': 1.5, 'b': 1.5}, 'b': {'a': 1.5, 'b': 1.5}},
    }
    coloring, cost = color_graph_to_targets(
        g, ['a', 'b'], ts, seed=1, return_cost=True,
        max_iterations=1000, verbose=False
    )
    assert isinstance(coloring, dict) and len(coloring) == 40
    assert isinstance(cost, float) and cost >= 0


# ---------------------------------------------------------------------------
# Differential-evolution radius optimizer (opt-in)
# ---------------------------------------------------------------------------

def test_mcp_setup_replicate_generator_method():
    """The setup_replicate_generator MCP handler wires method/n_restarts through."""
    import asyncio
    import json
    pytest.importorskip("mcp")
    from tissue_simulator.mcp.server import TissueSimulatorMCPServer

    server = TissueSimulatorMCPServer()
    server.target_stats = _reference_target()
    server.network_mode = "radius"
    server.network_radius = 30.0

    result = asyncio.run(server._handle_setup_replicate_generator({
        "height": 180, "width": 180, "thickness": 40,
        "cell_radii": {k: list(v) for k, v in RADII.items()},
        "seed": 7, "method": "graph_coloring", "n_restarts": 2,
    }))
    data = json.loads(result[0].text)
    assert data["status"] == "success", data
    assert data["method"] == "graph_coloring"
    assert server.replicate_generator.method == "graph_coloring"
    assert server.replicate_generator.n_restarts == 2


def test_differential_evolution_radius_optimizer_runs():
    target = _reference_target()
    gen = _make_generator(
        target, method="radius_tuning", radius_optimizer="differential_evolution",
        de_params={'maxiter': 3, 'popsize': 4, 'polish': False},
    )

    def run():
        reps = gen.generate_replicates(num_replicates=1, max_attempts=400)
        t, s = reps[0]
        return len(t.cells), round(s.divergence_score, 5)

    a = run()
    b = run()
    assert a[0] > 0
    assert a == b  # deterministic (fixed-seed objective)
