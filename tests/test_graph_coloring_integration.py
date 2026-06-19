"""
Integration tests for graph coloring and tissue network workflow.
"""

import random

import pytest
import numpy as np
from tissue_simulator import (
    TissueSection,
    SpherePacker,
    TissueSlicer,
    SpatialNetworkAnalyzer,
    GraphColorizer,
    evaluate_graph_coloring,
    TissueNetworkWorkflow,
    calculate_graph_statistics
)

@pytest.fixture
def base_tissue():
    """Fixture to provide a generated tissue section."""
    tissue = TissueSection(
        height=200, 
        width=200, 
        thickness=50,
        cell_radii={'placeholder': (8, 12)}
    )
    tissue.generate_cells(max_attempts=500, min_spacing=0.5)
    return tissue

@pytest.fixture
def tissue_slice(base_tissue):
    """Fixture to provide a 2D slice."""
    slicer = TissueSlicer(base_tissue)
    slicer.slice_plane(z_position=base_tissue.thickness / 2)
    return slicer

@pytest.fixture
def spatial_graph(tissue_slice):
    """Fixture to provide a network graph from a slice."""
    analyzer = SpatialNetworkAnalyzer()
    graph = analyzer.build_network_from_slice(tissue_slice, mode="radius", radius=50.0)
    return graph

def test_imports():
    """Verify that all components can be imported."""
    # This is mostly covered by the fact that this test runs,
    # but explicit check for the main classes.
    assert TissueSection is not None
    assert GraphColorizer is not None
    assert TissueNetworkWorkflow is not None

def test_tissue_generation(base_tissue):
    """Test tissue generation in the integration context."""
    assert len(base_tissue.cells) > 0
    stats = base_tissue.get_cell_statistics()
    assert stats['total_cells'] == len(base_tissue.cells)

def test_graph_colorizer_integration(spatial_graph):
    """Test the graph colorizer with synthetic target statistics."""
    num_nodes = spatial_graph.number_of_nodes()
    num_edges = spatial_graph.number_of_edges()
    
    assert num_nodes > 0
    
    target_stats = {
        'node_counts': {
            'cancer': max(1, num_nodes // 3),
            'immune': max(1, num_nodes // 3),
            'stroma': max(1, num_nodes - 2 * (num_nodes // 3))
        },
        'edge_counts': {
            'cancer-cancer': max(1, num_edges // 6),
            'cancer-immune': max(1, num_edges // 6),
            'cancer-stroma': max(1, num_edges // 6),
            'immune-immune': max(1, num_edges // 6),
            'immune-stroma': max(1, num_edges // 6),
            'stroma-stroma': max(1, num_edges // 6)
        },
        'neighbor_dist': {
            'cancer': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.0},
            'immune': {'cancer': 1.5, 'immune': 1.5, 'stroma': 1.0},
            'stroma': {'cancer': 1.0, 'immune': 1.0, 'stroma': 1.5}
        }
    }
    
    colorizer = GraphColorizer(
        target_graph=spatial_graph,
        colors=['cancer', 'immune', 'stroma'],
        target_statistics=target_stats
    )
    
    cell_type_assignment = colorizer.colorize(
        initial_temp=10.0,
        final_temp=0.1,
        cooling_rate=0.95,
        max_iterations=100,
        verbose=False
    )
    
    assert len(cell_type_assignment) == num_nodes
    assert set(cell_type_assignment.values()).issubset({'cancer', 'immune', 'stroma'})

def test_workflow_manager(base_tissue):
    """Test the TissueNetworkWorkflow orchestrator."""
    workflow = TissueNetworkWorkflow()
    workflow.set_tissue(base_tissue)
    
    # Run steps
    num_slice_cells = workflow.create_slice(z_position=base_tissue.thickness / 2)
    assert num_slice_cells > 0
    
    graph = workflow.build_network(mode="radius", radius=50.0)
    assert graph.number_of_nodes() == num_slice_cells
    
    target_stats = {
        'node_counts': {'cancer': num_slice_cells // 2, 'immune': num_slice_cells - (num_slice_cells // 2)},
        'edge_counts': {'cancer-cancer': 10, 'cancer-immune': 10, 'immune-immune': 10},
        'neighbor_dist': {
            'cancer': {'cancer': 1.0, 'immune': 1.0},
            'immune': {'cancer': 1.0, 'immune': 1.0}
        }
    }
    
    workflow.load_target_statistics(
        statistics=target_stats,
        cell_types=['cancer', 'immune']
    )
    
    assignment = workflow.assign_cell_types(max_iterations=100, verbose=False)
    assert len(assignment) == num_slice_cells
    
    evaluation = workflow.evaluate(print_report=False)
    assert 'js_divergence' in evaluation
    assert 'cosine_similarity' in evaluation


# ---------------------------------------------------------------------------
# Seed-reproducibility tests for GraphColorizer and workflow forwarding.
# ---------------------------------------------------------------------------


def _small_target_stats(num_nodes, num_edges):
    """Build a small, consistent target-statistics dict for three cell types."""
    return {
        'node_counts': {
            'cancer': max(1, num_nodes // 3),
            'immune': max(1, num_nodes // 3),
            'stroma': max(1, num_nodes - 2 * (num_nodes // 3)),
        },
        'edge_counts': {
            'cancer-cancer': max(1, num_edges // 6),
            'cancer-immune': max(1, num_edges // 6),
            'cancer-stroma': max(1, num_edges // 6),
            'immune-immune': max(1, num_edges // 6),
            'immune-stroma': max(1, num_edges // 6),
            'stroma-stroma': max(1, num_edges // 6),
        },
        'neighbor_dist': {
            'cancer': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.0},
            'immune': {'cancer': 1.5, 'immune': 1.5, 'stroma': 1.0},
            'stroma': {'cancer': 1.0, 'immune': 1.0, 'stroma': 1.5},
        },
    }


def test_graph_colorizer_seed_reproducibility(spatial_graph):
    """Two GraphColorizers with the same seed produce identical colorings."""
    num_nodes = spatial_graph.number_of_nodes()
    num_edges = spatial_graph.number_of_edges()
    assert num_nodes > 0

    target_stats = _small_target_stats(num_nodes, num_edges)
    colors = ['cancer', 'immune', 'stroma']

    colorizer1 = GraphColorizer(
        target_graph=spatial_graph,
        colors=colors,
        target_statistics=target_stats,
        seed=12345,
    )
    coloring1 = colorizer1.colorize(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=200,
        verbose=False,
    )

    colorizer2 = GraphColorizer(
        target_graph=spatial_graph,
        colors=colors,
        target_statistics=target_stats,
        seed=12345,
    )
    coloring2 = colorizer2.colorize(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=200,
        verbose=False,
    )

    # Identical node set and bit-for-bit identical assignment.
    assert set(coloring1.keys()) == set(coloring2.keys())
    assert len(coloring1) == num_nodes
    for node in coloring1:
        assert coloring1[node] == coloring2[node], (
            f"Mismatch at node {node}: {coloring1[node]!r} vs {coloring2[node]!r}"
        )


def test_graph_colorizer_seed_none_uses_global_random(spatial_graph):
    """seed=None must store the literal stdlib random module (backwards compat).

    seed=N must store a fresh random.Random instance that is *not* the module.
    """
    num_nodes = spatial_graph.number_of_nodes()
    num_edges = spatial_graph.number_of_edges()
    target_stats = _small_target_stats(num_nodes, num_edges)
    colors = ['cancer', 'immune', 'stroma']

    unseeded = GraphColorizer(
        target_graph=spatial_graph,
        colors=colors,
        target_statistics=target_stats,
    )
    # Byte-identity guarantee: the module itself, not a wrapper.
    assert unseeded._rng is random
    assert unseeded.seed is None

    seeded = GraphColorizer(
        target_graph=spatial_graph,
        colors=colors,
        target_statistics=target_stats,
        seed=42,
    )
    assert isinstance(seeded._rng, random.Random)
    assert seeded._rng is not random
    assert seeded.seed == 42


def test_graph_colorizer_different_seeds_diverge(spatial_graph):
    """Different seeds with the same inputs should produce different colorings."""
    num_nodes = spatial_graph.number_of_nodes()
    num_edges = spatial_graph.number_of_edges()
    # Skip cleanly if the graph happens to be trivially small.
    if num_nodes < 3:
        pytest.skip("Spatial graph too small to meaningfully diverge.")

    target_stats = _small_target_stats(num_nodes, num_edges)
    colors = ['cancer', 'immune', 'stroma']

    coloring_a = GraphColorizer(
        target_graph=spatial_graph,
        colors=colors,
        target_statistics=target_stats,
        seed=1,
    ).colorize(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=200,
        verbose=False,
    )

    coloring_b = GraphColorizer(
        target_graph=spatial_graph,
        colors=colors,
        target_statistics=target_stats,
        seed=2,
    ).colorize(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=200,
        verbose=False,
    )

    assert coloring_a != coloring_b, (
        "Different seeds produced identical colorings; the seed is not "
        "actually steering the simulated annealing."
    )


def test_assign_cell_types_seed_reproducibility(base_tissue):
    """TissueNetworkWorkflow.assign_cell_types(seed=N) must be reproducible."""
    # Set up two parallel workflows on the same tissue.
    def _build_workflow():
        wf = TissueNetworkWorkflow()
        wf.set_tissue(base_tissue)
        num_slice_cells = wf.create_slice(z_position=base_tissue.thickness / 2)
        wf.build_network(mode="radius", radius=50.0)

        target_stats = {
            'node_counts': {
                'cancer': num_slice_cells // 2,
                'immune': num_slice_cells - (num_slice_cells // 2),
            },
            'edge_counts': {
                'cancer-cancer': 10,
                'cancer-immune': 10,
                'immune-immune': 10,
            },
            'neighbor_dist': {
                'cancer': {'cancer': 1.0, 'immune': 1.0},
                'immune': {'cancer': 1.0, 'immune': 1.0},
            },
        }
        wf.load_target_statistics(
            statistics=target_stats,
            cell_types=['cancer', 'immune'],
        )
        return wf, num_slice_cells

    wf1, n1 = _build_workflow()
    assignment1 = wf1.assign_cell_types(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=200,
        verbose=False,
        seed=2024,
    )

    wf2, n2 = _build_workflow()
    assignment2 = wf2.assign_cell_types(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=200,
        verbose=False,
        seed=2024,
    )

    assert n1 == n2 > 0
    assert set(assignment1.keys()) == set(assignment2.keys())
    for node in assignment1:
        assert assignment1[node] == assignment2[node], (
            f"assign_cell_types not reproducible at node {node}: "
            f"{assignment1[node]!r} vs {assignment2[node]!r}"
        )


# ---------------------------------------------------------------------------
# Warm-start (caller-supplied initial coloring) and max_iterations=0 robustness.
# ---------------------------------------------------------------------------

import networkx as nx


def _ring_with_hub_graph():
    """Build a small, fully-connected graph with a known edge set.

    Nodes 0-5 form a cycle (each node has exactly two neighbors); node 6 is a
    hub connected to all cycle nodes. Every node has at least one neighbor.
    """
    graph = nx.Graph()
    cycle = [0, 1, 2, 3, 4, 5]
    nx.add_cycle(graph, cycle)
    for node in cycle:
        graph.add_edge(6, node)
    return graph


def _three_color_stats_for_graph(graph, colors=('cancer', 'immune', 'stroma')):
    """Build a target-statistics dict sized to ``graph`` for the given colors."""
    num_nodes = graph.number_of_nodes()
    num_edges = graph.number_of_edges()
    return _small_target_stats(num_nodes, num_edges)


def _make_colorizer(graph, seed=None, colors=('cancer', 'immune', 'stroma')):
    """Construct a GraphColorizer over ``graph`` with three-color stats."""
    return GraphColorizer(
        target_graph=graph,
        colors=list(colors),
        target_statistics=_three_color_stats_for_graph(graph, colors),
        seed=seed,
    )


def test_colorize_seeded_run_reproducible_and_deterministic():
    """The same seed yields identical colorings across two full runs."""
    graph = _ring_with_hub_graph()

    def _run():
        return _make_colorizer(graph, seed=555).colorize(
            initial_temp=10.0,
            final_temp=0.5,
            cooling_rate=0.95,
            max_iterations=300,
            verbose=False,
        )

    coloring_a = _run()
    coloring_b = _run()

    assert set(coloring_a.keys()) == set(graph.nodes())
    assert coloring_a == coloring_b, (
        "Seeded colorize is not reproducible across identical runs."
    )
    # Determinism check: every node assigned a valid color.
    assert all(c in ('cancer', 'immune', 'stroma') for c in coloring_a.values())


def test_colorize_warm_start_round_trip_max_iterations_zero(capsys):
    """colorize(initial_coloring=COLORING, max_iterations=0) returns COLORING.

    This both confirms warm-start is honored and that max_iterations=0 no
    longer raises NameError (the i = -1 fix). ``verbose=True`` is required so
    the post-loop ``print(f"...{i+1} iterations.")`` line actually executes --
    that is the only code path that references ``i`` after a zero-iteration
    loop, so it is the path the i = -1 guard protects.
    """
    graph = _ring_with_hub_graph()
    colorizer = _make_colorizer(graph, seed=3)

    # Supply a complete, valid coloring for every node.
    palette = ['cancer', 'immune', 'stroma']
    warm = {node: palette[i % len(palette)] for i, node in enumerate(graph.nodes())}

    result = colorizer.colorize(
        max_iterations=0,
        verbose=True,
        initial_coloring=warm,
    )

    assert result == warm, (
        "With max_iterations=0 and a complete warm-start, the returned coloring "
        "must equal the supplied initial_coloring exactly."
    )
    # The verbose post-loop print must report i+1 == 0 iterations, locking the
    # i = -1 contract (without the guard this print raises NameError).
    captured = capsys.readouterr()
    assert "Finished after 0 iterations." in captured.out


def test_colorize_max_iterations_zero_no_warm_start_no_error():
    """max_iterations=0 without warm-start returns the random-shuffle coloring.

    Regression guard for the latent NameError when the loop never executes.
    """
    graph = _ring_with_hub_graph()
    colorizer = _make_colorizer(graph, seed=11)

    result = colorizer.colorize(max_iterations=0, verbose=False)

    assert set(result.keys()) == set(graph.nodes())
    assert all(c in ('cancer', 'immune', 'stroma') for c in result.values())


def test_colorize_initial_coloring_invalid_color_raises():
    """A warm-start assigning a color not in self.colors must raise ValueError."""
    graph = _ring_with_hub_graph()
    colorizer = _make_colorizer(graph, seed=4)

    bad = {node: 'cancer' for node in graph.nodes()}
    # Corrupt one node with a color not in the palette.
    first_node = next(iter(graph.nodes()))
    bad[first_node] = 'not_a_real_color'

    with pytest.raises(ValueError) as exc_info:
        colorizer.colorize(max_iterations=10, verbose=False, initial_coloring=bad)
    # Error message should name the offending color.
    assert 'not_a_real_color' in str(exc_info.value)


def test_colorize_initial_coloring_missing_nodes_filled():
    """Nodes missing from initial_coloring are filled with a valid color.

    Every node in the target graph must receive a valid color even when the
    supplied warm-start covers only a subset of nodes.
    """
    graph = _ring_with_hub_graph()
    colorizer = _make_colorizer(graph, seed=8)

    # Provide colors for only a couple of nodes.
    partial = {0: 'cancer', 1: 'immune'}

    result = colorizer.colorize(
        max_iterations=0,
        verbose=False,
        initial_coloring=partial,
    )

    assert set(result.keys()) == set(graph.nodes()), (
        "Returned coloring must cover every node in the target graph."
    )
    assert all(c in ('cancer', 'immune', 'stroma') for c in result.values())
    # The explicitly supplied nodes keep their colors (max_iterations=0, no moves).
    assert result[0] == 'cancer'
    assert result[1] == 'immune'


def test_colorize_initial_coloring_ignores_foreign_keys():
    """Keys not in the target graph are ignored by warm-start."""
    graph = _ring_with_hub_graph()
    colorizer = _make_colorizer(graph, seed=9)

    warm = {node: 'cancer' for node in graph.nodes()}
    warm[12345] = 'immune'  # foreign key not present in the graph

    result = colorizer.colorize(
        max_iterations=0,
        verbose=False,
        initial_coloring=warm,
    )

    assert 12345 not in result
    assert set(result.keys()) == set(graph.nodes())


def test_workflow_assign_cell_types_forwards_warm_start(base_tissue):
    """TissueNetworkWorkflow.assign_cell_types forwards initial_coloring.

    Runs the workflow once normally, then re-runs with the prior assignment as
    initial_coloring (warm-start) to confirm the parameter is accepted and
    produces a complete, valid assignment over the same node set.
    """
    wf = TissueNetworkWorkflow()
    wf.set_tissue(base_tissue)
    num_slice_cells = wf.create_slice(z_position=base_tissue.thickness / 2)
    wf.build_network(mode="radius", radius=50.0)

    target_stats = {
        'node_counts': {
            'cancer': num_slice_cells // 2,
            'immune': num_slice_cells - (num_slice_cells // 2),
        },
        'edge_counts': {
            'cancer-cancer': 10,
            'cancer-immune': 10,
            'immune-immune': 10,
        },
        'neighbor_dist': {
            'cancer': {'cancer': 1.0, 'immune': 1.0},
            'immune': {'cancer': 1.0, 'immune': 1.0},
        },
    }
    wf.load_target_statistics(statistics=target_stats, cell_types=['cancer', 'immune'])

    first = wf.assign_cell_types(
        initial_temp=10.0,
        final_temp=0.5,
        cooling_rate=0.95,
        max_iterations=100,
        verbose=False,
        seed=77,
    )
    assert len(first) == num_slice_cells

    # Warm-start from the prior assignment with zero iterations so the returned
    # coloring must equal the supplied warm-start exactly.
    second = wf.assign_cell_types(
        max_iterations=0,
        verbose=False,
        initial_coloring=dict(first),
    )
    assert second == first
    assert set(second.keys()) == set(first.keys())


# ---------------------------------------------------------------------------
# generate_colored_replicates: cold-default diverse replicates, warm opt-in.
# ---------------------------------------------------------------------------


def _build_loaded_workflow(base_tissue):
    """Build a workflow with a slice, network, and two-color target loaded."""
    wf = TissueNetworkWorkflow()
    wf.set_tissue(base_tissue)
    num_slice_cells = wf.create_slice(z_position=base_tissue.thickness / 2)
    wf.build_network(mode="radius", radius=50.0)

    target_stats = {
        'node_counts': {
            'cancer': num_slice_cells // 2,
            'immune': num_slice_cells - (num_slice_cells // 2),
        },
        'edge_counts': {
            'cancer-cancer': 10,
            'cancer-immune': 10,
            'immune-immune': 10,
        },
        'neighbor_dist': {
            'cancer': {'cancer': 1.0, 'immune': 1.0},
            'immune': {'cancer': 1.0, 'immune': 1.0},
        },
    }
    wf.load_target_statistics(statistics=target_stats, cell_types=['cancer', 'immune'])
    return wf, num_slice_cells


def test_generate_colored_replicates_shape_and_validity(base_tissue):
    """Returns N colorings, each covering every node with valid colors."""
    wf, num_slice_cells = _build_loaded_workflow(base_tissue)

    replicates = wf.generate_colored_replicates(
        num_replicates=3, seed=5, max_iterations=100, verbose=False,
    )

    assert len(replicates) == 3
    for coloring in replicates:
        assert len(coloring) == num_slice_cells
        assert set(coloring.values()).issubset({'cancer', 'immune'})
    # The first replicate is exposed as the active assignment.
    assert wf.cell_type_assignment == replicates[0]


def test_generate_colored_replicates_seed_reproducible(base_tissue):
    """The same base seed yields a bit-identical set of replicates."""
    wf1, _ = _build_loaded_workflow(base_tissue)
    wf2, _ = _build_loaded_workflow(base_tissue)

    reps1 = wf1.generate_colored_replicates(num_replicates=3, seed=2024,
                                            max_iterations=150, verbose=False)
    reps2 = wf2.generate_colored_replicates(num_replicates=3, seed=2024,
                                            max_iterations=150, verbose=False)

    assert reps1 == reps2


def test_generate_colored_replicates_cold_default_is_diverse(base_tissue):
    """Cold-start replicates (the default) are not all identical.

    Uses max_iterations=0 so each replicate is just its own seeded
    node-counts shuffle; independent derived seeds must yield >1 distinct
    coloring.
    """
    wf, num_slice_cells = _build_loaded_workflow(base_tissue)
    if num_slice_cells < 4:
        pytest.skip("Slice too small to meaningfully diverge.")

    replicates = wf.generate_colored_replicates(
        num_replicates=4, seed=7, max_iterations=0, verbose=False,
    )
    distinct = {tuple(sorted(c.items())) for c in replicates}
    assert len(distinct) > 1, (
        "Cold-start replicates collapsed to a single coloring; they should be "
        "independent."
    )


def test_generate_colored_replicates_warm_start_chains_and_collapses(base_tissue):
    """warm_start=True chains each replicate from the previous one.

    With max_iterations=0 every warm-started replicate after the first must
    equal the first exactly (no annealing happens), demonstrating the
    diversity collapse that makes warm-start a refinement tool rather than a
    replicate generator.
    """
    wf, _ = _build_loaded_workflow(base_tissue)

    replicates = wf.generate_colored_replicates(
        num_replicates=3, seed=7, warm_start=True, max_iterations=0, verbose=False,
    )
    assert replicates[1] == replicates[0]
    assert replicates[2] == replicates[0]


def test_generate_colored_replicates_validation(base_tissue):
    """num_replicates < 1 and unbuilt prerequisites raise ValueError."""
    wf, _ = _build_loaded_workflow(base_tissue)
    with pytest.raises(ValueError):
        wf.generate_colored_replicates(num_replicates=0)

    # Missing network / target statistics also raise.
    empty = TissueNetworkWorkflow()
    with pytest.raises(ValueError):
        empty.generate_colored_replicates(num_replicates=2)
