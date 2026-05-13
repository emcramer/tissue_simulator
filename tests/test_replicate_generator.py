"""
Integration tests for the replicate generator.
"""

import pytest
import numpy as np
from tissue_simulator import (
    ReplicateGenerator,
    TargetStatistics,
    ReplicateStatistics,
    load_replicate_stats_csv,
    load_target_statistics_from_tissue,
    TissueSection,
    SpatialNetworkAnalyzer,
    InteractionStatistics
)

def test_imports():
    """Verify that all replicate generator components can be imported."""
    assert ReplicateGenerator is not None
    assert TargetStatistics is not None
    assert load_replicate_stats_csv is not None

def test_replicate_generation_workflow():
    """Test the full replicate generation workflow."""
    # Create small reference tissue
    reference = TissueSection(
        height=200,
        width=200,
        thickness=50,
        cell_radii={'type_a': (5, 8), 'type_b': (6, 10)}
    )
    
    num_cells = reference.generate_cells(max_attempts=500)
    assert num_cells > 0
    
    # Extract statistics. Use radius mode so the test exercises a
    # non-empty interaction network even when min_spacing prevents
    # surface-touching contacts (divergence would otherwise be nan).
    target_stats = load_target_statistics_from_tissue(
        reference,
        network_mode="radius",
        network_radius=20.0,
    )
    assert len(target_stats.interaction_stats) > 0

    # Setup generator
    generator = ReplicateGenerator(
        target_stats=target_stats,
        tissue_dimensions=(200, 200, 50),
        base_cell_radii={'type_a': (5, 8), 'type_b': (6, 10)},
        network_mode="radius",
        network_radius=20.0,
        seed=42
    )
    
    # Generate one replicate
    tissue, stats = generator.generate_single_replicate(
        replicate_id=0,
        max_attempts=500,
        max_iterations=2,
        tolerance=0.20
    )
    
    assert isinstance(tissue, TissueSection)
    assert isinstance(stats, ReplicateStatistics)
    assert len(tissue.cells) > 0
    assert stats.num_cells == len(tissue.cells)
    assert stats.divergence_score >= 0

def test_replicate_generator_seed_reproducibility():
    """Two ReplicateGenerator(seed=N) instances produce identical replicates."""
    # Build a small reference and extract target stats once.
    reference = TissueSection(
        height=150, width=150, thickness=50,
        cell_radii={'type_a': (5, 8), 'type_b': (6, 10)},
        seed=7,
    )
    reference.generate_cells(max_attempts=300)
    target_stats = load_target_statistics_from_tissue(
        reference, network_mode="radius", network_radius=20.0
    )

    common_kwargs = dict(
        target_stats=target_stats,
        tissue_dimensions=(150, 150, 50),
        base_cell_radii={'type_a': (5, 8), 'type_b': (6, 10)},
        network_mode="radius",
        network_radius=20.0,
        seed=2026,
    )

    gen1 = ReplicateGenerator(**common_kwargs)
    tissue1, stats1 = gen1.generate_single_replicate(
        replicate_id=3, max_attempts=300, max_iterations=2, tolerance=0.20
    )

    gen2 = ReplicateGenerator(**common_kwargs)
    tissue2, stats2 = gen2.generate_single_replicate(
        replicate_id=3, max_attempts=300, max_iterations=2, tolerance=0.20
    )

    assert len(tissue1.cells) == len(tissue2.cells)
    assert len(tissue1.cells) > 0
    for c1, c2 in zip(tissue1.cells, tissue2.cells):
        np.testing.assert_array_equal(c1.center, c2.center)
        assert c1.radius == c2.radius
        assert c1.cell_type == c2.cell_type

    # Divergence score must also be bit-identical for two identical runs.
    # nan != nan under IEEE 754, so treat both-nan as a match explicitly.
    d1, d2 = stats1.divergence_score, stats2.divergence_score
    assert (np.isnan(d1) and np.isnan(d2)) or d1 == d2
    # With radius-mode networking the target is non-empty, so the score
    # should be a finite non-negative number, not nan.
    assert not np.isnan(d1)
    assert d1 >= 0


def test_divergence_is_nan_when_all_zero():
    """A target+measured with zero interactions yields nan divergence."""
    # Build target stats that have zero signal for the single pair.
    zero_target = TargetStatistics(
        interaction_stats=[
            InteractionStatistics(
                type_a='type_a', type_b='type_a',
                num_interactions=0,
                normalized_interactions=0.0,
                avg_distance=0.0, median_distance=0.0,
            )
        ],
        cell_type_proportions={'type_a': 1.0},
    )

    gen = ReplicateGenerator(
        target_stats=zero_target,
        tissue_dimensions=(100, 100, 50),
        base_cell_radii={'type_a': (5, 8)},
        network_mode="contact",
        seed=1,
    )

    # Measured stats that are also zero everywhere for the same pair.
    measured = [
        InteractionStatistics(
            type_a='type_a', type_b='type_a',
            num_interactions=0,
            normalized_interactions=0.0,
            avg_distance=0.0, median_distance=0.0,
        )
    ]

    div = gen._compute_interaction_divergence(measured, zero_target.interaction_stats)
    assert np.isnan(div), f"Expected nan, got {div!r}"


def test_cell_types_ordering_is_deterministic():
    """ReplicateGenerator.cell_types must be a sorted tuple, not a set.

    A set's iteration order depends on PYTHONHASHSEED, and downstream
    code feeds cell_types order into the per-replicate RNG via dict
    insertion. Storing it as a sorted tuple keeps cross-process runs
    bit-reproducible.
    """
    target = TargetStatistics(
        interaction_stats=[
            InteractionStatistics(type_a="zebra", type_b="alpha",
                                  num_interactions=1,
                                  normalized_interactions=0.1,
                                  avg_distance=1.0, median_distance=1.0),
            InteractionStatistics(type_a="middle", type_b="middle",
                                  num_interactions=1,
                                  normalized_interactions=0.1,
                                  avg_distance=1.0, median_distance=1.0),
        ],
        cell_type_proportions={"zebra": 0.33, "alpha": 0.34, "middle": 0.33},
    )
    gen = ReplicateGenerator(
        target_stats=target,
        tissue_dimensions=(100, 100, 50),
        base_cell_radii={"zebra": (5, 8), "alpha": (5, 8), "middle": (5, 8)},
        network_mode="contact",
        seed=1,
    )
    assert isinstance(gen.cell_types, tuple)
    assert gen.cell_types == ("alpha", "middle", "zebra")


def test_mcp_server_integration():
    """Test that MCP server includes replicate generator tools."""
    try:
        from tissue_simulator.mcp.server import TissueSimulatorMCPServer
        
        server = TissueSimulatorMCPServer()
        
        # Verify new attributes exist in the server implementation
        assert hasattr(server, 'replicate_generator')
        assert hasattr(server, 'generated_replicates')
    except ImportError:
        pytest.skip("MCP not installed, skipping integration test")
