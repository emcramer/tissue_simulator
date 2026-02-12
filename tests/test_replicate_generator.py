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
    
    # Extract statistics
    target_stats = load_target_statistics_from_tissue(
        reference,
        network_mode="contact"
    )
    assert len(target_stats.interaction_stats) > 0
    
    # Setup generator
    generator = ReplicateGenerator(
        target_stats=target_stats,
        tissue_dimensions=(200, 200, 50),
        base_cell_radii={'type_a': (5, 8), 'type_b': (6, 10)},
        network_mode="contact",
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
