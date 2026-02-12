"""
Integration tests for graph coloring and tissue network workflow.
"""

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
