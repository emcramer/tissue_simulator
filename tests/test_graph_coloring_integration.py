"""
Test script to verify graph coloring integration is working correctly.

This script runs a minimal test of the complete workflow to ensure
all components are properly integrated.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def test_basic_integration():
    """Test basic integration of all components."""
    print("\n" + "="*60)
    print("Testing Graph Coloring Integration")
    print("="*60)
    
    # Test imports
    print("\n[1/7] Testing imports...")
    try:
        from tissue_simulator import (
            TissueSection,
            SpherePacker,
            TissueSlicer,
            SpatialNetworkAnalyzer,
            GraphColorizer,
            evaluate_graph_coloring,
            TissueNetworkWorkflow
        )
        print("✓ All imports successful")
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False
    
    # Test tissue generation
    print("\n[2/7] Testing tissue generation...")
    try:
        tissue = TissueSection(
            height=200, 
            width=200, 
            thickness=50,
            cell_radii={'placeholder': (8, 12)}
        )
        num_cells = tissue.generate_cells(
            max_attempts=500,
            min_spacing=0.5
        )
        print(f"✓ Generated {num_cells} cells")
    except Exception as e:
        print(f"✗ Tissue generation failed: {e}")
        return False
    
    # Test slicing
    print("\n[3/7] Testing slicing...")
    try:
        slicer = TissueSlicer(tissue)
        slice_cells = slicer.slice_plane(z_position=tissue.thickness / 2)
        print(f"✓ Created slice with {len(slice_cells)} cells")
    except Exception as e:
        print(f"✗ Slicing failed: {e}")
        return False
    
    # Test network building
    print("\n[4/7] Testing network building...")
    try:
        analyzer = SpatialNetworkAnalyzer()
        graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=50.0)
        print(f"✓ Built network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    except Exception as e:
        print(f"✗ Network building failed: {e}")
        return False
    
    # Test graph colorizer
    print("\n[5/7] Testing graph colorizer...")
    try:
        target_stats = {
            'node_counts': {
                'cancer': max(1, len(slice_cells) // 3),
                'immune': max(1, len(slice_cells) // 3),
                'stroma': max(1, len(slice_cells) - 2 * (len(slice_cells) // 3))
            },
            'edge_counts': {
                'cancer-cancer': max(1, graph.number_of_edges() // 6),
                'cancer-immune': max(1, graph.number_of_edges() // 6),
                'cancer-stroma': max(1, graph.number_of_edges() // 6),
                'immune-immune': max(1, graph.number_of_edges() // 6),
                'immune-stroma': max(1, graph.number_of_edges() // 6),
                'stroma-stroma': max(1, graph.number_of_edges() // 6)
            },
            'neighbor_dist': {
                'cancer': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.0},
                'immune': {'cancer': 1.5, 'immune': 1.5, 'stroma': 1.0},
                'stroma': {'cancer': 1.0, 'immune': 1.0, 'stroma': 1.5}
            }
        }
        
        colorizer = GraphColorizer(
            target_graph=graph,
            colors=['cancer', 'immune', 'stroma'],
            target_statistics=target_stats
        )
        
        cell_type_assignment = colorizer.colorize(
            initial_temp=100.0,
            final_temp=0.1,
            cooling_rate=0.995,
            max_iterations=1000,
            verbose=False
        )
        print(f"✓ Graph coloring completed, assigned {len(cell_type_assignment)} cells")
    except Exception as e:
        print(f"✗ Graph coloring failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test evaluation
    print("\n[6/7] Testing evaluation...")
    try:
        from tissue_simulator import calculate_graph_statistics
        
        final_stats = calculate_graph_statistics(
            graph,
            cell_type_assignment,
            ['cancer', 'immune', 'stroma']
        )
        
        # Convert target stats format
        target_stats_formatted = {}
        for color in ['cancer', 'immune', 'stroma']:
            target_stats_formatted[f'nodes_{color}'] = target_stats['node_counts'][color]
        for key, value in target_stats['edge_counts'].items():
            target_stats_formatted[f'edges_{key}'] = value
        
        evaluation = evaluate_graph_coloring(target_stats_formatted, final_stats)
        
        print("✓ Evaluation completed:")
        print(f"  JS Divergence: {evaluation['js_divergence']:.4f}")
        print(f"  Cosine Similarity: {evaluation['cosine_similarity']:.4f}")
    except Exception as e:
        print(f"✗ Evaluation failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Test workflow manager
    print("\n[7/7] Testing TissueNetworkWorkflow...")
    try:
        workflow = TissueNetworkWorkflow()
        workflow.set_tissue(tissue)
        workflow.create_slice(z_position=tissue.thickness / 2)
        workflow.build_network(mode="radius", radius=50.0)
        workflow.load_target_statistics(
            statistics=target_stats,
            cell_types=['cancer', 'immune', 'stroma']
        )
        workflow.assign_cell_types(
            initial_temp=100.0,
            max_iterations=500,
            verbose=False
        )
        workflow_eval = workflow.evaluate(print_report=False)
        print("✓ Workflow manager completed:")
        print(f"  JS Divergence: {workflow_eval['js_divergence']:.4f}")
    except Exception as e:
        print(f"✗ Workflow manager failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n" + "="*60)
    print("✓ ALL TESTS PASSED!")
    print("="*60)
    print("\nThe graph coloring integration is working correctly.")
    print("You can now use the complete workflow for your simulations.")
    return True


def main():
    """Run the test."""
    success = test_basic_integration()
    
    if success:
        print("\n" + "="*60)
        print("Next Steps:")
        print("="*60)
        print("1. Try the complete example:")
        print("   python examples/complete_graph_coloring_workflow.py")
        print("\n2. Read the comprehensive guide:")
        print("   docs/COMPLETE_WORKFLOW_GUIDE.md")
        print("\n3. Check the API documentation:")
        print("   docs/GRAPH_COLORING_GUIDE.md")
        sys.exit(0)
    else:
        print("\n" + "="*60)
        print("TESTS FAILED")
        print("="*60)
        print("Please check the error messages above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
