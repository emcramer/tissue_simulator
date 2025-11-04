"""
Complete workflow example: Tissue simulation with network-based cell type assignment.

This example demonstrates the full pipeline:
1. Generate a 3D tissue with placeholder cell types
2. Slice a 2D section from the tissue
3. Build a network graph from the slice
4. Assign cell types based on target network statistics
5. Visualize and export results
6. Evaluate the quality of cell type assignment
"""

import numpy as np
from tissue_simulator import (
    TissueSection, 
    SpherePacker,
    TissueNetworkWorkflow
)


def create_target_statistics():
    """
    Create example target statistics for a tissue with three cell types:
    - cancer: 40% of cells
    - immune: 30% of cells
    - stroma: 30% of cells
    """
    # For a network with ~100 cells
    target_stats = {
        'node_counts': {
            'cancer': 40,
            'immune': 30,
            'stroma': 30
        },
        'edge_counts': {
            # Expected edge counts between cell types
            'cancer-cancer': 45,  # Cancer cells cluster together
            'cancer-immune': 35,  # Cancer-immune interactions
            'cancer-stroma': 25,  # Cancer-stroma interactions
            'immune-immune': 20,  # Immune cells moderately cluster
            'immune-stroma': 15,  # Immune-stroma interactions
            'stroma-stroma': 30   # Stroma cells form support structure
        },
        'neighbor_dist': {
            # Average neighbors of each type for each cell type
            'cancer': {'cancer': 2.5, 'immune': 1.8, 'stroma': 1.2},
            'immune': {'cancer': 2.0, 'immune': 1.3, 'stroma': 1.0},
            'stroma': {'cancer': 1.5, 'immune': 1.0, 'stroma': 2.0}
        }
    }
    return target_stats


def main():
    """Run the complete workflow."""
    
    print("="*70)
    print("TISSUE NETWORK WORKFLOW EXAMPLE")
    print("="*70)
    
    # =========================================================================
    # Step 1: Create a 3D tissue with random cell types
    # =========================================================================
    print("\n[Step 1] Creating 3D tissue...")
    
    # Define tissue dimensions
    tissue = TissueSection(
        height=200,   # micrometers
        width=200,    # micrometers
        thickness=50  # micrometers
    )
    
    # Define cell types with size ranges (we'll reassign types later)
    cell_types = {
        'placeholder': [8, 12]  # Cell radius range in micrometers
    }
    
    # Generate cells using sphere packing
    packer = SpherePacker(tissue)
    num_cells = packer.pack_cells(
        cell_types=cell_types,
        max_attempts=1000,
        min_spacing=0.5,
        allow_boundary_cells=True
    )
    
    print(f"Generated {num_cells} cells")
    
    # =========================================================================
    # Step 2: Set up target statistics for cell type assignment
    # =========================================================================
    print("\n[Step 2] Setting up target statistics...")
    
    # Define the cell types we want to assign
    cell_types_to_assign = ['cancer', 'immune', 'stroma']
    
    # Create target statistics
    target_stats = create_target_statistics()
    
    print("Target cell type distribution:")
    for cell_type, count in target_stats['node_counts'].items():
        print(f"  {cell_type}: {count} cells")
    
    # =========================================================================
    # Step 3: Run the complete workflow
    # =========================================================================
    print("\n[Step 3] Running tissue network workflow...")
    
    # Create workflow manager
    workflow = TissueNetworkWorkflow()
    
    # Run complete workflow with custom annealing parameters
    evaluation = workflow.run_complete_workflow(
        tissue=tissue,
        z_position=tissue.thickness / 2,  # Middle of tissue
        network_radius=50.0,  # 50 micrometer radius for connections
        target_stats_dict=target_stats,
        cell_types=cell_types_to_assign,
        annealing_params={
            'initial_temp': 1000.0,
            'final_temp': 0.01,
            'cooling_rate': 0.998,
            'max_iterations': 15000,
            'verbose': True
        },
        export_dir="workflow_results",
        visualize=True
    )
    
    # =========================================================================
    # Step 4: Additional analysis and visualization
    # =========================================================================
    print("\n[Step 4] Additional analysis...")
    
    # Get final statistics
    final_stats = workflow.get_statistics()
    print("\nFinal cell type distribution:")
    for key, value in final_stats.items():
        if key.startswith('nodes_'):
            print(f"  {key}: {value}")
    
    # Compare with target
    print("\nStatistical comparison:")
    workflow.compare_statistics(verbose=True)
    
    # =========================================================================
    # Step 5: Summary
    # =========================================================================
    print("\n" + "="*70)
    print("WORKFLOW SUMMARY")
    print("="*70)
    print(f"\nTissue dimensions: {tissue.width} x {tissue.height} x {tissue.thickness} μm")
    print(f"Total cells generated: {len(tissue.cells)}")
    print(f"Cells in slice: {len(workflow.slicer.slice_cells)}")
    print(f"Network nodes: {workflow.graph.number_of_nodes()}")
    print(f"Network edges: {workflow.graph.number_of_edges()}")
    
    print("\nEvaluation metrics:")
    for metric, value in evaluation.items():
        if isinstance(value, float):
            print(f"  {metric}: {value:.4f}")
    
    print("\nResults exported to: workflow_results/")
    print("  - tissue_slice.csv: Slice cell data")
    print("  - tissue_network.graphml: Network graph")
    print("  - tissue_statistics.csv: Network statistics")
    print("  - tissue_slice.png: Slice visualization")
    print("  - tissue_network.png: Network visualization")
    
    print("\n" + "="*70)
    print("WORKFLOW COMPLETE!")
    print("="*70)


def quick_example():
    """
    Simplified example using the quick_workflow function.
    """
    print("\n" + "="*70)
    print("QUICK WORKFLOW EXAMPLE")
    print("="*70)
    
    # Create tissue
    tissue = TissueSection(height=200, width=200, thickness=50)
    cell_types = {'placeholder': [8, 12]}
    packer = SpherePacker(tissue)
    packer.pack_cells(cell_types=cell_types, max_attempts=1000)
    
    # Create target statistics (could also load from CSV)
    target_stats = create_target_statistics()
    
    # Save target statistics to CSV for demonstration
    import pandas as pd
    stats_df = {}
    for cell_type, count in target_stats['node_counts'].items():
        stats_df[f'nodes_{cell_type}'] = [count]
    for edge_key, count in target_stats['edge_counts'].items():
        stats_df[f'edges_{edge_key}'] = [count]
    pd.DataFrame(stats_df).to_csv("target_statistics.csv", index=False)
    
    # Run quick workflow
    from tissue_simulator import quick_workflow
    
    workflow = quick_workflow(
        tissue=tissue,
        cell_types=['cancer', 'immune', 'stroma'],
        target_stats_file="target_statistics.csv",
        network_radius=50.0,
        output_dir="quick_results"
    )
    
    print("\nQuick workflow complete! Results in quick_results/")


if __name__ == "__main__":
    # Run the detailed workflow
    main()
    
    # Optionally run the quick workflow
    # quick_example()
