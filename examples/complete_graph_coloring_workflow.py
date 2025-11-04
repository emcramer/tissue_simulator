"""
Complete workflow example for tissue simulation with graph-based cell type assignment.

This example demonstrates the full workflow:
1. Create a 3D tissue
2. Extract a 2D slice
3. Build a network graph from spatial relationships
4. Assign cell types using simulated annealing based on target statistics
5. Visualize results
6. Export data and graphs
7. Evaluate against target statistics

Author: Eric Cramer
Date: 2025
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tissue_simulator import (
    TissueSection,
    SpherePacker,
    TissueNetworkWorkflow,
    quick_workflow
)

# Create a global named tuple with the tissue bounds
Bounds = namedtuple('BOUNDS', ['width', 'height', 'thickness'])
BOUNDS = Bounds(width=300, height=300, thickness=80)

# Not initializing from cell type dict, so make a global shared radius
CELL_RADII = (8, 12)
CELL_RADII_CONFIG = {'placeholder': CELL_RADII}

def create_sample_target_statistics():
    """
    Create sample target statistics for demonstration.
    
    These represent the desired spatial organization of cell types:
    - 40% cancer cells
    - 35% immune cells  
    - 25% stromal cells
    
    With specific interaction patterns between cell types.
    """
    target_stats = {
        'node_counts': {
            'cancer': 40,
            'immune': 35,
            'stroma': 25
        },
        'edge_counts': {
            'cancer-cancer': 60,    # Cancer cells cluster together
            'cancer-immune': 45,    # Moderate cancer-immune interaction
            'cancer-stroma': 30,    # Some cancer-stroma interaction
            'immune-immune': 25,    # Immune cells moderately clustered
            'immune-stroma': 20,    # Immune-stroma interaction
            'stroma-stroma': 15     # Stroma less clustered
        },
        'neighbor_dist': {
            'cancer': {'cancer': 3.0, 'immune': 2.2, 'stroma': 1.5},
            'immune': {'cancer': 2.5, 'immune': 1.8, 'stroma': 1.3},
            'stroma': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.2}
        }
    }
    return target_stats


def example_quick_workflow():
    """
    Demonstrate the quick workflow function for rapid prototyping.
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Quick Workflow")
    print("="*80)
    
    # Create tissue
    print("\nCreating 3D tissue...")
    tissue = TissueSection(
        height=BOUNDS.height,
        width=BOUNDS.width,
        thickness=BOUNDS.thickness,
        cell_radii=CELL_RADII
    )
    
    # Generate cells with uniform size
    packer = SpherePacker(tissue)
    num_cells = packer.pack_cells(
        cell_types={'placeholder': [8, 12]},
        max_attempts=1500,
        min_spacing=0.5,
        allow_boundary_cells=True
    )
    print(f"Generated {num_cells} cells")
    
    # Create target statistics
    target_stats = create_sample_target_statistics()
    
    # Run quick workflow
    print("\nRunning quick workflow...")
    workflow = quick_workflow(
        tissue=tissue,
        cell_types=['cancer', 'immune', 'stroma'],
        target_stats_file=None,  # Using programmatic stats
        network_radius=50.0,
        output_dir="results/quick_workflow"
    )
    
    # Note: The quick_workflow doesn't accept target_stats_dict directly
    # We'll show how to do that in the detailed workflow
    
    print("\n✅ Quick workflow complete! Check results/quick_workflow/ for outputs.")


def example_detailed_workflow():
    """
    Demonstrate the detailed step-by-step workflow with full control.
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Detailed Step-by-Step Workflow")
    print("="*80)
    
    # Step 1: Create 3D tissue
    print("\n[Step 1/8] Creating 3D tissue...")

    # Create custom packer with specific parameters
    packer = SpherePacker(
        bounds=(BOUNDS.height, BOUNDS.width, BOUNDS.thickness),
        cell_radii_config=CELL_RADII_CONFIG,
        min_spacing=0.1,  # packing within 1/10 of 1 micron
        allow_boundary_cells=True
    )

    # Pack cells with progress tracking
    print("Packing cells...")

    def progress_callback(cells_placed, total_attempts):
        if cells_placed % 50 == 0:
            print(f"  Placed {cells_placed} cells (attempts: {total_attempts})")

    cells = packer.pack_with_progress(
        max_attempts=3000,
        callback=progress_callback
    )

    tissue = TissueSection(
        height=BOUNDS.height,
        width=BOUNDS.width,
        thickness=BOUNDS.thickness,
        cell_radii=CELL_RADII
    )
    tissue.cells = cells

    print(f"✓ Generated and placed {len(cells)} cells")
    
    # Step 2: Initialize workflow manager
    print("\n[Step 2/8] Initializing workflow manager...")
    workflow = TissueNetworkWorkflow()
    workflow.set_tissue(tissue)
    print("✓ Workflow initialized")
    
    # Step 3: Create 2D slice
    print("\n[Step 3/8] Creating 2D slice...")
    z_position = tissue.thickness / 2  # Middle of tissue
    num_slice_cells = workflow.create_slice(z_position=z_position)
    print(f"✓ Slice created with {num_slice_cells} cells")
    
    # Step 4: Build network graph
    print("\n[Step 4/8] Building network graph...")
    network_radius = 50.0  # Connect cells within 50 micrometers
    graph = workflow.build_network(mode="radius", radius=network_radius)
    print(f"✓ Network built: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    
    # Step 5: Load target statistics
    print("\n[Step 5/8] Loading target statistics...")
    target_stats = create_sample_target_statistics()
    cell_types = ['cancer', 'immune', 'stroma']
    workflow.load_target_statistics(
        statistics=target_stats,
        cell_types=cell_types
    )
    print("✓ Target statistics loaded")
    
    # Step 6: Assign cell types using simulated annealing
    print("\n[Step 6/8] Assigning cell types (this may take a minute)...")
    cell_type_assignment = workflow.assign_cell_types(
        initial_temp=1000.0,
        final_temp=0.01,
        cooling_rate=0.998,
        max_iterations=10000,
        verbose=True
    )
    workflow.apply_cell_types_to_slice()
    print("✓ Cell types assigned")
    
    # Step 7: Visualize results
    print("\n[Step 7/8] Creating visualizations...")
    output_dir = "results/detailed_workflow"
    os.makedirs(output_dir, exist_ok=True)
    
    # Visualize slice
    workflow.visualize_slice(
        save_path=os.path.join(output_dir, "tissue_slice.png"),
        figsize=(12, 12)
    )
    
    # Visualize network
    workflow.visualize_network(
        layout="spring",
        save_path=os.path.join(output_dir, "tissue_network.png"),
        figsize=(14, 12)
    )
    print("✓ Visualizations created")
    
    # Step 8: Evaluate and export
    print("\n[Step 8/8] Evaluating and exporting results...")
    
    # Compare statistics
    differences = workflow.compare_statistics(verbose=True)
    
    # Comprehensive evaluation
    evaluation = workflow.evaluate(print_report=True)
    
    # Export all results
    workflow.export_all(base_dir=output_dir, prefix="tissue")
    
    print(f"\n✅ Detailed workflow complete! Check {output_dir}/ for outputs.")
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY OF RESULTS")
    print("="*80)
    print(f"JS Divergence:          {evaluation['js_divergence']:.4f} (lower is better)")
    print(f"Cosine Similarity:      {evaluation['cosine_similarity']:.4f} (higher is better)")
    print(f"Avg % Difference:       {evaluation['avg_percent_difference']:.2f}%")
    print(f"RMSE:                   {evaluation['root_mean_squared_error']:.2f}")
    print("="*80)


def example_custom_annealing_parameters():
    """
    Demonstrate different annealing parameter settings for various use cases.
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Custom Annealing Parameters")
    print("="*80)
    
    # Create tissue
    print("\nCreating tissue...")
    tissue = TissueSection(height=250, width=250, thickness=60)
    packer = SpherePacker(tissue)
    packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1000)
    
    # Three different parameter sets
    param_sets = {
        'quick_test': {
            'initial_temp': 100.0,
            'final_temp': 0.1,
            'cooling_rate': 0.995,
            'max_iterations': 3000,
            'description': "Fast testing - lower quality results"
        },
        'balanced': {
            'initial_temp': 500.0,
            'final_temp': 0.01,
            'cooling_rate': 0.998,
            'max_iterations': 10000,
            'description': "Balanced speed and quality"
        },
        'high_quality': {
            'initial_temp': 2000.0,
            'final_temp': 0.001,
            'cooling_rate': 0.999,
            'max_iterations': 20000,
            'description': "Slow but highest quality"
        }
    }
    
    # Target statistics
    target_stats = create_sample_target_statistics()
    cell_types = ['cancer', 'immune', 'stroma']
    
    # Test each parameter set
    results = {}
    for name, params in param_sets.items():
        print(f"\n--- Testing: {name} ---")
        print(f"Description: {params['description']}")
        
        workflow = TissueNetworkWorkflow()
        workflow.set_tissue(tissue)
        workflow.create_slice(z_position=tissue.thickness / 2)
        workflow.build_network(mode="radius", radius=50.0)
        workflow.load_target_statistics(statistics=target_stats, cell_types=cell_types)
        
        # Extract annealing params
        annealing_params = {k: v for k, v in params.items() if k != 'description'}
        
        # Run assignment
        workflow.assign_cell_types(**annealing_params)
        
        # Evaluate
        evaluation = workflow.evaluate(print_report=False)
        results[name] = evaluation
        
        print(f"JS Divergence: {evaluation['js_divergence']:.4f}")
        print(f"Cosine Similarity: {evaluation['cosine_similarity']:.4f}")
    
    # Compare results
    print("\n" + "="*80)
    print("COMPARISON OF PARAMETER SETS")
    print("="*80)
    print(f"{'Parameter Set':<20} {'JS Divergence':<15} {'Cosine Sim':<15}")
    print("-" * 50)
    for name, eval_data in results.items():
        print(f"{name:<20} {eval_data['js_divergence']:<15.4f} {eval_data['cosine_similarity']:<15.4f}")
    print("="*80)


def example_network_modes():
    """
    Demonstrate different network building modes (contact vs radius).
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Different Network Modes")
    print("="*80)
    
    # Create tissue
    tissue = TissueSection(height=250, width=250, thickness=60)
    packer = SpherePacker(tissue)
    packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1000)
    
    target_stats = create_sample_target_statistics()
    cell_types = ['cancer', 'immune', 'stroma']
    
    # Test contact mode
    print("\n--- Contact Mode (cells must touch) ---")
    workflow_contact = TissueNetworkWorkflow()
    workflow_contact.set_tissue(tissue)
    workflow_contact.create_slice(z_position=tissue.thickness / 2)
    graph_contact = workflow_contact.build_network(mode="contact")
    print(f"Contact network: {graph_contact.number_of_edges()} edges")
    
    # Test radius mode with different radii
    radii = [30.0, 50.0, 70.0]
    for radius in radii:
        print(f"\n--- Radius Mode (radius={radius} μm) ---")
        workflow_radius = TissueNetworkWorkflow()
        workflow_radius.set_tissue(tissue)
        workflow_radius.create_slice(z_position=tissue.thickness / 2)
        graph_radius = workflow_radius.build_network(mode="radius", radius=radius)
        print(f"Radius network: {graph_radius.number_of_edges()} edges")
        
        # Run coloring
        workflow_radius.load_target_statistics(statistics=target_stats, cell_types=cell_types)
        workflow_radius.assign_cell_types(max_iterations=5000, verbose=False)
        evaluation = workflow_radius.evaluate(print_report=False)
        print(f"JS Divergence: {evaluation['js_divergence']:.4f}")


def example_batch_processing():
    """
    Demonstrate batch processing of multiple tissue slices.
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: Batch Processing Multiple Slices")
    print("="*80)
    
    # Create tissue
    tissue = TissueSection(height=300, width=300, thickness=100)
    packer = SpherePacker(tissue)
    num_cells = packer.pack_cells(
        cell_types={'placeholder': [8, 12]},
        max_attempts=1500
    )
    print(f"Created tissue with {num_cells} cells")
    
    target_stats = create_sample_target_statistics()
    cell_types = ['cancer', 'immune', 'stroma']
    
    # Process multiple slices at different z-positions
    z_positions = [20, 40, 60, 80]
    
    print(f"\nProcessing {len(z_positions)} slices...")
    
    results = []
    for i, z_pos in enumerate(z_positions):
        print(f"\n--- Processing slice {i+1}/{len(z_positions)} at z={z_pos} μm ---")
        
        workflow = TissueNetworkWorkflow()
        workflow.set_tissue(tissue)
        workflow.create_slice(z_position=z_pos)
        workflow.build_network(mode="radius", radius=50.0)
        workflow.load_target_statistics(statistics=target_stats, cell_types=cell_types)
        
        # Quick assignment for batch processing
        workflow.assign_cell_types(
            initial_temp=500.0,
            max_iterations=5000,
            verbose=False
        )
        
        evaluation = workflow.evaluate(print_report=False)
        results.append({
            'z_position': z_pos,
            'num_cells': len(workflow.slicer.slice_cells),
            'num_edges': workflow.graph.number_of_edges(),
            'js_divergence': evaluation['js_divergence'],
            'cosine_similarity': evaluation['cosine_similarity']
        })
        
        # Export
        output_dir = f"results/batch_processing/slice_{i+1}"
        workflow.export_all(base_dir=output_dir, prefix=f"slice_{i+1}")
    
    # Summary
    print("\n" + "="*80)
    print("BATCH PROCESSING SUMMARY")
    print("="*80)
    print(f"{'Z-Pos (μm)':<12} {'Cells':<8} {'Edges':<8} {'JS Div':<10} {'Cos Sim':<10}")
    print("-" * 58)
    for r in results:
        print(f"{r['z_position']:<12.1f} {r['num_cells']:<8} {r['num_edges']:<8} "
              f"{r['js_divergence']:<10.4f} {r['cosine_similarity']:<10.4f}")
    print("="*80)


def main():
    """Run all examples."""
    print("\n" + "="*80)
    print("TISSUE SIMULATOR: Graph-Based Cell Type Assignment Examples")
    print("="*80)
    print("\nThis script demonstrates the complete workflow for:")
    print("1. Generating 3D tissue structures")
    print("2. Extracting 2D slices")
    print("3. Building network graphs from spatial relationships")
    print("4. Assigning cell types using simulated annealing")
    print("5. Visualizing and evaluating results")
    
    # Run examples
    example_detailed_workflow()
    # example_quick_workflow()  # Commented out as it needs CSV file
    
    # uncomment these to test
    #example_custom_annealing_parameters()
    #example_network_modes()
    #example_batch_processing()
    
    print("\n" + "="*80)
    print("ALL EXAMPLES COMPLETE!")
    print("="*80)
    print("\nCheck the results/ directory for all output files.")
    print("\nGenerated files include:")
    print("  - tissue_slice.csv: 2D slice cell data")
    print("  - tissue_network.graphml: Network graph")
    print("  - tissue_statistics.csv: Network statistics")
    print("  - tissue_slice.png: Slice visualization")
    print("  - tissue_network.png: Network visualization")


if __name__ == "__main__":
    main()
