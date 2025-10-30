"""
Example: Generate tissue replicates matching spatial statistics.

This example demonstrates how to:
1. Create an initial tissue with specific spatial patterns
2. Extract spatial statistics from that tissue
3. Generate multiple replicates that match those statistics
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from tissue_simulator import (
    TissueSection,
    ReplicateGenerator,
    load_target_statistics_from_tissue
)


def main():
    print("=" * 80)
    print("Tissue Replicate Generation Example")
    print("=" * 80)
    
    # Step 1: Create an initial "reference" tissue
    print("\n1. Creating reference tissue...")
    reference_tissue = TissueSection(
        height=400,
        width=400,
        thickness=100,
        cell_radii={
            'cancer': (8, 12),
            'immune': (5, 8),
            'fibroblast': (6, 10)
        }
    )
    
    num_cells = reference_tissue.generate_cells(
        max_attempts=2000,
        min_spacing=0.5,
        allow_boundary_cells=True
    )
    
    print(f"   Generated {num_cells} cells in reference tissue")
    
    stats = reference_tissue.get_cell_statistics()
    print(f"   Cell types: {stats['cell_types']}")
    print(f"   Packing fraction: {stats['packing_fraction']:.4f}")
    
    # Step 2: Extract spatial statistics
    print("\n2. Extracting spatial statistics from reference tissue...")
    target_stats = load_target_statistics_from_tissue(
        reference_tissue,
        network_mode="contact"
    )
    
    print(f"   Found {len(target_stats.interaction_stats)} interaction types")
    print("   Interaction patterns:")
    for interaction in target_stats.interaction_stats:
        print(f"      {interaction.type_a} <-> {interaction.type_b}: "
              f"{interaction.normalized_interactions:.4f} (normalized)")
    
    print("\n   Cell type proportions:")
    for cell_type, proportion in target_stats.cell_type_proportions.items():
        print(f"      {cell_type}: {proportion:.3f}")
    
    # Step 3: Setup replicate generator
    print("\n3. Setting up replicate generator...")
    generator = ReplicateGenerator(
        target_stats=target_stats,
        tissue_dimensions=(400, 400, 100),
        base_cell_radii={
            'cancer': (8, 12),
            'immune': (5, 8),
            'fibroblast': (6, 10)
        },
        network_mode="contact",
        seed=42  # For reproducibility
    )
    
    # Step 4: Generate replicates
    print("\n4. Generating replicates...")
    print("   This may take a minute or two...")
    
    replicates = generator.generate_replicates(
        num_replicates=3,
        max_attempts=2000,
        min_spacing=0.5,
        allow_boundary=True,
        max_iterations=5,
        tolerance=0.15
    )
    
    print(f"\n   Successfully generated {len(replicates)} replicates!")
    
    # Step 5: Analyze replicates
    print("\n5. Replicate statistics:")
    print("   " + "-" * 70)
    print(f"   {'ID':<4} {'Cells':<8} {'Cancer':<8} {'Immune':<8} {'Fib':<8} {'Pack':<8} {'Div':<8}")
    print("   " + "-" * 70)
    
    for tissue, rep_stats in replicates:
        cancer_count = rep_stats.cell_type_counts.get('cancer', 0)
        immune_count = rep_stats.cell_type_counts.get('immune', 0)
        fibro_count = rep_stats.cell_type_counts.get('fibroblast', 0)
        
        print(f"   {rep_stats.replicate_id:<4} "
              f"{rep_stats.num_cells:<8} "
              f"{cancer_count:<8} "
              f"{immune_count:<8} "
              f"{fibro_count:<8} "
              f"{rep_stats.packing_fraction:<8.4f} "
              f"{rep_stats.divergence_score:<8.4f}")
    
    # Calculate averages
    avg_divergence = sum(s.divergence_score for _, s in replicates) / len(replicates)
    print("   " + "-" * 70)
    print(f"   Average divergence: {avg_divergence:.4f}")
    
    # Step 6: Export results
    print("\n6. Exporting results...")
    
    # Export statistics
    generator.export_replicate_statistics(replicates, "tmp/output/replicates")
    print("   ✓ Exported statistics to tmp/output/replicates_summary.csv")
    print("   ✓ Exported interactions to tmp/output/replicates_interactions.csv")
    
    # Export tissue files
    generator.export_replicate_tissues(replicates, "tmp/output/replicate_tissues")
    print("   ✓ Exported tissue files to tmp/output/replicate_tissues/")
    
    print("\n" + "=" * 80)
    print("Complete! Check the 'tmp/output' directory for exported files.")
    print("=" * 80)


if __name__ == "__main__":
    main()
