"""
Example: Generate replicates from CSV statistics file.

This example shows how to load target statistics from a CSV file
and generate replicates matching those statistics.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
from tissue_simulator import (
    ReplicateGenerator,
    load_target_statistics_from_csv,
    InteractionStatistics
)


def create_example_statistics_csv(filepath: str):
    """Create an example statistics CSV file."""
    # Example interaction statistics for cancer-immune-fibroblast tissue
    data = {
        'type_a': ['cancer', 'cancer', 'cancer', 'immune', 'immune', 'fibroblast'],
        'type_b': ['cancer', 'immune', 'fibroblast', 'immune', 'fibroblast', 'fibroblast'],
        'num_interactions': [45, 38, 42, 28, 35, 30],
        'normalized_interactions': [0.12, 0.15, 0.14, 0.18, 0.16, 0.11],
        'avg_distance': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        'median_distance': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    }
    
    df = pd.DataFrame(data)
    
    # Create output directory if needed
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(filepath, index=False)
    print(f"Created example statistics file: {filepath}")


def main():
    print("=" * 80)
    print("Generate Replicates from CSV Statistics")
    print("=" * 80)
    
    # Step 1: Create example CSV (in real use, you'd have this already)
    stats_file = "output/target_statistics.csv"
    print(f"\n1. Creating example statistics file...")
    create_example_statistics_csv(stats_file)
    
    # Display the CSV contents
    print("\n   CSV contents:")
    df = pd.read_csv(stats_file)
    print(df.to_string(index=False))
    
    # Step 2: Load statistics from CSV
    print(f"\n2. Loading target statistics from {stats_file}...")
    target_stats = load_target_statistics_from_csv(stats_file)
    
    print(f"   Loaded {len(target_stats.interaction_stats)} interaction types")
    
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
        seed=123
    )
    
    print("   Generator configured:")
    print(f"   - Tissue: 400 x 400 x 100 μm")
    print(f"   - Cell types: cancer, immune, fibroblast")
    print(f"   - Network mode: contact")
    
    # Step 4: Generate replicates
    print("\n4. Generating 3 replicates...")
    print("   (This will take a moment...)")
    
    replicates = generator.generate_replicates(
        num_replicates=3,
        max_attempts=1500,
        tolerance=0.15
    )
    
    print(f"\n   Generated {len(replicates)} replicates successfully!")
    
    # Step 5: Display results
    print("\n5. Replicate Summary:")
    print("   " + "-" * 60)
    
    for tissue, stats in replicates:
        print(f"\n   Replicate {stats.replicate_id}:")
        print(f"      Total cells: {stats.num_cells}")
        print(f"      Cell types: {stats.cell_type_counts}")
        print(f"      Packing fraction: {stats.packing_fraction:.4f}")
        print(f"      Divergence from target: {stats.divergence_score:.4f}")
    
    # Step 6: Export
    print("\n6. Exporting results...")
    generator.export_replicate_statistics(replicates, "output/csv_replicates")
    generator.export_replicate_tissues(replicates, "output/csv_replicate_tissues")
    
    print("   ✓ Statistics exported to output/csv_replicates_*.csv")
    print("   ✓ Tissues exported to output/csv_replicate_tissues/")
    
    print("\n" + "=" * 80)
    print("Complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
