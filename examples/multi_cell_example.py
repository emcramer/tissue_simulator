"""
Advanced example: Generate tissue with multiple cell types.
"""

from tissue_simulator import TissueSection

# Create a tissue section with multiple cell types
tissue = TissueSection(
    height=500,
    width=500,
    thickness=100,
    cell_radii={
        'epithelial': (6, 10),   # Epithelial cells: 6-10 μm
        'stromal': (8, 15),      # Stromal cells: 8-15 μm
        'immune': (3, 6),        # Immune cells: 3-6 μm
        'endothelial': (5, 8)    # Endothelial cells: 5-8 μm
    }
)

# Generate cells with custom parameters
print("Generating multi-cell-type tissue...")
num_cells = tissue.generate_cells(
    max_attempts=2000,
    min_spacing=1.0,
    allow_boundary_cells=True
)
print(f"Generated {num_cells} cells")

# Get detailed statistics
stats = tissue.get_cell_statistics()
print(f"\n=== Tissue Statistics ===")
print(f"Total cells: {stats['total_cells']}")
print(f"Interior cells: {stats['interior_cells']}")
print(f"Boundary cells: {stats['boundary_cells']}")
print(f"Packing fraction: {stats['packing_fraction']:.3f}")

print(f"\n=== Cell Type Breakdown ===")
for cell_type, count in stats['cell_types'].items():
    avg_radius = stats['avg_radii'][cell_type]
    percentage = (count / stats['total_cells']) * 100
    print(f"{cell_type}:")
    print(f"  Count: {count} ({percentage:.1f}%)")
    print(f"  Average radius: {avg_radius:.2f} μm")

# Visualize with custom viewing angle
print("\nDisplaying 3D visualization...")
tissue.visualize(elevation=30, azimuth=60)

# Export to CSV
output_file = "multi_cell_tissue.csv"
tissue.export_to_csv(output_file)
print(f"\nData exported to {output_file}")
