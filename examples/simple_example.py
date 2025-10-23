"""
Simple example: Generate a basic tissue section with uniform cell sizes.
"""

from tissue_simulator import TissueSection

# Create a tissue section with uniform cell sizes
tissue = TissueSection(
    height=400,
    width=400,
    thickness=80,
    cell_radii=(5, 12)  # All cells between 5-12 μm radius
)

# Generate cells
print("Generating tissue...")
num_cells = tissue.generate_cells(max_attempts=1000)
print(f"Generated {num_cells} cells")

# Get statistics
stats = tissue.get_cell_statistics()
print(f"\nStatistics:")
print(f"  Total cells: {stats['total_cells']}")
print(f"  Boundary cells: {stats['boundary_cells']}")
print(f"  Packing fraction: {stats['packing_fraction']:.3f}")

# Visualize
print("\nDisplaying 3D visualization...")
tissue.visualize()

# Export to CSV
output_file = "simple_tissue.csv"
tissue.export_to_csv(output_file)
print(f"\nData exported to {output_file}")
