"""
Simple 2D slicing example: Extract and visualize a horizontal slice.
"""

from tissue_simulator import TissueSection, TissueSlicer

# Generate a 3D tissue
print("Generating 3D tissue...")
tissue = TissueSection(
    height=400,
    width=400,
    thickness=100,
    cell_radii={'epithelial': (6, 10), 'stromal': (8, 14)}
)

num_cells = tissue.generate_cells(max_attempts=1500)
print(f"Generated {num_cells} cells in 3D tissue")

# Create a horizontal slice at z=50 (middle of tissue)
print("\nCreating 2D slice at z=50...")
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(z_position=50)

print(f"Captured {len(slice_cells)} cells in slice")

# Get slice statistics
stats = slicer.get_slice_statistics()
print(f"\nSlice Statistics:")
print(f"  Cell types: {stats['cell_types']}")
print(f"  Mean distance from plane: {stats['mean_distance_from_plane']:.2f} μm")

# Visualize the 2D slice
print("\nVisualizing 2D slice...")
slicer.visualize_slice_2d(title="Horizontal Slice at z=50 μm")

# Visualize the slice in 3D context
print("\nVisualizing slice in 3D context...")
slicer.visualize_slice_in_3d()

# Export slice data
output_file = "horizontal_slice.csv"
slicer.export_slice_csv(output_file)
print(f"\nSlice data exported to {output_file}")
