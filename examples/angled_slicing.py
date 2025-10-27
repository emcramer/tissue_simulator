"""
Angled slicing example: Create slices at various angles through the tissue.
"""

from tissue_simulator import TissueSection, TissueSlicer
import numpy as np

# Generate a 3D tissue
print("Generating 3D tissue...")
tissue = TissueSection(
    height=500,
    width=500,
    thickness=150,
    cell_radii={
        'type_a': (5, 10),
        'type_b': (8, 15),
        'type_c': (6, 12)
    }
)

num_cells = tissue.generate_cells(max_attempts=2000)
print(f"Generated {num_cells} cells")

# Example 1: Horizontal slice (0° rotation)
print("\n=== Example 1: Horizontal Slice ===")
slicer1 = TissueSlicer(tissue)
slice_cells1 = slicer1.slice_plane(z_position=75)
print(f"Captured {len(slice_cells1)} cells")
slicer1.visualize_slice_2d(title="Horizontal Slice (0°)")

# Example 2: Tilted slice (45° around X-axis)
print("\n=== Example 2: Tilted Slice (45° around X-axis) ===")
slicer2 = TissueSlicer(tissue)
slice_cells2 = slicer2.slice_plane(
    point=(250, 250, 75),
    angle_x=45,
    angle_y=0
)
print(f"Captured {len(slice_cells2)} cells")
slicer2.visualize_slice_2d(title="Tilted Slice (45° X-rotation)")

# Example 3: Diagonal slice (30° X, 30° Y)
print("\n=== Example 3: Diagonal Slice (30° both axes) ===")
slicer3 = TissueSlicer(tissue)
slice_cells3 = slicer3.slice_plane(
    point=(250, 250, 75),
    angle_x=30,
    angle_y=30
)
print(f"Captured {len(slice_cells3)} cells")
slicer3.visualize_slice_2d(title="Diagonal Slice (30° X, 30° Y)")

# Example 4: Custom normal vector
print("\n=== Example 4: Custom Normal Vector ===")
slicer4 = TissueSlicer(tissue)
# Normal vector pointing diagonally
normal = np.array([1, 1, 1])
slice_cells4 = slicer4.slice_plane(
    point=(250, 250, 75),
    normal=normal
)
print(f"Captured {len(slice_cells4)} cells")
stats4 = slicer4.get_slice_statistics()
print(f"Plane normal: {stats4['plane_normal']}")
slicer4.visualize_slice_2d(title="Custom Normal Vector (1,1,1)")

# Visualize all slices in 3D context
print("\nVisualizing slice planes in 3D...")
slicer2.visualize_slice_in_3d()

# Export angled slice
print("\nExporting angled slice data...")
slicer2.export_slice_csv("angled_slice_45deg.csv")
print("Data exported to angled_slice_45deg.csv")

# Compare slice statistics
print("\n=== Comparison of Slices ===")
for i, slicer in enumerate([slicer1, slicer2, slicer3, slicer4], 1):
    stats = slicer.get_slice_statistics()
    print(f"\nSlice {i}:")
    print(f"  Cells captured: {stats['num_cells']}")
    print(f"  Cell types: {stats['cell_types']}")
    print(f"  Mean distance: {stats['mean_distance_from_plane']:.2f} μm")
