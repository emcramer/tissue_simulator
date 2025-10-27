"""
Multiple parallel slices: Create serial sections through tissue.
"""

from tissue_simulator import TissueSection, create_standard_slices
import matplotlib.pyplot as plt
import numpy as np

# Generate a 3D tissue
print("Generating 3D tissue...")
tissue = TissueSection(
    height=400,
    width=400,
    thickness=120,
    cell_radii={
        'epithelial': (6, 10),
        'stromal': (9, 15),
        'immune': (4, 7)
    }
)

num_cells = tissue.generate_cells(max_attempts=1800)
print(f"Generated {num_cells} cells")

# Create 5 evenly-spaced parallel slices
print("\nCreating 5 parallel slices...")
num_slices = 5
slicers = create_standard_slices(tissue, num_slices=num_slices)

# Print statistics for each slice
print("\n=== Slice Statistics ===")
for i, slicer in enumerate(slicers, 1):
    stats = slicer.get_slice_statistics()
    z_pos = stats['plane_point'][2]
    print(f"\nSlice {i} (z = {z_pos:.1f} μm):")
    print(f"  Cells: {stats['num_cells']}")
    print(f"  Cell types: {stats['cell_types']}")

# Create a grid visualization of all slices
print("\nCreating grid visualization...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

# Plot each slice
for i, slicer in enumerate(slicers):
    ax = axes[i]
    
    # Get slice data
    stats = slicer.get_slice_statistics()
    z_pos = stats['plane_point'][2]
    
    # Color map for cell types
    cell_types = list(set(c.cell_type for c in slicer.slice_cells))
    colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
    color_map = dict(zip(cell_types, colors))
    
    # Plot cells
    for slice_cell in slicer.slice_cells:
        color = color_map[slice_cell.cell_type]
        max_dist = stats['max_distance_from_plane'] if stats['max_distance_from_plane'] > 0 else 1
        alpha = 1.0 - (slice_cell.distance_from_plane / max_dist) * 0.5
        
        circle = plt.Circle(
            slice_cell.center_2d,
            slice_cell.intersection_radius,
            color=color,
            alpha=alpha,
            edgecolor='black',
            linewidth=0.5
        )
        ax.add_patch(circle)
    
    # Set limits and labels
    if slicer.slice_cells:
        x_coords = [c.center_2d[0] for c in slicer.slice_cells]
        y_coords = [c.center_2d[1] for c in slicer.slice_cells]
        margin = 20
        ax.set_xlim(min(x_coords) - margin, max(x_coords) + margin)
        ax.set_ylim(min(y_coords) - margin, max(y_coords) + margin)
    
    ax.set_aspect('equal')
    ax.set_title(f'Slice {i+1} (z={z_pos:.1f} μm)\n{len(slicer.slice_cells)} cells')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('U (μm)')
    ax.set_ylabel('V (μm)')

# Add legend in the last subplot
ax_legend = axes[-1]
ax_legend.axis('off')
if slicers and slicers[0].slice_cells:
    cell_types = list(set(c.cell_type for c in slicers[0].slice_cells))
    colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
    color_map = dict(zip(cell_types, colors))
    
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w',
                  markerfacecolor=color_map[ct], markersize=15, label=ct)
        for ct in cell_types
    ]
    ax_legend.legend(handles=legend_elements, loc='center', 
                    fontsize=12, title='Cell Types')

plt.suptitle(f'Serial Tissue Sections ({num_slices} slices)', fontsize=16, y=0.98)
plt.tight_layout()
plt.savefig('serial_slices.png', dpi=150, bbox_inches='tight')
print("\nGrid visualization saved as 'serial_slices.png'")
plt.show()

# Export each slice to separate CSV files
print("\nExporting slice data...")
for i, slicer in enumerate(slicers, 1):
    filename = f"slice_{i}.csv"
    slicer.export_slice_csv(filename, include_3d=True)
    print(f"  Exported {filename}")

# Analyze cell count variation across slices
print("\n=== Cell Count Analysis ===")
cell_counts = [len(s.slice_cells) for s in slicers]
z_positions = [s.slice_plane_point[2] for s in slicers]

print(f"Mean cells per slice: {np.mean(cell_counts):.1f}")
print(f"Std dev: {np.std(cell_counts):.1f}")
print(f"Range: {min(cell_counts)} - {max(cell_counts)}")

# Plot cell count vs z-position
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(z_positions, cell_counts, 'o-', linewidth=2, markersize=8)
ax.set_xlabel('Z Position (μm)', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_title('Cell Count Variation Across Tissue Depth', fontsize=14)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('cell_count_vs_depth.png', dpi=150)
print("\nCell count plot saved as 'cell_count_vs_depth.png'")
plt.show()

print("\n✓ Analysis complete!")
