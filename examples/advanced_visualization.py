"""
Advanced visualization example: Create multiple views of the tissue.
"""

from tissue_simulator import TissueSection
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Generate tissue with multiple cell types
print("Generating tissue with multiple cell types...")
tissue = TissueSection(
    height=400,
    width=400,
    thickness=100,
    cell_radii={
        'epithelial': (7, 11),
        'stromal': (9, 16),
        'immune': (4, 7),
        'endothelial': (6, 9)
    }
)

num_cells = tissue.generate_cells(max_attempts=2000, min_spacing=0.5)
print(f"Generated {num_cells} cells")

# Create figure with multiple subplots
fig = plt.figure(figsize=(16, 12))

# Get cell data
cell_types = list(set(c.cell_type for c in tissue.cells))
colors_map = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
color_dict = dict(zip(cell_types, colors_map))

# 3D view 1: Isometric view
ax1 = fig.add_subplot(2, 3, 1, projection='3d')
for cell in tissue.cells:
    u = np.linspace(0, 2 * np.pi, 12)
    v = np.linspace(0, np.pi, 12)
    x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
    y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
    z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
    
    color = color_dict[cell.cell_type]
    alpha = 0.3 if cell.is_boundary else 0.6
    ax1.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)), alpha=alpha, linewidth=0, antialiased=True, shade=False)

ax1.set_xlabel('Width (μm)')
ax1.set_ylabel('Height (μm)')
ax1.set_zlabel('Thickness (μm)')
ax1.set_title('3D View (Isometric)')
ax1.view_init(elev=30, azim=45)

# 3D view 2: Top view
ax2 = fig.add_subplot(2, 3, 2, projection='3d')
for cell in tissue.cells:
    u = np.linspace(0, 2 * np.pi, 12)
    v = np.linspace(0, np.pi, 12)
    x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
    y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
    z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
    
    color = color_dict[cell.cell_type]
    alpha = 0.3 if cell.is_boundary else 0.6
    ax2.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)), alpha=alpha, linewidth=0, antialiased=True, shade=False)

ax2.set_xlabel('Width (μm)')
ax2.set_ylabel('Height (μm)')
ax2.set_zlabel('Thickness (μm)')
ax2.set_title('3D View (Top)')
ax2.view_init(elev=90, azim=0)

# 3D view 3: Side view
ax3 = fig.add_subplot(2, 3, 3, projection='3d')
for cell in tissue.cells:
    u = np.linspace(0, 2 * np.pi, 12)
    v = np.linspace(0, np.pi, 12)
    x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
    y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
    z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
    
    color = color_dict[cell.cell_type]
    alpha = 0.3 if cell.is_boundary else 0.6
    ax3.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)), alpha=alpha, linewidth=0, antialiased=True, shade=False)

ax3.set_xlabel('Width (μm)')
ax3.set_ylabel('Height (μm)')
ax3.set_zlabel('Thickness (μm)')
ax3.set_title('3D View (Side)')
ax3.view_init(elev=0, azim=0)

# 2D projection: XY plane
ax4 = fig.add_subplot(2, 3, 4)
for cell in tissue.cells:
    circle = plt.Circle(
        (cell.center[0], cell.center[1]),
        cell.radius,
        color=color_dict[cell.cell_type],
        alpha=0.4 if cell.is_boundary else 0.7,
        edgecolor='black',
        linewidth=0.5
    )
    ax4.add_patch(circle)

ax4.set_xlim(0, tissue.width)
ax4.set_ylim(0, tissue.height)
ax4.set_aspect('equal')
ax4.set_xlabel('Width (μm)')
ax4.set_ylabel('Height (μm)')
ax4.set_title('2D Projection (XY plane)')
ax4.grid(True, alpha=0.3)

# Cell size distribution
ax5 = fig.add_subplot(2, 3, 5)
for cell_type in cell_types:
    radii = [c.radius for c in tissue.cells if c.cell_type == cell_type]
    ax5.hist(radii, bins=15, alpha=0.6, label=cell_type, 
             color=color_dict[cell_type], edgecolor='black')

ax5.set_xlabel('Cell Radius (μm)')
ax5.set_ylabel('Frequency')
ax5.set_title('Cell Size Distribution')
ax5.legend()
ax5.grid(True, alpha=0.3)

# Cell type statistics
ax6 = fig.add_subplot(2, 3, 6)
stats = tissue.get_cell_statistics()

cell_type_counts = [stats['cell_types'][ct] for ct in cell_types]
colors_for_bar = [color_dict[ct] for ct in cell_types]

bars = ax6.bar(range(len(cell_types)), cell_type_counts, 
               color=colors_for_bar, alpha=0.7, edgecolor='black')
ax6.set_xticks(range(len(cell_types)))
ax6.set_xticklabels(cell_types, rotation=45, ha='right')
ax6.set_ylabel('Cell Count')
ax6.set_title('Cell Type Distribution')
ax6.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for i, (bar, count) in enumerate(zip(bars, cell_type_counts)):
    height = bar.get_height()
    ax6.text(bar.get_x() + bar.get_width()/2., height,
            f'{count}',
            ha='center', va='bottom')

plt.tight_layout()
plt.savefig('advanced_visualization.png', dpi=150, bbox_inches='tight')
print("\nVisualization saved as 'advanced_visualization.png'")
plt.show()

# Print detailed statistics
print("\n=== Detailed Statistics ===")
print(f"Tissue volume: {tissue.width * tissue.height * tissue.thickness:.0f} μm³")
print(f"Total cells: {stats['total_cells']}")
print(f"Packing fraction: {stats['packing_fraction']:.3f}")
print(f"Interior cells: {stats['interior_cells']}")
print(f"Boundary cells: {stats['boundary_cells']}")

print("\n=== Per Cell Type ===")
for cell_type in cell_types:
    count = stats['cell_types'][cell_type]
    avg_radius = stats['avg_radii'][cell_type]
    percentage = (count / stats['total_cells']) * 100
    print(f"\n{cell_type}:")
    print(f"  Count: {count} ({percentage:.1f}%)")
    print(f"  Avg radius: {avg_radius:.2f} μm")
    print(f"  Total volume: {count * (4/3) * np.pi * avg_radius**3:.0f} μm³")
