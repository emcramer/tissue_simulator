"""
Spatial analysis of 2D tissue slices.
"""

from tissue_simulator import TissueSection, TissueSlicer, SpatialNetworkAnalyzer
import matplotlib.pyplot as plt
import numpy as np

# Generate 3D tissue
print("Generating 3D tissue...")
tissue = TissueSection(
    height=500,
    width=500,
    thickness=150,
    cell_radii={
        'epithelial': (7, 11),
        'stromal': (10, 16),
        'immune': (4, 7)
    }
)

num_cells = tissue.generate_cells(max_attempts=2200)
print(f"Generated {num_cells} cells in 3D tissue\n")

# Create multiple slices
print("Creating slices at different depths...")
z_positions = [30, 60, 90, 120]
slice_results = []

for z_pos in z_positions:
    print(f"\nAnalyzing slice at z={z_pos} μm")
    print("-" * 40)
    
    # Create slice
    slicer = TissueSlicer(tissue)
    slicer.slice_plane(z_position=z_pos)
    
    print(f"  Cells in slice: {len(slicer.slice_cells)}")
    
    # Build network
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_slice(slicer, mode="contact")
    
    # Get statistics
    global_stats = analyzer.compute_global_statistics()
    cell_type_stats = analyzer.compute_cell_type_statistics()
    
    print(f"  Edges: {global_stats.total_edges}")
    print(f"  Avg degree: {global_stats.avg_degree:.2f}")
    print(f"  Avg clustering: {global_stats.avg_clustering:.4f}")
    
    slice_results.append({
        'z_position': z_pos,
        'num_cells': len(slicer.slice_cells),
        'global': global_stats,
        'cell_types': cell_type_stats,
        'analyzer': analyzer
    })

# Compare slices
print("\n" + "="*60)
print("COMPARISON ACROSS SLICES")
print("="*60)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Number of cells and edges
ax = axes[0, 0]
z_pos = [r['z_position'] for r in slice_results]
num_cells = [r['num_cells'] for r in slice_results]
num_edges = [r['global'].total_edges for r in slice_results]

ax2 = ax.twinx()
ax.plot(z_pos, num_cells, 'o-', color='steelblue', linewidth=2, markersize=8, label='Cells')
ax2.plot(z_pos, num_edges, 's-', color='coral', linewidth=2, markersize=8, label='Edges')

ax.set_xlabel('Z Position (μm)')
ax.set_ylabel('Number of Cells', color='steelblue')
ax2.set_ylabel('Number of Edges', color='coral')
ax.tick_params(axis='y', labelcolor='steelblue')
ax2.tick_params(axis='y', labelcolor='coral')
ax.set_title('Cells and Edges vs Depth')
ax.grid(True, alpha=0.3)

# Plot 2: Average degree
ax = axes[0, 1]
avg_degrees = [r['global'].avg_degree for r in slice_results]
ax.plot(z_pos, avg_degrees, 'o-', color='green', linewidth=2, markersize=8)
ax.set_xlabel('Z Position (μm)')
ax.set_ylabel('Average Degree')
ax.set_title('Average Connections vs Depth')
ax.grid(True, alpha=0.3)

# Plot 3: Clustering coefficient
ax = axes[1, 0]
clusterings = [r['global'].avg_clustering for r in slice_results]
ax.plot(z_pos, clusterings, 'o-', color='purple', linewidth=2, markersize=8)
ax.set_xlabel('Z Position (μm)')
ax.set_ylabel('Average Clustering')
ax.set_title('Clustering Coefficient vs Depth')
ax.grid(True, alpha=0.3)

# Plot 4: Network density
ax = axes[1, 1]
densities = [r['global'].network_density for r in slice_results]
ax.plot(z_pos, densities, 'o-', color='orange', linewidth=2, markersize=8)
ax.set_xlabel('Z Position (μm)')
ax.set_ylabel('Network Density')
ax.set_title('Network Density vs Depth')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('slice_spatial_comparison.png', dpi=150, bbox_inches='tight')
print("\nSaved: slice_spatial_comparison.png")
plt.show()

# Detailed analysis of middle slice
print("\n" + "="*60)
print("DETAILED ANALYSIS - MIDDLE SLICE (z=90 μm)")
print("="*60)

middle_result = slice_results[2]  # z=90

print(f"\nGlobal Statistics:")
print(f"  Total cells: {middle_result['global'].total_nodes}")
print(f"  Total edges: {middle_result['global'].total_edges}")
print(f"  Avg degree: {middle_result['global'].avg_degree:.2f}")
print(f"  Network density: {middle_result['global'].network_density:.4f}")
print(f"  Avg clustering: {middle_result['global'].avg_clustering:.4f}")
print(f"  Transitivity: {middle_result['global'].transitivity:.4f}")

print(f"\nCell Type Statistics:")
for cell_type, stats in middle_result['cell_types'].items():
    print(f"\n  {cell_type}:")
    print(f"    Count: {stats.count}")
    print(f"    Avg degree: {stats.avg_degree:.2f}")
    print(f"    Avg clustering: {stats.avg_clustering:.4f}")
    print(f"    Degree centrality: {stats.degree_centrality:.4f}")

# Visualize middle slice network
print("\nVisualizing middle slice network...")
middle_result['analyzer'].visualize_network(
    layout="spatial",
    save_path="slice_network_z90.png"
)
print("Saved: slice_network_z90.png")

# Export analysis
print("\nExporting slice analyses...")
for i, result in enumerate(slice_results):
    z = result['z_position']
    result['analyzer'].export_statistics_csv(f"slice_z{z}_analysis")
    print(f"  Exported: slice_z{z}_analysis_*.csv")

print("\n" + "="*60)
print("Slice analysis complete!")
print("="*60)
