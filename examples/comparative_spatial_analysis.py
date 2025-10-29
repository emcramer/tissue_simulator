"""
Comparative spatial analysis: contact mode vs. radius mode.
"""

from tissue_simulator import TissueSection, SpatialNetworkAnalyzer
import matplotlib.pyplot as plt
import numpy as np

# Generate tissue
print("Generating tissue...")
tissue = TissueSection(
    height=500,
    width=500,
    thickness=120,
    cell_radii={
        'epithelial': (7, 11),
        'stromal': (10, 16),
        'immune': (4, 7),
        'endothelial': (6, 9)
    }
)

num_cells = tissue.generate_cells(max_attempts=2500)
print(f"Generated {num_cells} cells\n")

# Compare different modes
modes = [
    ("contact", None),
    ("radius", 20),
    ("radius", 40),
    ("radius", 60)
]

results = []

for mode, radius in modes:
    print("="*60)
    if mode == "contact":
        print(f"Analyzing: CONTACT mode")
    else:
        print(f"Analyzing: RADIUS mode (r={radius} μm)")
    print("="*60)
    
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_tissue(tissue, mode=mode, radius=radius)
    
    # Get statistics
    global_stats = analyzer.compute_global_statistics()
    cell_type_stats = analyzer.compute_cell_type_statistics()
    interaction_stats = analyzer.compute_interaction_statistics()
    
    print(f"\nGlobal statistics:")
    print(f"  Edges: {global_stats.total_edges}")
    print(f"  Avg degree: {global_stats.avg_degree:.2f}")
    print(f"  Network density: {global_stats.network_density:.4f}")
    print(f"  Avg clustering: {global_stats.avg_clustering:.4f}")
    print(f"  Connected: {global_stats.is_connected}")
    
    # Store for comparison
    mode_label = "Contact" if mode == "contact" else f"Radius {radius}μm"
    results.append({
        'label': mode_label,
        'mode': mode,
        'radius': radius,
        'global': global_stats,
        'cell_types': cell_type_stats,
        'interactions': interaction_stats
    })
    
    print()

# Create comparison plots
print("="*60)
print("Creating comparison plots...")
print("="*60)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: Number of edges
ax = axes[0, 0]
labels = [r['label'] for r in results]
edges = [r['global'].total_edges for r in results]
ax.bar(range(len(labels)), edges, color='steelblue', alpha=0.7)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Number of Edges')
ax.set_title('Total Interactions')
ax.grid(True, alpha=0.3, axis='y')

# Plot 2: Average degree
ax = axes[0, 1]
avg_degrees = [r['global'].avg_degree for r in results]
ax.bar(range(len(labels)), avg_degrees, color='coral', alpha=0.7)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Average Degree')
ax.set_title('Average Connections per Cell')
ax.grid(True, alpha=0.3, axis='y')

# Plot 3: Network density
ax = axes[0, 2]
densities = [r['global'].network_density for r in results]
ax.bar(range(len(labels)), densities, color='lightgreen', alpha=0.7)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Network Density')
ax.set_title('Network Density')
ax.grid(True, alpha=0.3, axis='y')

# Plot 4: Average clustering
ax = axes[1, 0]
clusterings = [r['global'].avg_clustering for r in results]
ax.bar(range(len(labels)), clusterings, color='plum', alpha=0.7)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Average Clustering')
ax.set_title('Clustering Coefficient')
ax.grid(True, alpha=0.3, axis='y')

# Plot 5: Degree by cell type (contact mode only)
ax = axes[1, 1]
contact_result = results[0]
cell_types = list(contact_result['cell_types'].keys())
degrees_by_type = [contact_result['cell_types'][ct].avg_degree for ct in cell_types]
colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
ax.bar(range(len(cell_types)), degrees_by_type, color=colors, alpha=0.7)
ax.set_xticks(range(len(cell_types)))
ax.set_xticklabels(cell_types, rotation=45, ha='right')
ax.set_ylabel('Average Degree')
ax.set_title('Avg Degree by Cell Type (Contact)')
ax.grid(True, alpha=0.3, axis='y')

# Plot 6: Number of components
ax = axes[1, 2]
components = [r['global'].num_components for r in results]
ax.bar(range(len(labels)), components, color='orange', alpha=0.7)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_ylabel('Number of Components')
ax.set_title('Network Components')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('spatial_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: spatial_comparison.png")
plt.show()

# Print interaction analysis for contact mode
print("\n" + "="*60)
print("INTERACTION ANALYSIS (Contact Mode)")
print("="*60)

contact_interactions = results[0]['interactions']
for interaction in contact_interactions:
    print(f"\n{interaction.type_a} <-> {interaction.type_b}:")
    print(f"  Raw interactions: {interaction.num_interactions}")
    print(f"  Normalized: {interaction.normalized_interactions:.4f}")
    if interaction.num_interactions > 0:
        print(f"  Avg distance: {interaction.avg_distance:.2f} μm")

# Compare epithelial-stromal interactions across modes
print("\n" + "="*60)
print("EPITHELIAL-STROMAL INTERACTIONS ACROSS MODES")
print("="*60)

for result in results:
    for interaction in result['interactions']:
        if (interaction.type_a == 'epithelial' and interaction.type_b == 'stromal') or \
           (interaction.type_a == 'stromal' and interaction.type_b == 'epithelial'):
            print(f"\n{result['label']}:")
            print(f"  Interactions: {interaction.num_interactions}")
            print(f"  Normalized: {interaction.normalized_interactions:.4f}")
            break

print("\n" + "="*60)
print("Analysis complete!")
print("="*60)
