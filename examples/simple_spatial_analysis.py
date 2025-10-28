"""
Simple spatial network analysis example.
"""

from tissue_simulator import TissueSection, SpatialNetworkAnalyzer
import json

# Generate tissue with multiple cell types
print("Generating tissue...")
tissue = TissueSection(
    height=400,
    width=400,
    thickness=100,
    cell_radii={
        'epithelial': (6, 10),
        'stromal': (9, 15),
        'immune': (4, 7)
    }
)

num_cells = tissue.generate_cells(max_attempts=1800)
print(f"Generated {num_cells} cells\n")

# Create spatial network analyzer
print("="*60)
print("Spatial Network Analysis - Contact Mode")
print("="*60)

analyzer = SpatialNetworkAnalyzer()

# Build network based on cell contact
print("\nBuilding network (contact mode)...")
graph = analyzer.build_network_from_tissue(tissue, mode="contact")

print(f"Network created:")
print(f"  Nodes (cells): {graph.number_of_nodes()}")
print(f"  Edges (contacts): {graph.number_of_edges()}")

# Get comprehensive analysis
print("\n" + "="*60)
print("GLOBAL NETWORK STATISTICS")
print("="*60)

global_stats = analyzer.compute_global_statistics()
print(f"\nTotal nodes: {global_stats.total_nodes}")
print(f"Total edges: {global_stats.total_edges}")
print(f"Average degree: {global_stats.avg_degree:.2f}")
print(f"Network density: {global_stats.network_density:.4f}")
print(f"Average clustering: {global_stats.avg_clustering:.4f}")
print(f"Transitivity: {global_stats.transitivity:.4f}")
print(f"Is connected: {global_stats.is_connected}")
print(f"Number of components: {global_stats.num_components}")

if global_stats.avg_path_length:
    print(f"Average path length: {global_stats.avg_path_length:.2f}")
    print(f"Diameter: {global_stats.diameter}")

# Cell type statistics
print("\n" + "="*60)
print("CELL TYPE STATISTICS")
print("="*60)

cell_type_stats = analyzer.compute_cell_type_statistics()
for cell_type, stats in cell_type_stats.items():
    print(f"\n{cell_type}:")
    print(f"  Count: {stats.count}")
    print(f"  Avg degree: {stats.avg_degree:.2f}")
    print(f"  Avg clustering: {stats.avg_clustering:.4f}")
    print(f"  Degree centrality: {stats.degree_centrality:.4f}")
    print(f"  Betweenness centrality: {stats.betweenness_centrality:.4f}")
    print(f"  Closeness centrality: {stats.closeness_centrality:.4f}")

# Interaction statistics
print("\n" + "="*60)
print("INTERACTION STATISTICS")
print("="*60)

interaction_stats = analyzer.compute_interaction_statistics()
for stats in interaction_stats:
    print(f"\n{stats.type_a} <-> {stats.type_b}:")
    print(f"  Interactions: {stats.num_interactions}")
    print(f"  Normalized: {stats.normalized_interactions:.4f}")
    print(f"  Avg distance: {stats.avg_distance:.2f} μm")
    print(f"  Median distance: {stats.median_distance:.2f} μm")

# Export analysis
print("\n" + "="*60)
print("EXPORTING RESULTS")
print("="*60)

analyzer.export_statistics_csv("contact_analysis")
print("  Exported: contact_analysis_global.csv")
print("  Exported: contact_analysis_cell_types.csv")
print("  Exported: contact_analysis_interactions.csv")

analyzer.export_network("contact_network.graphml", format="graphml")
print("  Exported: contact_network.graphml")

# Visualize
print("\nVisualizing network...")
analyzer.visualize_network(layout="spring", save_path="contact_network.png")
print("  Saved: contact_network.png")

print("\n" + "="*60)
print("Analysis complete!")
print("="*60)
