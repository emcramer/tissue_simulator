# Spatial Network Analysis Documentation

## Overview

The spatial analysis module provides network-based analysis of cell-cell interactions using NetworkX. It creates networks where nodes represent cells and edges represent spatial relationships (contact or proximity).

## Key Features

- ✅ **Contact-based networks**: Connect cells that are touching
- ✅ **Radius-based networks**: Connect cells within a distance threshold
- ✅ **3D tissue analysis**: Analyze spatial relationships in 3D tissues
- ✅ **2D slice analysis**: Analyze spatial relationships in tissue slices
- ✅ **Comprehensive statistics**: Global, cell-type-specific, and interaction metrics
- ✅ **Network visualization**: Visualize cell interaction networks
- ✅ **Export capabilities**: Save networks and statistics to files

## Installation

NetworkX is required:

```bash
pip install networkx
```

## Quick Start

### Analyze 3D Tissue

```python
from tissue_simulator import TissueSection, SpatialNetworkAnalyzer

# Generate tissue
tissue = TissueSection(400, 400, 100, cell_radii={'type_a': (5, 10)})
tissue.generate_cells(max_attempts=1000)

# Create analyzer
analyzer = SpatialNetworkAnalyzer()

# Build network (contact mode)
analyzer.build_network_from_tissue(tissue, mode="contact")

# Get statistics
global_stats = analyzer.compute_global_statistics()
print(f"Average degree: {global_stats.avg_degree:.2f}")
print(f"Clustering: {global_stats.avg_clustering:.4f}")
```

### Analyze 2D Slice

```python
from tissue_simulator import TissueSlicer

# Create slice
slicer = TissueSlicer(tissue)
slicer.slice_plane(z_position=50)

# Build network from slice
analyzer.build_network_from_slice(slicer, mode="contact")

# Analyze
stats = analyzer.compute_global_statistics()
```

## Network Modes

### Contact Mode

Connects cells that are touching (surfaces in contact):

```python
analyzer.build_network_from_tissue(tissue, mode="contact")
```

**How it works:**
- Calculates distance between cell centers
- Compares to sum of radii
- Creates edge if `distance <= radius_a + radius_b`
- Uses 1% tolerance for numerical stability

**Use cases:**
- Direct cell-cell contact
- Adhesion-based interactions
- Physical touching relationships

### Radius Mode

Connects cells within a specified distance:

```python
analyzer.build_network_from_tissue(tissue, mode="radius", radius=30)
```

**Parameters:**
- `radius`: Distance threshold in micrometers

**How it works:**
- Calculates distance between cell centers
- Creates edge if `distance <= radius`

**Use cases:**
- Paracrine signaling
- Proximity-based interactions
- Neighborhood analysis

## Statistics

### Global Network Statistics

```python
stats = analyzer.compute_global_statistics()
```

**Returns NetworkStatistics with:**
- `total_nodes`: Number of cells
- `total_edges`: Number of connections
- `avg_degree`: Average connections per cell
- `network_density`: Ratio of actual/possible connections
- `avg_clustering`: Average clustering coefficient
- `transitivity`: Global clustering coefficient
- `avg_path_length`: Average shortest path (if connected)
- `diameter`: Maximum shortest path (if connected)
- `is_connected`: Whether network is fully connected
- `num_components`: Number of connected components

**Example:**
```python
print(f"Cells: {stats.total_nodes}")
print(f"Contacts: {stats.total_edges}")
print(f"Avg contacts per cell: {stats.avg_degree:.2f}")
print(f"Clustering: {stats.avg_clustering:.4f}")
```

### Cell Type Statistics

```python
type_stats = analyzer.compute_cell_type_statistics()
```

**Returns dict mapping cell type to CellTypeStatistics:**
- `cell_type`: Cell type name
- `count`: Number of cells
- `avg_degree`: Average connections
- `avg_clustering`: Average clustering
- `degree_centrality`: Degree centrality measure
- `betweenness_centrality`: Betweenness centrality
- `closeness_centrality`: Closeness centrality

**Example:**
```python
for cell_type, stats in type_stats.items():
    print(f"{cell_type}:")
    print(f"  Count: {stats.count}")
    print(f"  Avg degree: {stats.avg_degree:.2f}")
    print(f"  Clustering: {stats.avg_clustering:.4f}")
```

### Interaction Statistics

```python
interactions = analyzer.compute_interaction_statistics()
```

**Returns list of InteractionStatistics:**
- `type_a`, `type_b`: Cell types
- `num_interactions`: Number of connections
- `normalized_interactions`: Normalized by cell counts
- `avg_distance`: Average connection distance
- `median_distance`: Median connection distance

**Normalization:**
- Self-interactions: divided by n*(n-1)/2
- Cross-interactions: divided by n*m

**Example:**
```python
for interaction in interactions:
    print(f"{interaction.type_a} <-> {interaction.type_b}:")
    print(f"  Raw: {interaction.num_interactions}")
    print(f"  Normalized: {interaction.normalized_interactions:.4f}")
    print(f"  Avg distance: {interaction.avg_distance:.2f} μm")
```

## Comprehensive Analysis

Get all statistics at once:

```python
analysis = analyzer.get_comprehensive_analysis()
```

**Returns dict with:**
- `global`: Global network statistics
- `by_cell_type`: Cell type statistics
- `interactions`: Interaction statistics

## Export Functions

### Export Statistics to CSV

```python
analyzer.export_statistics_csv("analysis")
```

Creates three files:
- `analysis_global.csv`: Global statistics
- `analysis_cell_types.csv`: Per-cell-type statistics
- `analysis_interactions.csv`: Pairwise interactions

### Export Network

```python
analyzer.export_network("network.graphml", format="graphml")
```

**Supported formats:**
- `graphml`: GraphML (recommended)
- `gexf`: GEXF format
- `gml`: GML format
- `edgelist`: Simple edge list

**Usage with other tools:**
- Cytoscape: Use GraphML or GEXF
- Gephi: Use GEXF
- igraph: Any format
- NetworkX: Any format

## Visualization

```python
analyzer.visualize_network(
    figsize=(12, 10),
    layout="spring",
    save_path="network.png"
)
```

**Parameters:**
- `figsize`: Figure size tuple
- `layout`: "spring", "kamada_kawai", or "spatial"
- `save_path`: Optional path to save figure

**Layout options:**
- `spring`: Force-directed layout
- `kamada_kawai`: Energy-minimizing layout
- `spatial`: Use actual cell positions (2D only)

## Advanced Usage

### Compare Contact vs Radius

```python
# Contact network
analyzer1 = SpatialNetworkAnalyzer()
analyzer1.build_network_from_tissue(tissue, mode="contact")
contact_stats = analyzer1.compute_global_statistics()

# Radius network
analyzer2 = SpatialNetworkAnalyzer()
analyzer2.build_network_from_tissue(tissue, mode="radius", radius=30)
radius_stats = analyzer2.compute_global_statistics()

print(f"Contact edges: {contact_stats.total_edges}")
print(f"Radius edges: {radius_stats.total_edges}")
```

### Analyze Multiple Slices

```python
from tissue_simulator import create_standard_slices

slicers = create_standard_slices(tissue, num_slices=5)

for i, slicer in enumerate(slicers):
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_slice(slicer, mode="contact")
    stats = analyzer.compute_global_statistics()
    print(f"Slice {i+1}: {stats.avg_degree:.2f} avg degree")
```

### Custom Network Analysis

Access the NetworkX graph directly:

```python
graph = analyzer.graph

# Custom analysis
import networkx as nx

# Find most connected cell
degrees = dict(graph.degree())
max_node = max(degrees, key=degrees.get)
print(f"Most connected cell: {max_node} with {degrees[max_node]} connections")

# Find communities
from networkx.algorithms import community
communities = community.greedy_modularity_communities(graph)
print(f"Found {len(communities)} communities")

# Calculate additional metrics
pagerank = nx.pagerank(graph)
eigenvector = nx.eigenvector_centrality(graph)
```

## Interpretation Guide

### Network Density

**Range:** 0 to 1

- **Low (< 0.1)**: Sparse network, few interactions
- **Medium (0.1-0.3)**: Moderate connectivity
- **High (> 0.3)**: Dense network, many interactions

**Biological interpretation:**
- Contact networks typically: 0.01-0.05
- Radius networks (large r): 0.1-0.5

### Clustering Coefficient

**Range:** 0 to 1

- **Low (< 0.3)**: Random-like structure
- **Medium (0.3-0.6)**: Some local structure
- **High (> 0.6)**: Highly clustered, organized

**Biological interpretation:**
- High clustering suggests organized tissue structure
- Low clustering suggests random cell distribution

### Average Degree

**Range:** 0 to n-1 (n = number of cells)

- **Contact mode**: Typically 4-8 neighbors
- **Radius mode**: Depends on radius and cell density

**Biological interpretation:**
- Higher degree = more interactions
- Compare across cell types for differential connectivity

### Normalized Interactions

**Range:** 0 to 1

- **0**: No interactions between types
- **Low (< 0.1)**: Rare interactions
- **Medium (0.1-0.3)**: Moderate mixing
- **High (> 0.3)**: Strong mixing

**Biological interpretation:**
- Low values: Segregated cell types
- High values: Well-mixed cell types

## Examples

See the `examples/` directory:
- `simple_spatial_analysis.py`: Basic network analysis
- `comparative_spatial_analysis.py`: Compare modes
- `slice_spatial_analysis.py`: Analyze multiple slices

## Troubleshooting

### Issue: "NetworkX not installed"

```bash
pip install networkx
```

### Issue: Network visualization doesn't show

Make sure matplotlib is installed and use appropriate backend:

```python
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg'
```

### Issue: "Network not built"

Call `build_network_from_tissue()` or `build_network_from_slice()` first:

```python
analyzer = SpatialNetworkAnalyzer()
analyzer.build_network_from_tissue(tissue, mode="contact")
# Now you can compute statistics
```

### Issue: Large networks are slow

For tissues with > 1000 cells:
- Use contact mode (fewer edges than radius mode)
- Reduce visualization resolution
- Export to file and analyze in specialized tools

## Performance Notes

- Network building: O(n²) where n = number of cells
- Contact mode: Faster than radius mode
- Statistics computation: O(n + m) where m = number of edges
- Large networks (> 1000 cells): May take several minutes

## References

### Network Metrics

- **Degree**: Number of connections
- **Clustering**: Fraction of neighbors that are connected
- **Centrality**: Importance in network
- **Path length**: Distance between nodes

### Biological Context

- Cell-cell contact is fundamental to tissue organization
- Network analysis reveals spatial patterns
- Clustering indicates organized vs random structure
- Degree distribution shows connectivity patterns

## See Also

- [Tissue Simulator Guide](GUIDE.md)
- [Slicing Documentation](SLICING.md)
- [NetworkX Documentation](https://networkx.org/)
