"""
Spatial analysis module for tissue simulator.

This module provides network-based spatial analysis of cell-cell interactions
using NetworkX. It can analyze both 3D tissues and 2D slices.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
import csv
from dataclasses import dataclass, asdict

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    print("Warning: NetworkX not installed. Install with: pip install networkx")

from .tissue import TissueSection, Cell
from .slicing import TissueSlicer, SliceCell


@dataclass
class NetworkStatistics:
    """
    Container for network analysis statistics.
    
    Attributes:
        total_nodes: Total number of nodes (cells) in network
        total_edges: Total number of edges (connections)
        avg_degree: Average number of connections per cell
        network_density: Ratio of actual to possible connections
        avg_clustering: Average clustering coefficient
        transitivity: Global clustering coefficient
        avg_path_length: Average shortest path length (if connected)
        diameter: Maximum shortest path length (if connected)
        is_connected: Whether the network is fully connected
        num_components: Number of connected components
    """
    total_nodes: int
    total_edges: int
    avg_degree: float
    network_density: float
    avg_clustering: float
    transitivity: float
    avg_path_length: Optional[float]
    diameter: Optional[int]
    is_connected: bool
    num_components: int
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class CellTypeStatistics:
    """
    Statistics for a specific cell type in the network.
    
    Attributes:
        cell_type: Name of the cell type
        count: Number of cells of this type
        avg_degree: Average number of connections
        avg_clustering: Average clustering coefficient
        degree_centrality: Average degree centrality
        betweenness_centrality: Average betweenness centrality
        closeness_centrality: Average closeness centrality
    """
    cell_type: str
    count: int
    avg_degree: float
    avg_clustering: float
    degree_centrality: float
    betweenness_centrality: float
    closeness_centrality: float
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class InteractionStatistics:
    """
    Statistics for interactions between two cell types.
    
    Attributes:
        type_a: First cell type
        type_b: Second cell type
        num_interactions: Number of connections between types
        normalized_interactions: Interactions normalized by cell counts
        avg_distance: Average distance of connections
        median_distance: Median distance of connections
    """
    type_a: str
    type_b: str
    num_interactions: int
    normalized_interactions: float
    avg_distance: float
    median_distance: float
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return asdict(self)


class SpatialNetworkAnalyzer:
    """
    Analyze spatial relationships between cells using network analysis.
    
    Creates a network where nodes are cells and edges represent spatial
    relationships (contact or proximity).
    """
    
    def __init__(self):
        """Initialize the analyzer."""
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX is required. Install with: pip install networkx")
        
        self.graph: Optional[nx.Graph] = None
        self.cell_positions: Dict[int, np.ndarray] = {}
        self.cell_types: Dict[int, str] = {}
        self.cell_radii: Dict[int, float] = {}
    
    def build_network_from_tissue(self, 
                                  tissue: TissueSection,
                                  mode: str = "contact",
                                  radius: Optional[float] = None) -> nx.Graph:
        """
        Build a spatial network from a 3D tissue.
        
        Args:
            tissue: TissueSection to analyze
            mode: "contact" for touching cells, "radius" for proximity
            radius: Distance threshold for "radius" mode (in micrometers)
        
        Returns:
            NetworkX graph with cells as nodes and spatial relationships as edges
        """
        if not tissue.cells:
            raise ValueError("Tissue has no cells. Generate cells first.")
        
        # Create graph
        self.graph = nx.Graph()
        
        # Add nodes for each cell
        for i, cell in enumerate(tissue.cells):
            self.graph.add_node(
                i,
                position=cell.center,
                cell_type=cell.cell_type,
                radius=cell.radius,
                is_boundary=cell.is_boundary
            )
            self.cell_positions[i] = cell.center
            self.cell_types[i] = cell.cell_type
            self.cell_radii[i] = cell.radius
        
        # Add edges based on spatial relationships
        if mode == "contact":
            self._add_contact_edges(tissue.cells)
        elif mode == "radius":
            if radius is None:
                raise ValueError("Radius must be specified for 'radius' mode")
            self._add_radius_edges(tissue.cells, radius)
        else:
            raise ValueError(f"Unknown mode: {mode}. Use 'contact' or 'radius'")
        
        return self.graph
    
    def build_network_from_slice(self,
                                 slicer: TissueSlicer,
                                 mode: str = "contact",
                                 radius: Optional[float] = None) -> nx.Graph:
        """
        Build a spatial network from a 2D slice.
        
        Args:
            slicer: TissueSlicer with computed slice
            mode: "contact" for touching cells, "radius" for proximity
            radius: Distance threshold for "radius" mode (in micrometers)
        
        Returns:
            NetworkX graph with cells as nodes and spatial relationships as edges
        """
        if not slicer.slice_cells:
            raise ValueError("Slicer has no slice cells. Create slice first.")
        
        # Create graph
        self.graph = nx.Graph()
        
        # Add nodes for each cell in slice
        for i, slice_cell in enumerate(slicer.slice_cells):
            self.graph.add_node(
                i,
                position=slice_cell.center_2d,
                position_3d=slice_cell.center_3d,
                cell_type=slice_cell.cell_type,
                radius=slice_cell.intersection_radius,
                radius_3d=slice_cell.radius,
                is_boundary=slice_cell.is_boundary,
                distance_from_plane=slice_cell.distance_from_plane
            )
            self.cell_positions[i] = slice_cell.center_2d
            self.cell_types[i] = slice_cell.cell_type
            self.cell_radii[i] = slice_cell.intersection_radius
        
        # Add edges based on spatial relationships
        if mode == "contact":
            self._add_contact_edges_2d(slicer.slice_cells)
        elif mode == "radius":
            if radius is None:
                raise ValueError("Radius must be specified for 'radius' mode")
            self._add_radius_edges_2d(slicer.slice_cells, radius)
        else:
            raise ValueError(f"Unknown mode: {mode}. Use 'contact' or 'radius'")
        
        return self.graph
    
    def _add_contact_edges(self, cells: List[Cell]):
        """Add edges for cells in contact (3D)."""
        for i, cell_a in enumerate(cells):
            for j, cell_b in enumerate(cells[i+1:], start=i+1):
                distance = np.linalg.norm(cell_a.center - cell_b.center)
                contact_distance = cell_a.radius + cell_b.radius
                
                # Check if cells are in contact (with small tolerance)
                if distance <= contact_distance * 1.01:  # 1% tolerance
                    self.graph.add_edge(i, j, weight=distance, distance=distance)
    
    def _add_radius_edges(self, cells: List[Cell], radius: float):
        """Add edges for cells within radius (3D)."""
        for i, cell_a in enumerate(cells):
            for j, cell_b in enumerate(cells[i+1:], start=i+1):
                distance = np.linalg.norm(cell_a.center - cell_b.center)
                
                if distance <= radius:
                    self.graph.add_edge(i, j, weight=distance, distance=distance)
    
    def _add_contact_edges_2d(self, slice_cells: List[SliceCell]):
        """Add edges for cells in contact (2D)."""
        for i, cell_a in enumerate(slice_cells):
            for j, cell_b in enumerate(slice_cells[i+1:], start=i+1):
                distance = np.linalg.norm(cell_a.center_2d - cell_b.center_2d)
                contact_distance = cell_a.intersection_radius + cell_b.intersection_radius
                
                # Check if cells are in contact (with small tolerance)
                if distance <= contact_distance * 1.01:  # 1% tolerance
                    self.graph.add_edge(i, j, weight=distance, distance=distance)
    
    def _add_radius_edges_2d(self, slice_cells: List[SliceCell], radius: float):
        """Add edges for cells within radius (2D)."""
        for i, cell_a in enumerate(slice_cells):
            for j, cell_b in enumerate(slice_cells[i+1:], start=i+1):
                distance = np.linalg.norm(cell_a.center_2d - cell_b.center_2d)
                
                if distance <= radius:
                    self.graph.add_edge(i, j, weight=distance, distance=distance)
    
    def compute_global_statistics(self) -> NetworkStatistics:
        """
        Compute global network statistics.
        
        Returns:
            NetworkStatistics object with global metrics
        """
        if self.graph is None:
            raise ValueError("Network not built. Call build_network_from_tissue/slice first.")
        
        # Basic metrics
        n_nodes = self.graph.number_of_nodes()
        n_edges = self.graph.number_of_edges()
        
        # Degree statistics
        degrees = [d for n, d in self.graph.degree()]
        avg_degree = np.mean(degrees) if degrees else 0.0
        
        # Density
        density = nx.density(self.graph)
        
        # Clustering
        avg_clustering = nx.average_clustering(self.graph)
        transitivity = nx.transitivity(self.graph)
        
        # Connectivity
        is_connected = nx.is_connected(self.graph)
        num_components = nx.number_connected_components(self.graph)
        
        # Path lengths (only if connected)
        avg_path_length = None
        diameter = None
        if is_connected:
            avg_path_length = nx.average_shortest_path_length(self.graph)
            diameter = nx.diameter(self.graph)
        
        return NetworkStatistics(
            total_nodes=n_nodes,
            total_edges=n_edges,
            avg_degree=avg_degree,
            network_density=density,
            avg_clustering=avg_clustering,
            transitivity=transitivity,
            avg_path_length=avg_path_length,
            diameter=diameter,
            is_connected=is_connected,
            num_components=num_components
        )
    
    def compute_cell_type_statistics(self) -> Dict[str, CellTypeStatistics]:
        """
        Compute statistics for each cell type.
        
        Returns:
            Dictionary mapping cell type to CellTypeStatistics
        """
        if self.graph is None:
            raise ValueError("Network not built. Call build_network_from_tissue/slice first.")
        
        # Group nodes by cell type
        type_nodes = {}
        for node, data in self.graph.nodes(data=True):
            cell_type = data['cell_type']
            if cell_type not in type_nodes:
                type_nodes[cell_type] = []
            type_nodes[cell_type].append(node)
        
        # Compute centralities once
        degree_centrality = nx.degree_centrality(self.graph)
        betweenness_centrality = nx.betweenness_centrality(self.graph)
        closeness_centrality = nx.closeness_centrality(self.graph)
        clustering = nx.clustering(self.graph)
        
        # Compute statistics for each type
        results = {}
        for cell_type, nodes in type_nodes.items():
            # Degree statistics
            degrees = [self.graph.degree(n) for n in nodes]
            avg_degree = np.mean(degrees)
            
            # Clustering
            type_clustering = [clustering[n] for n in nodes]
            avg_clustering = np.mean(type_clustering)
            
            # Centralities
            type_degree_cent = [degree_centrality[n] for n in nodes]
            avg_degree_cent = np.mean(type_degree_cent)
            
            type_between_cent = [betweenness_centrality[n] for n in nodes]
            avg_between_cent = np.mean(type_between_cent)
            
            type_close_cent = [closeness_centrality[n] for n in nodes]
            avg_close_cent = np.mean(type_close_cent)
            
            results[cell_type] = CellTypeStatistics(
                cell_type=cell_type,
                count=len(nodes),
                avg_degree=avg_degree,
                avg_clustering=avg_clustering,
                degree_centrality=avg_degree_cent,
                betweenness_centrality=avg_between_cent,
                closeness_centrality=avg_close_cent
            )
        
        return results
    
    def compute_interaction_statistics(self) -> List[InteractionStatistics]:
        """
        Compute pairwise interaction statistics between cell types.
        
        Returns:
            List of InteractionStatistics for each cell type pair
        """
        if self.graph is None:
            raise ValueError("Network not built. Call build_network_from_tissue/slice first.")
        
        # Count cells by type
        type_counts = {}
        for node, data in self.graph.nodes(data=True):
            cell_type = data['cell_type']
            type_counts[cell_type] = type_counts.get(cell_type, 0) + 1
        
        # Get all unique cell types
        cell_types = sorted(type_counts.keys())
        
        # Compute interactions for each pair
        results = []
        for i, type_a in enumerate(cell_types):
            for type_b in cell_types[i:]:  # Include self-interactions
                interactions = []
                distances = []
                
                # Find all edges between these types
                for u, v, data in self.graph.edges(data=True):
                    u_type = self.graph.nodes[u]['cell_type']
                    v_type = self.graph.nodes[v]['cell_type']
                    
                    if (u_type == type_a and v_type == type_b) or \
                       (u_type == type_b and v_type == type_a):
                        interactions.append((u, v))
                        distances.append(data['distance'])
                
                # Compute statistics
                num_interactions = len(interactions)
                
                # Normalize by cell counts
                if type_a == type_b:
                    # Self-interactions: normalize by n*(n-1)/2
                    n = type_counts[type_a]
                    max_possible = n * (n - 1) / 2 if n > 1 else 1
                else:
                    # Cross-interactions: normalize by n*m
                    max_possible = type_counts[type_a] * type_counts[type_b]
                
                normalized = num_interactions / max_possible if max_possible > 0 else 0.0
                
                # Distance statistics
                avg_dist = np.mean(distances) if distances else 0.0
                median_dist = np.median(distances) if distances else 0.0
                
                results.append(InteractionStatistics(
                    type_a=type_a,
                    type_b=type_b,
                    num_interactions=num_interactions,
                    normalized_interactions=normalized,
                    avg_distance=avg_dist,
                    median_distance=median_dist
                ))
        
        return results
    
    def get_comprehensive_analysis(self) -> Dict:
        """
        Get a comprehensive analysis of the spatial network.
        
        Returns:
            Dictionary containing all statistics
        """
        return {
            'global': self.compute_global_statistics().to_dict(),
            'by_cell_type': {
                ct: stats.to_dict() 
                for ct, stats in self.compute_cell_type_statistics().items()
            },
            'interactions': [
                stats.to_dict() 
                for stats in self.compute_interaction_statistics()
            ]
        }
    
    def export_network(self, filename: str, format: str = "graphml"):
        """
        Export the network to a file.
        
        Args:
            filename: Output filename
            format: "graphml", "gexf", "gml", or "edgelist"
        """
        if self.graph is None:
            raise ValueError("Network not built.")
        
        if format == "graphml":
            nx.write_graphml(self.graph, filename)
        elif format == "gexf":
            nx.write_gexf(self.graph, filename)
        elif format == "gml":
            nx.write_gml(self.graph, filename)
        elif format == "edgelist":
            nx.write_edgelist(self.graph, filename, data=['weight', 'distance'])
        else:
            raise ValueError(f"Unknown format: {format}")
    
    def export_statistics_csv(self, base_filename: str):
        """
        Export all statistics to CSV files.
        
        Args:
            base_filename: Base name for output files (without extension)
        """
        analysis = self.get_comprehensive_analysis()
        
        # Export global statistics
        with open(f"{base_filename}_global.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Metric', 'Value'])
            for key, value in analysis['global'].items():
                writer.writerow([key, value])
        
        # Export cell type statistics
        with open(f"{base_filename}_cell_types.csv", 'w', newline='') as f:
            if analysis['by_cell_type']:
                first_type = list(analysis['by_cell_type'].values())[0]
                writer = csv.DictWriter(f, fieldnames=first_type.keys())
                writer.writeheader()
                for stats in analysis['by_cell_type'].values():
                    writer.writerow(stats)
        
        # Export interaction statistics
        with open(f"{base_filename}_interactions.csv", 'w', newline='') as f:
            if analysis['interactions']:
                writer = csv.DictWriter(f, fieldnames=analysis['interactions'][0].keys())
                writer.writeheader()
                for stats in analysis['interactions']:
                    writer.writerow(stats)
    
    def visualize_network(self, 
                         figsize: Tuple[int, int] = (12, 10),
                         layout: str = "spring",
                         save_path: Optional[str] = None):
        """
        Visualize the spatial network.
        
        Args:
            figsize: Figure size
            layout: "spring", "kamada_kawai", or "spatial" (use actual positions)
            save_path: If provided, save figure to this path
        """
        if self.graph is None:
            raise ValueError("Network not built.")
        
        import matplotlib.pyplot as plt
        from ._viz_utils import make_color_map

        fig, ax = plt.subplots(figsize=figsize)
        
        # Compute layout
        if layout == "spatial":
            # Use actual positions (works for 2D)
            pos = {node: data['position'][:2] for node, data in self.graph.nodes(data=True)}
        elif layout == "spring":
            pos = nx.spring_layout(self.graph, seed=42)
        elif layout == "kamada_kawai":
            pos = nx.kamada_kawai_layout(self.graph)
        else:
            raise ValueError(f"Unknown layout: {layout}")
        
        # Color nodes by cell type (sorted for deterministic assignment)
        color_map = make_color_map(self.cell_types.values())
        cell_types = list(color_map.keys())
        
        node_colors = [color_map[self.cell_types[node]] for node in self.graph.nodes()]
        
        # Draw network
        nx.draw_networkx_nodes(
            self.graph, pos,
            node_color=node_colors,
            node_size=100,
            alpha=0.8,
            ax=ax
        )
        
        nx.draw_networkx_edges(
            self.graph, pos,
            alpha=0.3,
            width=0.5,
            ax=ax
        )
        
        # Legend
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=color_map[ct], markersize=10, label=ct)
            for ct in cell_types
        ]
        ax.legend(handles=legend_elements, loc='best')
        
        ax.set_title(f'Spatial Network: {self.graph.number_of_nodes()} nodes, '
                    f'{self.graph.number_of_edges()} edges')
        ax.axis('off')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        
        plt.show()


def analyze_tissue_network(tissue: TissueSection,
                           mode: str = "contact",
                           radius: Optional[float] = None) -> Dict:
    """
    Convenience function to analyze a tissue's spatial network.
    
    Args:
        tissue: TissueSection to analyze
        mode: "contact" or "radius"
        radius: Distance threshold for "radius" mode
    
    Returns:
        Comprehensive analysis dictionary
    """
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_tissue(tissue, mode=mode, radius=radius)
    return analyzer.get_comprehensive_analysis()


def analyze_slice_network(slicer: TissueSlicer,
                          mode: str = "contact",
                          radius: Optional[float] = None) -> Dict:
    """
    Convenience function to analyze a slice's spatial network.
    
    Args:
        slicer: TissueSlicer with computed slice
        mode: "contact" or "radius"
        radius: Distance threshold for "radius" mode
    
    Returns:
        Comprehensive analysis dictionary
    """
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_slice(slicer, mode=mode, radius=radius)
    return analyzer.get_comprehensive_analysis()
