"""
Graph coloring module for assigning cell types based on network statistics.

This module integrates simulated annealing-based graph coloring to assign
cell types to tissue networks based on target statistical properties.
"""

import numpy as np
import random
import math
from collections import defaultdict, Counter
from itertools import combinations_with_replacement
from typing import Dict, List, Tuple, Optional
import csv

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    print("Warning: NetworkX not installed. Install with: pip install networkx")


class GraphColorizer:
    """
    Imposes network statistics from target data onto a graph by
    finding an optimal node coloring using Simulated Annealing.

    The statistics matched are:
    1.  Node color counts.
    2.  Pairwise edge counts between all combinations of colors.
    3.  Neighbor color distribution (average number of neighbors of color X for a node of color Y).
    """

    def __init__(self, source_graph: nx.Graph = None, 
                 target_graph: nx.Graph = None, 
                 colors: list = None,
                 target_statistics: Dict = None):
        """
        Initializes the GraphColorizer.

        Args:
            source_graph: Graph with 'color' attribute (optional if target_statistics provided)
            target_graph: Graph to be colored
            colors: List of possible color strings (e.g., ['cancer', 'immune', 'stroma'])
            target_statistics: Pre-calculated target statistics (optional)
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX is required. Install with: pip install networkx")
        
        if target_graph is None:
            raise ValueError("Target graph must be provided.")
        
        self.target_graph = target_graph
        # Ensure target_graph has no self-loops as they complicate edge counting
        self.target_graph.remove_edges_from(nx.selfloop_edges(self.target_graph))
        
        self.nodes = list(self.target_graph.nodes())
        self.colors = colors if colors is not None else []
        self.color_map = {color: i for i, color in enumerate(self.colors)}

        # Get target statistics either from source graph or directly
        if target_statistics is not None:
            print("Using provided target statistics...")
            self.target_stats = target_statistics
        elif source_graph is not None:
            if not nx.get_node_attributes(source_graph, 'color'):
                raise ValueError("Source graph nodes must have a 'color' attribute.")
            print("Calculating target statistics from source graph...")
            source_coloring = nx.get_node_attributes(source_graph, 'color')
            self.target_stats = self._calculate_statistics(source_graph, source_coloring)
        else:
            raise ValueError("Either source_graph or target_statistics must be provided.")
        
        print("Target Statistics:", self.target_stats)

    def _calculate_statistics(self, graph: nx.Graph, coloring: dict):
        """Calculates the statistics for a given graph and coloring."""
        stats = {
            'node_counts': defaultdict(int),
            'edge_counts': defaultdict(int),
            'neighbor_dist': defaultdict(lambda: defaultdict(float))
        }

        # Node counts
        for node in graph.nodes():
            color = coloring.get(node)
            if color:
                stats['node_counts'][color] += 1

        # Edge counts
        color_pairs = list(combinations_with_replacement(self.colors, 2))
        for c1, c2 in color_pairs:
            # Generate sorted key to handle ('blue', 'red') vs ('red', 'blue')
            key = '-'.join(sorted([c1, c2]))
            stats['edge_counts'][key] = 0

        for u, v in graph.edges():
            c1 = coloring.get(u)
            c2 = coloring.get(v)
            if c1 and c2:
                key = '-'.join(sorted([c1, c2]))
                stats['edge_counts'][key] += 1
        
        # Neighbor distribution
        neighbor_counts = defaultdict(lambda: defaultdict(int))
        for node in graph.nodes():
            node_color = coloring.get(node)
            if not node_color:
                continue
            for neighbor in graph.neighbors(node):
                neighbor_color = coloring.get(neighbor)
                if neighbor_color:
                    neighbor_counts[node_color][neighbor_color] += 1
        
        for c1 in self.colors:
            total_nodes_of_c1 = stats['node_counts'][c1]
            if total_nodes_of_c1 > 0:
                for c2 in self.colors:
                    # Avg number of c2 neighbors for a c1 node
                    stats['neighbor_dist'][c1][c2] = neighbor_counts[c1][c2] / total_nodes_of_c1

        return stats

    def _calculate_cost(self, current_stats: dict):
        """Calculates the cost (error) between current and target stats."""
        cost = 0.0
        
        # Weighted sum of squared errors
        # Edge count error
        for key, target_val in self.target_stats['edge_counts'].items():
            current_val = current_stats['edge_counts'].get(key, 0)
            cost += ((current_val - target_val) ** 2) * 1.0  # Weight 1.0

        # Neighbor distribution error
        for c1, dist in self.target_stats['neighbor_dist'].items():
            for c2, target_val in dist.items():
                current_val = current_stats['neighbor_dist'].get(c1, {}).get(c2, 0)
                cost += ((current_val - target_val) ** 2) * 2.0  # Weight 2.0 (often more sensitive)

        return cost
        
    def colorize(self, initial_temp=100.0, final_temp=0.1, cooling_rate=0.995, 
                 max_iterations=100000, verbose=True):
        """
        Performs the simulated annealing process to find the optimal coloring.

        Args:
            initial_temp: Starting temperature
            final_temp: Temperature at which to stop
            cooling_rate: Rate at which temperature decreases (e.g., 0.99 -> slow, 0.9 -> fast)
            max_iterations: Maximum number of iterations
            verbose: If True, prints progress updates

        Returns:
            dict: Best coloring found for the target graph
        """
        # Initial coloring: match the node counts from the target statistics
        initial_coloring = {}
        color_list = []
        for color, count in self.target_stats['node_counts'].items():
            color_list.extend([color] * count)
        
        # Ensure we have a color for every node in the target graph
        if len(color_list) < len(self.nodes):
            # Fill remaining nodes with the most frequent color
            most_frequent_color = max(self.target_stats['node_counts'], 
                                     key=self.target_stats['node_counts'].get)
            color_list.extend([most_frequent_color] * (len(self.nodes) - len(color_list)))
        
        random.shuffle(color_list)
        current_coloring = {node: color for node, color in zip(self.nodes, color_list[:len(self.nodes)])}
        
        best_coloring = current_coloring.copy()
        
        current_stats = self._calculate_statistics(self.target_graph, current_coloring)
        current_cost = self._calculate_cost(current_stats)
        best_cost = current_cost
        
        temperature = initial_temp
        
        for i in range(max_iterations):
            if temperature < final_temp:
                if verbose:
                    print("Temperature fell below final threshold. Stopping.")
                break

            # Propose a new state by swapping colors of two random nodes
            node1, node2 = random.sample(self.nodes, 2)
            
            new_coloring = current_coloring.copy()
            new_coloring[node1], new_coloring[node2] = new_coloring[node2], new_coloring[node1]

            # This is the performance bottleneck. For large graphs, this should be
            # an incremental update, not a full recalculation.
            new_stats = self._calculate_statistics(self.target_graph, new_coloring)
            new_cost = self._calculate_cost(new_stats)
            
            delta_cost = new_cost - current_cost
            
            # Acceptance criteria
            if delta_cost < 0 or random.random() < math.exp(-delta_cost / temperature):
                current_coloring = new_coloring
                current_cost = new_cost
                current_stats = new_stats
                
                if current_cost < best_cost:
                    best_coloring = current_coloring
                    best_cost = current_cost
            
            # Cool down
            temperature *= cooling_rate
            
            if verbose and i % 500 == 0:
                print(f"Iter {i}: Temp={temperature:.2f}, Cost={current_cost:.4f}, Best Cost={best_cost:.4f}")
        
        if verbose:
            print(f"\nFinished after {i+1} iterations.")
            print(f"Final best cost: {best_cost}")
        return best_coloring


def calculate_graph_statistics(graph: nx.Graph, colors_map: Dict, color_names: List[str]) -> Dict:
    """
    Calculate node and edge count statistics for a colored graph.
    
    Args:
        graph: NetworkX graph
        colors_map: Dictionary mapping node IDs to colors
        color_names: List of all possible color names
        
    Returns:
        Dictionary of statistics
    """
    stats = {}
    
    # Node counts
    node_counts = Counter(colors_map.values())
    for color in color_names:
        stats[f'nodes_{color}'] = node_counts.get(color, 0)
        
    # Edge counts
    edge_counts = Counter()
    for u, v in graph.edges():
        c1, c2 = sorted((colors_map[u], colors_map[v]))
        edge_counts[(c1, c2)] += 1
        
    for i in range(len(color_names)):
        for j in range(i, len(color_names)):
            c1, c2 = color_names[i], color_names[j]
            key = tuple(sorted((c1, c2)))
            stats[f'edges_{c1}-{c2}'] = edge_counts.get(key, 0)
            
    return stats


def compare_graph_statistics(source_stats: Dict, target_stats: Dict, verbose: bool = True) -> Dict:
    """
    Calculate and optionally print the percent difference between two sets of statistics.
    
    Args:
        source_stats: Statistics from source/target graph
        target_stats: Statistics from generated graph
        verbose: Whether to print comparison details
        
    Returns:
        Dictionary of differences for each statistic
    """
    if verbose:
        print("--- Statistics Comparison ---")
    differences = {}
    
    edge_keys = [k for k in source_stats if k.startswith('edges_')]
    
    for key in edge_keys:
        source_val = source_stats.get(key, 0)
        target_val = target_stats.get(key, 0)
        
        if source_val == 0 and target_val == 0:
            diff = 0.0
        elif source_val == 0:
            diff = float('inf')  # Or handle as a special case
        else:
            diff = (abs(target_val - source_val) / source_val) * 100
        
        differences[key] = diff
        if verbose:
            print(f"{key}:")
            print(f"  Source: {source_val}, Target: {target_val}")
            print(f"  Percent Difference: {diff:.2f}%")
        
    avg_diff = sum(differences.values()) / len(differences) if differences else 0
    if verbose:
        print(f"\nAverage Percent Difference (Edges): {avg_diff:.2f}%")
    return differences


def load_target_statistics_from_csv(filepath: str, color_names: List[str]) -> Dict:
    """
    Load target statistics from a CSV file.
    
    Args:
        filepath: Path to CSV file with statistics
        color_names: List of color names to use
        
    Returns:
        Dictionary with target statistics in GraphColorizer format
    """
    import pandas as pd
    
    # Read statistics
    df = pd.read_csv(filepath)
    
    if len(df) == 0:
        raise ValueError("CSV file is empty")
    
    # Extract statistics from first row
    row = df.iloc[0]
    
    # Build target statistics structure
    target_stats = {
        'node_counts': {},
        'edge_counts': {},
        'neighbor_dist': defaultdict(lambda: defaultdict(float))
    }
    
    # Extract node counts
    for color in color_names:
        node_key = f'nodes_{color}'
        if node_key in row:
            target_stats['node_counts'][color] = int(row[node_key])
    
    # Extract edge counts
    for i in range(len(color_names)):
        for j in range(i, len(color_names)):
            c1, c2 = color_names[i], color_names[j]
            edge_key = f'edges_{c1}-{c2}'
            if edge_key in row:
                key = '-'.join(sorted([c1, c2]))
                target_stats['edge_counts'][key] = int(row[edge_key])
    
    return target_stats


def export_colored_graph_statistics(graph: nx.Graph, colors_map: Dict, 
                                    color_names: List[str], filename: str):
    """
    Export graph statistics to CSV file.
    
    Args:
        graph: NetworkX graph
        colors_map: Dictionary mapping node IDs to colors
        color_names: List of color names
        filename: Output CSV filename
    """
    import pandas as pd
    
    stats = calculate_graph_statistics(graph, colors_map, color_names)
    pd.DataFrame([stats]).to_csv(filename, index=False)
    print(f"Saved statistics to {filename}")


def visualize_colored_graph(graph: nx.Graph, 
                           colors_map: Dict,
                           color_palette: Dict = None,
                           layout: str = "spring",
                           title: str = "Colored Graph",
                           save_path: str = None,
                           figsize: Tuple[int, int] = (12, 10)):
    """
    Visualize a colored graph.
    
    Args:
        graph: NetworkX graph
        colors_map: Dictionary mapping node IDs to colors
        color_palette: Dictionary mapping color names to RGB values
        layout: Layout algorithm ("spring", "kamada_kawai", "circular")
        title: Plot title
        save_path: If provided, save figure to this path
        figsize: Figure size
    """
    import matplotlib.pyplot as plt
    
    # Default color palette if not provided
    if color_palette is None:
        unique_colors = list(set(colors_map.values()))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_colors)))
        color_palette = dict(zip(unique_colors, colors))
    
    # Compute layout
    if layout == "spring":
        pos = nx.spring_layout(graph, seed=42)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(graph)
    elif layout == "circular":
        pos = nx.circular_layout(graph)
    else:
        raise ValueError(f"Unknown layout: {layout}")
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Draw graph
    node_colors = [color_palette.get(colors_map[node], 'gray') for node in graph.nodes()]
    
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, 
                           node_size=500, alpha=0.8)
    nx.draw_networkx_edges(graph, pos, alpha=0.3, width=1.0)
    nx.draw_networkx_labels(graph, pos, font_size=8, font_color='white')
    
    # Add legend
    unique_colors = list(set(colors_map.values()))
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w',
                  markerfacecolor=color_palette.get(ct, 'gray'), 
                  markersize=10, label=ct)
        for ct in unique_colors
    ]
    plt.legend(handles=legend_elements, loc='best')
    
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved visualization to {save_path}")
    
    plt.show()


def visualize_graph_comparison(source_graph: nx.Graph,
                               source_colors: Dict,
                               target_graph: nx.Graph,
                               target_colors: Dict,
                               color_palette: Dict = None,
                               save_path: str = "graph_comparison.png",
                               figsize: Tuple[int, int] = (16, 8)):
    """
    Visualize source and target graphs side by side.
    
    Args:
        source_graph: Source graph with target statistics
        source_colors: Color mapping for source graph
        target_graph: Target graph after coloring
        target_colors: Color mapping for target graph
        color_palette: Dictionary mapping color names to RGB values
        save_path: Path to save visualization
        figsize: Figure size
    """
    import matplotlib.pyplot as plt
    
    # Default color palette if not provided
    if color_palette is None:
        unique_colors = list(set(list(source_colors.values()) + list(target_colors.values())))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_colors)))
        color_palette = dict(zip(unique_colors, colors))
    
    plt.figure(figsize=figsize)
    
    # Source Graph
    plt.subplot(1, 2, 1)
    source_node_colors = [color_palette.get(source_colors[node], 'gray') 
                         for node in source_graph.nodes()]
    pos_source = nx.spring_layout(source_graph, seed=42)
    nx.draw(source_graph, pos_source, with_labels=True, 
           node_color=source_node_colors, node_size=500, font_color='white')
    plt.title("Source Graph (Target Statistics)")

    # Target Graph
    plt.subplot(1, 2, 2)
    target_node_colors = [color_palette.get(target_colors[node], 'gray') 
                         for node in target_graph.nodes()]
    pos_target = nx.spring_layout(target_graph, seed=42)
    nx.draw(target_graph, pos_target, with_labels=True, 
           node_color=target_node_colors, node_size=500, font_color='white')
    plt.title("Target Graph with Imposed Statistics")

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison visualization to {save_path}")
    plt.show()
