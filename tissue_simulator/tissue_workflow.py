"""
Complete workflow for tissue simulation with graph-based cell type assignment.

This module provides a unified interface for the complete workflow:
1. Generate 3D tissue
2. Slice tissue section
3. Build network graph
4. Assign cell types using simulated annealing
5. Visualize and export results
6. Evaluate against target statistics
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import csv
import os

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False

from .tissue import TissueSection
from .slicing import TissueSlicer
from .spatial_analysis import SpatialNetworkAnalyzer
from .graph_coloring import (
    GraphColorizer,
    color_graph_to_targets,
    calculate_graph_statistics,
    compare_graph_statistics,
    load_target_statistics_from_csv,
    export_colored_graph_statistics,
    visualize_colored_graph,
    visualize_graph_comparison
)
from .evaluation import evaluate_graph_coloring, print_evaluation_report


class TissueNetworkWorkflow:
    """
    Complete workflow manager for tissue simulation with network-based cell type assignment.
    """
    
    def __init__(self):
        """Initialize the workflow manager."""
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX is required. Install with: pip install networkx")
        
        self.tissue: Optional[TissueSection] = None
        self.slicer: Optional[TissueSlicer] = None
        self.analyzer: Optional[SpatialNetworkAnalyzer] = None
        self.graph: Optional[nx.Graph] = None
        self.cell_type_assignment: Optional[Dict] = None
        self.target_statistics: Optional[Dict] = None
        self.color_names: Optional[List[str]] = None
        
    def set_tissue(self, tissue: TissueSection):
        """
        Set the tissue for the workflow.
        
        Args:
            tissue: TissueSection object
        """
        self.tissue = tissue
        print(f"Tissue set: {len(tissue.cells)} cells")
    
    def create_slice(self, z_position: float = None, 
                    point: Tuple[float, float, float] = None,
                    normal: Tuple[float, float, float] = None,
                    angle_x: float = 0.0,
                    angle_y: float = 0.0) -> int:
        """
        Create a 2D slice from the tissue.
        
        Args:
            z_position: Z-position for horizontal slice (simplified)
            point: Point on the slice plane
            normal: Normal vector to the plane
            angle_x: Rotation angle around X-axis (degrees)
            angle_y: Rotation angle around Y-axis (degrees)
            
        Returns:
            Number of cells in the slice
        """
        if self.tissue is None:
            raise ValueError("Tissue not set. Call set_tissue() first.")
        
        self.slicer = TissueSlicer(self.tissue)
        slice_cells = self.slicer.slice_plane(
            z_position=z_position,
            point=point,
            normal=normal,
            angle_x=angle_x,
            angle_y=angle_y
        )
        
        print(f"Slice created: {len(slice_cells)} cells captured")
        return len(slice_cells)
    
    def build_network(self, mode: str = "radius", radius: float = None) -> nx.Graph:
        """
        Build a network graph from the tissue slice.
        
        Args:
            mode: "contact" for touching cells, "radius" for proximity
            radius: Distance threshold for "radius" mode (micrometers)
            
        Returns:
            NetworkX graph
        """
        if self.slicer is None:
            raise ValueError("Slice not created. Call create_slice() first.")
        
        self.analyzer = SpatialNetworkAnalyzer()
        self.graph = self.analyzer.build_network_from_slice(
            self.slicer, 
            mode=mode, 
            radius=radius
        )
        
        print(f"Network built: {self.graph.number_of_nodes()} nodes, "
              f"{self.graph.number_of_edges()} edges")
        return self.graph
    
    def load_target_statistics(self, filepath: str = None, 
                              statistics: Dict = None,
                              cell_types: List[str] = None):
        """
        Load target statistics for cell type assignment.
        
        Args:
            filepath: Path to CSV file with target statistics
            statistics: Dictionary with pre-computed statistics
            cell_types: List of cell type names
        """
        if cell_types is None:
            raise ValueError("cell_types must be provided")
        
        self.color_names = cell_types
        
        if filepath is not None:
            self.target_statistics = load_target_statistics_from_csv(filepath, cell_types)
            print(f"Loaded target statistics from {filepath}")
        elif statistics is not None:
            self.target_statistics = statistics
            print("Target statistics set")
        else:
            raise ValueError("Either filepath or statistics must be provided")
    
    def assign_cell_types(self,
                         initial_temp: float = 100.0,
                         final_temp: float = 0.1,
                         cooling_rate: float = 0.995,
                         max_iterations: int = 100000,
                         verbose: bool = True,
                         seed: Optional[int] = None,
                         initial_coloring: Optional[Dict] = None) -> Dict:
        """
        Assign cell types to network nodes using simulated annealing.

        Args:
            initial_temp: Starting temperature for annealing
            final_temp: Final temperature for annealing
            cooling_rate: Rate of temperature decrease
            max_iterations: Maximum iterations
            verbose: Whether to print progress
            seed: Optional integer seed forwarded to GraphColorizer for
                bit-reproducible simulated annealing. When None (default),
                behavior is unchanged.
            initial_coloring: Optional dict mapping node IDs to cell types to
                use as the warm-start coloring instead of a random assignment.
                When None (default), a random initial coloring is used.

        Returns:
            Dictionary mapping node IDs to cell types
        """
        if self.graph is None:
            raise ValueError("Network not built. Call build_network() first.")
        if self.target_statistics is None:
            raise ValueError("Target statistics not loaded. Call load_target_statistics() first.")

        print("\nInitializing GraphColorizer...")
        colorizer = GraphColorizer(
            target_graph=self.graph,
            colors=self.color_names,
            target_statistics=self.target_statistics,
            seed=seed
        )
        
        print("Running simulated annealing to assign cell types...")
        self.cell_type_assignment = colorizer.colorize(
            initial_temp=initial_temp,
            final_temp=final_temp,
            cooling_rate=cooling_rate,
            max_iterations=max_iterations,
            verbose=verbose,
            initial_coloring=initial_coloring
        )
        
        print("Cell type assignment complete!")
        return self.cell_type_assignment

    def generate_colored_replicates(self,
                                    num_replicates: int,
                                    seed: Optional[int] = None,
                                    warm_start: bool = False,
                                    initial_temp: float = 100.0,
                                    final_temp: float = 0.1,
                                    cooling_rate: float = 0.995,
                                    max_iterations: int = 100000,
                                    verbose: bool = True) -> List[Dict]:
        """
        Generate multiple cell-type colorings of the current network.

        Each replicate is an independent simulated-annealing colorization of
        the same target graph against the loaded target statistics, producing
        diverse cell-type assignments that all match the target. This is the
        cell-type-assignment analogue of replicate generation: use it to draw
        several statistically-equivalent-but-distinct labelings of one tissue
        slice.

        By default each replicate **cold-starts** from its own deterministic
        per-replicate seed (derived from ``seed`` via ``numpy.random.SeedSequence``,
        matching ``ReplicateGenerator``), so the replicates are genuinely
        independent. Benchmarks show cold replicates both converge quickly and
        remain diverse.

        Args:
            num_replicates: Number of colorings to generate (>= 1).
            seed: Optional base seed. When provided, replicate ``k`` uses a
                deterministic seed derived from ``(seed, k)``, so the whole set
                of replicates is bit-reproducible. When None, each replicate is
                unseeded.
            warm_start: When True, each replicate after the first warm-starts
                from the previous replicate's coloring instead of cold-starting.
                This speeds convergence on strongly-structured targets but
                **collapses replicate diversity** — subsequent replicates become
                near-identical to the first. Intended for *refining* a coloring,
                not for generating independent replicates. Default False.
            initial_temp: Starting temperature for annealing.
            final_temp: Final temperature for annealing.
            cooling_rate: Rate of temperature decrease.
            max_iterations: Maximum iterations per replicate.
            verbose: Whether to print per-replicate progress.

        Returns:
            List of ``num_replicates`` colorings, each a dict mapping node IDs
            to cell types. ``self.cell_type_assignment`` is also set to the
            first replicate so the downstream ``apply_cell_types_to_slice`` /
            ``evaluate`` helpers operate on a concrete coloring.
        """
        if self.graph is None:
            raise ValueError("Network not built. Call build_network() first.")
        if self.target_statistics is None:
            raise ValueError("Target statistics not loaded. Call load_target_statistics() first.")
        if num_replicates < 1:
            raise ValueError(f"num_replicates must be >= 1, got {num_replicates}.")

        replicates: List[Dict] = []
        previous_coloring: Optional[Dict] = None

        for k in range(num_replicates):
            # Derive a deterministic per-replicate seed from (seed, k) using
            # SeedSequence, matching the convention in ReplicateGenerator. When
            # seed is None each replicate stays unseeded (backward-compatible).
            if seed is not None:
                replicate_seed = int(
                    np.random.SeedSequence([seed, k]).generate_state(1)[0]
                )
            else:
                replicate_seed = None

            if verbose:
                mode = "warm-start" if (warm_start and previous_coloring is not None) else "cold-start"
                print(f"\nColoring replicate {k + 1}/{num_replicates} ({mode})...")

            initial_coloring = previous_coloring if (warm_start and previous_coloring is not None) else None
            coloring = color_graph_to_targets(
                self.graph,
                self.color_names,
                self.target_statistics,
                seed=replicate_seed,
                initial_coloring=initial_coloring,
                initial_temp=initial_temp,
                final_temp=final_temp,
                cooling_rate=cooling_rate,
                max_iterations=max_iterations,
                verbose=verbose,
            )

            replicates.append(coloring)
            previous_coloring = coloring

        # Expose the first replicate as the active assignment so the existing
        # apply/evaluate helpers have a concrete coloring to work with.
        self.cell_type_assignment = replicates[0]
        return replicates

    def apply_cell_types_to_slice(self):
        """
        Apply the assigned cell types back to the slice cells.
        Updates the cell_type attribute of each SliceCell.
        """
        if self.cell_type_assignment is None:
            raise ValueError("Cell types not assigned. Call assign_cell_types() first.")
        if self.slicer is None:
            raise ValueError("Slice not created.")
        
        # Update slice cells with new cell types
        for i, slice_cell in enumerate(self.slicer.slice_cells):
            if i in self.cell_type_assignment:
                slice_cell.cell_type = self.cell_type_assignment[i]
        
        print("Applied cell types to slice cells")
    
    def get_statistics(self) -> Dict:
        """
        Get statistics for the colored graph.
        
        Returns:
            Dictionary of graph statistics
        """
        if self.cell_type_assignment is None:
            raise ValueError("Cell types not assigned. Call assign_cell_types() first.")
        
        return calculate_graph_statistics(
            self.graph, 
            self.cell_type_assignment, 
            self.color_names
        )
    
    def compare_statistics(self, verbose: bool = True) -> Dict:
        """
        Compare generated statistics with target statistics.
        
        Args:
            verbose: Whether to print comparison details
            
        Returns:
            Dictionary of differences
        """
        if self.target_statistics is None:
            raise ValueError("Target statistics not loaded.")
        
        current_stats = self.get_statistics()
        
        # Convert target_statistics format to match current_stats format
        target_stats_formatted = {}
        for color in self.color_names:
            target_stats_formatted[f'nodes_{color}'] = \
                self.target_statistics['node_counts'].get(color, 0)
        
        for key, value in self.target_statistics['edge_counts'].items():
            # key is in format 'color1-color2'
            target_stats_formatted[f'edges_{key}'] = value
        
        return compare_graph_statistics(target_stats_formatted, current_stats, verbose)
    
    def evaluate(self, print_report: bool = True) -> Dict:
        """
        Comprehensive evaluation of the cell type assignment.
        
        Args:
            print_report: Whether to print detailed evaluation report
            
        Returns:
            Dictionary of evaluation metrics
        """
        if self.target_statistics is None:
            raise ValueError("Target statistics not loaded.")
        
        current_stats = self.get_statistics()
        
        # Convert target_statistics format
        target_stats_formatted = {}
        for color in self.color_names:
            target_stats_formatted[f'nodes_{color}'] = \
                self.target_statistics['node_counts'].get(color, 0)
        
        for key, value in self.target_statistics['edge_counts'].items():
            target_stats_formatted[f'edges_{key}'] = value
        
        evaluation = evaluate_graph_coloring(target_stats_formatted, current_stats)
        
        if print_report:
            print_evaluation_report(target_stats_formatted, current_stats, evaluation)
        
        return evaluation
    
    def visualize_slice(self, save_path: str = None, figsize: Tuple[int, int] = (10, 10)):
        """
        Visualize the 2D tissue slice.
        
        Args:
            save_path: Path to save figure (optional)
            figsize: Figure size
        """
        if self.slicer is None:
            raise ValueError("Slice not created.")
        
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Color map for cell types
        if self.color_names:
            cell_types = self.color_names
        else:
            cell_types = list(set(c.cell_type for c in self.slicer.slice_cells))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        # Plot cells as circles
        for slice_cell in self.slicer.slice_cells:
            color = color_map.get(slice_cell.cell_type, 'gray')
            
            circle = plt.Circle(
                slice_cell.center_2d,
                slice_cell.intersection_radius,
                color=color,
                alpha=0.7,
                edgecolor='black',
                linewidth=0.5
            )
            ax.add_patch(circle)
        
        # Set equal aspect ratio
        ax.set_aspect('equal')
        
        # Set limits
        if self.slicer.slice_cells:
            x_coords = [c.center_2d[0] for c in self.slicer.slice_cells]
            y_coords = [c.center_2d[1] for c in self.slicer.slice_cells]
            margin = 20
            ax.set_xlim(min(x_coords) - margin, max(x_coords) + margin)
            ax.set_ylim(min(y_coords) - margin, max(y_coords) + margin)
        
        ax.set_xlabel('U coordinate (μm)')
        ax.set_ylabel('V coordinate (μm)')
        ax.grid(True, alpha=0.3)
        
        # Add legend
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=color_map[ct], markersize=10, label=ct)
            for ct in cell_types
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        ax.set_title(f'2D Tissue Slice: {len(self.slicer.slice_cells)} cells')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved slice visualization to {save_path}")
        
        plt.show()
    
    def visualize_network(self, 
                         layout: str = "spring",
                         save_path: str = None,
                         figsize: Tuple[int, int] = (12, 10)):
        """
        Visualize the network graph.
        
        Args:
            layout: Layout algorithm ("spring", "kamada_kawai", "circular")
            save_path: Path to save figure (optional)
            figsize: Figure size
        """
        if self.graph is None:
            raise ValueError("Network not built.")
        if self.cell_type_assignment is None:
            raise ValueError("Cell types not assigned.")
        
        visualize_colored_graph(
            self.graph,
            self.cell_type_assignment,
            layout=layout,
            title=f"Tissue Network: {self.graph.number_of_nodes()} cells",
            save_path=save_path,
            figsize=figsize
        )
    
    def export_slice_csv(self, filename: str = "tissue_slice.csv", include_3d: bool = True):
        """
        Export slice data to CSV.
        
        Args:
            filename: Output filename
            include_3d: Whether to include 3D coordinates
        """
        if self.slicer is None:
            raise ValueError("Slice not created.")
        
        self.slicer.export_slice_csv(filename, include_3d)
        print(f"Exported slice to {filename}")
    
    def export_network(self, filename: str = "tissue_network.graphml", 
                      format: str = "graphml"):
        """
        Export network graph to file.
        
        Args:
            filename: Output filename
            format: File format ("graphml", "gexf", "gml", "edgelist")
        """
        if self.analyzer is None:
            raise ValueError("Network not built.")
        
        # Set node colors as attributes
        if self.cell_type_assignment:
            nx.set_node_attributes(self.graph, self.cell_type_assignment, 'cell_type')
        
        self.analyzer.export_network(filename, format)
        print(f"Exported network to {filename}")
    
    def export_statistics_csv(self, filename: str = "tissue_statistics.csv"):
        """
        Export statistics to CSV.
        
        Args:
            filename: Output filename
        """
        if self.cell_type_assignment is None:
            raise ValueError("Cell types not assigned.")
        
        export_colored_graph_statistics(
            self.graph,
            self.cell_type_assignment,
            self.color_names,
            filename
        )
    
    def export_all(self, base_dir: str = ".", prefix: str = "tissue"):
        """
        Export all results (slice, network, statistics).
        
        Args:
            base_dir: Base directory for outputs
            prefix: Prefix for filenames
        """
        os.makedirs(base_dir, exist_ok=True)
        
        if self.slicer:
            slice_path = os.path.join(base_dir, f"{prefix}_slice.csv")
            self.export_slice_csv(slice_path)
        
        if self.analyzer:
            network_path = os.path.join(base_dir, f"{prefix}_network.graphml")
            self.export_network(network_path)
        
        if self.cell_type_assignment:
            stats_path = os.path.join(base_dir, f"{prefix}_statistics.csv")
            self.export_statistics_csv(stats_path)
        
        print(f"\nAll results exported to {base_dir}/")
    
    def run_complete_workflow(self,
                            tissue: TissueSection,
                            z_position: float = None,
                            network_radius: float = 50.0,
                            target_stats_file: str = None,
                            target_stats_dict: Dict = None,
                            cell_types: List[str] = None,
                            annealing_params: Dict = None,
                            export_dir: str = "results",
                            visualize: bool = True,
                            seed: Optional[int] = None) -> Dict:
        """
        Run the complete workflow from tissue to cell type assignment.

        Args:
            tissue: TissueSection object
            z_position: Z-position for slice (default: middle of tissue)
            network_radius: Radius for network building (micrometers)
            target_stats_file: Path to target statistics CSV
            target_stats_dict: Pre-computed target statistics
            cell_types: List of cell type names
            annealing_params: Parameters for simulated annealing (optional)
            export_dir: Directory for exports
            visualize: Whether to create visualizations
            seed: Optional integer seed forwarded to ``assign_cell_types`` so
                the graph-coloring step is bit-reproducible. When None
                (default), behavior is unchanged.

        Returns:
            Dictionary with evaluation metrics
        """
        print("\n" + "="*60)
        print("STARTING COMPLETE TISSUE NETWORK WORKFLOW")
        print("="*60)
        
        # Step 1: Set tissue
        print("\n[1/7] Setting tissue...")
        self.set_tissue(tissue)
        
        # Step 2: Create slice
        print("\n[2/7] Creating slice...")
        if z_position is None:
            z_position = tissue.thickness / 2
        self.create_slice(z_position=z_position)
        
        # Step 3: Build network
        print("\n[3/7] Building network...")
        self.build_network(mode="radius", radius=network_radius)
        
        # Step 4: Load target statistics
        print("\n[4/7] Loading target statistics...")
        self.load_target_statistics(
            filepath=target_stats_file,
            statistics=target_stats_dict,
            cell_types=cell_types
        )
        
        # Step 5: Assign cell types
        print("\n[5/7] Assigning cell types...")
        if annealing_params is None:
            annealing_params = {
                'initial_temp': 100.0,
                'final_temp': 0.1,
                'cooling_rate': 0.995,
                'max_iterations': 10000
            }
        self.assign_cell_types(**annealing_params, seed=seed)
        self.apply_cell_types_to_slice()
        
        # Step 6: Visualize
        if visualize:
            print("\n[6/7] Creating visualizations...")
            os.makedirs(export_dir, exist_ok=True)
            
            self.visualize_slice(
                save_path=os.path.join(export_dir, "tissue_slice.png")
            )
            self.visualize_network(
                save_path=os.path.join(export_dir, "tissue_network.png")
            )
        
        # Step 7: Export and evaluate
        print("\n[7/7] Exporting results and evaluating...")
        self.export_all(base_dir=export_dir, prefix="tissue")
        evaluation = self.evaluate(print_report=True)
        
        print("\n" + "="*60)
        print("WORKFLOW COMPLETE!")
        print("="*60)
        
        return evaluation


def quick_workflow(tissue: TissueSection,
                  cell_types: List[str],
                  target_stats_file: str = None,
                  network_radius: float = 50.0,
                  z_position: float = None,
                  output_dir: str = "results",
                  seed: Optional[int] = None) -> TissueNetworkWorkflow:
    """
    Convenience function to run the complete workflow with minimal configuration.

    Args:
        tissue: TissueSection object
        cell_types: List of cell type names
        target_stats_file: Path to target statistics CSV
        network_radius: Radius for network building
        z_position: Z-position for slice (default: middle)
        output_dir: Directory for outputs
        seed: Optional integer seed forwarded to ``run_complete_workflow`` so
            the graph-coloring step is bit-reproducible. When None (default),
            behavior is unchanged.

    Returns:
        TissueNetworkWorkflow object with completed workflow
    """
    workflow = TissueNetworkWorkflow()

    workflow.run_complete_workflow(
        tissue=tissue,
        z_position=z_position,
        network_radius=network_radius,
        target_stats_file=target_stats_file,
        cell_types=cell_types,
        export_dir=output_dir,
        visualize=True,
        seed=seed,
    )

    return workflow
