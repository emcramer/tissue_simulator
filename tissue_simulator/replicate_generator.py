"""
Replicate tissue generator based on spatial statistics.

This module allows generation of random tissue replicates that match specified
spatial interaction statistics. It uses an optimization-based approach to tune
tissue parameters to achieve desired cell-cell interaction patterns.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
import csv
from dataclasses import dataclass, asdict
import warnings
from pathlib import Path

from .tissue import TissueSection, Cell
from .packing import SpherePacker
from .spatial_analysis import SpatialNetworkAnalyzer, InteractionStatistics

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    warnings.warn("NetworkX not installed. Install with: pip install networkx")


@dataclass
class TargetStatistics:
    """
    Target spatial statistics for replicate generation.
    
    Attributes:
        interaction_stats: List of InteractionStatistics defining target patterns
        cell_type_proportions: Dict mapping cell types to target proportions (0-1)
        target_cell_count: Optional target for total cell count
        target_density: Optional target for packing fraction
    """
    interaction_stats: List[InteractionStatistics]
    cell_type_proportions: Optional[Dict[str, float]] = None
    target_cell_count: Optional[int] = None
    target_density: Optional[float] = None
    
    def validate(self):
        """Validate that statistics are consistent."""
        if self.cell_type_proportions:
            total_prop = sum(self.cell_type_proportions.values())
            if not (0.99 <= total_prop <= 1.01):
                raise ValueError(f"Cell type proportions must sum to 1.0, got {total_prop}")
        
        if self.target_density and not (0 < self.target_density < 1):
            raise ValueError(f"Target density must be between 0 and 1, got {self.target_density}")


@dataclass
class ReplicateStatistics:
    """
    Statistics for a generated replicate.
    
    Attributes:
        replicate_id: Unique identifier
        num_cells: Total cell count
        cell_type_counts: Dict of counts per cell type
        packing_fraction: Volume fraction occupied by cells
        interaction_stats: Measured interaction statistics
        divergence_score: Overall divergence from target statistics
    """
    replicate_id: int
    num_cells: int
    cell_type_counts: Dict[str, int]
    packing_fraction: float
    interaction_stats: List[InteractionStatistics]
    divergence_score: float
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        result = asdict(self)
        result['interaction_stats'] = [s.to_dict() for s in self.interaction_stats]
        return result


class ReplicateGenerator:
    """
    Generate tissue replicates matching target spatial statistics.
    
    This class uses an iterative approach to generate tissue samples that
    match specified spatial interaction patterns. It adjusts tissue parameters
    to achieve the desired statistics.
    """
    
    def __init__(self, 
                 target_stats: TargetStatistics,
                 tissue_dimensions: Tuple[float, float, float],
                 base_cell_radii: Dict[str, Tuple[float, float]],
                 network_mode: str = "contact",
                 network_radius: Optional[float] = None,
                 seed: Optional[int] = None):
        """
        Initialize replicate generator.
        
        Args:
            target_stats: Target spatial statistics to match
            tissue_dimensions: (height, width, thickness) in micrometers
            base_cell_radii: Dict mapping cell types to (min_radius, max_radius)
            network_mode: "contact" or "radius" for spatial analysis
            network_radius: Distance threshold if using "radius" mode
            seed: Random seed for reproducibility
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX required for replicate generation")
        
        target_stats.validate()
        
        self.target_stats = target_stats
        self.tissue_dimensions = tissue_dimensions
        self.base_cell_radii = base_cell_radii
        self.network_mode = network_mode
        self.network_radius = network_radius
        self.seed = seed
        
        if seed is not None:
            np.random.seed(seed)
        
        # Extract cell types from target stats
        self.cell_types = set()
        for stat in target_stats.interaction_stats:
            self.cell_types.add(stat.type_a)
            self.cell_types.add(stat.type_b)
        
        # Set default proportions if not provided
        if target_stats.cell_type_proportions is None:
            n_types = len(self.cell_types)
            self.target_stats.cell_type_proportions = {
                ct: 1.0 / n_types for ct in self.cell_types
            }
        
        # Validate cell types match
        config_types = set(base_cell_radii.keys())
        if not self.cell_types.issubset(config_types):
            missing = self.cell_types - config_types
            raise ValueError(f"Cell types in target stats not in radii config: {missing}")
    
    def _compute_interaction_divergence(self, 
                                       measured: List[InteractionStatistics],
                                       target: List[InteractionStatistics]) -> float:
        """
        Compute divergence between measured and target interaction statistics.
        
        Uses relative difference in normalized interactions as the metric.
        
        Args:
            measured: Measured interaction statistics
            target: Target interaction statistics
        
        Returns:
            Average relative divergence (0 = perfect match)
        """
        # Create lookup for measured stats
        measured_dict = {}
        for stat in measured:
            key = tuple(sorted([stat.type_a, stat.type_b]))
            measured_dict[key] = stat
        
        # Compute divergences
        divergences = []
        for target_stat in target:
            key = tuple(sorted([target_stat.type_a, target_stat.type_b]))
            
            if key not in measured_dict:
                # Missing interaction type - large penalty
                divergences.append(1.0)
                continue
            
            measured_stat = measured_dict[key]
            
            # Use normalized interactions as primary metric
            target_val = target_stat.normalized_interactions
            measured_val = measured_stat.normalized_interactions
            
            # Relative difference
            if target_val > 0:
                rel_diff = abs(measured_val - target_val) / target_val
            else:
                rel_diff = abs(measured_val - target_val)
            
            divergences.append(min(rel_diff, 2.0))  # Cap at 2.0
        
        return np.mean(divergences)
    
    def _adjust_cell_type_proportions(self, 
                                     tissue: TissueSection,
                                     iteration: int) -> Dict[str, Tuple[float, float]]:
        """
        Adjust cell type proportions for next iteration.
        
        Args:
            tissue: Current tissue
            iteration: Current iteration number
        
        Returns:
            Updated cell radii configuration
        """
        # Get current proportions
        stats = tissue.get_cell_statistics()
        current_counts = stats['cell_types']
        total_cells = sum(current_counts.values())
        
        current_props = {
            ct: current_counts.get(ct, 0) / total_cells 
            for ct in self.cell_types
        }
        
        # Compute adjustment factors
        adjusted_radii = {}
        for cell_type in self.cell_types:
            target_prop = self.target_stats.cell_type_proportions[cell_type]
            current_prop = current_props[cell_type]
            
            # Adjust radii to change proportions
            # More cells -> smaller radii, fewer cells -> larger radii
            if current_prop > 0:
                adjustment = np.sqrt(target_prop / current_prop)
                adjustment = np.clip(adjustment, 0.7, 1.3)  # Limit adjustments
            else:
                adjustment = 1.1  # Slightly increase if missing
            
            min_r, max_r = self.base_cell_radii[cell_type]
            adjusted_radii[cell_type] = (
                min_r * adjustment,
                max_r * adjustment
            )
        
        return adjusted_radii
    
    def generate_single_replicate(self,
                                 replicate_id: int,
                                 max_attempts: int = 1000,
                                 min_spacing: float = 0.5,
                                 allow_boundary: bool = True,
                                 max_iterations: int = 5,
                                 tolerance: float = 0.15) -> Tuple[TissueSection, ReplicateStatistics]:
        """
        Generate a single tissue replicate.
        
        Iteratively generates tissues and adjusts parameters to approach
        target statistics.
        
        Args:
            replicate_id: Unique identifier for this replicate
            max_attempts: Max attempts for cell packing
            min_spacing: Minimum spacing between cells
            allow_boundary: Allow cells extending beyond bounds
            max_iterations: Max parameter adjustment iterations
            tolerance: Acceptable divergence threshold
        
        Returns:
            Tuple of (TissueSection, ReplicateStatistics)
        """
        best_tissue = None
        best_divergence = float('inf')
        best_stats = None
        
        current_radii = self.base_cell_radii.copy()
        
        for iteration in range(max_iterations):
            # Generate tissue
            tissue = TissueSection(
                height=self.tissue_dimensions[0],
                width=self.tissue_dimensions[1],
                thickness=self.tissue_dimensions[2],
                cell_radii=current_radii
            )
            
            num_cells = tissue.generate_cells(
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary_cells=allow_boundary
            )
            
            if num_cells == 0:
                warnings.warn(f"Replicate {replicate_id}, iteration {iteration}: No cells generated")
                continue
            
            # Analyze spatial interactions
            analyzer = SpatialNetworkAnalyzer()
            analyzer.build_network_from_tissue(
                tissue, 
                mode=self.network_mode,
                radius=self.network_radius
            )
            
            measured_interactions = analyzer.compute_interaction_statistics()
            
            # Compute divergence
            divergence = self._compute_interaction_divergence(
                measured_interactions,
                self.target_stats.interaction_stats
            )
            
            # Track best result
            if divergence < best_divergence:
                best_divergence = divergence
                best_tissue = tissue
                best_stats = measured_interactions
            
            # Check if we've met tolerance
            if divergence <= tolerance:
                break
            
            # Adjust parameters for next iteration
            if iteration < max_iterations - 1:
                current_radii = self._adjust_cell_type_proportions(tissue, iteration)
        
        # Create statistics object
        if best_tissue is None:
            raise RuntimeError(f"Failed to generate replicate {replicate_id}")
        
        tissue_stats = best_tissue.get_cell_statistics()
        
        replicate_stats = ReplicateStatistics(
            replicate_id=replicate_id,
            num_cells=tissue_stats['total_cells'],
            cell_type_counts=tissue_stats['cell_types'],
            packing_fraction=tissue_stats['packing_fraction'],
            interaction_stats=best_stats,
            divergence_score=best_divergence
        )
        
        return best_tissue, replicate_stats
    
    def generate_replicates(self,
                           num_replicates: int,
                           max_attempts: int = 1000,
                           min_spacing: float = 0.5,
                           allow_boundary: bool = True,
                           max_iterations: int = 5,
                           tolerance: float = 0.15,
                           parallel: bool = False) -> List[Tuple[TissueSection, ReplicateStatistics]]:
        """
        Generate multiple tissue replicates.
        
        Args:
            num_replicates: Number of replicates to generate
            max_attempts: Max attempts for cell packing per replicate
            min_spacing: Minimum spacing between cells
            allow_boundary: Allow cells extending beyond bounds
            max_iterations: Max parameter adjustment iterations
            tolerance: Acceptable divergence threshold
            parallel: Whether to use parallel processing (not yet implemented)
        
        Returns:
            List of (TissueSection, ReplicateStatistics) tuples
        """
        replicates = []
        
        for i in range(num_replicates):
            print(f"Generating replicate {i+1}/{num_replicates}...")
            
            tissue, stats = self.generate_single_replicate(
                replicate_id=i,
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary=allow_boundary,
                max_iterations=max_iterations,
                tolerance=tolerance
            )
            
            replicates.append((tissue, stats))
            
            print(f"  Cells: {stats.num_cells}, Divergence: {stats.divergence_score:.4f}")
        
        return replicates
    
    def export_replicate_statistics(self,
                                   replicates: List[Tuple[TissueSection, ReplicateStatistics]],
                                   base_filename: str):
        """
        Export statistics for all replicates to CSV files.
        
        Args:
            replicates: List of (tissue, stats) tuples
            base_filename: Base name for output files
        """
        # Summary statistics
        summary_data = []
        for tissue, stats in replicates:
            row = {
                'replicate_id': stats.replicate_id,
                'num_cells': stats.num_cells,
                'packing_fraction': stats.packing_fraction,
                'divergence_score': stats.divergence_score
            }
            # Add cell type counts
            for ct, count in stats.cell_type_counts.items():
                row[f'count_{ct}'] = count
                row[f'proportion_{ct}'] = count / stats.num_cells
            
            summary_data.append(row)
        
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv(f"{base_filename}_summary.csv", index=False)
        
        # Interaction statistics
        interaction_data = []
        for tissue, stats in replicates:
            for interaction in stats.interaction_stats:
                row = {
                    'replicate_id': stats.replicate_id,
                    'type_a': interaction.type_a,
                    'type_b': interaction.type_b,
                    'num_interactions': interaction.num_interactions,
                    'normalized_interactions': interaction.normalized_interactions,
                    'avg_distance': interaction.avg_distance,
                    'median_distance': interaction.median_distance
                }
                interaction_data.append(row)
        
        df_interactions = pd.DataFrame(interaction_data)
        df_interactions.to_csv(f"{base_filename}_interactions.csv", index=False)
        
        print(f"Exported statistics to {base_filename}_summary.csv and {base_filename}_interactions.csv")
    
    def export_replicate_tissues(self,
                                replicates: List[Tuple[TissueSection, ReplicateStatistics]],
                                output_dir: str):
        """
        Export each replicate tissue to a separate CSV file.
        
        Args:
            replicates: List of (tissue, stats) tuples
            output_dir: Directory for output files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for tissue, stats in replicates:
            filename = output_path / f"replicate_{stats.replicate_id:03d}_tissue.csv"
            tissue.export_to_csv(str(filename))
        
        print(f"Exported {len(replicates)} tissue files to {output_dir}")


def load_target_statistics_from_csv(filepath: str) -> TargetStatistics:
    """
    Load target statistics from a CSV file.
    
    Expected CSV format:
    type_a, type_b, num_interactions, normalized_interactions, avg_distance, median_distance
    
    Args:
        filepath: Path to CSV file
    
    Returns:
        TargetStatistics object
    """
    df = pd.read_csv(filepath)
    
    required_cols = ['type_a', 'type_b', 'normalized_interactions']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    
    # Create InteractionStatistics objects
    interaction_stats = []
    for _, row in df.iterrows():
        stat = InteractionStatistics(
            type_a=str(row['type_a']),
            type_b=str(row['type_b']),
            num_interactions=int(row.get('num_interactions', 0)),
            normalized_interactions=float(row['normalized_interactions']),
            avg_distance=float(row.get('avg_distance', 0.0)),
            median_distance=float(row.get('median_distance', 0.0))
        )
        interaction_stats.append(stat)
    
    return TargetStatistics(interaction_stats=interaction_stats)


def load_target_statistics_from_tissue(tissue: TissueSection,
                                      network_mode: str = "contact",
                                      network_radius: Optional[float] = None) -> TargetStatistics:
    """
    Extract target statistics from an existing tissue.
    
    Args:
        tissue: TissueSection to analyze
        network_mode: "contact" or "radius"
        network_radius: Distance threshold for "radius" mode
    
    Returns:
        TargetStatistics object
    """
    # Analyze the tissue
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_tissue(
        tissue,
        mode=network_mode,
        radius=network_radius
    )
    
    # Get interaction statistics
    interaction_stats = analyzer.compute_interaction_statistics()
    
    # Get cell type proportions
    tissue_stats = tissue.get_cell_statistics()
    total_cells = tissue_stats['total_cells']
    cell_type_proportions = {
        ct: count / total_cells
        for ct, count in tissue_stats['cell_types'].items()
    }
    
    return TargetStatistics(
        interaction_stats=interaction_stats,
        cell_type_proportions=cell_type_proportions,
        target_cell_count=total_cells,
        target_density=tissue_stats['packing_fraction']
    )
