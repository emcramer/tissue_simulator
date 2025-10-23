"""
Core tissue section and cell classes.
"""

import numpy as np
import csv
from typing import Dict, Tuple, List, Union, Optional


class Cell:
    """
    Represents a single cell in the tissue section.
    
    Attributes:
        center: (x, y, z) coordinates of cell center
        radius: Cell radius
        cell_type: Type/classification of the cell
        is_boundary: Whether cell extends beyond tissue bounds
    """
    
    def __init__(self, center: Tuple[float, float, float], 
                 radius: float, cell_type: str = "default",
                 is_boundary: bool = False):
        self.center = np.array(center)
        self.radius = radius
        self.cell_type = cell_type
        self.is_boundary = is_boundary
        
    def intersects(self, other: 'Cell') -> bool:
        """Check if this cell intersects with another cell."""
        distance = np.linalg.norm(self.center - other.center)
        return distance < (self.radius + other.radius)
    
    def is_within_bounds(self, bounds: Tuple[float, float, float]) -> bool:
        """
        Check if cell is completely within tissue bounds.
        
        Args:
            bounds: (height, width, thickness) of tissue section
        """
        height, width, thickness = bounds
        x, y, z = self.center
        
        within_x = (x - self.radius >= 0) and (x + self.radius <= width)
        within_y = (y - self.radius >= 0) and (y + self.radius <= height)
        within_z = (z - self.radius >= 0) and (z + self.radius <= thickness)
        
        return within_x and within_y and within_z
    
    def intersects_bounds(self, bounds: Tuple[float, float, float]) -> bool:
        """
        Check if cell center is within bounds (allowing partial overlap).
        
        Args:
            bounds: (height, width, thickness) of tissue section
        """
        height, width, thickness = bounds
        x, y, z = self.center
        
        within_x = 0 <= x <= width
        within_y = 0 <= y <= height
        within_z = 0 <= z <= thickness
        
        return within_x and within_y and within_z
    
    def __repr__(self) -> str:
        return (f"Cell(center={self.center}, radius={self.radius}, "
                f"type={self.cell_type}, boundary={self.is_boundary})")


class TissueSection:
    """
    Represents a 3D tissue section with packed cells.
    
    Attributes:
        height: Y-dimension of tissue (micrometers)
        width: X-dimension of tissue (micrometers)
        thickness: Z-dimension of tissue (micrometers)
        cell_radii: Dictionary mapping cell types to (min_radius, max_radius)
        cells: List of Cell objects in the tissue
    """
    
    def __init__(self, height: float, width: float, thickness: float,
                 cell_radii: Union[Tuple[float, float], 
                                  Dict[str, Tuple[float, float]]]):
        """
        Initialize a tissue section.
        
        Args:
            height: Y-dimension in micrometers
            width: X-dimension in micrometers
            thickness: Z-dimension in micrometers
            cell_radii: Either a tuple (min, max) for uniform cells,
                       or dict mapping cell types to (min, max) radii
        """
        self.height = height
        self.width = width
        self.thickness = thickness
        
        # Convert single tuple to dictionary format
        if isinstance(cell_radii, tuple):
            self.cell_radii = {"default": cell_radii}
        else:
            self.cell_radii = cell_radii
            
        self.cells: List[Cell] = []
        self._rng = np.random.default_rng()
        
    def get_bounds(self) -> Tuple[float, float, float]:
        """Return tissue dimensions as (height, width, thickness)."""
        return (self.height, self.width, self.thickness)
    
    def _select_cell_type_and_radius(self) -> Tuple[str, float]:
        """Randomly select a cell type and radius based on configuration."""
        # Randomly select a cell type
        cell_type = self._rng.choice(list(self.cell_radii.keys()))
        min_radius, max_radius = self.cell_radii[cell_type]
        
        # Generate random radius within range
        radius = self._rng.uniform(min_radius, max_radius)
        
        return cell_type, radius
    
    def generate_cells(self, max_attempts: int = 1000, 
                      min_spacing: float = 0.5,
                      allow_boundary_cells: bool = True) -> int:
        """
        Generate cells using random sphere packing.
        
        Args:
            max_attempts: Maximum placement attempts before stopping
            min_spacing: Minimum spacing between cell surfaces
            allow_boundary_cells: If True, allow cells that extend beyond bounds
            
        Returns:
            Number of cells successfully placed
        """
        from .packing import SpherePacker
        
        packer = SpherePacker(
            bounds=self.get_bounds(),
            cell_radii_config=self.cell_radii,
            min_spacing=min_spacing,
            allow_boundary_cells=allow_boundary_cells
        )
        
        self.cells = packer.pack(max_attempts=max_attempts)
        
        return len(self.cells)
    
    def get_cell_statistics(self) -> Dict:
        """Calculate statistics about the packed cells."""
        if not self.cells:
            return {"total_cells": 0}
        
        stats = {
            "total_cells": len(self.cells),
            "boundary_cells": sum(1 for c in self.cells if c.is_boundary),
            "interior_cells": sum(1 for c in self.cells if not c.is_boundary),
        }
        
        # Cell type breakdown
        type_counts = {}
        type_radii = {}
        for cell in self.cells:
            if cell.cell_type not in type_counts:
                type_counts[cell.cell_type] = 0
                type_radii[cell.cell_type] = []
            type_counts[cell.cell_type] += 1
            type_radii[cell.cell_type].append(cell.radius)
        
        stats["cell_types"] = type_counts
        
        # Average radii per type
        stats["avg_radii"] = {
            cell_type: np.mean(radii) 
            for cell_type, radii in type_radii.items()
        }
        
        # Volume calculations
        tissue_volume = self.height * self.width * self.thickness
        cell_volume = sum((4/3) * np.pi * c.radius**3 for c in self.cells)
        stats["packing_fraction"] = cell_volume / tissue_volume
        
        return stats
    
    def export_to_csv(self, filename: str):
        """
        Export cell data to CSV file.
        
        Args:
            filename: Output CSV file path
        """
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['x', 'y', 'z', 'radius', 'cell_type', 'is_boundary'])
            
            for cell in self.cells:
                writer.writerow([
                    cell.center[0],
                    cell.center[1],
                    cell.center[2],
                    cell.radius,
                    cell.cell_type,
                    cell.is_boundary
                ])
    
    def visualize(self, show_boundary: bool = True, 
                 elevation: float = 20, azimuth: float = 45):
        """
        Create a 3D visualization of the tissue section.
        
        Args:
            show_boundary: Whether to show boundary box
            elevation: Viewing elevation angle
            azimuth: Viewing azimuth angle
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Color map for cell types
        cell_types = list(set(c.cell_type for c in self.cells))
        colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        # Plot cells
        for cell in self.cells:
            # Create sphere surface
            u = np.linspace(0, 2 * np.pi, 20)
            v = np.linspace(0, np.pi, 20)
            x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
            y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
            z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
            
            color = color_map[cell.cell_type]
            alpha = 0.3 if cell.is_boundary else 0.6
            
            ax.plot_surface(x, y, z, color=color, alpha=alpha, 
                          edgecolors='none', shade=True)
        
        # Draw boundary box
        if show_boundary:
            # Define box edges
            edges = [
                [[0, self.width], [0, 0], [0, 0]],
                [[0, self.width], [self.height, self.height], [0, 0]],
                [[0, self.width], [0, 0], [self.thickness, self.thickness]],
                [[0, self.width], [self.height, self.height], [self.thickness, self.thickness]],
                [[0, 0], [0, self.height], [0, 0]],
                [[self.width, self.width], [0, self.height], [0, 0]],
                [[0, 0], [0, self.height], [self.thickness, self.thickness]],
                [[self.width, self.width], [0, self.height], [self.thickness, self.thickness]],
                [[0, 0], [0, 0], [0, self.thickness]],
                [[self.width, self.width], [0, 0], [0, self.thickness]],
                [[0, 0], [self.height, self.height], [0, self.thickness]],
                [[self.width, self.width], [self.height, self.height], [0, self.thickness]],
            ]
            
            for edge in edges:
                ax.plot3D(*edge, 'k-', linewidth=1, alpha=0.3)
        
        # Set labels and limits
        ax.set_xlabel('Width (μm)')
        ax.set_ylabel('Height (μm)')
        ax.set_zlabel('Thickness (μm)')
        ax.set_xlim(0, self.width)
        ax.set_ylim(0, self.height)
        ax.set_zlim(0, self.thickness)
        
        # Set viewing angle
        ax.view_init(elev=elevation, azim=azimuth)
        
        # Add legend
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor=color_map[ct], markersize=10, label=ct)
            for ct in cell_types
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        ax.set_title(f'3D Tissue Section: {len(self.cells)} cells')
        plt.tight_layout()
        plt.show()
    
    def clear_cells(self):
        """Remove all cells from the tissue section."""
        self.cells = []
