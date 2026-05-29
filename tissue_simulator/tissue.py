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
                                  Dict[str, Tuple[float, float]]],
                 seed: Optional[int] = None):
        """
        Initialize a tissue section.

        Args:
            height: Y-dimension in micrometers
            width: X-dimension in micrometers
            thickness: Z-dimension in micrometers
            cell_radii: Either a tuple (min, max) for uniform cells,
                       or dict mapping cell types to (min, max) radii
            seed: Optional integer seed controlling randomness for this
                tissue (both the internal sampler and the SpherePacker
                created by ``generate_cells``). When None (default), the
                RNG is seeded from system entropy as before.
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
        self.seed = seed

    def get_bounds(self) -> Tuple[float, float, float]:
        """Return tissue dimensions as (height, width, thickness)."""
        return (self.height, self.width, self.thickness)

    def generate_cells(self, max_attempts: int = 1000,
                      min_spacing: float = 0.5,
                      allow_boundary_cells: bool = True,
                      seed: Optional[int] = None) -> int:
        """
        Generate cells using random sphere packing.

        Args:
            max_attempts: Maximum placement attempts before stopping
            min_spacing: Minimum spacing between cell surfaces
            allow_boundary_cells: If True, allow cells that extend beyond bounds
            seed: Optional integer seed that overrides ``self.seed`` for this
                call. When provided, the TissueSection's internal RNG is
                re-constructed with this seed and the seed is propagated into
                the internal SpherePacker, making generation bit-reproducible.
                When None (default), ``self.seed`` is used (which itself may
                be None for unseeded behavior).

        Returns:
            Number of cells successfully placed
        """
        from .packing import SpherePacker

        # Resolve the effective seed for this generation run. An explicit
        # ``seed=`` argument overrides ``self.seed``; otherwise we reuse the
        # seed supplied at construction.
        if seed is not None:
            self.seed = seed
        effective_seed = self.seed

        packer = SpherePacker(
            bounds=self.get_bounds(),
            cell_radii_config=self.cell_radii,
            min_spacing=min_spacing,
            allow_boundary_cells=allow_boundary_cells,
            seed=effective_seed,
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
    
    @classmethod
    def from_cells(cls, cells: List["Cell"], height: Optional[float] = None,
                   width: Optional[float] = None, thickness: Optional[float] = None,
                   cell_radii: Optional[Dict[str, Tuple[float, float]]] = None) -> "TissueSection":
        """
        Wrap a collection of pre-positioned cells in a TissueSection.

        This is the inverse complement to cell generation: instead of packing
        new cells into an empty tissue, it adopts cells that already have
        positions (e.g. imported from an external source or a CSV file) so the
        spatial-analysis and replicate-generation API can run on externally
        sourced tissue.

        Dimension inference:
            Any dimension left as None is inferred from the bounding box of the
            cell CENTERS along its axis (span = max - min): width from axis 0
            (x), height from axis 1 (y), thickness from axis 2 (z). If a span is
            0 -- as happens for a 2D slice where every z is equal -- the
            dimension falls back to the largest cell diameter
            ``2 * max(c.radius for c in cells)`` so it stays strictly positive
            (this keeps ``get_cell_statistics()["packing_fraction"]`` within
            (0, 1)). A dimension passed in explicitly is used as-is.

        Args:
            cells: Non-empty list of Cell objects with positions already set.
            height: Optional Y-dimension; inferred from cell centers if None.
            width: Optional X-dimension; inferred from cell centers if None.
            thickness: Optional Z-dimension; inferred from cell centers if None.
            cell_radii: Optional mapping of cell type to (min, max) radius. When
                None, it is derived by grouping the cells on cell_type and
                mapping each type to its observed (min_radius, max_radius).

        Returns:
            A TissueSection whose ``cells`` are the provided cells.

        Raises:
            ValueError: If ``cells`` is empty.
        """
        if not cells:
            raise ValueError("Cannot create a TissueSection from an empty cell list.")

        # Derive cell_radii from observed radii per cell type when not provided.
        if cell_radii is None:
            type_radii: Dict[str, List[float]] = {}
            for cell in cells:
                type_radii.setdefault(cell.cell_type, []).append(cell.radius)
            cell_radii = {
                cell_type: (min(radii), max(radii))
                for cell_type, radii in type_radii.items()
            }

        # Fall-back dimension for zero-span axes keeps the tissue volume
        # strictly positive (largest cell diameter).
        fallback = 2 * max(c.radius for c in cells)

        def _infer(axis: int) -> float:
            coords = [c.center[axis] for c in cells]
            span = max(coords) - min(coords)
            return span if span > 0 else fallback

        if width is None:
            width = _infer(0)
        if height is None:
            height = _infer(1)
        if thickness is None:
            thickness = _infer(2)

        tissue = cls(height=height, width=width, thickness=thickness,
                     cell_radii=cell_radii, seed=None)
        tissue.cells = list(cells)
        return tissue

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
        from ._viz_utils import make_color_map

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Color map for cell types (sorted for deterministic assignment)
        color_map = make_color_map(c.cell_type for c in self.cells)
        cell_types = list(color_map.keys())
        
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
            
            # Convert color to proper format for plot_surface
            ax.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)), 
                          alpha=alpha, linewidth=0, antialiased=True, shade=False)
        
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


def load_tissue_from_csv(filepath: str, height: Optional[float] = None,
                         width: Optional[float] = None, thickness: Optional[float] = None,
                         default_radius: float = 10.0) -> TissueSection:
    """
    Load a tissue section from a CSV file.

    This is the exact inverse of ``TissueSection.export_to_csv``: it reads the
    rows that ``export_to_csv`` writes back into Cell objects and wraps them in
    a TissueSection (via ``TissueSection.from_cells``), inferring any omitted
    dimensions from the cell positions.

    The CSV is read with ``csv.DictReader``. The expected columns are
    ``x``, ``y``, ``z``, ``radius``, ``cell_type`` and ``is_boundary``. Only the
    coordinate columns are required: ``radius``, ``cell_type`` and
    ``is_boundary`` are optional. When ``radius`` is missing or blank for a row,
    ``default_radius`` is used; a missing ``cell_type`` defaults to "default";
    and ``is_boundary`` is treated as True only when its value is the string
    "true" (case-insensitive), otherwise False.

    ``cell_type`` is the canonical column name.  As a convenience, the column
    ``type`` (a common alias used by tools such as PhysiCell) is also accepted
    when ``cell_type`` is absent or blank.  When both columns are present,
    ``cell_type`` takes precedence.

    Args:
        filepath: Path to the CSV file to read.
        height: Optional Y-dimension; inferred from cell centers if None.
        width: Optional X-dimension; inferred from cell centers if None.
        thickness: Optional Z-dimension; inferred from cell centers if None.
        default_radius: Radius assigned to rows whose ``radius`` column is
            missing or empty.

    Returns:
        A TissueSection containing the loaded cells.
    """
    cells: List[Cell] = []
    with open(filepath, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            center = (float(row['x']), float(row['y']), float(row['z']))

            radius_str = row.get('radius')
            if radius_str is not None and str(radius_str).strip() != '':
                radius = float(radius_str)
            else:
                radius = default_radius

            cell_type = row.get('cell_type') or row.get('type') or 'default'
            is_boundary = str(row.get('is_boundary')).strip().lower() == 'true'

            cells.append(Cell(center=center, radius=radius,
                              cell_type=cell_type, is_boundary=is_boundary))

    return TissueSection.from_cells(cells, height=height, width=width,
                                    thickness=thickness)
