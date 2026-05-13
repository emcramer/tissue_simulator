"""
Sphere packing algorithm for cell placement.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from .tissue import Cell


class SpherePacker:
    """
    Random sphere packing algorithm for placing cells in tissue.
    """

    def __init__(self, bounds: Tuple[float, float, float],
                 cell_radii_config: Dict[str, Tuple[float, float]],
                 min_spacing: float = 0.5,
                 allow_boundary_cells: bool = True,
                 seed: Optional[int] = None):
        """
        Initialize sphere packer.

        Args:
            bounds: (height, width, thickness) of tissue
            cell_radii_config: Dict mapping cell types to (min, max) radii
            min_spacing: Minimum spacing between cell surfaces
            allow_boundary_cells: Allow cells extending beyond bounds
            seed: Optional integer seed for the instance RNG. When provided,
                the packing process is deterministic; when None the RNG is
                seeded from system entropy (previous default behavior).
        """
        self.bounds = bounds
        self.cell_radii_config = cell_radii_config
        self.min_spacing = min_spacing
        self.allow_boundary_cells = allow_boundary_cells
        self.seed = seed
        self._rng = np.random.default_rng(seed)
        
    def _select_cell_type_and_radius(self) -> Tuple[str, float]:
        """Randomly select a cell type and radius."""
        cell_type = self._rng.choice(list(self.cell_radii_config.keys()))
        min_radius, max_radius = self.cell_radii_config[cell_type]
        radius = self._rng.uniform(min_radius, max_radius)
        return cell_type, radius
    
    def _generate_random_position(self, radius: float) -> np.ndarray:
        """
        Generate a random position for a cell.
        
        Args:
            radius: Cell radius
            
        Returns:
            3D position array [x, y, z]
        """
        height, width, thickness = self.bounds
        
        if self.allow_boundary_cells:
            # Allow centers anywhere in the tissue (cells can extend beyond)
            x = self._rng.uniform(0, width)
            y = self._rng.uniform(0, height)
            z = self._rng.uniform(0, thickness)
        else:
            # Constrain centers so cells stay within bounds
            x = self._rng.uniform(radius, width - radius)
            y = self._rng.uniform(radius, height - radius)
            z = self._rng.uniform(radius, thickness - radius)
        
        return np.array([x, y, z])
    
    def _check_collision(self, candidate: Cell, existing_cells: List[Cell]) -> bool:
        """
        Check if candidate cell collides with existing cells.
        
        Args:
            candidate: Cell to test
            existing_cells: List of already placed cells
            
        Returns:
            True if collision detected, False otherwise
        """
        for cell in existing_cells:
            distance = np.linalg.norm(candidate.center - cell.center)
            min_distance = candidate.radius + cell.radius + self.min_spacing
            
            if distance < min_distance:
                return True
        
        return False
    
    def _is_valid_placement(self, cell: Cell) -> bool:
        """
        Check if cell placement is valid given boundary constraints.
        
        Args:
            cell: Cell to validate
            
        Returns:
            True if valid, False otherwise
        """
        if self.allow_boundary_cells:
            # Only require center to be within bounds
            return cell.intersects_bounds(self.bounds)
        else:
            # Require entire cell to be within bounds
            return cell.is_within_bounds(self.bounds)
    
    def pack(self, max_attempts: int = 1000) -> List[Cell]:
        """
        Pack cells into tissue using random placement.
        
        Args:
            max_attempts: Maximum placement attempts before stopping
            
        Returns:
            List of successfully placed Cell objects
        """
        cells = []
        failed_attempts = 0
        total_attempts = 0
        
        while failed_attempts < max_attempts:
            total_attempts += 1
            
            # Select cell type and radius
            cell_type, radius = self._select_cell_type_and_radius()
            
            # Generate random position
            position = self._generate_random_position(radius)
            
            # Create candidate cell
            candidate = Cell(
                center=position,
                radius=radius,
                cell_type=cell_type
            )
            
            # Check if placement is valid
            if not self._is_valid_placement(candidate):
                failed_attempts += 1
                continue
            
            # Check for collisions
            if self._check_collision(candidate, cells):
                failed_attempts += 1
                continue
            
            # Mark as boundary cell if it extends beyond bounds
            candidate.is_boundary = not candidate.is_within_bounds(self.bounds)
            
            # Add cell and reset failure counter
            cells.append(candidate)
            failed_attempts = 0
        
        return cells
    
    def pack_with_progress(self, max_attempts: int = 1000, 
                          callback=None) -> List[Cell]:
        """
        Pack cells with progress callback for GUI updates.
        
        Args:
            max_attempts: Maximum placement attempts before stopping
            callback: Function called with (cells_placed, total_attempts)
            
        Returns:
            List of successfully placed Cell objects
        """
        cells = []
        failed_attempts = 0
        total_attempts = 0
        
        while failed_attempts < max_attempts:
            total_attempts += 1
            
            # Select cell type and radius
            cell_type, radius = self._select_cell_type_and_radius()
            
            # Generate random position
            position = self._generate_random_position(radius)
            
            # Create candidate cell
            candidate = Cell(
                center=position,
                radius=radius,
                cell_type=cell_type
            )
            
            # Check if placement is valid
            if not self._is_valid_placement(candidate):
                failed_attempts += 1
                continue
            
            # Check for collisions
            if self._check_collision(candidate, cells):
                failed_attempts += 1
                continue
            
            # Mark as boundary cell if it extends beyond bounds
            candidate.is_boundary = not candidate.is_within_bounds(self.bounds)
            
            # Add cell and reset failure counter
            cells.append(candidate)
            failed_attempts = 0
            
            # Call progress callback
            if callback and len(cells) % 10 == 0:
                callback(len(cells), total_attempts)
        
        return cells
