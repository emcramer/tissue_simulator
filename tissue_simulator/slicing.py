"""
2D tissue slicing module for extracting planar sections through 3D tissue.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import csv

from .tissue import Cell, TissueSection


@dataclass
class SliceCell:
    """
    Represents a cell captured in a 2D slice.
    
    Attributes:
        center_3d: Original 3D center position
        center_2d: Position in 2D slice coordinates
        radius: Cell radius
        cell_type: Cell classification
        is_boundary: Whether cell was a boundary cell in 3D tissue
        distance_from_plane: Perpendicular distance from slice plane
        intersection_radius: Radius of circular intersection with plane
    """
    center_3d: np.ndarray
    center_2d: np.ndarray
    radius: float
    cell_type: str
    is_boundary: bool
    distance_from_plane: float
    intersection_radius: float


class TissueSlicer:
    """
    Handles 2D slicing operations on 3D tissue sections.
    """
    
    def __init__(self, tissue: TissueSection):
        """
        Initialize slicer with a tissue section.
        
        Args:
            tissue: TissueSection object to slice
        """
        self.tissue = tissue
        self.slice_cells: List[SliceCell] = []
        self.slice_plane_point: Optional[np.ndarray] = None
        self.slice_plane_normal: Optional[np.ndarray] = None
        self.slice_basis_u: Optional[np.ndarray] = None
        self.slice_basis_v: Optional[np.ndarray] = None
    
    def slice_plane(self, 
                   point: Tuple[float, float, float] = None,
                   normal: Tuple[float, float, float] = None,
                   angle_x: float = 0.0,
                   angle_y: float = 0.0,
                   z_position: float = None) -> List[SliceCell]:
        """
        Create a 2D slice through the tissue.
        
        Args:
            point: Point on the plane (default: center of tissue)
            normal: Normal vector to the plane (default: [0, 0, 1] for XY plane)
            angle_x: Rotation angle around X-axis in degrees (alternative to normal)
            angle_y: Rotation angle around Y-axis in degrees (alternative to normal)
            z_position: Z-position for horizontal slice (simplified interface)
            
        Returns:
            List of SliceCell objects intersecting the plane
        """
        # Handle simplified interface for horizontal slices
        if z_position is not None:
            point = (self.tissue.width / 2, self.tissue.height / 2, z_position)
            normal = (0, 0, 1)
        
        # Set default point (center of tissue)
        if point is None:
            point = (self.tissue.width / 2, 
                    self.tissue.height / 2, 
                    self.tissue.thickness / 2)
        
        self.slice_plane_point = np.array(point)
        
        # Set or calculate normal vector
        if normal is not None:
            self.slice_plane_normal = np.array(normal)
        else:
            # Calculate normal from rotation angles
            self.slice_plane_normal = self._normal_from_angles(angle_x, angle_y)
        
        # Normalize the normal vector
        self.slice_plane_normal = self.slice_plane_normal / np.linalg.norm(self.slice_plane_normal)
        
        # Create orthonormal basis for the plane
        self._create_plane_basis()
        
        # Find all cells intersecting the plane
        self.slice_cells = []
        for cell in self.tissue.cells:
            slice_cell = self._intersect_cell_with_plane(cell)
            if slice_cell is not None:
                self.slice_cells.append(slice_cell)
        
        return self.slice_cells
    
    def _normal_from_angles(self, angle_x: float, angle_y: float) -> np.ndarray:
        """
        Calculate normal vector from rotation angles.
        
        Args:
            angle_x: Rotation around X-axis (degrees)
            angle_y: Rotation around Y-axis (degrees)
            
        Returns:
            Normal vector
        """
        # Convert to radians
        ax = np.radians(angle_x)
        ay = np.radians(angle_y)
        
        # Start with Z-axis pointing up
        normal = np.array([0.0, 0.0, 1.0])
        
        # Rotation matrix around X-axis
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(ax), -np.sin(ax)],
            [0, np.sin(ax), np.cos(ax)]
        ])
        
        # Rotation matrix around Y-axis
        Ry = np.array([
            [np.cos(ay), 0, np.sin(ay)],
            [0, 1, 0],
            [-np.sin(ay), 0, np.cos(ay)]
        ])
        
        # Apply rotations
        normal = Ry @ Rx @ normal
        
        return normal
    
    def _create_plane_basis(self):
        """Create orthonormal basis vectors for the slice plane."""
        # Find a vector not parallel to normal
        if abs(self.slice_plane_normal[2]) < 0.9:
            arbitrary = np.array([0, 0, 1])
        else:
            arbitrary = np.array([1, 0, 0])
        
        # Create orthonormal basis using Gram-Schmidt
        self.slice_basis_u = np.cross(self.slice_plane_normal, arbitrary)
        self.slice_basis_u = self.slice_basis_u / np.linalg.norm(self.slice_basis_u)
        
        self.slice_basis_v = np.cross(self.slice_plane_normal, self.slice_basis_u)
        self.slice_basis_v = self.slice_basis_v / np.linalg.norm(self.slice_basis_v)
    
    def _intersect_cell_with_plane(self, cell: Cell) -> Optional[SliceCell]:
        """
        Check if a cell intersects with the slice plane.
        
        Args:
            cell: Cell to check
            
        Returns:
            SliceCell if intersection exists, None otherwise
        """
        # Calculate distance from cell center to plane
        vec_to_point = cell.center - self.slice_plane_point
        distance = np.dot(vec_to_point, self.slice_plane_normal)
        
        # Check if plane intersects the sphere
        if abs(distance) > cell.radius:
            return None
        
        # Calculate radius of circular intersection
        intersection_radius = np.sqrt(cell.radius**2 - distance**2)
        
        # Project cell center onto plane
        projected_center = cell.center - distance * self.slice_plane_normal
        
        # Convert to 2D plane coordinates
        center_2d = self._project_to_2d(projected_center)
        
        return SliceCell(
            center_3d=cell.center.copy(),
            center_2d=center_2d,
            radius=cell.radius,
            cell_type=cell.cell_type,
            is_boundary=cell.is_boundary,
            distance_from_plane=abs(distance),
            intersection_radius=intersection_radius
        )
    
    def _project_to_2d(self, point_3d: np.ndarray) -> np.ndarray:
        """
        Project a 3D point onto the 2D slice plane.
        
        Args:
            point_3d: 3D point to project
            
        Returns:
            2D coordinates in plane basis
        """
        vec_from_origin = point_3d - self.slice_plane_point
        u_coord = np.dot(vec_from_origin, self.slice_basis_u)
        v_coord = np.dot(vec_from_origin, self.slice_basis_v)
        
        return np.array([u_coord, v_coord])
    
    def export_slice_csv(self, filename: str, include_3d: bool = True):
        """
        Export slice data to CSV file.
        
        Args:
            filename: Output CSV file path
            include_3d: Whether to include original 3D coordinates
        """
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            if include_3d:
                writer.writerow([
                    'x_2d', 'y_2d', 'intersection_radius',
                    'x_3d', 'y_3d', 'z_3d', 'radius_3d',
                    'cell_type', 'is_boundary', 'distance_from_plane'
                ])
            else:
                writer.writerow([
                    'x_2d', 'y_2d', 'intersection_radius',
                    'cell_type', 'distance_from_plane'
                ])
            
            # Write data
            for slice_cell in self.slice_cells:
                if include_3d:
                    writer.writerow([
                        slice_cell.center_2d[0],
                        slice_cell.center_2d[1],
                        slice_cell.intersection_radius,
                        slice_cell.center_3d[0],
                        slice_cell.center_3d[1],
                        slice_cell.center_3d[2],
                        slice_cell.radius,
                        slice_cell.cell_type,
                        slice_cell.is_boundary,
                        slice_cell.distance_from_plane
                    ])
                else:
                    writer.writerow([
                        slice_cell.center_2d[0],
                        slice_cell.center_2d[1],
                        slice_cell.intersection_radius,
                        slice_cell.cell_type,
                        slice_cell.distance_from_plane
                    ])
    
    def visualize_slice_2d(self, show_radii: bool = True, 
                          figsize: Tuple[int, int] = (10, 10),
                          title: str = None):
        """
        Visualize the 2D slice showing cell cross-sections.
        
        Args:
            show_radii: Whether to show cells as circles with intersection radii
            figsize: Figure size (width, height)
            title: Plot title
        """
        import matplotlib.pyplot as plt
        
        if not self.slice_cells:
            print("No cells in slice. Run slice_plane() first.")
            return
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Color map for cell types
        cell_types = list(set(c.cell_type for c in self.slice_cells))
        colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        # Plot cells as circles
        for slice_cell in self.slice_cells:
            if show_radii:
                radius = slice_cell.intersection_radius
            else:
                # Show original 3D radius for reference
                radius = slice_cell.radius
            
            color = color_map[slice_cell.cell_type]
            
            # Adjust alpha based on distance from plane
            # Cells farther from plane are more transparent
            max_dist = max(c.distance_from_plane for c in self.slice_cells) if self.slice_cells else 1
            alpha = 1.0 - (slice_cell.distance_from_plane / max_dist) * 0.5
            
            circle = plt.Circle(
                slice_cell.center_2d,
                radius,
                color=color,
                alpha=alpha,
                edgecolor='black',
                linewidth=0.5
            )
            ax.add_patch(circle)
        
        # Set equal aspect ratio
        ax.set_aspect('equal')
        
        # Set limits
        if self.slice_cells:
            x_coords = [c.center_2d[0] for c in self.slice_cells]
            y_coords = [c.center_2d[1] for c in self.slice_cells]
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
        
        if title:
            ax.set_title(title)
        else:
            ax.set_title(f'2D Tissue Slice: {len(self.slice_cells)} cells')
        
        plt.tight_layout()
        plt.show()
    
    def visualize_slice_in_3d(self, show_plane: bool = True,
                             plane_alpha: float = 0.3,
                             plane_size: float = None):
        """
        Visualize the slice plane within the 3D tissue context.
        
        Args:
            show_plane: Whether to show the slice plane
            plane_alpha: Transparency of the slice plane
            plane_size: Size of plane to display (default: tissue dimensions)
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        if self.slice_plane_point is None:
            print("No slice defined. Run slice_plane() first.")
            return
        
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot all cells in tissue
        cell_types = list(set(c.cell_type for c in self.tissue.cells))
        colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        # Plot cells (reduced resolution for performance)
        slice_cell_centers = {tuple(sc.center_3d) for sc in self.slice_cells}
        
        for cell in self.tissue.cells:
            u = np.linspace(0, 2 * np.pi, 12)
            v = np.linspace(0, np.pi, 12)
            x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
            y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
            z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
            
            color = color_map[cell.cell_type]
            
            # Highlight cells in the slice
            if tuple(cell.center) in slice_cell_centers:
                alpha = 0.8
                linewidth = 1
            else:
                alpha = 0.15
                linewidth = 0
            
            ax.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)),
                          alpha=alpha, linewidth=linewidth, antialiased=True, shade=False)
        
        # Draw slice plane
        if show_plane:
            if plane_size is None:
                plane_size = max(self.tissue.width, self.tissue.height, self.tissue.thickness)
            
            # Create plane mesh
            u_range = np.linspace(-plane_size/2, plane_size/2, 10)
            v_range = np.linspace(-plane_size/2, plane_size/2, 10)
            u_grid, v_grid = np.meshgrid(u_range, v_range)
            
            # Convert to 3D coordinates
            x_plane = (self.slice_plane_point[0] + 
                      u_grid * self.slice_basis_u[0] + 
                      v_grid * self.slice_basis_v[0])
            y_plane = (self.slice_plane_point[1] + 
                      u_grid * self.slice_basis_u[1] + 
                      v_grid * self.slice_basis_v[1])
            z_plane = (self.slice_plane_point[2] + 
                      u_grid * self.slice_basis_u[2] + 
                      v_grid * self.slice_basis_v[2])
            
            ax.plot_surface(x_plane, y_plane, z_plane, 
                          color='red', alpha=plane_alpha, 
                          edgecolor='darkred', linewidth=1)
        
        # Set labels and limits
        ax.set_xlabel('Width (μm)')
        ax.set_ylabel('Height (μm)')
        ax.set_zlabel('Thickness (μm)')
        ax.set_xlim(0, self.tissue.width)
        ax.set_ylim(0, self.tissue.height)
        ax.set_zlim(0, self.tissue.thickness)
        
        ax.set_title(f'3D Tissue with Slice Plane ({len(self.slice_cells)} cells captured)')
        plt.tight_layout()
        plt.show()
    
    def get_slice_statistics(self) -> Dict:
        """
        Calculate statistics about the slice.
        
        Returns:
            Dictionary with slice statistics
        """
        if not self.slice_cells:
            return {"num_cells": 0}
        
        stats = {
            "num_cells": len(self.slice_cells),
            "plane_point": self.slice_plane_point.tolist(),
            "plane_normal": self.slice_plane_normal.tolist(),
        }
        
        # Cell type distribution
        type_counts = {}
        type_radii = {}
        for slice_cell in self.slice_cells:
            if slice_cell.cell_type not in type_counts:
                type_counts[slice_cell.cell_type] = 0
                type_radii[slice_cell.cell_type] = []
            type_counts[slice_cell.cell_type] += 1
            type_radii[slice_cell.cell_type].append(slice_cell.intersection_radius)
        
        stats["cell_types"] = type_counts
        stats["avg_intersection_radii"] = {
            cell_type: np.mean(radii)
            for cell_type, radii in type_radii.items()
        }
        
        # Distance statistics
        distances = [sc.distance_from_plane for sc in self.slice_cells]
        stats["mean_distance_from_plane"] = np.mean(distances)
        stats["max_distance_from_plane"] = np.max(distances)
        
        return stats


def create_standard_slices(tissue: TissueSection, 
                          num_slices: int = 5) -> List[TissueSlicer]:
    """
    Create multiple parallel slices through the tissue.
    
    Args:
        tissue: TissueSection to slice
        num_slices: Number of evenly-spaced slices
        
    Returns:
        List of TissueSlicer objects with computed slices
    """
    slicers = []
    
    for i in range(num_slices):
        z_position = tissue.thickness * (i + 1) / (num_slices + 1)
        slicer = TissueSlicer(tissue)
        slicer.slice_plane(z_position=z_position)
        slicers.append(slicer)
    
    return slicers
