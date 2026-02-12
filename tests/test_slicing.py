"""
Unit tests for tissue slicing functionality.
"""

import unittest
import numpy as np
import tempfile
import os

from tissue_simulator import TissueSection, TissueSlicer, SliceCell, create_standard_slices


class TestTissueSlicer(unittest.TestCase):
    """Test TissueSlicer class functionality."""
    
    def setUp(self):
        """Set up test tissue."""
        self.tissue = TissueSection(
            height=100,
            width=100,
            thickness=50,
            cell_radii={'test': (5, 10)}
        )
        self.tissue.generate_cells(max_attempts=200)
    
    def test_slicer_creation(self):
        """Test slicer initialization."""
        slicer = TissueSlicer(self.tissue)
        
        self.assertEqual(slicer.tissue, self.tissue)
        self.assertEqual(len(slicer.slice_cells), 0)
        self.assertIsNone(slicer.slice_plane_point)
    
    def test_horizontal_slice(self):
        """Test creating a horizontal slice."""
        slicer = TissueSlicer(self.tissue)
        slice_cells = slicer.slice_plane(z_position=25)
        
        self.assertGreater(len(slice_cells), 0)
        self.assertEqual(len(slicer.slice_cells), len(slice_cells))
        
        # Check plane is horizontal
        np.testing.assert_array_almost_equal(
            slicer.slice_plane_normal,
            np.array([0, 0, 1])
        )
    
    def test_angled_slice(self):
        """Test creating an angled slice."""
        slicer = TissueSlicer(self.tissue)
        slice_cells = slicer.slice_plane(
            point=(50, 50, 25),
            angle_x=45,
            angle_y=0
        )
        
        self.assertIsNotNone(slicer.slice_plane_normal)
        self.assertGreater(len(slice_cells), 0)
        
        # Normal should not be vertical
        self.assertNotAlmostEqual(slicer.slice_plane_normal[2], 1.0)
    
    def test_custom_normal(self):
        """Test slice with custom normal vector."""
        slicer = TissueSlicer(self.tissue)
        custom_normal = np.array([1, 1, 1])
        
        slice_cells = slicer.slice_plane(
            point=(50, 50, 25),
            normal=custom_normal
        )
        
        # Normal should be normalized
        self.assertAlmostEqual(np.linalg.norm(slicer.slice_plane_normal), 1.0)
        
        # Direction should be correct
        expected_normal = custom_normal / np.linalg.norm(custom_normal)
        np.testing.assert_array_almost_equal(
            slicer.slice_plane_normal,
            expected_normal
        )
    
    def test_slice_cell_properties(self):
        """Test properties of sliced cells."""
        slicer = TissueSlicer(self.tissue)
        slice_cells = slicer.slice_plane(z_position=25)
        
        if len(slice_cells) > 0:
            slice_cell = slice_cells[0]
            
            # Check types
            self.assertIsInstance(slice_cell, SliceCell)
            self.assertIsInstance(slice_cell.center_3d, np.ndarray)
            self.assertIsInstance(slice_cell.center_2d, np.ndarray)
            self.assertIsInstance(slice_cell.radius, float)
            
            # Check dimensions
            self.assertEqual(len(slice_cell.center_3d), 3)
            self.assertEqual(len(slice_cell.center_2d), 2)
            
            # Check intersection radius is less than or equal to cell radius
            self.assertLessEqual(slice_cell.intersection_radius, slice_cell.radius)
            
            # Check distance is non-negative
            self.assertGreaterEqual(slice_cell.distance_from_plane, 0)
    
    def test_plane_basis_orthogonal(self):
        """Test that plane basis vectors are orthonormal."""
        slicer = TissueSlicer(self.tissue)
        slicer.slice_plane(z_position=25)
        
        # Check normalization
        self.assertAlmostEqual(np.linalg.norm(slicer.slice_basis_u), 1.0)
        self.assertAlmostEqual(np.linalg.norm(slicer.slice_basis_v), 1.0)
        self.assertAlmostEqual(np.linalg.norm(slicer.slice_plane_normal), 1.0)
        
        # Check orthogonality
        self.assertAlmostEqual(np.dot(slicer.slice_basis_u, slicer.slice_basis_v), 0.0)
        self.assertAlmostEqual(np.dot(slicer.slice_basis_u, slicer.slice_plane_normal), 0.0)
        self.assertAlmostEqual(np.dot(slicer.slice_basis_v, slicer.slice_plane_normal), 0.0)
    
    def test_slice_statistics(self):
        """Test slice statistics calculation."""
        slicer = TissueSlicer(self.tissue)
        slicer.slice_plane(z_position=25)
        
        stats = slicer.get_slice_statistics()
        
        self.assertIn('num_cells', stats)
        self.assertIn('plane_point', stats)
        self.assertIn('plane_normal', stats)
        self.assertIn('cell_types', stats)
        
        self.assertEqual(stats['num_cells'], len(slicer.slice_cells))
    
    def test_export_slice_csv(self):
        """Test CSV export of slice data."""
        slicer = TissueSlicer(self.tissue)
        slicer.slice_plane(z_position=25)
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            temp_path = f.name
        
        try:
            # Export with 3D coordinates
            slicer.export_slice_csv(temp_path, include_3d=True)
            
            # Check file exists
            self.assertTrue(os.path.exists(temp_path))
            
            # Read and verify content
            with open(temp_path, 'r') as f:
                lines = f.readlines()
                self.assertGreater(len(lines), 1)  # Header + data
                
                # Check header
                header = lines[0].strip()
                self.assertIn('x_2d', header)
                self.assertIn('y_2d', header)
                self.assertIn('x_3d', header)
                self.assertIn('cell_type', header)
            
            # Export without 3D coordinates
            slicer.export_slice_csv(temp_path, include_3d=False)
            
            with open(temp_path, 'r') as f:
                lines = f.readlines()
                header = lines[0].strip()
                self.assertIn('x_2d', header)
                self.assertNotIn('x_3d', header)
                
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
    
    def test_empty_slice(self):
        """Test slicing when no cells are captured."""
        slicer = TissueSlicer(self.tissue)
        
        # Slice far outside tissue
        slice_cells = slicer.slice_plane(z_position=500)
        
        self.assertEqual(len(slice_cells), 0)
        
        stats = slicer.get_slice_statistics()
        self.assertEqual(stats['num_cells'], 0)
    
    def test_all_cells_captured(self):
        """Test that slice through tissue center captures many cells."""
        # Create tissue with known cells
        tissue = TissueSection(100, 100, 50, (5, 8))
        tissue.generate_cells(max_attempts=500)
        
        slicer = TissueSlicer(tissue)
        slice_cells = slicer.slice_plane(z_position=25)  # Middle
        
        # Should capture a reasonable number of cells
        self.assertGreater(len(slice_cells), 0)
        self.assertLessEqual(len(slice_cells), len(tissue.cells))


class TestStandardSlices(unittest.TestCase):
    """Test create_standard_slices functionality."""
    
    def test_create_multiple_slices(self):
        """Test creating multiple parallel slices."""
        tissue = TissueSection(100, 100, 50, (5, 10))
        tissue.generate_cells(max_attempts=200)
        
        num_slices = 5
        slicers = create_standard_slices(tissue, num_slices=num_slices)
        
        self.assertEqual(len(slicers), num_slices)
        
        # Check all are TissueSlicer instances
        for slicer in slicers:
            self.assertIsInstance(slicer, TissueSlicer)
            self.assertGreater(len(slicer.slice_cells), 0)
    
    def test_slice_spacing(self):
        """Test that slices are evenly spaced."""
        tissue = TissueSection(100, 100, 60, (5, 10))
        tissue.generate_cells(max_attempts=200)
        
        num_slices = 4
        slicers = create_standard_slices(tissue, num_slices=num_slices)
        
        # Get z-positions
        z_positions = [s.slice_plane_point[2] for s in slicers]
        
        # Check they're in order
        self.assertEqual(z_positions, sorted(z_positions))
        
        # Check spacing is roughly equal
        spacings = [z_positions[i+1] - z_positions[i] 
                   for i in range(len(z_positions)-1)]
        
        mean_spacing = np.mean(spacings)
        for spacing in spacings:
            self.assertAlmostEqual(spacing, mean_spacing, places=1)
    
    def test_slices_within_tissue(self):
        """Test that all slice planes are within tissue bounds."""
        tissue = TissueSection(100, 100, 80, (5, 10))
        tissue.generate_cells(max_attempts=200)
        
        slicers = create_standard_slices(tissue, num_slices=10)
        
        for slicer in slicers:
            z_pos = slicer.slice_plane_point[2]
            self.assertGreater(z_pos, 0)
            self.assertLess(z_pos, tissue.thickness)


class TestSlicingGeometry(unittest.TestCase):
    """Test geometric calculations in slicing."""
    
    def test_intersection_radius_calculation(self):
        """Test intersection radius for known geometry."""
        tissue = TissueSection(100, 100, 100, (10, 10))
        
        # Add a cell at known position
        from tissue_simulator.tissue import Cell
        cell = Cell(center=(50, 50, 50), radius=10, cell_type='test')
        tissue.cells = [cell]
        
        slicer = TissueSlicer(tissue)
        
        # Slice through center of cell
        slice_cells = slicer.slice_plane(z_position=50)
        
        self.assertEqual(len(slice_cells), 1)
        
        # Intersection radius should equal cell radius (plane through center)
        self.assertAlmostEqual(
            slice_cells[0].intersection_radius,
            10.0,
            places=5
        )
        
        # Distance from plane should be ~0
        self.assertAlmostEqual(
            slice_cells[0].distance_from_plane,
            0.0,
            places=5
        )
    
    def test_off_center_intersection(self):
        """Test intersection radius when plane doesn't pass through center."""
        tissue = TissueSection(100, 100, 100, (10, 10))
        
        from tissue_simulator.tissue import Cell
        cell = Cell(center=(50, 50, 50), radius=10, cell_type='test')
        tissue.cells = [cell]
        
        slicer = TissueSlicer(tissue)
        
        # Slice 5 units above center
        slice_cells = slicer.slice_plane(z_position=55)
        
        self.assertEqual(len(slice_cells), 1)
        
        # Distance should be 5
        self.assertAlmostEqual(
            slice_cells[0].distance_from_plane,
            5.0,
            places=5
        )
        
        # Intersection radius should be sqrt(10^2 - 5^2) = sqrt(75)
        expected_radius = np.sqrt(10**2 - 5**2)
        self.assertAlmostEqual(
            slice_cells[0].intersection_radius,
            expected_radius,
            places=5
        )
    
    def test_grazing_intersection(self):
        """Test intersection at edge of sphere."""
        tissue = TissueSection(100, 100, 100, (10, 10))
        
        from tissue_simulator.tissue import Cell
        cell = Cell(center=(50, 50, 50), radius=10, cell_type='test')
        tissue.cells = [cell]
        
        slicer = TissueSlicer(tissue)
        
        # Slice at edge of sphere
        slice_cells = slicer.slice_plane(z_position=59.9)
        
        self.assertEqual(len(slice_cells), 1)
        
        # Intersection radius should be small
        self.assertLess(slice_cells[0].intersection_radius, 2.0)
        
        # Distance should be close to radius
        self.assertAlmostEqual(
            slice_cells[0].distance_from_plane,
            9.9,
            places=1
        )
    
    def test_no_intersection(self):
        """Test that cell outside slice plane is not captured."""
        tissue = TissueSection(100, 100, 100, (10, 10))
        
        from tissue_simulator.tissue import Cell
        cell = Cell(center=(50, 50, 50), radius=10, cell_type='test')
        tissue.cells = [cell]
        
        slicer = TissueSlicer(tissue)
        
        # Slice beyond sphere
        slice_cells = slicer.slice_plane(z_position=70)
        
        self.assertEqual(len(slice_cells), 0)


if __name__ == '__main__':
    unittest.main()
