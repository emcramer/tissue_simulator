"""
Unit tests for tissue simulator package.
"""

import unittest
import numpy as np
from tissue_simulator import TissueSection, Cell, SpherePacker


class TestCell(unittest.TestCase):
    """Test Cell class functionality."""
    
    def test_cell_creation(self):
        """Test basic cell creation."""
        cell = Cell(center=(10, 20, 30), radius=5, cell_type="test")
        
        self.assertEqual(cell.radius, 5)
        self.assertEqual(cell.cell_type, "test")
        np.testing.assert_array_equal(cell.center, np.array([10, 20, 30]))
    
    def test_cell_intersection(self):
        """Test cell intersection detection."""
        cell1 = Cell(center=(0, 0, 0), radius=5)
        cell2 = Cell(center=(8, 0, 0), radius=5)  # Touching
        cell3 = Cell(center=(20, 0, 0), radius=5)  # Far away
        
        self.assertTrue(cell1.intersects(cell2))
        self.assertFalse(cell1.intersects(cell3))
    
    def test_within_bounds(self):
        """Test bounds checking."""
        bounds = (100, 100, 100)
        
        # Cell completely within bounds
        cell1 = Cell(center=(50, 50, 50), radius=10)
        self.assertTrue(cell1.is_within_bounds(bounds))
        
        # Cell partially outside
        cell2 = Cell(center=(95, 50, 50), radius=10)
        self.assertFalse(cell2.is_within_bounds(bounds))
    
    def test_intersects_bounds(self):
        """Test if cell center is within bounds."""
        bounds = (100, 100, 100)
        
        # Center inside
        cell1 = Cell(center=(50, 50, 50), radius=10)
        self.assertTrue(cell1.intersects_bounds(bounds))
        
        # Center outside
        cell2 = Cell(center=(150, 50, 50), radius=10)
        self.assertFalse(cell2.intersects_bounds(bounds))


class TestTissueSection(unittest.TestCase):
    """Test TissueSection class functionality."""
    
    def test_tissue_creation_simple(self):
        """Test tissue creation with simple radii."""
        tissue = TissueSection(
            height=100,
            width=100,
            thickness=50,
            cell_radii=(5, 10)
        )
        
        self.assertEqual(tissue.height, 100)
        self.assertEqual(tissue.width, 100)
        self.assertEqual(tissue.thickness, 50)
        self.assertIn("default", tissue.cell_radii)
    
    def test_tissue_creation_dict(self):
        """Test tissue creation with cell type dictionary."""
        radii_dict = {
            'type_a': (5, 10),
            'type_b': (8, 15)
        }
        tissue = TissueSection(
            height=100,
            width=100,
            thickness=50,
            cell_radii=radii_dict
        )
        
        self.assertEqual(len(tissue.cell_radii), 2)
        self.assertIn('type_a', tissue.cell_radii)
        self.assertIn('type_b', tissue.cell_radii)
    
    def test_get_bounds(self):
        """Test bounds getter."""
        tissue = TissueSection(100, 200, 50, (5, 10))
        bounds = tissue.get_bounds()
        
        self.assertEqual(bounds, (100, 200, 50))
    
    def test_cell_generation(self):
        """Test cell generation."""
        tissue = TissueSection(100, 100, 50, (5, 10))
        num_cells = tissue.generate_cells(max_attempts=100)
        
        self.assertGreater(num_cells, 0)
        self.assertEqual(len(tissue.cells), num_cells)
    
    def test_statistics(self):
        """Test statistics calculation."""
        tissue = TissueSection(100, 100, 50, {'type_a': (5, 10)})
        tissue.generate_cells(max_attempts=100)
        
        stats = tissue.get_cell_statistics()
        
        self.assertIn('total_cells', stats)
        self.assertIn('packing_fraction', stats)
        self.assertEqual(stats['total_cells'], len(tissue.cells))
        self.assertGreater(stats['packing_fraction'], 0)
        self.assertLess(stats['packing_fraction'], 1)
    
    def test_clear_cells(self):
        """Test clearing cells."""
        tissue = TissueSection(100, 100, 50, (5, 10))
        tissue.generate_cells(max_attempts=100)
        
        self.assertGreater(len(tissue.cells), 0)
        
        tissue.clear_cells()
        self.assertEqual(len(tissue.cells), 0)


class TestSpherePacker(unittest.TestCase):
    """Test SpherePacker class functionality."""
    
    def test_packer_creation(self):
        """Test packer initialization."""
        packer = SpherePacker(
            bounds=(100, 100, 50),
            cell_radii_config={'type_a': (5, 10)},
            min_spacing=0.5,
            allow_boundary_cells=True
        )
        
        self.assertEqual(packer.bounds, (100, 100, 50))
        self.assertEqual(packer.min_spacing, 0.5)
        self.assertTrue(packer.allow_boundary_cells)
    
    def test_packing(self):
        """Test basic packing."""
        packer = SpherePacker(
            bounds=(100, 100, 50),
            cell_radii_config={'type_a': (5, 10)},
            min_spacing=1.0,
            allow_boundary_cells=True
        )
        
        cells = packer.pack(max_attempts=200)
        
        self.assertGreater(len(cells), 0)
        self.assertIsInstance(cells[0], Cell)
    
    def test_no_collisions(self):
        """Test that packed cells don't collide."""
        packer = SpherePacker(
            bounds=(100, 100, 50),
            cell_radii_config={'type_a': (5, 8)},
            min_spacing=1.0,
            allow_boundary_cells=False
        )
        
        cells = packer.pack(max_attempts=200)
        
        # Check all pairs for collisions
        for i, cell1 in enumerate(cells):
            for cell2 in cells[i+1:]:
                distance = np.linalg.norm(cell1.center - cell2.center)
                min_distance = cell1.radius + cell2.radius + packer.min_spacing
                self.assertGreaterEqual(distance, min_distance - 1e-6)
    
    def test_boundary_enforcement(self):
        """Test boundary cell restrictions."""
        # With boundary cells allowed
        packer1 = SpherePacker(
            bounds=(50, 50, 50),
            cell_radii_config={'type_a': (10, 15)},
            min_spacing=0.5,
            allow_boundary_cells=True
        )
        cells1 = packer1.pack(max_attempts=100)
        
        # Without boundary cells allowed
        packer2 = SpherePacker(
            bounds=(50, 50, 50),
            cell_radii_config={'type_a': (10, 15)},
            min_spacing=0.5,
            allow_boundary_cells=False
        )
        cells2 = packer2.pack(max_attempts=100)
        
        # All cells in packer2 should be completely within bounds
        for cell in cells2:
            self.assertTrue(cell.is_within_bounds(packer2.bounds))


class TestIntegration(unittest.TestCase):
    """Integration tests for complete workflows."""
    
    def test_full_workflow(self):
        """Test complete tissue generation workflow."""
        # Create tissue
        tissue = TissueSection(
            height=200,
            width=200,
            thickness=100,
            cell_radii={
                'type_a': (5, 10),
                'type_b': (7, 12)
            }
        )
        
        # Generate cells
        num_cells = tissue.generate_cells(
            max_attempts=500,
            min_spacing=0.5,
            allow_boundary_cells=True
        )
        
        self.assertGreater(num_cells, 0)
        
        # Get statistics
        stats = tissue.get_cell_statistics()
        self.assertEqual(stats['total_cells'], num_cells)
        
        # Check cell types are present
        cell_types = [c.cell_type for c in tissue.cells]
        self.assertTrue(any(ct == 'type_a' for ct in cell_types) or
                       any(ct == 'type_b' for ct in cell_types))
    
    def test_export_import(self):
        """Test CSV export functionality."""
        import tempfile
        import os
        
        tissue = TissueSection(100, 100, 50, (5, 10))
        tissue.generate_cells(max_attempts=100)
        
        # Export to temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            temp_path = f.name
        
        try:
            tissue.export_to_csv(temp_path)
            
            # Check file exists and has content
            self.assertTrue(os.path.exists(temp_path))
            
            with open(temp_path, 'r') as f:
                lines = f.readlines()
                self.assertGreater(len(lines), 1)  # Header + data
                self.assertIn('x', lines[0])  # Check header
        finally:
            os.unlink(temp_path)


if __name__ == '__main__':
    unittest.main()
