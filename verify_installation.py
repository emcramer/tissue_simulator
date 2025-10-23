#!/usr/bin/env python3
"""
Installation verification and test script for tissue_simulator package.
Run this after installation to verify everything works correctly.
"""

import sys
import os

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    try:
        import numpy
        print("  ✓ numpy imported successfully")
    except ImportError as e:
        print(f"  ✗ Failed to import numpy: {e}")
        return False
    
    try:
        import matplotlib
        print("  ✓ matplotlib imported successfully")
    except ImportError as e:
        print(f"  ✗ Failed to import matplotlib: {e}")
        return False
    
    try:
        import PyQt5
        print("  ✓ PyQt5 imported successfully")
    except ImportError as e:
        print(f"  ✗ Failed to import PyQt5: {e}")
        print("    Note: PyQt5 is only needed for GUI functionality")
    
    try:
        from tissue_simulator import TissueSection, Cell, SpherePacker
        print("  ✓ tissue_simulator package imported successfully")
    except ImportError as e:
        print(f"  ✗ Failed to import tissue_simulator: {e}")
        return False
    
    return True

def test_basic_functionality():
    """Test basic tissue generation."""
    print("\nTesting basic functionality...")
    
    try:
        from tissue_simulator import TissueSection
        
        # Create tissue
        tissue = TissueSection(
            height=100,
            width=100,
            thickness=50,
            cell_radii=(5, 10)
        )
        print("  ✓ TissueSection created successfully")
        
        # Generate cells
        num_cells = tissue.generate_cells(max_attempts=200)
        print(f"  ✓ Generated {num_cells} cells")
        
        if num_cells == 0:
            print("  ⚠ Warning: No cells generated (this is unusual)")
        
        # Get statistics
        stats = tissue.get_cell_statistics()
        print(f"  ✓ Statistics calculated: {stats['total_cells']} cells, "
              f"packing fraction = {stats['packing_fraction']:.3f}")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error during basic functionality test: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_multi_cell_type():
    """Test multiple cell types."""
    print("\nTesting multiple cell types...")
    
    try:
        from tissue_simulator import TissueSection
        
        tissue = TissueSection(
            height=150,
            width=150,
            thickness=75,
            cell_radii={
                'type_a': (5, 10),
                'type_b': (7, 12)
            }
        )
        print("  ✓ Multi-cell-type tissue created")
        
        num_cells = tissue.generate_cells(max_attempts=300)
        print(f"  ✓ Generated {num_cells} cells")
        
        stats = tissue.get_cell_statistics()
        print(f"  ✓ Cell types: {list(stats['cell_types'].keys())}")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error during multi-cell-type test: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_export():
    """Test CSV export functionality."""
    print("\nTesting CSV export...")
    
    try:
        from tissue_simulator import TissueSection
        import tempfile
        
        tissue = TissueSection(100, 100, 50, (5, 10))
        tissue.generate_cells(max_attempts=200)
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            temp_path = f.name
        
        try:
            tissue.export_to_csv(temp_path)
            print(f"  ✓ Exported to {temp_path}")
            
            # Verify file exists and has content
            if os.path.exists(temp_path):
                with open(temp_path, 'r') as f:
                    lines = f.readlines()
                    print(f"  ✓ CSV file contains {len(lines)} lines")
            else:
                print("  ✗ Export file not found")
                return False
                
        finally:
            # Clean up
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error during export test: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_packing_algorithm():
    """Test SpherePacker directly."""
    print("\nTesting SpherePacker...")
    
    try:
        from tissue_simulator import SpherePacker
        
        packer = SpherePacker(
            bounds=(100, 100, 50),
            cell_radii_config={'type_a': (5, 10)},
            min_spacing=0.5,
            allow_boundary_cells=True
        )
        print("  ✓ SpherePacker created")
        
        cells = packer.pack(max_attempts=200)
        print(f"  ✓ Packed {len(cells)} cells")
        
        # Verify no collisions
        for i, cell1 in enumerate(cells):
            for cell2 in cells[i+1:]:
                import numpy as np
                distance = np.linalg.norm(cell1.center - cell2.center)
                min_distance = cell1.radius + cell2.radius + packer.min_spacing
                if distance < min_distance - 1e-6:
                    print(f"  ✗ Collision detected between cells!")
                    return False
        
        print("  ✓ No collisions detected")
        return True
        
    except Exception as e:
        print(f"  ✗ Error during packing test: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("Tissue Simulator Installation Verification")
    print("=" * 60)
    
    results = []
    
    # Run tests
    results.append(("Imports", test_imports()))
    results.append(("Basic Functionality", test_basic_functionality()))
    results.append(("Multiple Cell Types", test_multi_cell_type()))
    results.append(("CSV Export", test_export()))
    results.append(("Packing Algorithm", test_packing_algorithm()))
    
    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{test_name:.<40} {status}")
    
    print("=" * 60)
    print(f"\nResults: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ All tests passed! Installation is working correctly.")
        print("\nYou can now:")
        print("  1. Run the GUI: python -m tissue_simulator.gui")
        print("  2. Try examples: python examples/simple_example.py")
        print("  3. Read the guide: cat GUIDE.md")
        return 0
    else:
        print("\n✗ Some tests failed. Please check the errors above.")
        print("\nTroubleshooting:")
        print("  1. Ensure all dependencies are installed: pip install -r requirements.txt")
        print("  2. Reinstall the package: pip install -e .")
        print("  3. Check Python version: python --version (need 3.8+)")
        return 1

if __name__ == '__main__':
    sys.exit(main())
