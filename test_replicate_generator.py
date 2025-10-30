#!/usr/bin/env python3
"""
Quick test to verify replicate generator imports and basic functionality.
"""

import sys
from pathlib import Path

# Add package to path
sys.path.insert(0, str(Path(__file__).parent))

def test_imports():
    """Test that all new modules import correctly."""
    print("Testing imports...")
    
    try:
        from tissue_simulator import (
            ReplicateGenerator,
            TargetStatistics,
            ReplicateStatistics,
            load_target_statistics_from_csv,
            load_target_statistics_from_tissue,
            TissueSection,
            SpatialNetworkAnalyzer,
            InteractionStatistics
        )
        print("✓ All imports successful")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False


def test_basic_functionality():
    """Test basic replicate generator functionality."""
    print("\nTesting basic functionality...")
    
    try:
        from tissue_simulator import (
            TissueSection,
            load_target_statistics_from_tissue,
            ReplicateGenerator
        )
        
        # Create small reference tissue
        print("  Creating reference tissue...")
        reference = TissueSection(
            height=200,
            width=200,
            thickness=50,
            cell_radii={'type_a': (5, 8), 'type_b': (6, 10)}
        )
        
        num_cells = reference.generate_cells(max_attempts=500)
        print(f"  Generated {num_cells} cells")
        
        # Extract statistics
        print("  Extracting spatial statistics...")
        target_stats = load_target_statistics_from_tissue(
            reference,
            network_mode="contact"
        )
        print(f"  Found {len(target_stats.interaction_stats)} interaction types")
        
        # Setup generator
        print("  Setting up replicate generator...")
        generator = ReplicateGenerator(
            target_stats=target_stats,
            tissue_dimensions=(200, 200, 50),
            base_cell_radii={'type_a': (5, 8), 'type_b': (6, 10)},
            network_mode="contact",
            seed=42
        )
        
        # Generate one replicate
        print("  Generating test replicate...")
        tissue, stats = generator.generate_single_replicate(
            replicate_id=0,
            max_attempts=500,
            max_iterations=2,
            tolerance=0.20
        )
        
        print(f"  Replicate cells: {stats.num_cells}")
        print(f"  Divergence: {stats.divergence_score:.4f}")
        
        print("✓ Basic functionality test passed")
        return True
        
    except Exception as e:
        print(f"✗ Functionality test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mcp_imports():
    """Test that MCP server includes new tools."""
    print("\nTesting MCP integration...")
    
    try:
        from tissue_simulator.mcp.server import TissueSimulatorMCPServer
        
        # Check that new attributes exist
        server = TissueSimulatorMCPServer()
        
        # Verify new attributes
        assert hasattr(server, 'replicate_generator'), "Missing replicate_generator attribute"
        assert hasattr(server, 'generated_replicates'), "Missing generated_replicates attribute"
        
        print("✓ MCP server integration verified")
        return True
        
    except Exception as e:
        print(f"✗ MCP integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("=" * 70)
    print("Replicate Generator Verification")
    print("=" * 70)
    
    results = []
    
    # Run tests
    results.append(("Imports", test_imports()))
    results.append(("Basic Functionality", test_basic_functionality()))
    results.append(("MCP Integration", test_mcp_imports()))
    
    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    for test_name, passed in results:
        status = "PASSED" if passed else "FAILED"
        symbol = "✓" if passed else "✗"
        print(f"{symbol} {test_name}: {status}")
    
    all_passed = all(result[1] for result in results)
    
    print("=" * 70)
    if all_passed:
        print("All tests passed! ✓")
        print("\nReplicate generation functionality is ready to use.")
        print("See examples/replicate_generation_example.py for usage examples.")
    else:
        print("Some tests failed. ✗")
        print("Please check the errors above.")
    print("=" * 70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
