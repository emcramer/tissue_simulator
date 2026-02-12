"""
Simple test client for the Tissue Simulator MCP server.

This script demonstrates how to interact with the MCP server programmatically.
"""

import asyncio
import json
from typing import Dict, Any


async def test_mcp_server():
    """Test the MCP server with a complete workflow."""
    
    print("=" * 60)
    print("Tissue Simulator MCP Server Test")
    print("=" * 60)
    print()
    
    # Note: This is a demonstration of the tool calls.
    # In practice, these would be called through the MCP protocol.
    
    test_workflow = [
        {
            "step": 1,
            "tool": "create_tissue",
            "description": "Create a 3D tissue",
            "args": {
                "height": 400,
                "width": 400,
                "thickness": 100,
                "cell_types": {
                    "epithelial": [6, 10],
                    "stromal": [8, 15],
                    "immune": [3, 6]
                }
            }
        },
        {
            "step": 2,
            "tool": "generate_cells",
            "description": "Populate with cells",
            "args": {
                "max_attempts": 1500,
                "min_spacing": 0.5,
                "allow_boundary_cells": True
            }
        },
        {
            "step": 3,
            "tool": "get_tissue_statistics",
            "description": "Get tissue statistics",
            "args": {}
        },
        {
            "step": 4,
            "tool": "create_slice",
            "description": "Create horizontal slice",
            "args": {
                "z_position": 50
            }
        },
        {
            "step": 5,
            "tool": "get_slice_statistics",
            "description": "Get slice statistics",
            "args": {}
        },
        {
            "step": 6,
            "tool": "export_tissue_csv",
            "description": "Export tissue data",
            "args": {
                "filename": "test_tissue.csv"
            }
        },
        {
            "step": 7,
            "tool": "export_slice_csv",
            "description": "Export slice data",
            "args": {
                "filename": "test_slice.csv",
                "include_3d": True
            }
        }
    ]
    
    print("Test Workflow:")
    print()
    
    for test in test_workflow:
        print(f"Step {test['step']}: {test['description']}")
        print(f"  Tool: {test['tool']}")
        print(f"  Args: {json.dumps(test['args'], indent=4)}")
        print()
    
    print("=" * 60)
    print("Expected Results:")
    print("=" * 60)
    print()
    
    print("Step 1 - create_tissue:")
    print(json.dumps({
        "status": "success",
        "message": "Created tissue: 400x400x100 μm",
        "cell_types": ["epithelial", "stromal", "immune"]
    }, indent=2))
    print()
    
    print("Step 2 - generate_cells:")
    print(json.dumps({
        "status": "success",
        "num_cells_generated": 234,
        "interior_cells": 198,
        "boundary_cells": 36,
        "packing_fraction": 0.312
    }, indent=2))
    print()
    
    print("Step 3 - get_tissue_statistics:")
    print(json.dumps({
        "total_cells": 234,
        "packing_fraction": 0.312,
        "cell_types": {"epithelial": 152, "stromal": 68, "immune": 14}
    }, indent=2))
    print()
    
    print("Step 4 - create_slice:")
    print(json.dumps({
        "status": "success",
        "num_cells_in_slice": 89,
        "plane_point": [200, 200, 50],
        "plane_normal": [0, 0, 1]
    }, indent=2))
    print()
    
    print("Step 5 - get_slice_statistics:")
    print(json.dumps({
        "num_cells": 89,
        "cell_types": {"epithelial": 58, "stromal": 26, "immune": 5}
    }, indent=2))
    print()
    
    print("Step 6 & 7 - export:")
    print(json.dumps({
        "status": "success",
        "filepath": "/tmp/tissue_sim_xyz/test_tissue.csv"
    }, indent=2))
    print()
    
    print("=" * 60)
    print("To actually run the MCP server:")
    print("  python run_mcp_server.py")
    print()
    print("Then configure it in Claude Desktop to use these tools!")
    print("=" * 60)


async def demonstrate_tool_schemas():
    """Show the tool schemas that would be available."""
    
    print("\n" + "=" * 60)
    print("Available MCP Tools")
    print("=" * 60)
    print()
    
    tools = [
        "create_tissue - Create 3D tissue structure",
        "generate_cells - Populate tissue with cells",
        "get_tissue_statistics - Get tissue metrics",
        "create_slice - Create 2D slice",
        "get_slice_statistics - Get slice metrics",
        "create_serial_slices - Create multiple slices",
        "export_tissue_csv - Export 3D data",
        "export_slice_csv - Export 2D data",
        "visualize_tissue - Create 3D visualization",
        "visualize_slice_2d - Create 2D visualization",
        "reset_tissue - Clear and start fresh"
    ]
    
    for i, tool in enumerate(tools, 1):
        print(f"{i:2d}. {tool}")
    
    print()
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(test_mcp_server())
    asyncio.run(demonstrate_tool_schemas())
