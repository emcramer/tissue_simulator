# Model Context Protocol (MCP) Integration

## Overview

The Tissue Simulator now includes a Model Context Protocol (MCP) server that exposes tissue generation and slicing functionality to Large Language Models (LLMs). This allows LLMs like Claude to generate, analyze, and visualize tissue simulations through structured tool calls.

## What is MCP?

Model Context Protocol (MCP) is a standardized protocol that enables LLMs to interact with external tools and data sources. It provides:

- **Structured tool definitions**: Clear schemas for what each tool does
- **Type-safe interfaces**: Input/output validation
- **Stateful sessions**: Maintain context across multiple tool calls
- **Standard communication**: Works with multiple LLM providers

## Features

The MCP server exposes 11 tools for tissue simulation:

### Tissue Generation
- `create_tissue`: Define tissue dimensions and cell types
- `generate_cells`: Populate tissue with cells using sphere packing
- `get_tissue_statistics`: Get comprehensive tissue statistics
- `reset_tissue`: Clear and start fresh

### 2D Slicing
- `create_slice`: Create a 2D slice at any angle
- `get_slice_statistics`: Get slice statistics
- `create_serial_slices`: Create multiple parallel slices

### Data Export
- `export_tissue_csv`: Export 3D tissue data
- `export_slice_csv`: Export 2D slice data

### Visualization
- `visualize_tissue`: Generate 3D tissue visualization
- `visualize_slice_2d`: Generate 2D slice visualization

## Installation

### 1. Install MCP Library

```bash
pip install mcp
```

### 2. Update Package

```bash
cd /Users/cramere/tissue_simulator
pip install -e .
```

### 3. Test the Server

```bash
python run_mcp_server.py
```

## Configuration

### For Claude Desktop

1. Locate your Claude Desktop config file:
   - macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - Windows: `%APPDATA%\Claude\claude_desktop_config.json`
   - Linux: `~/.config/Claude/claude_desktop_config.json`

2. Add the tissue simulator server:

```json
{
  "mcpServers": {
    "tissue-simulator": {
      "command": "python",
      "args": [
        "/Users/cramere/tissue_simulator/run_mcp_server.py"
      ],
      "description": "3D tissue generation and 2D slicing simulator"
    }
  }
}
```

3. Restart Claude Desktop

### For Other MCP Clients

Use the configuration template in `mcp_config_claude_desktop.json` and adapt for your client.

## Usage Examples

### Example 1: Basic Tissue Generation

```
LLM: I'll create a tissue simulation for you.

[Calls create_tissue]
{
  "height": 400,
  "width": 400,
  "thickness": 100,
  "cell_types": {
    "epithelial": [6, 10],
    "stromal": [8, 15]
  }
}

[Calls generate_cells]
{
  "max_attempts": 1500,
  "min_spacing": 0.5
}

[Calls get_tissue_statistics]

The tissue has been generated with 234 cells:
- 152 epithelial cells (avg radius: 7.8 μm)
- 82 stromal cells (avg radius: 11.2 μm)
- Packing fraction: 0.312
```

### Example 2: Create and Analyze Slice

```
LLM: I'll create a horizontal slice through the middle.

[Calls create_slice]
{
  "z_position": 50
}

[Calls get_slice_statistics]

The slice captured 89 cells at z=50 μm:
- 58 epithelial cells
- 31 stromal cells
- Mean distance from plane: 3.2 μm
```

### Example 3: Serial Sections

```
LLM: I'll create 5 serial sections for histology analysis.

[Calls create_serial_slices]
{
  "num_slices": 5
}

Created 5 slices:
- Slice 1 (z=16.7): 45 cells
- Slice 2 (z=33.3): 68 cells
- Slice 3 (z=50.0): 89 cells
- Slice 4 (z=66.7): 72 cells
- Slice 5 (z=83.3): 51 cells
```

### Example 4: Export and Visualize

```
LLM: I'll export the data and create visualizations.

[Calls export_tissue_csv]
{
  "filename": "my_tissue.csv"
}

[Calls visualize_tissue]
{
  "elevation": 30,
  "azimuth": 60,
  "filename": "tissue_view.png"
}

Data exported to: /tmp/tissue_sim_xyz/my_tissue.csv
Visualization saved to: /tmp/tissue_sim_xyz/tissue_view.png
```

## Tool Reference

### create_tissue

**Description**: Create a 3D tissue section with specified dimensions and cell types.

**Parameters**:
- `height` (number, required): Height in micrometers (50-2000)
- `width` (number, required): Width in micrometers (50-2000)
- `thickness` (number, required): Thickness in micrometers (20-500)
- `cell_types` (object, required): Dict mapping cell type names to [min_radius, max_radius]

**Example**:
```json
{
  "height": 500,
  "width": 500,
  "thickness": 100,
  "cell_types": {
    "epithelial": [6, 10],
    "stromal": [8, 15],
    "immune": [3, 6]
  }
}
```

**Returns**:
```json
{
  "status": "success",
  "message": "Created tissue: 500x500x100 μm",
  "cell_types": ["epithelial", "stromal", "immune"],
  "cell_type_radii": {...}
}
```

### generate_cells

**Description**: Populate tissue with cells using random sphere packing.

**Parameters**:
- `max_attempts` (integer, optional): Max failed attempts (100-10000, default: 1000)
- `min_spacing` (number, optional): Min spacing between cells (0-10, default: 0.5)
- `allow_boundary_cells` (boolean, optional): Allow boundary cells (default: true)

**Returns**:
```json
{
  "status": "success",
  "num_cells_generated": 234,
  "interior_cells": 198,
  "boundary_cells": 36,
  "packing_fraction": 0.312,
  "cell_type_counts": {"epithelial": 152, "stromal": 82}
}
```

### get_tissue_statistics

**Description**: Get comprehensive statistics about the tissue.

**Parameters**: None

**Returns**:
```json
{
  "total_cells": 234,
  "interior_cells": 198,
  "boundary_cells": 36,
  "packing_fraction": 0.312,
  "cell_types": {"epithelial": 152, "stromal": 82},
  "average_radii": {"epithelial": 7.8, "stromal": 11.2},
  "tissue_dimensions": {"width": 500, "height": 500, "thickness": 100}
}
```

### create_slice

**Description**: Create a 2D slice through the tissue.

**Parameters** (one method required):
- Method 1: `z_position` (number): Z-position for horizontal slice
- Method 2: `angle_x` (number), `angle_y` (number): Rotation angles
- Method 3: `point` (array), `normal` (array): Custom plane definition

**Example (horizontal)**:
```json
{
  "z_position": 50
}
```

**Example (angled)**:
```json
{
  "angle_x": 45,
  "angle_y": 30,
  "point": [250, 250, 50]
}
```

**Returns**:
```json
{
  "status": "success",
  "num_cells_in_slice": 89,
  "plane_point": [250, 250, 50],
  "plane_normal": [0, 0, 1],
  "cell_type_counts": {"epithelial": 58, "stromal": 31},
  "mean_distance_from_plane": 3.2
}
```

### get_slice_statistics

**Description**: Get statistics about the current slice.

**Parameters**: None

**Returns**:
```json
{
  "num_cells": 89,
  "plane_point": [250, 250, 50],
  "plane_normal": [0, 0, 1],
  "cell_types": {"epithelial": 58, "stromal": 31},
  "avg_intersection_radii": {"epithelial": 5.4, "stromal": 7.8},
  "mean_distance_from_plane": 3.2,
  "max_distance_from_plane": 8.7
}
```

### create_serial_slices

**Description**: Create multiple evenly-spaced parallel slices.

**Parameters**:
- `num_slices` (integer, required): Number of slices (2-20)

**Returns**:
```json
{
  "status": "success",
  "num_slices_created": 5,
  "slices": [
    {"slice_number": 1, "z_position": 16.7, "num_cells": 45, "cell_types": {...}},
    {"slice_number": 2, "z_position": 33.3, "num_cells": 68, "cell_types": {...}},
    ...
  ]
}
```

### export_tissue_csv

**Description**: Export 3D tissue data to CSV.

**Parameters**:
- `filename` (string, optional): CSV filename (default: "tissue_data.csv")

**Returns**:
```json
{
  "status": "success",
  "filepath": "/tmp/tissue_sim_xyz/tissue_data.csv",
  "num_cells_exported": 234
}
```

### export_slice_csv

**Description**: Export 2D slice data to CSV.

**Parameters**:
- `filename` (string, optional): CSV filename (default: "slice_data.csv")
- `include_3d` (boolean, optional): Include 3D coordinates (default: true)

**Returns**:
```json
{
  "status": "success",
  "filepath": "/tmp/tissue_sim_xyz/slice_data.csv",
  "num_cells_exported": 89,
  "include_3d_coordinates": true
}
```

### visualize_tissue

**Description**: Create 3D visualization and save as PNG.

**Parameters**:
- `elevation` (number, optional): Viewing elevation (default: 20)
- `azimuth` (number, optional): Viewing azimuth (default: 45)
- `filename` (string, optional): PNG filename (default: "tissue_3d.png")

**Returns**:
```json
{
  "status": "success",
  "filepath": "/tmp/tissue_sim_xyz/tissue_3d.png",
  "cells_visualized": 100
}
```

### visualize_slice_2d

**Description**: Create 2D slice visualization and save as PNG.

**Parameters**:
- `filename` (string, optional): PNG filename (default: "slice_2d.png")

**Returns**:
```json
{
  "status": "success",
  "filepath": "/tmp/tissue_sim_xyz/slice_2d.png",
  "cells_visualized": 89
}
```

### reset_tissue

**Description**: Clear current tissue and start fresh.

**Parameters**: None

**Returns**:
```json
{
  "status": "success",
  "message": "Tissue and slice data cleared. Ready for new simulation."
}
```

## Workflow Examples

### Complete Analysis Workflow

1. Create tissue: `create_tissue`
2. Generate cells: `generate_cells`
3. Get statistics: `get_tissue_statistics`
4. Create slice: `create_slice`
5. Analyze slice: `get_slice_statistics`
6. Export data: `export_tissue_csv`, `export_slice_csv`
7. Create visualizations: `visualize_tissue`, `visualize_slice_2d`

### Serial Section Analysis

1. Create tissue: `create_tissue`
2. Generate cells: `generate_cells`
3. Create serial sections: `create_serial_slices`
4. Analyze each section's statistics
5. Export and visualize

### Parameter Exploration

1. Create tissue with parameters
2. Generate cells multiple times with different `max_attempts`
3. Compare packing fractions
4. Reset and try different cell type configurations

## Error Handling

All tools return JSON with either:
- `"status": "success"` with results
- `"error": "message"` if something went wrong

Common errors:
- "No tissue created. Call create_tissue first."
- "No tissue with cells available."
- "No slice created. Call create_slice first."

## Performance Notes

- Tissue generation is memory-efficient
- Slicing operations are very fast (O(n))
- Visualizations limited to 100 cells for performance
- CSV exports handle thousands of cells efficiently

## Temporary Files

The MCP server stores generated files in a temporary directory:
- Location: `/tmp/tissue_sim_*`
- Includes: CSV exports, PNG visualizations
- Cleaned up on server restart

## Debugging

### Test the server manually:

```bash
# Run server
python run_mcp_server.py

# Server will wait for MCP protocol messages on stdin/stdout
```

### Check if MCP library is installed:

```bash
python -c "import mcp; print('MCP installed')"
```

### Verify package installation:

```bash
python -c "from tissue_simulator.mcp import TissueSimulatorMCPServer; print('Server module loaded')"
```

## Limitations

1. **Stateful**: One tissue/slice at a time per session
2. **No persistence**: Data cleared on server restart
3. **Visualization**: Limited to 100 cells for 3D renders
4. **Single-threaded**: Sequential tool calls only

## Future Enhancements

- [ ] Multiple tissue sessions
- [ ] Persistent storage
- [ ] Streaming progress updates
- [ ] Resource cleanup controls
- [ ] Advanced visualization options
- [ ] Batch operations

## Support

For issues with MCP integration:
1. Check MCP library installation
2. Verify configuration file syntax
3. Test server independently
4. Check Claude Desktop logs

## See Also

- [MCP Documentation](https://modelcontextprotocol.io/)
- [Tissue Simulator Guide](../GUIDE.md)
- [Slicing Documentation](../SLICING.md)
