# 🤖 Using Tissue Simulator with Large Language Models

## Quick Start

The Tissue Simulator can be used as a tool by Large Language Models (LLMs) through the Model Context Protocol (MCP).

### 1. Install MCP Support

```bash
pip install mcp
```

### 2. Configure Claude Desktop

Edit your Claude Desktop config file and add:

**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "tissue-simulator": {
      "command": "python",
      "args": ["/Users/cramere/tissue_simulator/run_mcp_server.py"]
    }
  }
}
```

### 3. Restart Claude Desktop

The tissue simulator tools will now be available!

## Example Conversations

### Generate and Analyze Tissue

```
You: Can you create a simulated tissue section with epithelial and stromal cells?

Claude: I'll create a tissue simulation for you.

[Uses create_tissue tool]
[Uses generate_cells tool]
[Uses get_tissue_statistics tool]

I've created a 400x400x100 μm tissue with:
- 152 epithelial cells (radius: 6-10 μm)
- 82 stromal cells (radius: 8-15 μm)
- Packing fraction: 0.312
- 36 boundary cells

Would you like me to create visualizations or export the data?
```

### Create Serial Sections

```
You: Create 5 serial sections through the tissue and analyze the cell distribution.

Claude: I'll create serial sections and analyze them.

[Uses create_serial_slices tool]
[Uses get_slice_statistics for each]

Created 5 sections:
- Section 1 (z=16.7 μm): 45 cells
- Section 2 (z=33.3 μm): 68 cells
- Section 3 (z=50.0 μm): 89 cells (peak)
- Section 4 (z=66.7 μm): 72 cells
- Section 5 (z=83.3 μm): 51 cells

The cell density peaks in the middle section, which is typical
for uniform cell distributions.
```

### Export and Visualize

```
You: Export the data and create visualizations.

Claude: I'll export the data and generate visualizations.

[Uses export_tissue_csv tool]
[Uses export_slice_csv tool]
[Uses visualize_tissue tool]
[Uses visualize_slice_2d tool]

Data exported to:
- /tmp/tissue_sim_xyz/tissue_data.csv (234 cells)
- /tmp/tissue_sim_xyz/slice_data.csv (89 cells)

Visualizations created:
- /tmp/tissue_sim_xyz/tissue_3d.png
- /tmp/tissue_sim_xyz/slice_2d.png

You can download these files from the paths shown.
```

## Available Tools

1. **create_tissue** - Define tissue dimensions and cell types
2. **generate_cells** - Populate with cells using sphere packing
3. **get_tissue_statistics** - Analyze tissue composition
4. **create_slice** - Extract 2D slice at any angle
5. **get_slice_statistics** - Analyze slice composition
6. **create_serial_slices** - Create multiple parallel slices
7. **export_tissue_csv** - Export 3D tissue data
8. **export_slice_csv** - Export 2D slice data
9. **visualize_tissue** - Generate 3D visualization
10. **visualize_slice_2d** - Generate 2D visualization
11. **reset_tissue** - Start fresh simulation

## Use Cases with LLMs

### 1. Exploratory Analysis

```
"Generate tissues with different cell type ratios and compare 
their packing efficiencies."
```

The LLM will:
- Create multiple tissues
- Generate cells with different parameters
- Compare statistics
- Provide insights

### 2. Histology Simulation

```
"Create a tissue and generate serial sections like in real 
histology. Analyze how cell counts vary across sections."
```

The LLM will:
- Create tissue
- Generate serial sections
- Analyze each section
- Identify patterns

### 3. Parameter Optimization

```
"Find the optimal max_attempts parameter to achieve a packing 
fraction of at least 0.35."
```

The LLM will:
- Try different parameters
- Generate and analyze tissues
- Converge on optimal value
- Explain the results

### 4. Educational Demonstrations

```
"Explain how tissue slicing works and demonstrate it with 
visualizations at different angles."
```

The LLM will:
- Create tissue
- Make slices at various angles
- Generate visualizations
- Explain the geometry

## Advanced Usage

### Custom Workflows

The LLM can create sophisticated workflows like:

```python
# Pseudo-code of what the LLM might do:
1. create_tissue(specific parameters)
2. generate_cells(optimized parameters)
3. create_serial_slices(10 slices)
4. for each slice:
     analyze statistics
     export data
5. compare all slices
6. generate visualizations
7. provide comprehensive report
```

### Comparative Studies

```
"Compare tissues with 2 cell types vs 4 cell types. Which has 
better packing efficiency?"
```

The LLM will design and execute a comparative study.

### Data Analysis

```
"Analyze how cell density varies with depth in the tissue."
```

The LLM will:
- Create serial sections
- Extract statistics
- Identify trends
- Provide analysis

## Tips for LLM Interactions

### Be Specific

Avoid vague prompts like "Make a tissue".

Prefer concrete prompts like "Create a 500x500x100 μm tissue with epithelial cells (7-11 μm) and stromal cells (9-15 μm)".

### Sequential Requests

The tools maintain state, so you can build on previous actions:

```
1. "Create a tissue with 3 cell types"
2. "Now create a slice at z=50"
3. "Export that slice data"
4. "Create another slice at 45 degrees"
```

### Ask for Analysis

The LLM can interpret the results:

```
"What does a packing fraction of 0.35 mean? Is that good?"
"Why are there fewer cells in the edge slices?"
"How does the min_spacing parameter affect the results?"
```

## Troubleshooting

### Tools Not Showing Up

1. Check MCP is installed: `pip install mcp`
2. Verify config file location
3. Check JSON syntax in config
4. Restart Claude Desktop
5. Look for error messages

### Server Not Running

1. Test manually: `python run_mcp_server.py`
2. Check Python path in config
3. Verify package is installed: `pip install -e .`

### Unexpected Results

1. Use `reset_tissue` to start fresh
2. Check parameter ranges
3. Ask the LLM to explain what happened

## Performance Notes

- Tissue generation: ~10-30 seconds for medium tissues
- Slicing: Near-instant (milliseconds)
- Visualization: ~2-5 seconds
- CSV export: Fast for thousands of cells

## Security Notes

- MCP server runs locally
- No network access required
- Files stored in temp directory
- Cleaned on server restart

## Learn More

- Full documentation: [MCP API reference](../api/mcp.md)
- MCP Protocol: https://modelcontextprotocol.io/
- Tissue Simulator: [README.md](README.md)

## Examples Repository

More example conversations and workflows: [examples/mcp_examples/](examples/mcp_examples/)

---

**Now you can use natural language to create and analyze tissue simulations!** 🧬🤖
