# Example MCP Conversations

This directory contains example conversations showing how LLMs can use the Tissue Simulator MCP tools.

## Example 1: Basic Tissue Generation

**User**: "Create a simple tissue with two cell types and show me the statistics."

**LLM Response**:
```
I'll create a tissue simulation with two cell types.

[Calls create_tissue with:
- height: 400, width: 400, thickness: 100
- cell_types: {"epithelial": [6, 10], "stromal": [8, 15]}]

[Calls generate_cells with:
- max_attempts: 1500
- min_spacing: 0.5]

[Calls get_tissue_statistics]

I've created a 400x400x100 μm tissue with:
- 152 epithelial cells (avg radius: 7.8 μm)
- 82 stromal cells (avg radius: 11.2 μm)
- Total: 234 cells
- Packing fraction: 0.312
- 36 cells extend beyond boundaries

This packing fraction of 0.312 indicates moderate packing density, 
which is realistic for biological tissues.
```

## Example 2: Serial Sectioning Analysis

**User**: "Create serial sections through the tissue and analyze how cell counts vary with depth."

**LLM Response**:
```
I'll create 7 serial sections and analyze the distribution.

[Calls create_serial_slices with num_slices: 7]

[Calls get_slice_statistics for each slice]

Here's the cell count profile through the tissue:

Slice 1 (z=12.5 μm): 38 cells
Slice 2 (z=25.0 μm): 56 cells
Slice 3 (z=37.5 μm): 73 cells
Slice 4 (z=50.0 μm): 89 cells ← Peak
Slice 5 (z=62.5 μm): 78 cells
Slice 6 (z=75.0 μm): 61 cells
Slice 7 (z=87.5 μm): 42 cells

The distribution shows:
- Peak cell count in the middle (z=50 μm)
- Gradual decrease toward edges
- This bell-shaped curve is expected for uniform 3D cell distribution
- Edge slices capture fewer cells due to boundary effects
```

## Example 3: Parameter Optimization

**User**: "Find the best max_attempts value to get at least 250 cells in a 400x400x100 tissue."

**LLM Response**:
```
I'll test different max_attempts values to optimize cell count.

Test 1: max_attempts = 1000
[Calls create_tissue, generate_cells]
Result: 234 cells (below target)

Test 2: max_attempts = 2000
[Calls reset_tissue, create_tissue, generate_cells]
Result: 267 cells (above target!)

Test 3: max_attempts = 1500
[Calls reset_tissue, create_tissue, generate_cells]
Result: 248 cells (close but below)

Test 4: max_attempts = 1700
[Calls reset_tissue, create_tissue, generate_cells]
Result: 255 cells (target achieved!)

Optimal value: max_attempts = 1700
This gives 255 cells with a packing fraction of 0.335.

The relationship isn't linear because:
- Early cells place easily
- Later cells require many attempts due to limited space
- Packing efficiency plateaus around 0.35-0.38 for this algorithm
```

## Example 4: Angled Slicing Exploration

**User**: "Show me how slicing at different angles affects the view of the tissue."

**LLM Response**:
```
I'll create slices at different angles and visualize them.

[Calls create_tissue, generate_cells]

Horizontal Slice (0°):
[Calls create_slice with z_position: 50]
[Calls visualize_slice_2d with filename: "slice_0deg.png"]
Result: 89 cells, plane normal: [0, 0, 1]

Tilted 30°:
[Calls create_slice with angle_x: 30, angle_y: 0]
[Calls visualize_slice_2d with filename: "slice_30deg.png"]
Result: 112 cells, plane normal: [0, -0.5, 0.866]

Tilted 45°:
[Calls create_slice with angle_x: 45, angle_y: 0]
[Calls visualize_slice_2d with filename: "slice_45deg.png"]
Result: 127 cells, plane normal: [0, -0.707, 0.707]

Diagonal (45° both axes):
[Calls create_slice with angle_x: 45, angle_y: 45]
[Calls visualize_slice_2d with filename: "slice_diagonal.png"]
Result: 143 cells, plane normal: [0.5, -0.5, 0.707]

Key findings:
- More tilted slices capture more cells (longer path through tissue)
- Horizontal slice: 89 cells
- 45° diagonal: 143 cells (60% more!)
- The diagonal slice travels ~141 μm through the tissue vs 100 μm for horizontal
```

## Example 5: Cell Type Distribution Analysis

**User**: "Create a tissue with 4 cell types and analyze their spatial distribution."

**LLM Response**:
```
I'll create a complex tissue and analyze cell type distribution.

[Calls create_tissue with cell_types:
- "epithelial": [6, 10]
- "stromal": [9, 15]
- "immune": [3, 6]
- "endothelial": [5, 8]
]

[Calls generate_cells with max_attempts: 2500]

Total tissue: 312 cells
- Epithelial: 142 cells (45.5%)
- Stromal: 98 cells (31.4%)
- Immune: 41 cells (13.1%)
- Endothelial: 31 cells (9.9%)

Now analyzing distribution across depth:

[Calls create_serial_slices with num_slices: 5]

Slice 1 (z=16.7):
- Epithelial: 28, Stromal: 19, Immune: 8, Endothelial: 6

Slice 2 (z=33.3):
- Epithelial: 42, Stromal: 28, Immune: 12, Endothelial: 9

Slice 3 (z=50.0):
- Epithelial: 53, Stromal: 37, Immune: 15, Endothelial: 12

Slice 4 (z=66.7):
- Epithelial: 45, Stromal: 31, Immune: 13, Endothelial: 10

Slice 5 (z=83.3):
- Epithelial: 31, Stromal: 21, Immune: 9, Endothelial: 7

Analysis:
- All cell types maintain similar proportions across slices
- Peak density at middle slice (z=50)
- Small cells (immune: 3-6 μm) pack more efficiently
- Large cells (stromal: 9-15 μm) dominate volume despite lower count
```

## Example 6: Export and Data Analysis

**User**: "Generate a tissue, create multiple slices, and export all data for external analysis."

**LLM Response**:
```
I'll create a complete dataset for you.

[Calls create_tissue, generate_cells]
Generated: 267 cells

[Calls export_tissue_csv with filename: "full_tissue.csv"]
Exported: /tmp/tissue_sim_xyz/full_tissue.csv

[Calls create_serial_slices with num_slices: 10]
Created 10 slices

Exporting each slice:
[Calls export_slice_csv for each slice with filenames: 
 "slice_01.csv", "slice_02.csv", ..., "slice_10.csv"]

Files created:
- full_tissue.csv: 267 cells with 3D coordinates
- slice_01.csv through slice_10.csv: 2D projections

Each slice CSV includes:
- 2D coordinates (x_2d, y_2d)
- 3D coordinates (x_3d, y_3d, z_3d)
- Intersection radius
- Cell type
- Distance from plane

You can now:
1. Load these CSVs into pandas/R for analysis
2. Visualize with custom plotting tools
3. Run spatial statistics
4. Compare slice-to-slice variation
```

## Example 7: Comparative Study

**User**: "Compare tissues with different cell size ranges. Which configuration achieves higher packing density?"

**LLM Response**:
```
I'll run a comparative study of different cell size configurations.

Configuration A: Small cells (5-8 μm)
[Calls create_tissue with cell_types: {"cells": [5, 8]}]
[Calls generate_cells with max_attempts: 2000]
Results: 342 cells, packing fraction: 0.298

Configuration B: Medium cells (7-12 μm)
[Calls reset_tissue, create_tissue, generate_cells]
Results: 234 cells, packing fraction: 0.312

Configuration C: Large cells (10-15 μm)
[Calls reset_tissue, create_tissue, generate_cells]
Results: 156 cells, packing fraction: 0.335

Configuration D: Mixed sizes (5-15 μm)
[Calls reset_tissue, create_tissue, generate_cells]
Results: 287 cells, packing fraction: 0.351

Winner: Configuration D (mixed sizes)
- Highest packing fraction: 0.351
- Mixed sizes allow small cells to fill gaps
- This mimics the Apollonian packing principle
- More biologically realistic

Key insight: Polydisperse (mixed size) distributions achieve
higher packing densities than monodisperse (uniform size) distributions.
This is well-established in packing theory!
```

## Example 8: Educational Demonstration

**User**: "Explain sphere packing and demonstrate it visually."

**LLM Response**:
```
I'll demonstrate sphere packing concepts with visualizations.

Concept: Random Sequential Addition (RSA)
The simulator uses RSA to pack spheres (cells) into a volume:

1. Pick a random position
2. Pick a random radius from allowed range
3. Check if it fits without overlapping
4. If yes, place it; if no, try again
5. Stop after too many failed attempts

Let me demonstrate:

Small tissue (200x200x50) for clarity:
[Calls create_tissue, generate_cells with max_attempts: 500]
Result: 45 cells

[Calls visualize_tissue with elevation: 30, azimuth: 45]
[Calls visualize_tissue with elevation: 90, azimuth: 0]

Created views:
- tissue_side.png: Side view showing 3D structure
- tissue_top.png: Top view showing packing pattern

Now a slice to see cross-sections:
[Calls create_slice with z_position: 25]
[Calls visualize_slice_2d]

Key observations:
1. No overlaps (enforced by collision detection)
2. Random positions (not crystalline)
3. Packing fraction ~32% (typical for RSA)
4. Boundary effects (fewer cells at edges)
5. Circle packing in 2D slice shows cell cross-sections

Compare to optimal packing:
- Random: ~32-38%
- FCC crystal: ~74%
- Biology typically: ~30-40%

Our random packing is biologically realistic!
```

## Running These Examples

These are示例 example conversations. To actually execute them:

1. Configure the MCP server in Claude Desktop
2. Start a conversation with Claude
3. Use natural language similar to the examples above
4. Claude will automatically use the appropriate tools

## Creating Your Own

You can ask Claude to:
- Generate specific tissue configurations
- Run comparative studies
- Create visualizations
- Export data for analysis
- Explore parameter spaces
- Explain biological concepts

The tools are flexible and can be combined in many ways!
