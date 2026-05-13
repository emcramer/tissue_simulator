# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tissue Simulator is a Python package for generating 3D simulated biological tissue sections using random sphere packing algorithms. It supports network-based spatial analysis, graph-based cell type assignment via simulated annealing, and tissue replicate generation matching target spatial statistics.

## Development Commands

### Installation and Setup
```bash
# Install in development mode
pip install -e .

# Install with all dependencies
pip install -r requirements.txt

# Verify installation
python verify_installation.py
```

### Running Tests
```bash
# Run all tests
pytest tests/ -v

# Run specific test files
python tests/test_tissue_simulator.py
python tests/test_slicing.py

# Quick integration test
python test_graph_coloring_integration.py
```

### Running Examples
```bash
# Basic tissue generation
python examples/simple_example.py
python examples/multi_cell_example.py

# Tissue slicing
python examples/simple_slicing.py
python examples/serial_slices.py

# Spatial analysis
python examples/simple_spatial_analysis.py
python examples/comparative_spatial_analysis.py

# Complete workflow with graph coloring
python examples/complete_graph_coloring_workflow.py

# Replicate generation
python examples/replicate_generation_example.py
```

### GUI Application
```bash
# Launch interactive GUI
python -m tissue_simulator.gui
# OR
tissue-simulator
```

### MCP Server (LLM Integration)
```bash
# Run MCP server for LLM integration
python run_mcp_server.py

# Test MCP client
python test_mcp_client.py
```

## Code Architecture

### Core Modules (tissue_simulator/)

**tissue.py** - Foundation classes
- `Cell`: Represents individual cells with position, radius, type, boundary status
- `TissueSection`: Container for 3D tissue with dimensions and cell collection
- Handles visualization, statistics, CSV export/import

**packing.py** - Cell placement algorithm
- `SpherePacker`: Random Sequential Addition (RSA) algorithm for sphere packing
- Collision detection with minimum spacing constraints
- Supports multiple cell types with different size distributions
- Boundary cell handling (cells can extend beyond tissue bounds)

**slicing.py** - 2D sectioning
- `TissueSlicer`: Extracts 2D planar sections from 3D tissue
- `SliceCell`: Represents cells captured in 2D slice with intersection geometry
- Supports arbitrary slice angles and orientations
- `create_standard_slices()`: Generate multiple parallel slices

**spatial_analysis.py** - Network-based analysis
- `SpatialNetworkAnalyzer`: Builds NetworkX graphs from tissues/slices
- Two modes: "contact" (touching cells) and "radius" (within distance threshold)
- Computes global metrics: degree, density, clustering, path lengths
- Per-cell-type statistics: centrality measures, degree distributions
- Pairwise interaction statistics: edge counts, distances
- Exports to CSV and network formats (GraphML, GEXF, GML)

**graph_coloring.py** - Cell type assignment
- `GraphColorizer`: Assigns cell types to networks using simulated annealing
- Matches target statistics: node counts, edge counts, neighbor distributions
- Optimizes via energy function with configurable weights
- `calculate_graph_statistics()`: Extract statistics from colored graphs
- `compare_graph_statistics()`: Compare source vs target statistics

**replicate_generator.py** - Batch generation
- `ReplicateGenerator`: Generates multiple tissues matching target statistics
- Iteratively tunes parameters to achieve spatial interaction patterns
- `load_target_statistics_from_tissue()`: Extract stats from existing tissue
- `load_target_statistics_from_csv()`: Load stats from CSV file
- Tracks divergence scores across replicates

**tissue_workflow.py** - Unified interface
- `TissueNetworkWorkflow`: Complete pipeline manager
- `quick_workflow()`: Convenience function for full workflow
- Integrates: tissue generation → slicing → network building → cell type assignment → evaluation

**evaluation.py** - Comparison metrics
- `js_divergence()`: Jensen-Shannon divergence (lower = better match)
- `cosine_similarity()`: Vector similarity (higher = better match)
- `evaluate_graph_coloring()`: Comprehensive comparison with multiple metrics
- `print_evaluation_report()`: Formatted output of evaluation results

**gui.py** - Interactive application
- `TissueSimulatorGUI`: PyQt5 main window with controls
- `TissueViewer3D`: Matplotlib 3D visualization canvas
- `PackingThread`: Background thread for non-blocking generation
- Real-time statistics, JSON cell type configuration

**mcp/server.py** - LLM integration
- Exposes tissue simulator as Model Context Protocol (MCP) server
- Tools for: tissue creation, cell generation, slicing, analysis, replicates
- Enables natural language control via LLMs (e.g., Claude Desktop)

### Key Design Patterns

**Coordinate System**:
- X-axis: width dimension
- Y-axis: height dimension
- Z-axis: thickness dimension
- All units in micrometers (μm)

**Cell Types**:
- Specified as dictionaries: `{'type_name': (min_radius, max_radius)}`
- Single type: Can use tuple `(min, max)` which becomes `{'default': (min, max)}`
- Random sampling from uniform distribution within radius range

**Network Building**:
- "contact" mode: Edges between cells whose surfaces touch (distance < r1 + r2)
- "radius" mode: Edges between cells within specified distance threshold
- Both modes work on 3D tissues or 2D slices

**Target Statistics Format**:
- `node_counts`: Dict of cell type counts `{'cancer': 40, 'immune': 35}`
- `edge_counts`: Dict of pairwise edge counts `{'cancer-cancer': 60, 'cancer-immune': 45}`
- `neighbor_dist`: Nested dict of average neighbor counts by type
  ```python
  {'cancer': {'cancer': 3.0, 'immune': 2.2, 'stroma': 1.5}}
  ```

**Simulated Annealing Parameters**:
- Temperature schedule: exponential cooling `T = T_initial * (cooling_rate)^step`
- Energy weights: balance node count, edge count, neighbor distribution matching
- Default iterations: 100-500 depending on graph size
- Random neighbor swaps to explore color assignment space

## Important Implementation Details

**Sphere Packing Algorithm**:
- Maximum theoretical packing fraction: ~38% for monodisperse spheres
- Polydisperse distributions achieve higher packing (40-50%)
- `max_attempts` controls thoroughness vs. speed tradeoff
- Progress tracking via callbacks for GUI integration

**Boundary Cells**:
- When `allow_boundary_cells=True`: Cell centers must be in bounds, but cells can extend beyond
- When `allow_boundary_cells=False`: Entire cell must fit within bounds
- Boundary cells marked with `is_boundary=True` attribute

**Slicing Geometry**:
- Plane defined by point and normal vector
- Cells intersecting plane contribute circular cross-sections
- Intersection radius calculated: `sqrt(r² - d²)` where d = distance to plane
- 2D coordinates computed using orthonormal basis in slice plane

**CSV Export Formats**:

3D Tissue:
```csv
x, y, z, radius, cell_type, is_boundary
250.3, 125.7, 50.2, 8.5, epithelial, False
```

2D Slice:
```csv
x, y, radius, cell_type, center_3d_x, center_3d_y, center_3d_z, distance_from_plane
125.4, 78.2, 7.3, cancer, 125.4, 78.2, 50.1, 2.1
```

Network Statistics (3 files):
- `*_global_stats.csv`: Network-level metrics
- `*_cell_type_stats.csv`: Per-type statistics
- `*_interaction_stats.csv`: Pairwise interactions

**Graph Coloring Workflow**:
1. Generate tissue with placeholder cell types
2. Extract 2D slice
3. Build network graph (nodes = cells, edges = interactions)
4. Load/define target statistics
5. Run simulated annealing to find optimal cell type assignment
6. Evaluate match quality with multiple metrics
7. Export colored graph and statistics

## Common Pitfalls

**NetworkX Dependency**:
- Required for spatial analysis, graph coloring, and replicate generation
- Not required for basic tissue generation and visualization
- Check `NETWORKX_AVAILABLE` flag before using network features

**MCP Configuration**:
- MCP server requires `mcp` package: `pip install mcp`
- Path in `claude_desktop_config.json` must be absolute
- Server runs as subprocess, check logs in Claude Desktop for errors

**Statistics Matching**:
- Target statistics must be achievable given network structure
- Highly constrained targets may not converge
- Increase simulated annealing iterations if poor match
- Edge count statistics are symmetric: 'A-B' = 'B-A'

**Memory Considerations**:
- Large tissues (>2000 cells) can be slow to visualize
- Network analysis scales as O(n²) for dense graphs
- Consider using slices for analysis of large 3D tissues

**Cell Type Names**:
- Must be valid Python dictionary keys (strings)
- Consistent naming required across target statistics and cell types list
- Case-sensitive matching

## Testing Strategy

**Unit Tests** (`tests/test_tissue_simulator.py`):
- Cell class: intersection, bounds checking
- TissueSection: initialization, generation, statistics
- SpherePacker: collision detection, packing algorithm

**Integration Tests**:
- Full workflow tests in root directory (e.g., `test_graph_coloring_integration.py`)
- Example scripts serve as functional tests
- MCP client tests for LLM integration

**Test Data**:
- Small tissues (100×100×50) for fast tests
- CSV exports for round-trip validation
- Known packing fractions for validation

## Documentation Structure

- `README.md`: Main project documentation with examples
- `CHANGELOG.md`: Release history (v0.1.0+)
- `docs/quickstart.md`: Installation and basic usage
- `docs/api/core.md`: Core tissue/cell/packer API reference
- `docs/api/slicing.md`: 2D slicing API
- `docs/api/spatial-analysis.md`: Network analysis API
- `docs/api/graph-coloring.md`: Cell-type-assignment API
- `docs/api/replicate-generation.md`: Replicate generator API
- `docs/api/mcp.md`: MCP server API (LLM integration)
- `docs/guides/complete-workflow.md`: End-to-end tutorial
- `docs/guides/mcp.md`: MCP 5-minute quickstart
- `docs/notes/`: Maintainer notes (e.g., known-issue regressions)
- `docs/design/`: Research / paper-track artifacts

## Special Considerations

**When modifying packing algorithm**: Test with various `max_attempts` values and verify packing fractions don't exceed theoretical limits

**When adding new statistics**: Update both calculation and comparison functions in `spatial_analysis.py` and `graph_coloring.py`

**When changing network building**: Ensure symmetry of edge creation and consistent distance calculations

**When modifying GUI**: Use `PackingThread` for long operations to prevent UI freezing

**When updating MCP tools**: Update schema in `mcp/server.py` and test with `test_mcp_client.py`
