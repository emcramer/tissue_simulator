# 🧬 Tissue Simulator

A comprehensive Python package for generating 3D simulated biological tissue sections using random sphere packing algorithms.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-emcramer.github.io%2Ftissue__simulator-blue)](https://emcramer.github.io/tissue_simulator/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17465675.svg)](https://doi.org/10.5281/zenodo.17465675)

<!-- get github_id at https://api.github.com/repos/emcramer/tissue_simulator}, but need to wait to make the repo public. use to replace doi button above. -->
<!--[![DOI](https://zenodo.org/badge/{github_id}.svg)](https://zenodo.org/badge/latestdoi/{github_id})-->

## ✨ Features

- **🎯 Flexible Cell Configuration**: Define single or multiple cell types with customizable size ranges
- **🔬 Realistic 3D Sphere Packing**: Random sequential addition algorithm with collision detection
- **🎨 Interactive GUI**: Real-time 3D visualization with PyQt5 interface
- **✂️ 2D Slicing**: Extract planar sections at any angle through the tissue
- **📊 Statistical Analysis**: Automatic calculation of packing fractions and cell distributions
- **💾 Data Export**: CSV output for further analysis and integration with other tools
- **🧪 Boundary Cell Support**: Handle cells that extend beyond tissue section boundaries
- **⚡ Performance Optimized**: Efficient algorithms for generating hundreds of cells

Full documentation, guides, and API reference: [https://emcramer.github.io/tissue_simulator/](https://emcramer.github.io/tissue_simulator/)

## 📦 Installation

### Quick Install

```bash
cd tissue_simulator
pip install -e .
```

### Manual Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Install package
pip install -e .
```

### Verify Installation

```bash
python verify_installation.py
```

## 🚀 Quick Start

### Using the GUI

Launch the interactive interface:

```bash
python -m tissue_simulator.gui
```

Or:

```bash
tissue-simulator
```

### Python API

```python
from tissue_simulator import TissueSection

# Create tissue with uniform cell sizes
tissue = TissueSection(
    height=500,      # μm (Y-axis)
    width=500,       # μm (X-axis)
    thickness=100,   # μm (Z-axis)
    cell_radii=(5, 15)  # min and max radius in μm
)

# Generate cells
num_cells = tissue.generate_cells(max_attempts=1000)
print(f"Generated {num_cells} cells")

# Visualize in 3D
tissue.visualize()

# Export data
tissue.export_to_csv('tissue_data.csv')

# Get statistics
stats = tissue.get_cell_statistics()
print(f"Packing fraction: {stats['packing_fraction']:.3f}")
```

### Multiple Cell Types

```python
tissue = TissueSection(
    height=500,
    width=500,
    thickness=100,
    cell_radii={
        'epithelial': (6, 10),   # Epithelial cells
        'stromal': (8, 15),      # Stromal cells
        'immune': (3, 6),        # Immune cells
        'endothelial': (5, 8)    # Endothelial cells
    }
)

tissue.generate_cells(max_attempts=2000)
tissue.visualize()

# Extract a 2D slice
from tissue_simulator import TissueSlicer
slicer = TissueSlicer(tissue)
slicer.slice_plane(z_position=50)  # Horizontal slice at z=50
slicer.visualize_slice_2d()        # View 2D cross-section
slicer.export_slice_csv('slice_data.csv')  # Export slice data
```

## 📚 Documentation

**Documentation site:** <https://emcramer.github.io/tissue_simulator/> (auto-deployed; per-version archives via the navbar switcher).

**📊 Code-driven tour:** <https://emcramer.github.io/tissue_simulator/latest/slides/tour.html> — a scrolling, end-to-end tour (code + matplotlib output side by side), rendered from `docs/slides/tour.py` (Marimo) on every release.

**Community wiki:** <https://github.com/emcramer/tissue_simulator/wiki> — FAQ, troubleshooting, roadmap (community-editable).

- **[Quick Start Guide](docs/quickstart.md)**: Get up and running in minutes
- **[Core API](docs/api/core.md)**: TissueSection, Cell, SpherePacker reference
- **[Spatial Analysis](docs/api/spatial-analysis.md)**: Network-based spatial analysis guide
- **[Replicate Generation](docs/api/replicate-generation.md)**: Generate tissues matching spatial statistics
- **[Complete Workflow](docs/guides/complete-workflow.md)**: End-to-end tutorial
- **[CHANGELOG](CHANGELOG.md)**: Release history and migration notes

## 💡 Examples

The `examples/` directory contains several demonstration scripts:

### Tissue Generation Examples

### 1. Simple Example
Basic tissue generation with uniform cell sizes:
```bash
python examples/simple_example.py
```

### 2. Multi-Cell Example
Multiple cell types with different size distributions:
```bash
python examples/multi_cell_example.py
```

### 3. Custom Packing
Advanced control using SpherePacker directly:
```bash
python examples/custom_packing.py
```

### 4. Batch Generation
Generate multiple tissue sections for statistical analysis:
```bash
python examples/batch_generation.py
```

### 5. Advanced Visualization
Multiple viewing angles and detailed plots:
```bash
python examples/advanced_visualization.py
```

### Tissue Slicing Examples

### 6. Simple Slicing
Extract a horizontal 2D slice:
```bash
python examples/simple_slicing.py
```

### 7. Angled Slicing
Create slices at various angles:
```bash
python examples/angled_slicing.py
```

### 8. Serial Slices
Generate multiple parallel slices through tissue:
```bash
python examples/serial_slices.py
```

## 🎮 GUI Features

The interactive GUI provides:

### Control Panel
- **Tissue Dimensions**: Sliders for height, width, and thickness (100-1000 μm)
- **Cell Radii Configuration**:
  - Simple mode: Min/max radius sliders
  - JSON mode: Define multiple cell types with custom ranges
- **Packing Parameters**:
  - Max attempts: Controls packing thoroughness
  - Min spacing: Gap between cell surfaces
  - Boundary cells: Allow/disallow cells extending beyond bounds
- **Actions**: Generate, clear, and export buttons with progress tracking

### Visualization Panel
- **3D Viewer**: Interactive matplotlib canvas with rotation controls
- **Statistics Tab**: Real-time metrics including:
  - Total cell count
  - Cell type distributions
  - Packing efficiency
  - Average radii per type

### JSON Configuration Format

Define multiple cell types in the text box:

```json
{
  "epithelial": [5, 10],
  "stromal": [8, 15],
  "immune": [3, 6],
  "endothelial": [5, 8]
}
```
## 🕸️ Spatial Analysis & Graph-Based Cell Type Assignment

### NEW: Graph-Based Cell Type Assignment

Assign cell types to tissues based on target spatial interaction patterns using simulated annealing optimization!

**Complete Workflow:**
1. Generate 3D tissue → 2. Extract 2D slice → 3. Build network graph → 4. Assign cell types → 5. Visualize → 6. Export → 7. Evaluate

```python
from tissue_simulator import TissueSection, SpherePacker, quick_workflow

# Create tissue
tissue = TissueSection(height=300, width=300, thickness=80)
packer = SpherePacker(tissue)
packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1500)

# Define target spatial statistics
target_stats = {
    'node_counts': {'cancer': 40, 'immune': 35, 'stroma': 25},
    'edge_counts': {
        'cancer-cancer': 60,    # Cancer cells cluster
        'cancer-immune': 45,    # Moderate interaction
        'cancer-stroma': 30,
        'immune-immune': 25,
        'immune-stroma': 20,
        'stroma-stroma': 15
    },
    'neighbor_dist': {
        'cancer': {'cancer': 3.0, 'immune': 2.2, 'stroma': 1.5},
        'immune': {'cancer': 2.5, 'immune': 1.8, 'stroma': 1.3},
        'stroma': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.2}
    }
}

# Run complete workflow
from tissue_simulator import TissueNetworkWorkflow
workflow = TissueNetworkWorkflow()
evaluation = workflow.run_complete_workflow(
    tissue=tissue,
    z_position=40,
    network_radius=50.0,
    target_stats_dict=target_stats,
    cell_types=['cancer', 'immune', 'stroma'],
    export_dir="results",
    visualize=True
)

print(f"JS Divergence: {evaluation['js_divergence']:.4f}")  # Lower is better
print(f"Cosine Similarity: {evaluation['cosine_similarity']:.4f}")  # Higher is better
```

**Key Features:**
- **Target-based assignment**: Match spatial statistics from real tissues or hypothetical patterns
- **Simulated annealing**: Optimization-based approach for realistic cell distributions
- **Comprehensive evaluation**: Multiple metrics (JS divergence, cosine similarity, MAE, RMSE)
- **Full workflow integration**: Seamless integration with tissue generation and analysis

**Documentation:**
- **[Complete Workflow Guide](docs/guides/complete-workflow.md)** - Step-by-step tutorial with examples
- **[Graph Coloring Guide](docs/api/graph-coloring.md)** - Detailed API documentation
- **[Example Script](examples/complete_graph_coloring_workflow.py)** - Comprehensive working example

### Network-Based Spatial Analysis

**Network Construction:**
- **Contact mode:** Connects touching cells
- **Radius mode:** Connects cells within distance threshold
- Works with 3D tissues and 2D slices

**Comprehensive Statistics:**
- **Global:** degree, density, clustering, path lengths
- **Per cell type:** degree, clustering, centrality measures
- **Pairwise interactions:** counts (normalized), distances

**Export & Visualization:**
- CSV exports (3 files per analysis)
- Network formats (GraphML, GEXF, GML)
- Network visualizations

### Example Usage

```python
from tissue_simulator import TissueSection, SpatialNetworkAnalyzer

# Generate tissue
tissue = TissueSection(400, 400, 100, cell_radii={'epithelial': (6, 10)})
tissue.generate_cells(max_attempts=1000)

# Analyze
analyzer = SpatialNetworkAnalyzer()
analyzer.build_network_from_tissue(tissue, mode="contact")

# Get statistics
stats = analyzer.compute_global_statistics()
print(f"Avg degree: {stats.avg_degree:.2f}")
print(f"Clustering: {stats.avg_clustering:.4f}")

# Export
analyzer.export_statistics_csv("analysis")
analyzer.visualize_network(save_path="network.png")
```

## 🔄 Replicate Generation

Generate multiple tissue samples matching specific spatial interaction patterns:

### Key Features:

- **Target-based generation**: Match spatial statistics from existing tissues or CSV
- **Iterative optimization**: Automatically tunes parameters to achieve targets
- **Batch processing**: Generate multiple replicates efficiently
- **Statistical validation**: Track divergence from target patterns
- **Full MCP support**: Accessible to LLM coding assistants

### Example Usage

```python
from tissue_simulator import (
    TissueSection,
    load_target_statistics_from_tissue,
    ReplicateGenerator
)

# Create reference tissue
reference = TissueSection(400, 400, 100, 
                         cell_radii={'cancer': (8, 12), 'immune': (5, 8)})
reference.generate_cells(max_attempts=1000)

# Extract spatial statistics
target_stats = load_target_statistics_from_tissue(
    reference, 
    network_mode="contact"
)

# Setup generator
generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'cancer': (8, 12), 'immune': (5, 8)},
    network_mode="contact"
)

# Generate replicates
replicates = generator.generate_replicates(num_replicates=10)

# Export results
generator.export_replicate_statistics(replicates, "output/stats")
generator.export_replicate_tissues(replicates, "output/tissues")
```

### Load from CSV

```python
from tissue_simulator import load_target_statistics_from_csv

# Load target statistics from CSV file
target_stats = load_target_statistics_from_csv("statistics.csv")

generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'cancer': (8, 12), 'immune': (5, 8)},
    network_mode="contact"
)

replicates = generator.generate_replicates(num_replicates=10)
```

## 🤖 LLM Integration (MCP)

The Tissue Simulator can be used as a tool by Large Language Models through the Model Context Protocol (MCP).

### Quick Setup

```bash
# Install MCP support
pip install mcp

# Configure Claude Desktop
# Edit: ~/Library/Application Support/Claude/claude_desktop_config.json
{
  "mcpServers": {
    "tissue-simulator": {
      "command": "python",
      "args": ["/Users/cramere/tissue_simulator/run_mcp_server.py"]
    }
  }
}
```

### Example Usage

Once configured, you can use natural language:

```
You: "Create a tissue with epithelial and stromal cells, 
      then create 5 serial sections and analyze them."

Claude: [Uses create_tissue, generate_cells, create_serial_slices tools]
        I've created a 400x400x100 μm tissue with 234 cells...
```

### Available Tools

#### Tissue Generation
- **create_tissue** - Define tissue dimensions and cell types
- **generate_cells** - Populate with sphere packing
- **get_tissue_statistics** - Analyze composition
- **reset_tissue** - Start fresh

#### 2D Slicing
- **create_slice** - Extract 2D slice at any angle
- **get_slice_statistics** - Analyze slice
- **create_serial_slices** - Create multiple slices

#### Export & Visualization
- **export_tissue_csv** - Export 3D data
- **export_slice_csv** - Export 2D data
- **visualize_tissue** - Generate 3D visualization
- **visualize_slice_2d** - Generate 2D visualization

#### Replicate Generation (New!)
- **load_target_statistics** - Load target spatial statistics from CSV or current tissue
- **setup_replicate_generator** - Configure replicate generator
- **generate_replicates** - Generate multiple replicates matching targets
- **get_replicate_summary** - Get statistics across all replicates
- **export_replicate_statistics** - Export replicate statistics to CSV
- **export_replicate_tissues** - Export each replicate tissue to CSV

### Documentation

- **[MCP Quick Start](docs/guides/mcp.md)** - Get started in 5 minutes
- **[MCP Complete Guide](docs/api/mcp.md)** - Full API reference
- **[Example Conversations](examples/mcp_examples/)** - Usage examples

## 🔬 Scientific Background

### Cell Sizes (typical ranges)
- **Red blood cells**: ~7-8 μm diameter
- **Lymphocytes**: ~6-10 μm
- **Epithelial cells**: ~10-30 μm
- **Fibroblasts**: ~10-20 μm
- **Hepatocytes**: ~20-30 μm

### Tissue Dimensions
- **Histological sections**: 5-10 μm thick
- **Simulated sections**: 50-200 μm thick (for 3D context)
- **Field of view**: 100-1000 μm typical

### Packing Algorithm
The Random Sequential Addition (RSA) algorithm:
- Maximum packing fraction ~38% for monodisperse spheres
- Higher fractions achievable with polydisperse distributions
- Biologically realistic cell distributions
- Efficient for up to thousands of cells

## 🧪 Testing

Run the test suite:

```bash
cd tests
python test_tissue_simulator.py
```

Or with pytest:

```bash
pytest tests/ -v
```

Tests cover:
- Cell class functionality
- TissueSection operations
- SpherePacker algorithm
- Collision detection
- Export/import workflows

## 📊 Output Format

### CSV Export Structure

```csv
x, y, z, radius, cell_type, is_boundary
250.3, 125.7, 50.2, 8.5, epithelial, False
378.1, 442.9, 75.8, 12.3, stromal, True
...
```

Fields:
- `x, y, z`: Cell center coordinates (μm)
- `radius`: Cell radius (μm)
- `cell_type`: Classification string
- `is_boundary`: Boolean indicating if cell extends beyond tissue bounds

## 🎯 Use Cases

1. **Algorithm Development**: Generate synthetic data for testing image analysis pipelines
2. **Machine Learning**: Create training datasets for cell segmentation models
3. **Education**: Visualize 3D tissue structure and spatial organization
4. **Research**: Study cell packing efficiency and spatial statistics
5. **Validation**: Compare simulated vs. real tissue morphology

## 🛠️ Requirements

- **Python**: 3.8 or higher
- **NumPy**: ≥1.20.0 (numerical operations)
- **Matplotlib**: ≥3.3.0 (3D visualization)
- **PyQt5**: ≥5.15.0 (GUI framework)

## 📈 Performance

Typical generation times (on modern hardware):

| Tissue Size | Cell Count | Time |
|------------|------------|------|
| 100×100×50 μm | ~20-50 | <1 second |
| 500×500×100 μm | ~200-400 | 10-30 seconds |
| 1000×1000×200 μm | ~800-1500 | 1-3 minutes |

Performance tips:
- Reduce `max_attempts` for faster generation (lower cell count)
- Smaller cell radii = more cells = longer generation time
- Enable boundary cells for higher packing efficiency

## 📚 Citation

If you use `tissue_simulator` in your research, please cite it. GitHub's
"Cite this repository" button (top right of the repo page) generates APA and
BibTeX from [CITATION.cff](CITATION.cff). BibTeX:

```bibtex
@software{cramer_tissue_simulator,
  author    = {Cramer, Eric},
  title     = {Tissue Simulator: 3D simulated biological tissue section
               generator with network-based spatial analysis},
  year      = {2026},
  version   = {0.1.15},
  doi       = {10.5281/zenodo.17465675},
  url       = {https://github.com/emcramer/tissue_simulator}
}
```

The DOI above always resolves to the latest release; Zenodo also mints a
version-specific DOI for each release if you need to cite an exact version.

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This package implements concepts from:
- Random sphere packing algorithms in computational geometry
- Histological tissue structure principles
- Spatial cell organization in biological systems

## 📧 Support

For issues and questions:

1. Check the [Core API reference](docs/api/core.md)
2. Review the [examples](examples/) directory
3. Run `verify_installation.py` to check setup
4. File an issue on the project repository

---

**Made with ❤️ for the computational biology community**
