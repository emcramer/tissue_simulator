# 🧬 Tissue Simulator

A comprehensive Python package for generating 3D simulated biological tissue sections using random sphere packing algorithms.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!--[![DOI]()](https://doi.org/10.5281/zenodo.17465675)-->
<!-- get github_id at https://api.github.com/repos/emcramer/tissue_simulator}, but need to wait to make the repo public-->
[![DOI](https://zenodo.org/badge/{github_id}.svg)](https://zenodo.org/badge/latestdoi/{github_id})

## ✨ Features

- **🎯 Flexible Cell Configuration**: Define single or multiple cell types with customizable size ranges
- **🔬 Realistic 3D Sphere Packing**: Random sequential addition algorithm with collision detection
- **🎨 Interactive GUI**: Real-time 3D visualization with PyQt5 interface
- **✂️ 2D Slicing**: Extract planar sections at any angle through the tissue
- **📊 Statistical Analysis**: Automatic calculation of packing fractions and cell distributions
- **💾 Data Export**: CSV output for further analysis and integration with other tools
- **🧪 Boundary Cell Support**: Handle cells that extend beyond tissue section boundaries
- **⚡ Performance Optimized**: Efficient algorithms for generating hundreds of cells

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

- **[Quick Start Guide](QUICKSTART.md)**: Get up and running in minutes
- **[Comprehensive Guide](GUIDE.md)**: Detailed documentation with examples
- **[Package Structure](STRUCTURE.md)**: Architecture and module descriptions
- **[API Reference](GUIDE.md#api-reference)**: Complete function and class documentation

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

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This package implements concepts from:
- Random sphere packing algorithms in computational geometry
- Histological tissue structure principles
- Spatial cell organization in biological systems

## 📧 Support

For issues and questions:

1. Check the [Comprehensive Guide](GUIDE.md)
2. Review the [examples](examples/) directory
3. Run `verify_installation.py` to check setup
4. File an issue on the project repository

---

**Made with ❤️ for the computational biology community**
