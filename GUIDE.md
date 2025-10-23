# Tissue Simulator User Guide

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Core Concepts](#core-concepts)
4. [API Reference](#api-reference)
5. [GUI Usage](#gui-usage)
6. [Advanced Topics](#advanced-topics)
7. [Troubleshooting](#troubleshooting)

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Install from source

```bash
cd tissue_simulator
pip install -e .
```

### Dependencies

The following packages will be automatically installed:
- numpy: Numerical computations
- matplotlib: 3D visualization
- PyQt5: GUI framework

## Quick Start

### Using the GUI

Launch the graphical interface:

```bash
python -m tissue_simulator.gui
```

Or if installed as a command-line tool:

```bash
tissue-simulator
```

### Basic Python Usage

```python
from tissue_simulator import TissueSection

# Create tissue
tissue = TissueSection(
    height=500,      # μm
    width=500,       # μm
    thickness=100,   # μm
    cell_radii=(5, 15)  # min and max radius in μm
)

# Generate cells
tissue.generate_cells(max_attempts=1000)

# Visualize
tissue.visualize()

# Export data
tissue.export_to_csv('tissue_data.csv')
```

## Core Concepts

### Tissue Dimensions

The tissue section is defined as a 3D rectangular prism with:
- **Height**: Y-axis dimension (micrometers)
- **Width**: X-axis dimension (micrometers)
- **Thickness**: Z-axis dimension (micrometers)

### Cell Types and Radii

Cells can be defined in two ways:

1. **Uniform cells** - Single radius range:
   ```python
   cell_radii=(5, 15)  # All cells between 5-15 μm
   ```

2. **Multiple cell types** - Dictionary mapping types to ranges:
   ```python
   cell_radii={
       'epithelial': (6, 10),
       'stromal': (8, 15),
       'immune': (3, 6)
   }
   ```

### Sphere Packing Algorithm

The package uses a random sphere packing algorithm:

1. Randomly select a cell type and radius
2. Generate a random position within the tissue
3. Check for collisions with existing cells
4. Place cell if valid, otherwise retry
5. Continue until max attempts reached

### Boundary Cells

Cells whose boundaries extend beyond the tissue section are called **boundary cells**. These represent cells that would be partially sectioned in a real tissue sample.

- When `allow_boundary_cells=True`: Cell centers can be anywhere in tissue
- When `allow_boundary_cells=False`: Entire cell must fit within bounds

## API Reference

### TissueSection Class

```python
TissueSection(height, width, thickness, cell_radii)
```

**Parameters:**
- `height` (float): Y-dimension in micrometers
- `width` (float): X-dimension in micrometers
- `thickness` (float): Z-dimension in micrometers
- `cell_radii`: Either `(min, max)` tuple or dict of cell types

**Methods:**

#### generate_cells()

```python
generate_cells(max_attempts=1000, min_spacing=0.5, allow_boundary_cells=True)
```

Generate cells using random sphere packing.

**Parameters:**
- `max_attempts` (int): Maximum failed placement attempts before stopping
- `min_spacing` (float): Minimum gap between cell surfaces (μm)
- `allow_boundary_cells` (bool): Allow cells extending beyond bounds

**Returns:** Number of cells placed (int)

#### get_cell_statistics()

```python
get_cell_statistics()
```

Calculate statistics about the packed cells.

**Returns:** Dictionary containing:
- `total_cells`: Total number of cells
- `boundary_cells`: Number of boundary cells
- `interior_cells`: Number of interior cells
- `cell_types`: Count per cell type
- `avg_radii`: Average radius per cell type
- `packing_fraction`: Volume fraction occupied by cells

#### visualize()

```python
visualize(show_boundary=True, elevation=20, azimuth=45)
```

Create 3D visualization of the tissue.

**Parameters:**
- `show_boundary` (bool): Show bounding box
- `elevation` (float): Viewing elevation angle (degrees)
- `azimuth` (float): Viewing azimuth angle (degrees)

#### export_to_csv()

```python
export_to_csv(filename)
```

Export cell data to CSV file.

**Parameters:**
- `filename` (str): Output file path

**CSV Format:**
```
x, y, z, radius, cell_type, is_boundary
250.3, 125.7, 50.2, 8.5, epithelial, False
...
```

### Cell Class

```python
Cell(center, radius, cell_type="default", is_boundary=False)
```

Represents a single cell.

**Attributes:**
- `center` (numpy.ndarray): 3D position [x, y, z]
- `radius` (float): Cell radius in micrometers
- `cell_type` (str): Classification of the cell
- `is_boundary` (bool): Whether cell extends beyond bounds

**Methods:**

#### intersects()

```python
intersects(other_cell)
```

Check if this cell intersects another cell.

#### is_within_bounds()

```python
is_within_bounds(bounds)
```

Check if cell is completely within tissue bounds.

### SpherePacker Class

```python
SpherePacker(bounds, cell_radii_config, min_spacing=0.5, allow_boundary_cells=True)
```

Low-level packing algorithm for advanced usage.

**Methods:**

#### pack()

```python
pack(max_attempts=1000)
```

Pack cells and return list of Cell objects.

#### pack_with_progress()

```python
pack_with_progress(max_attempts=1000, callback=None)
```

Pack cells with progress callback for GUI or monitoring.

**Parameters:**
- `callback`: Function called as `callback(cells_placed, total_attempts)`

## GUI Usage

### Main Interface

The GUI consists of two main panels:

#### Left Panel - Controls

1. **Tissue Dimensions**: Sliders for height, width, and thickness
2. **Cell Radii Configuration**:
   - Simple Range: Two sliders for min/max radius
   - JSON Cell Types: Text box for complex configurations
3. **Packing Parameters**:
   - Max Attempts: How many failed placements before stopping
   - Min Spacing: Gap between cell surfaces
   - Allow Boundary Cells: Checkbox to enable/disable
4. **Action Buttons**:
   - Generate Tissue: Start cell generation
   - Clear: Remove all cells
   - Export to CSV: Save data to file

#### Right Panel - Visualization

Two tabs:

1. **3D Viewer**:
   - Interactive 3D visualization
   - Elevation and Azimuth sliders to rotate view
   - Color-coded by cell type
   - Transparent boundary cells

2. **Statistics**:
   - Cell counts and distributions
   - Packing efficiency metrics
   - Per-cell-type statistics

### JSON Configuration Format

For multiple cell types, use JSON format in the text box:

```json
{
  "epithelial": [5, 10],
  "stromal": [8, 15],
  "immune": [3, 6],
  "endothelial": [5, 8]
}
```

Each entry is: `"cell_type_name": [min_radius, max_radius]`

### Workflow

1. Set tissue dimensions using sliders
2. Configure cell radii (simple or JSON)
3. Adjust packing parameters if needed
4. Click "Generate Tissue"
5. Wait for generation (progress bar shows status)
6. Explore in 3D viewer
7. View statistics in Statistics tab
8. Export to CSV for further analysis

## Advanced Topics

### Custom Packing Strategies

For more control, use `SpherePacker` directly:

```python
from tissue_simulator import SpherePacker

packer = SpherePacker(
    bounds=(500, 500, 100),
    cell_radii_config={'type_a': (5, 10)},
    min_spacing=0.5,
    allow_boundary_cells=True
)

cells = packer.pack(max_attempts=2000)
```

### Analyzing Cell Distributions

```python
import numpy as np

# Get cell positions
positions = np.array([cell.center for cell in tissue.cells])

# Spatial analysis
x_coords = positions[:, 0]
y_coords = positions[:, 1]
z_coords = positions[:, 2]

# Calculate nearest neighbor distances
from scipy.spatial.distance import cdist
distances = cdist(positions, positions)
np.fill_diagonal(distances, np.inf)
nearest_neighbors = distances.min(axis=1)

print(f"Mean nearest neighbor: {nearest_neighbors.mean():.2f} μm")
```

### Batch Processing

Generate multiple sections for statistical analysis:

```python
sections = []
for i in range(10):
    tissue = TissueSection(500, 500, 100, cell_radii=(5, 15))
    tissue.generate_cells(max_attempts=1500)
    sections.append(tissue)
    
# Compare packing fractions
fractions = [t.get_cell_statistics()['packing_fraction'] 
             for t in sections]
print(f"Mean packing: {np.mean(fractions):.3f} ± {np.std(fractions):.3f}")
```

### Optimizing Performance

For large tissues or many cells:

1. **Increase max_attempts**: More attempts = more cells, but slower
2. **Adjust min_spacing**: Smaller spacing = more cells, but harder to pack
3. **Use appropriate dimensions**: Smaller tissues generate faster
4. **Batch with threading**: Generate multiple sections in parallel

```python
from concurrent.futures import ThreadPoolExecutor

def generate_section(params):
    tissue = TissueSection(**params)
    tissue.generate_cells()
    return tissue

params_list = [{'height': 500, 'width': 500, 'thickness': 100, 
                'cell_radii': (5, 15)} for _ in range(5)]

with ThreadPoolExecutor(max_workers=4) as executor:
    tissues = list(executor.map(generate_section, params_list))
```

## Troubleshooting

### Common Issues

**1. Few cells generated**

Problem: Only a small number of cells are placed.

Solutions:
- Increase `max_attempts`
- Reduce `min_spacing`
- Check that cell radii are appropriate for tissue size
- Enable boundary cells

**2. GUI freezes during generation**

Problem: Application becomes unresponsive.

Solutions:
- This is normal for large generations (30+ seconds)
- Reduce tissue size or max_attempts for faster results
- The progress bar updates every 10 cells

**3. JSON parse error**

Problem: "Invalid JSON format" error in GUI.

Solutions:
- Check JSON syntax (commas, brackets, quotes)
- Use double quotes, not single quotes
- Ensure array format: `[min, max]`
- Validate JSON at jsonlint.com

**4. Import errors**

Problem: `ModuleNotFoundError` when running examples.

Solutions:
- Install package: `pip install -e .`
- Check Python version: `python --version` (need 3.8+)
- Verify dependencies: `pip list`

**5. Visualization doesn't show**

Problem: `visualize()` produces no output.

Solutions:
- Check matplotlib backend
- Try: `import matplotlib; matplotlib.use('TkAgg')`
- On macOS, may need: `export MPLBACKEND=TkAgg`
- Ensure cells were generated: check `len(tissue.cells)`

### Performance Tips

- **Small tissues**: < 500 μm sides generate in seconds
- **Medium tissues**: 500-1000 μm sides generate in 10-30 seconds
- **Large tissues**: > 1000 μm sides may take minutes
- Use lower `max_attempts` for quick tests (500-1000)
- Use higher `max_attempts` for final simulations (2000-5000)

### Getting Help

For issues not covered here:

1. Check the examples in the `examples/` directory
2. Review the source code in `tissue_simulator/`
3. File an issue on the project repository

## References

### Biological Context

Cell sizes vary by type:
- Red blood cells: ~7-8 μm diameter
- Lymphocytes: ~6-10 μm
- Epithelial cells: ~10-30 μm
- Fibroblasts: ~10-20 μm
- Hepatocytes: ~20-30 μm

Typical tissue dimensions:
- Histological sections: 5-10 μm thick
- Simulated sections: 50-200 μm thick (for 3D context)
- Field of view: 100-1000 μm wide

### Algorithm Notes

The random sequential addition (RSA) algorithm used here:
- Simple and robust
- Maximum packing fraction ~38% for monodisperse spheres
- Higher fractions possible with polydisperse distributions
- Not optimal but biologically realistic
