# Tissue Simulator Package Structure

```
tissue_simulator/
│
├── README.md                          # Project overview and basic usage
├── GUIDE.md                           # Comprehensive user guide
├── QUICKSTART.md                      # Quick installation and usage
├── LICENSE                            # MIT License
├── MANIFEST.in                        # Packaging manifest
├── requirements.txt                   # Python dependencies
├── setup.py                           # Package installation configuration
├── .gitignore                         # Git ignore patterns
│
├── tissue_simulator/                  # Main package directory
│   ├── __init__.py                    # Package initialization
│   ├── tissue.py                      # Core TissueSection and Cell classes
│   ├── packing.py                     # SpherePacker algorithm
│   └── gui.py                         # PyQt5 GUI application
│
├── examples/                          # Example scripts
│   ├── simple_example.py              # Basic usage example
│   ├── multi_cell_example.py          # Multiple cell types
│   ├── custom_packing.py              # Advanced packing control
│   ├── batch_generation.py            # Batch processing
│   └── advanced_visualization.py      # Multiple visualization views
│
└── tests/                             # Unit tests
    ├── __init__.py                    # Test package initialization
    └── test_tissue_simulator.py       # Comprehensive test suite
```

## Module Descriptions

### Core Modules

**tissue.py**
- `Cell` class: Represents individual cells with position, radius, and type
- `TissueSection` class: Main container for tissue simulation
- Methods for generation, statistics, visualization, and export

**packing.py**
- `SpherePacker` class: Random sphere packing algorithm
- Collision detection and boundary handling
- Progress tracking for GUI integration

**gui.py**
- `TissueSimulatorGUI`: Main PyQt5 window
- `TissueViewer3D`: 3D matplotlib visualization canvas
- `PackingThread`: Background thread for non-blocking generation
- Interactive sliders, JSON configuration, and real-time updates

### Example Scripts

**simple_example.py**
- Basic tissue generation with uniform cell sizes
- Visualization and CSV export
- Good starting point for new users

**multi_cell_example.py**
- Multiple cell types with different size ranges
- Detailed statistics output
- Custom viewing angles

**custom_packing.py**
- Direct use of SpherePacker for fine control
- Progress callbacks
- Density calculations

**batch_generation.py**
- Generate multiple tissue sections
- Statistical analysis across batches
- Comparison plots

**advanced_visualization.py**
- Multiple 3D viewing angles
- 2D projections
- Size distribution histograms
- Cell type breakdowns

### Tests

**test_tissue_simulator.py**
- Unit tests for Cell class
- Unit tests for TissueSection class
- Unit tests for SpherePacker class
- Integration tests for complete workflows
- Export/import validation

## Key Features

### 1. Flexible Cell Configuration
- Simple radius range: `(min, max)`
- Multiple cell types: `{'type': (min, max), ...}`

### 2. Realistic Sphere Packing
- Random sequential addition algorithm
- Configurable minimum spacing
- Boundary cell handling

### 3. Interactive GUI
- Real-time 3D visualization
- Adjustable tissue dimensions
- JSON cell type configuration
- Progress monitoring
- Statistics display

### 4. Data Export
- CSV format with cell positions and metadata
- Compatible with analysis tools
- Batch processing support

### 5. Comprehensive Testing
- Unit tests for all core functionality
- Integration tests for workflows
- Validation of packing algorithms

## Installation

```bash
cd tissue_simulator
pip install -e .
```

## Quick Usage

### Python API
```python
from tissue_simulator import TissueSection

tissue = TissueSection(500, 500, 100, cell_radii=(5, 15))
tissue.generate_cells(max_attempts=1000)
tissue.visualize()
tissue.export_to_csv('output.csv')
```

### GUI
```bash
python -m tissue_simulator.gui
```

## Dependencies

- **numpy**: Numerical computations and array operations
- **matplotlib**: 3D visualization and plotting
- **PyQt5**: GUI framework for interactive application

## Use Cases

1. **Histology simulation**: Generate synthetic tissue sections for algorithm testing
2. **Cell packing studies**: Analyze spatial distributions and packing efficiency
3. **Educational tool**: Visualize 3D tissue structure and cell arrangements
4. **Data generation**: Create training data for machine learning models
5. **Parameter exploration**: Test effects of cell sizes and packing parameters

## Future Enhancements

Potential additions:
- Non-spherical cell shapes (ellipsoids)
- Tissue layers with different compositions
- Cell interaction forces (soft-sphere models)
- Time-dependent growth simulations
- Import from microscopy data
- VTK export for advanced 3D rendering
- Multi-threading for large-scale simulations
