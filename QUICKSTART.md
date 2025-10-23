# Quick Start Guide

## Installation

1. Navigate to the package directory:
```bash
cd tissue_simulator
```

2. Install the package:
```bash
pip install -e .
```

Or install dependencies first:
```bash
pip install -r requirements.txt
pip install -e .
```

## Running the GUI

```bash
python -m tissue_simulator.gui
```

Or:
```bash
tissue-simulator
```

## Running Examples

Simple example:
```bash
cd examples
python simple_example.py
```

Multi-cell example:
```bash
python multi_cell_example.py
```

Batch generation:
```bash
python batch_generation.py
```

## Running Tests

```bash
cd tests
python test_tissue_simulator.py
```

Or with pytest (if installed):
```bash
pytest tests/
```

## Python API Quick Example

```python
from tissue_simulator import TissueSection

# Create and generate
tissue = TissueSection(500, 500, 100, cell_radii=(5, 15))
tissue.generate_cells(max_attempts=1000)

# Visualize
tissue.visualize()

# Export
tissue.export_to_csv('my_tissue.csv')
```

## Troubleshooting

### macOS: GUI doesn't appear
```bash
export MPLBACKEND=TkAgg
python -m tissue_simulator.gui
```

### Linux: Missing Qt platform plugin
```bash
sudo apt-get install python3-pyqt5
```

### Windows: DLL load failed
Make sure Visual C++ Redistributable is installed.

## More Information

See GUIDE.md for complete documentation.
