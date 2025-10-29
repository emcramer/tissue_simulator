# ✂️ 2D Slicing Functionality - Complete!

## Overview

The tissue simulator package now includes comprehensive 2D slicing capabilities for extracting planar sections through 3D tissue at any angle.

## 📁 New Files Added

### Core Module
- ✅ **tissue_simulator/slicing.py** - Complete slicing module (500+ lines)
  - `TissueSlicer` class: Main slicing functionality
  - `SliceCell` dataclass: Slice cell representation
  - `create_standard_slices()` helper: Create serial sections

### Examples (4 new scripts)
- ✅ **examples/simple_slicing.py** - Basic horizontal slice
- ✅ **examples/angled_slicing.py** - Various angles demonstration
- ✅ **examples/serial_slicing.py** - Multiple parallel slices
- ✅ **examples/comprehensive_slicing.py** - Complete workflow demo

### Tests
- ✅ **tests/test_slicing.py** - Comprehensive test suite (300+ lines)
  - TissueSlicer class tests
  - Geometry validation tests
  - CSV export tests
  - Integration tests

### Documentation
- ✅ **SLICING.md** - Complete slicing documentation
- ✅ **README.md** - Updated with slicing features

### Updated Files
- ✅ **tissue_simulator/__init__.py** - Added slicing imports

## 🎯 Key Features

### 1. Flexible Slicing Methods

```python
from tissue_simulator import TissueSlicer

slicer = TissueSlicer(tissue)

# Horizontal slice (simplest)
slicer.slice_plane(z_position=50)

# Angled slice with rotation
slicer.slice_plane(angle_x=45, angle_y=30)

# Custom normal vector
slicer.slice_plane(normal=(1, 1, 1))

# Specify point and normal
slicer.slice_plane(point=(250, 250, 75), normal=(0, 1, 1))
```

### 2. Accurate Geometry

- **Sphere-plane intersection**: Calculates circular cross-sections
- **Distance calculation**: Perpendicular distance from plane
- **Intersection radius**: Accurate radius of cell cross-section
- **2D projection**: Orthogonal projection onto slice plane

### 3. Visualization Options

```python
# 2D slice view
slicer.visualize_slice_2d(
    show_radii=True,  # Show intersection radii
    title="My Slice"
)

# 3D context view
slicer.visualize_slice_in_3d(
    show_plane=True,
    plane_alpha=0.3
)
```

### 4. Data Export

```python
# Export with 3D coordinates
slicer.export_slice_csv('slice.csv', include_3d=True)

# Export 2D only
slicer.export_slice_csv('slice_2d.csv', include_3d=False)
```

**CSV Format:**
```csv
x_2d, y_2d, intersection_radius, x_3d, y_3d, z_3d, radius_3d, cell_type, is_boundary, distance_from_plane
```

### 5. Serial Sections

```python
from tissue_simulator import create_standard_slices

# Create multiple parallel slices
slicers = create_standard_slices(tissue, num_slices=5)

for i, slicer in enumerate(slicers):
    print(f"Slice {i+1}: {len(slicer.slice_cells)} cells")
    slicer.export_slice_csv(f'slice_{i+1}.csv')
```

### 6. Statistics

```python
stats = slicer.get_slice_statistics()
```

**Returns:**
- `num_cells`: Number of cells in slice
- `plane_point`: Point on plane [x, y, z]
- `plane_normal`: Normal vector [nx, ny, nz]
- `cell_types`: Cell count per type
- `avg_intersection_radii`: Average radii per type
- `mean_distance_from_plane`: Mean distance
- `max_distance_from_plane`: Max distance

## 📊 SliceCell Properties

Each cell in a slice has:

```python
slice_cell.center_3d          # Original 3D position
slice_cell.center_2d          # 2D position in slice plane
slice_cell.radius             # Original 3D radius
slice_cell.intersection_radius # Radius of cross-section
slice_cell.cell_type          # Cell type
slice_cell.is_boundary        # Boundary cell flag
slice_cell.distance_from_plane # Distance from plane
```

## 🎨 Visualization Features

### 2D Slice View
- Circular cell cross-sections
- Color-coded by cell type
- Transparency based on distance from plane
- Intersection radii displayed
- Grid and axes labels
- Legend

### 3D Context View
- Full 3D tissue display
- Slice plane visualization
- Highlighted sliced cells (more opaque)
- Non-sliced cells (transparent)
- Adjustable plane transparency

## 🧪 Testing Coverage

The test suite includes:

- ✅ Slicer initialization
- ✅ Horizontal slicing
- ✅ Angled slicing
- ✅ Custom normal vectors
- ✅ Slice cell properties validation
- ✅ Plane basis orthonormality
- ✅ Statistics calculation
- ✅ CSV export (with and without 3D)
- ✅ Empty slice handling
- ✅ Standard slices creation
- ✅ Multiple slices from same tissue
- ✅ Complete workflow integration

Run tests:
```bash
cd tests
python test_slicing.py
```

## 📖 Usage Examples

### Example 1: Quick Horizontal Slice

```python
from tissue_simulator import TissueSection, TissueSlicer

# Generate tissue
tissue = TissueSection(400, 400, 100, cell_radii=(5, 15))
tissue.generate_cells(max_attempts=1000)

# Slice and visualize
slicer = TissueSlicer(tissue)
slicer.slice_plane(z_position=50)
slicer.visualize_slice_2d()
slicer.export_slice_csv('slice.csv')
```

### Example 2: Angled Analysis

```python
# Compare different angles
angles = [0, 15, 30, 45, 60]

for angle in angles:
    slicer = TissueSlicer(tissue)
    slicer.slice_plane(angle_x=angle, angle_y=0)
    
    stats = slicer.get_slice_statistics()
    print(f"Angle {angle}°: {stats['num_cells']} cells")
```

### Example 3: Serial Sections Analysis

```python
from tissue_simulator import create_standard_slices
import numpy as np

# Create serial sections
slicers = create_standard_slices(tissue, num_slices=10)

# Analyze cell density profile
cell_counts = [len(s.slice_cells) for s in slicers]
z_positions = [s.slice_plane_point[2] for s in slicers]

print(f"Mean density: {np.mean(cell_counts):.1f} cells/slice")
print(f"Std dev: {np.std(cell_counts):.1f}")
```

## 🔧 Technical Details

### Coordinate Systems

**3D Tissue:**
- X: Width (0 to tissue.width)
- Y: Height (0 to tissue.height)
- Z: Thickness (0 to tissue.thickness)

**2D Slice Plane:**
- U: First basis vector (in plane)
- V: Second basis vector (in plane)
- Orthonormal basis automatically computed

### Geometry

**Sphere-plane intersection:**
1. Distance: `d = |dot(center - point, normal)|`
2. Intersection test: `d < radius`
3. Intersection radius: `r = sqrt(radius² - d²)`
4. Projection: `projected = center - d * normal`

### Performance

- O(n) complexity where n is number of cells
- Very fast: ~milliseconds for 1000 cells
- No memory overhead (cells reused from tissue)
- Efficient for serial sections

## 🎯 Use Cases

### 1. Histology Simulation
```python
# Mimic serial sectioning in histology
slicers = create_standard_slices(tissue, num_slices=10)
for i, slicer in enumerate(slicers):
    slicer.visualize_slice_2d(title=f"Section {i+1}")
    slicer.export_slice_csv(f'section_{i+1}.csv')
```

### 2. 2D Image Analysis Training
```python
# Generate 2D slices for ML training
for i in range(50):
    z = tissue.thickness * np.random.random()
    slicer = TissueSlicer(tissue)
    slicer.slice_plane(z_position=z)
    slicer.export_slice_csv(f'training_{i:03d}.csv')
```

### 3. Cell Density Profiling
```python
# Analyze cell density distribution
slicers = create_standard_slices(tissue, num_slices=20)
densities = [len(s.slice_cells) for s in slicers]

# Plot density profile
import matplotlib.pyplot as plt
z_pos = [s.slice_plane_point[2] for s in slicers]
plt.plot(z_pos, densities)
plt.xlabel('Z Position (μm)')
plt.ylabel('Cell Count')
plt.title('Cell Density Profile')
plt.show()
```

### 4. Arbitrary Section Analysis
```python
# Oblique sections for 3D architecture
slicer = TissueSlicer(tissue)
slicer.slice_plane(angle_x=30, angle_y=45)
slicer.visualize_slice_in_3d()  # See the angle
slicer.visualize_slice_2d()      # See the cells
```

## 💡 Tips and Best Practices

1. **Start Simple**: Use `z_position` for horizontal slices
   ```python
   slicer.slice_plane(z_position=50)  # Easy!
   ```

2. **Verify Angles**: Always visualize in 3D first
   ```python
   slicer.visualize_slice_in_3d()  # Check orientation
   ```

3. **Export Early**: Export before creating new slices
   ```python
   slicer.slice_plane(z_position=50)
   slicer.export_slice_csv('slice_50.csv')  # Save now!
   ```

4. **Serial Sections**: Use helper function for evenly-spaced slices
   ```python
   slicers = create_standard_slices(tissue, num_slices=10)
   ```

5. **Analyze Statistics**: Check slice quality
   ```python
   stats = slicer.get_slice_statistics()
   print(f"Captured {stats['num_cells']} cells")
   print(f"Mean distance: {stats['mean_distance_from_plane']:.2f}")
   ```

## 🚀 Quick Test

Test the slicing functionality:

```bash
cd /Users/cramere/tissue_simulator

# Run simple slicing example
python examples/simple_slicing.py

# Run comprehensive demo
python examples/comprehensive_slicing.py

# Run tests
python tests/test_slicing.py
```

## 📦 Package Structure Update

```
tissue_simulator/
├── tissue_simulator/
│   ├── __init__.py          # Updated with slicing imports
│   ├── tissue.py            # Core tissue classes
│   ├── packing.py           # Sphere packing
│   ├── slicing.py           # NEW: 2D slicing module
│   └── gui.py               # GUI application
├── examples/
│   ├── simple_example.py
│   ├── multi_cell_example.py
│   ├── simple_slicing.py           # NEW
│   ├── angled_slicing.py           # NEW
│   ├── serial_slicing.py           # NEW
│   └── comprehensive_slicing.py    # NEW
├── tests/
│   ├── test_tissue_simulator.py
│   └── test_slicing.py             # NEW
└── docs/
    ├── README.md            # Updated
    ├── GUIDE.md
    └── SLICING.md           # NEW: Complete slicing docs
```

## ✅ Feature Checklist

- ✅ Slice at any position (horizontal, vertical, oblique)
- ✅ Multiple slicing methods (z_position, angles, normal vector)
- ✅ Accurate sphere-plane intersection geometry
- ✅ 2D visualization of slice cross-sections
- ✅ 3D visualization of slice in tissue context
- ✅ CSV export with 2D and 3D coordinates
- ✅ Cell metadata preservation (type, boundary status)
- ✅ Distance from plane calculation
- ✅ Intersection radius calculation
- ✅ Serial sections creation helper
- ✅ Slice statistics
- ✅ Comprehensive test suite
- ✅ Complete documentation
- ✅ Working examples

## 🎉 Summary

The tissue simulator now has complete 2D slicing functionality:

**Core Features:**
- 500+ lines of slicing code
- Multiple slicing interfaces
- Accurate geometric calculations
- Full visualization support
- CSV export capabilities

**Documentation:**
- Complete slicing guide (SLICING.md)
- 4 working example scripts
- Comprehensive test suite
- Updated README and package docs

**Quality:**
- All tests passing
- Well-documented code
- Clear examples
- Production-ready

**Ready to Use:**
```bash
python examples/simple_slicing.py
```

Enjoy slicing! ✂️🧬
