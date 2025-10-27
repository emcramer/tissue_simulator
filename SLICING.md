# 2D Slicing Module Documentation

## Overview

The slicing module enables extraction of 2D planar sections through 3D tissue at any angle. This simulates histological sectioning and provides tools for visualization and analysis of tissue slices.

## Key Features

- ✅ **Flexible slicing angles**: Horizontal, vertical, or any arbitrary angle
- ✅ **Accurate geometry**: Calculates circular intersections of spheres with planes
- ✅ **2D visualization**: View cell cross-sections in the slice plane
- ✅ **3D context**: Visualize slice plane within the 3D tissue
- ✅ **CSV export**: Export slice data with metadata
- ✅ **Serial sections**: Create multiple parallel slices
- ✅ **Statistics**: Analyze cell distributions in slices

## Quick Start

```python
from tissue_simulator import TissueSection, TissueSlicer

# Generate 3D tissue
tissue = TissueSection(400, 400, 100, cell_radii=(5, 15))
tissue.generate_cells(max_attempts=1000)

# Create horizontal slice
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(z_position=50)

# Visualize
slicer.visualize_slice_2d()
slicer.visualize_slice_in_3d()

# Export
slicer.export_slice_csv('slice_data.csv')
```

## Classes

### TissueSlicer

Main class for slicing operations.

#### Initialization

```python
TissueSlicer(tissue: TissueSection)
```

**Parameters:**
- `tissue`: TissueSection object to slice

#### Methods

##### slice_plane()

```python
slice_plane(point=None, normal=None, angle_x=0.0, angle_y=0.0, z_position=None)
```

Create a 2D slice through the tissue.

**Parameters:**
- `point`: (x, y, z) point on the plane (default: tissue center)
- `normal`: (nx, ny, nz) normal vector to plane
- `angle_x`: Rotation angle around X-axis in degrees
- `angle_y`: Rotation angle around Y-axis in degrees
- `z_position`: Z-position for simplified horizontal slice

**Returns:** List of SliceCell objects

**Examples:**

```python
# Horizontal slice at z=50
slicer.slice_plane(z_position=50)

# Tilted slice (45° around X-axis)
slicer.slice_plane(point=(250, 250, 75), angle_x=45, angle_y=0)

# Custom normal vector
slicer.slice_plane(point=(250, 250, 75), normal=(1, 1, 1))
```

##### visualize_slice_2d()

```python
visualize_slice_2d(show_radii=True, figsize=(10, 10), title=None)
```

Visualize the 2D slice showing cell cross-sections.

**Parameters:**
- `show_radii`: Show cells with intersection radii (vs. full 3D radii)
- `figsize`: Figure size tuple
- `title`: Plot title

##### visualize_slice_in_3d()

```python
visualize_slice_in_3d(show_plane=True, plane_alpha=0.3, plane_size=None)
```

Visualize the slice plane within 3D tissue context.

**Parameters:**
- `show_plane`: Whether to show the slice plane
- `plane_alpha`: Transparency of plane (0-1)
- `plane_size`: Size of plane to display

##### export_slice_csv()

```python
export_slice_csv(filename, include_3d=True)
```

Export slice data to CSV file.

**Parameters:**
- `filename`: Output file path
- `include_3d`: Include original 3D coordinates

**CSV Format (include_3d=True):**
```csv
x_2d, y_2d, intersection_radius, x_3d, y_3d, z_3d, radius_3d, cell_type, is_boundary, distance_from_plane
```

**CSV Format (include_3d=False):**
```csv
x_2d, y_2d, intersection_radius, cell_type, distance_from_plane
```

##### get_slice_statistics()

```python
get_slice_statistics()
```

Calculate statistics about the slice.

**Returns:** Dictionary containing:
- `num_cells`: Number of cells in slice
- `plane_point`: Point on plane [x, y, z]
- `plane_normal`: Normal vector [nx, ny, nz]
- `cell_types`: Count per cell type
- `avg_intersection_radii`: Average intersection radius per type
- `mean_distance_from_plane`: Mean distance of cell centers from plane
- `max_distance_from_plane`: Maximum distance

### SliceCell

Dataclass representing a cell captured in a slice.

**Attributes:**
- `center_3d`: Original 3D center position (numpy array)
- `center_2d`: Position in 2D slice coordinates (numpy array)
- `radius`: Original 3D cell radius
- `cell_type`: Cell classification string
- `is_boundary`: Whether cell was boundary cell in 3D
- `distance_from_plane`: Perpendicular distance from slice plane
- `intersection_radius`: Radius of circular intersection with plane

## Helper Functions

### create_standard_slices()

```python
create_standard_slices(tissue, num_slices=5)
```

Create multiple evenly-spaced parallel slices through tissue.

**Parameters:**
- `tissue`: TissueSection to slice
- `num_slices`: Number of slices to create

**Returns:** List of TissueSlicer objects with computed slices

**Example:**
```python
slicers = create_standard_slices(tissue, num_slices=5)
for i, slicer in enumerate(slicers):
    print(f"Slice {i+1}: {len(slicer.slice_cells)} cells")
    slicer.export_slice_csv(f'slice_{i+1}.csv')
```

## Usage Examples

### Example 1: Simple Horizontal Slice

```python
from tissue_simulator import TissueSection, TissueSlicer

# Generate tissue
tissue = TissueSection(400, 400, 100, cell_radii=(5, 15))
tissue.generate_cells(max_attempts=1000)

# Create slice at middle
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(z_position=50)

print(f"Captured {len(slice_cells)} cells")

# Visualize
slicer.visualize_slice_2d()

# Export
slicer.export_slice_csv('horizontal_slice.csv')
```

### Example 2: Angled Slice

```python
# Create slice tilted 45° around X-axis
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(
    point=(200, 200, 50),  # Center point
    angle_x=45,
    angle_y=0
)

# Visualize in 3D to see the angle
slicer.visualize_slice_in_3d(show_plane=True)

# Get statistics
stats = slicer.get_slice_statistics()
print(f"Plane normal: {stats['plane_normal']}")
print(f"Cell types: {stats['cell_types']}")
```

### Example 3: Serial Sections

```python
from tissue_simulator import create_standard_slices

# Create 5 parallel slices
slicers = create_standard_slices(tissue, num_slices=5)

# Process each slice
for i, slicer in enumerate(slicers, 1):
    stats = slicer.get_slice_statistics()
    z_pos = stats['plane_point'][2]
    
    print(f"Slice {i} at z={z_pos:.1f}: {stats['num_cells']} cells")
    
    # Export each slice
    slicer.export_slice_csv(f'slice_{i}.csv')
```

### Example 4: Custom Analysis

```python
# Create slice
slicer = TissueSlicer(tissue)
slicer.slice_plane(z_position=50)

# Analyze slice cells
for slice_cell in slicer.slice_cells:
    print(f"Cell type: {slice_cell.cell_type}")
    print(f"  2D position: {slice_cell.center_2d}")
    print(f"  3D position: {slice_cell.center_3d}")
    print(f"  Intersection radius: {slice_cell.intersection_radius:.2f}")
    print(f"  Distance from plane: {slice_cell.distance_from_plane:.2f}")
    print()
```

## Coordinate Systems

### 3D Tissue Coordinates

- **X-axis**: Width (0 to tissue.width)
- **Y-axis**: Height (0 to tissue.height)
- **Z-axis**: Thickness (0 to tissue.thickness)

### 2D Slice Coordinates

The slice plane has its own 2D coordinate system:
- **U-axis**: First basis vector in the plane
- **V-axis**: Second basis vector in the plane
- Origin at the specified plane point

The basis vectors are automatically calculated to be:
1. Orthogonal to the plane normal
2. Orthogonal to each other
3. Unit length (normalized)

## Geometry Details

### Sphere-Plane Intersection

When a sphere (cell) intersects a plane:

1. **Distance calculation**: Distance from cell center to plane
   ```
   d = |dot(center - point, normal)|
   ```

2. **Intersection test**: Cell intersects if `d < radius`

3. **Intersection radius**: Radius of circular cross-section
   ```
   r_intersection = sqrt(radius² - d²)
   ```

4. **2D projection**: Cell center projected onto plane
   ```
   projected = center - d * normal
   ```

### Rotation Angles

- `angle_x`: Rotation around X-axis (pitch)
  - 0°: Horizontal plane (XY)
  - 90°: Vertical plane (YZ)
  
- `angle_y`: Rotation around Y-axis (yaw)
  - 0°: No rotation around Y
  - 90°: Vertical plane (XZ)

Combined rotations: Y-rotation applied first, then X-rotation

## Performance Considerations

- Slicing is O(n) where n is number of cells in tissue
- Each cell is tested for intersection with plane
- 2D visualization is faster than 3D
- Serial sections reuse tissue data efficiently

## Tips and Best Practices

1. **Horizontal slices**: Use `z_position` parameter for simplicity
   ```python
   slicer.slice_plane(z_position=50)
   ```

2. **Angled slices**: Use `angle_x` and `angle_y` for intuitive control
   ```python
   slicer.slice_plane(angle_x=45, angle_y=30)
   ```

3. **Custom orientations**: Use `normal` vector for precise control
   ```python
   slicer.slice_plane(normal=(1, 1, 1))  # Diagonal
   ```

4. **Visualizing angles**: Always use `visualize_slice_in_3d()` to verify slice orientation

5. **Export early**: Export slice data before creating new slices
   ```python
   slicer.slice_plane(z_position=50)
   slicer.export_slice_csv('slice_50.csv')  # Export before next slice
   ```

6. **Serial sections**: Use `create_standard_slices()` for evenly-spaced slices

## Common Use Cases

### Histology Simulation
```python
# Multiple parallel sections like in histology
slicers = create_standard_slices(tissue, num_slices=10)
for i, slicer in enumerate(slicers):
    slicer.visualize_slice_2d(title=f"Section {i+1}")
```

### Arbitrary Sectioning
```python
# Oblique section for examining tissue architecture
slicer = TissueSlicer(tissue)
slicer.slice_plane(angle_x=30, angle_y=45)
slicer.visualize_slice_in_3d()
```

### Cell Counting
```python
# Count cells in specific region via slicing
slicer = TissueSlicer(tissue)
slicer.slice_plane(z_position=50)
stats = slicer.get_slice_statistics()
print(f"Cell density in slice: {stats['num_cells']} cells")
```

### 2D Image Analysis Training
```python
# Generate 2D slices for training segmentation algorithms
for i in range(20):
    z = tissue.thickness * (i + 1) / 21
    slicer = TissueSlicer(tissue)
    slicer.slice_plane(z_position=z)
    slicer.export_slice_csv(f'training_slice_{i:02d}.csv')
```

## Troubleshooting

**Q: No cells in slice?**
- Check slice position is within tissue bounds
- Verify tissue has cells (`len(tissue.cells) > 0`)
- Try different z-positions

**Q: Unexpected slice orientation?**
- Use `visualize_slice_in_3d()` to verify
- Check normal vector is normalized
- Try simpler interface (`z_position` or `angle_x/y`)

**Q: 2D coordinates seem wrong?**
- Coordinates are in plane basis (U, V), not (X, Y)
- Use `center_3d` from SliceCell for original coordinates
- Check plane basis vectors with `slicer.slice_basis_u/v`

## See Also

- [Main Guide](GUIDE.md) - Complete package documentation
- [Examples](examples/) - Working code samples
- [API Reference](GUIDE.md#api-reference) - Full API documentation
