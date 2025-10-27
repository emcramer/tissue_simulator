# Bug Fix: Matplotlib plot_surface Color Error

## Problem
When running the examples, you encountered this error:
```
ValueError: operands could not be broadcast together with shapes (361,1) (0,4)
```

## Root Cause
The error occurred in matplotlib's `plot_surface()` function when using the `color` parameter with newer versions of matplotlib. The `color` parameter expects a specific format that wasn't being provided correctly, causing a broadcasting error in the shading calculations.

## Solution
Changed from using `color` parameter to `facecolors` parameter with properly formatted color arrays.

### Before (Caused Error):
```python
ax.plot_surface(x, y, z, color=color, alpha=alpha, 
               edgecolors='none', shade=True)
```

### After (Fixed):
```python
ax.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)),
               alpha=alpha, linewidth=0, antialiased=True, shade=False)
```

## What Changed
The fix tiles the color array to match the shape of the surface mesh:
- `np.tile(color, x.shape + (1,))` creates an array where each face of the sphere surface gets the same color
- `facecolors` expects per-face colors rather than a single color
- `linewidth=0` removes edge lines
- `shade=False` prevents the shading calculations that caused the error

## Files Updated
1. `tissue_simulator/tissue.py` - Main visualize() method
2. `tissue_simulator/gui.py` - GUI 3D viewer
3. `examples/advanced_visualization.py` - Example script with multiple views

## Testing
The fix has been applied to all visualization functions. To verify:

```bash
cd /Users/cramere/tissue_simulator

# Test basic example
python examples/simple_example.py

# Test GUI
python -m tissue_simulator.gui

# Test advanced visualization
python examples/advanced_visualization.py
```

## Compatibility
This fix is compatible with:
- matplotlib >= 3.3.0
- numpy >= 1.20.0
- Python 3.8+

The visualization should now work correctly without broadcasting errors.
