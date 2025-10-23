"""
Custom packing example: Using the SpherePacker directly for fine control.
"""

from tissue_simulator import TissueSection, SpherePacker
import numpy as np

# Define tissue dimensions
bounds = (500, 500, 100)  # height, width, thickness

# Define cell types
cell_radii_config = {
    'cancer': (10, 18),      # Larger cancer cells
    'fibroblast': (7, 12),   # Medium fibroblasts
    'lymphocyte': (4, 7)     # Small lymphocytes
}

# Create custom packer with specific parameters
packer = SpherePacker(
    bounds=bounds,
    cell_radii_config=cell_radii_config,
    min_spacing=0.5,  # Very tight packing
    allow_boundary_cells=True
)

# Pack cells with progress tracking
print("Packing cells...")

def progress_callback(cells_placed, total_attempts):
    if cells_placed % 50 == 0:
        print(f"  Placed {cells_placed} cells (attempts: {total_attempts})")

cells = packer.pack_with_progress(
    max_attempts=3000,
    callback=progress_callback
)

print(f"\nPacking complete: {len(cells)} cells placed")

# Create tissue from packed cells
tissue = TissueSection(
    height=bounds[0],
    width=bounds[1],
    thickness=bounds[2],
    cell_radii=cell_radii_config
)
tissue.cells = cells

# Analyze results
stats = tissue.get_cell_statistics()
print(f"\n=== Results ===")
print(f"Total cells: {stats['total_cells']}")
print(f"Packing fraction: {stats['packing_fraction']:.3f}")

# Calculate cell density
tissue_volume = bounds[0] * bounds[1] * bounds[2]
cell_density = len(cells) / tissue_volume
print(f"Cell density: {cell_density:.6f} cells/μm³")
print(f"Cell density: {cell_density * 1e9:.2f} cells/mm³")

# Visualize
tissue.visualize(elevation=25, azimuth=135)

# Export
tissue.export_to_csv("custom_packed_tissue.csv")
