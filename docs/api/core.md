# Core API: TissueSection, Cell, SpherePacker

Per-module API reference for the core tissue-generation primitives.
For installation, a runnable quick tour, and headline-workflow links,
see [`../quickstart.md`](../quickstart.md).

## Core Concepts

### Tissue Dimensions

The tissue section is a 3D rectangular prism specified in micrometers:

- **height**: Y-axis dimension
- **width**: X-axis dimension
- **thickness**: Z-axis dimension

### Cell Types and Radii

Cells can be specified in two ways:

1. **Uniform cells** - single radius range:

   ```python
   cell_radii=(5, 15)  # all cells between 5-15 um, type "default"
   ```

2. **Multiple cell types** - dict mapping type name to a `(min, max)`
   radius range:

   ```python
   cell_radii={
       "epithelial": (6, 10),
       "stromal": (8, 15),
       "immune": (3, 6),
   }
   ```

A bare tuple is internally normalized to `{"default": (min, max)}`.

### Sphere Packing Algorithm

The package uses Random Sequential Addition (RSA):

1. Randomly select a cell type and a radius from its range.
2. Generate a random position within the tissue bounds.
3. Check for collisions with existing cells (centers must be at least
   `r_i + r_j + min_spacing` apart).
4. Place the cell if the placement is valid; otherwise retry.
5. Continue until the failed-attempt counter exceeds `max_attempts`.

Theoretical packing fraction is ~38% for monodisperse spheres; polydisperse
mixtures can reach 40-50%.

### Boundary Cells

Cells whose surfaces extend beyond the tissue volume are flagged with
`is_boundary=True` and represent the partially sectioned cells seen in
real histological samples.

- `allow_boundary_cells=True`: cell centers must be in bounds, but cell
  surfaces may extend beyond.
- `allow_boundary_cells=False`: the entire cell (center plus radius)
  must fit within bounds.

## API Reference

### `TissueSection`

```python
TissueSection(height, width, thickness, cell_radii, seed=None)
```

Container for a 3D tissue section and its packed cells.

**Parameters:**

- `height` (float): Y-dimension in micrometers.
- `width` (float): X-dimension in micrometers.
- `thickness` (float): Z-dimension in micrometers.
- `cell_radii` (tuple or dict): a `(min, max)` tuple for a single
  cell type, or a dict mapping each cell-type name to a
  `(min, max)` radius range.
- `seed` (int, optional): RNG seed controlling randomness for this
  tissue (both the internal sampler and the `SpherePacker` created by
  `generate_cells`). When `None` (default), the RNG is seeded from
  system entropy. Added in v0.1.2; see the reproducibility note in
  [`../quickstart.md`](../quickstart.md).

**Attributes:**

- `height`, `width`, `thickness` (float)
- `cell_radii` (dict): always stored in dict form, even when constructed
  from a tuple.
- `cells` (list of `Cell`): populated by `generate_cells`.
- `seed` (int or None): the active seed, updated when an explicit
  `seed=` is passed to `generate_cells`.

#### `get_bounds()`

Return tissue dimensions as `(height, width, thickness)`.

#### `generate_cells(max_attempts=1000, min_spacing=0.5, allow_boundary_cells=True, seed=None)`

Generate cells via random sphere packing.

**Parameters:**

- `max_attempts` (int): maximum *consecutive* failed-placement attempts
  before stopping.
- `min_spacing` (float): minimum gap between cell surfaces, in
  micrometers.
- `allow_boundary_cells` (bool): allow cells whose surfaces extend
  beyond tissue bounds.
- `seed` (int, optional): per-call seed override. When provided, the
  internal RNG is reconstructed with this seed and the seed is
  propagated to the underlying `SpherePacker`, making generation
  bit-reproducible. When `None` (default), `self.seed` is used (which
  may itself be `None` for unseeded behavior). Added in v0.1.2.

**Returns:** number of cells placed (int).

#### `get_cell_statistics()`

Compute statistics about the currently packed cells.

**Returns:** dict containing:

- `total_cells`: total cell count.
- `boundary_cells`: number of cells with `is_boundary=True`.
- `interior_cells`: number of cells with `is_boundary=False`.
- `cell_types`: dict mapping each cell-type name to its count.
- `avg_radii`: dict mapping each cell-type name to its mean radius.
- `packing_fraction`: cell volume divided by tissue volume.

When the tissue has no cells, only `total_cells: 0` is returned.

#### `export_to_csv(filename)`

Write every cell to a CSV with the columns:

```text
x, y, z, radius, cell_type, is_boundary
```

#### `from_cells(cells, height=None, width=None, thickness=None, cell_radii=None)` (classmethod)

Build a `TissueSection` from already-positioned `Cell` objects instead of
packing new ones — the inverse of `export_to_csv`. This lets externally
sourced tissue (a measured sample, or a layout exported from another tool)
flow into the analysis and replicate APIs (`get_cell_statistics`,
`SpatialNetworkAnalyzer`, `load_target_statistics_from_tissue`).
Added in v0.1.7.

**Parameters:**

- `cells` (list of `Cell`): pre-positioned cells. Must be non-empty.
- `height`, `width`, `thickness` (float, optional): tissue dimensions in
  micrometers. Any omitted dimension is inferred from the cell-center
  bounding box (`max - min`) on its axis. When every coordinate on an axis
  is equal (e.g. a 2D slice with constant z), that dimension falls back to
  the largest cell diameter so it stays positive and `packing_fraction`
  remains in `(0, 1)`.
- `cell_radii` (dict, optional): `{cell_type: (min, max)}`. When `None`,
  derived per type from the observed radii in `cells`.

**Returns:** a `TissueSection` with `cells` populated and `seed=None`.

**Raises:** `ValueError` if `cells` is empty.

#### `visualize(show_boundary=True, elevation=20, azimuth=45)`

Render the tissue as a 3D matplotlib figure. Boundary cells are drawn
semi-transparent. Requires `matplotlib`.

**Parameters:**

- `show_boundary` (bool): draw the tissue bounding box.
- `elevation` (float): viewing elevation in degrees.
- `azimuth` (float): viewing azimuth in degrees.

Color assignment is deterministic — cell types are sorted alphabetically
before being mapped to the `tab10` palette, so the same set of types always
yields the same colors regardless of the order they appear or the Python
interpreter's `PYTHONHASHSEED`. Added in v0.1.12.

#### `clear_cells()`

Remove all cells from the tissue (does not reset `seed`).

### `Cell`

```python
Cell(center, radius, cell_type="default", is_boundary=False)
```

A single cell.

**Attributes:**

- `center` (numpy.ndarray): 3D position `[x, y, z]`.
- `radius` (float): cell radius in micrometers.
- `cell_type` (str): cell-type name.
- `is_boundary` (bool): whether the cell extends beyond tissue bounds.

#### `intersects(other)`

Return True if this cell overlaps another cell (distance between
centers less than the sum of radii).

#### `is_within_bounds(bounds)`

Return True if the cell fits entirely inside `bounds` (a
`(height, width, thickness)` triple). Used by `SpherePacker` when
`allow_boundary_cells=False`.

#### `intersects_bounds(bounds)`

Return True if the cell's *center* is inside `bounds`, regardless of
radius. Used by `SpherePacker` when `allow_boundary_cells=True`.

### `SpherePacker`

```python
SpherePacker(bounds, cell_radii_config, min_spacing=0.5,
             allow_boundary_cells=True, seed=None)
```

Lower-level packer for callers who want to drive packing without
constructing a full `TissueSection`. `TissueSection.generate_cells`
delegates to this class.

**Parameters:**

- `bounds` (tuple): `(height, width, thickness)` of the tissue.
- `cell_radii_config` (dict): mapping of cell-type name to
  `(min, max)` radii.
- `min_spacing` (float): minimum gap between cell surfaces, in
  micrometers.
- `allow_boundary_cells` (bool): allow cells whose surfaces extend
  beyond tissue bounds.
- `seed` (int, optional): RNG seed. When provided, packing is
  deterministic; when `None` the RNG is seeded from system entropy.
  Added in v0.1.2.

#### `pack(max_attempts=1000)`

Run packing and return the resulting list of `Cell` objects. The loop
exits after `max_attempts` *consecutive* failed placements.

#### `pack_with_progress(max_attempts=1000, callback=None)`

Identical to `pack` but invokes `callback(cells_placed, total_attempts)`
every ten successful placements. Used by the GUI for progress bars.

### `load_tissue_from_csv`

```python
from tissue_simulator import load_tissue_from_csv

tissue = load_tissue_from_csv("cells.csv")
```

Module-level function that loads a `TissueSection` from a CSV written by
`TissueSection.export_to_csv` — the exact inverse of that method, built on
`TissueSection.from_cells`. Added in v0.1.7.

**Parameters:**

- `filepath` (str): path to the CSV.
- `height`, `width`, `thickness` (float, optional): tissue dimensions;
  inferred from the coordinate bounds when omitted (see `from_cells`).
- `default_radius` (float): radius used for rows whose `radius` column is
  missing or blank. Defaults to `10.0`.

**Expected columns:** `x, y, z, radius, cell_type, is_boundary`. Only the
coordinate columns are required: a missing `radius` uses `default_radius`,
a missing `cell_type` becomes `"default"`, and a missing or unparseable
`is_boundary` defaults to `False` (the `"True"`/`"False"` strings written
by `export_to_csv` round-trip exactly).

**Returns:** a populated `TissueSection`.

For full target statistics (interactions plus proportions and density)
straight from a coordinate CSV, see `load_target_statistics_from_coordinates`
in [`replicate-generation.md`](replicate-generation.md).

## Examples

### Direct `SpherePacker` usage

```python
from tissue_simulator import SpherePacker

packer = SpherePacker(
    bounds=(500, 500, 100),
    cell_radii_config={"type_a": (5, 10), "type_b": (3, 7)},
    min_spacing=0.5,
    allow_boundary_cells=True,
    seed=42,
)
cells = packer.pack(max_attempts=2000)
```

### Cell position analysis

```python
import numpy as np
from scipy.spatial.distance import cdist

positions = np.array([cell.center for cell in tissue.cells])
distances = cdist(positions, positions)
np.fill_diagonal(distances, np.inf)
print(f"Mean nearest-neighbor distance: {distances.min(axis=1).mean():.2f} um")
```

### Batch generation for ensemble statistics

```python
import numpy as np
from tissue_simulator import TissueSection

sections = []
for i in range(10):
    tissue = TissueSection(500, 500, 100, cell_radii=(5, 15), seed=i)
    tissue.generate_cells(max_attempts=1500)
    sections.append(tissue)

fractions = [t.get_cell_statistics()["packing_fraction"] for t in sections]
print(f"Mean packing fraction: {np.mean(fractions):.3f} +/- {np.std(fractions):.3f}")
```

For target-statistics-driven batch generation, see
[`replicate-generation.md`](replicate-generation.md) rather than rolling
your own loop.

## See also

- [`../quickstart.md`](../quickstart.md) - five-minute tour, install
  notes, reproducibility guarantees.
- [`slicing.md`](slicing.md) - 2D planar sections of a `TissueSection`.
- [`spatial-analysis.md`](spatial-analysis.md) - network construction
  and contact / radius edge statistics.
- [`replicate-generation.md`](replicate-generation.md) - target-matched
  replicate batches via `ReplicateGenerator`.
- [`graph-coloring.md`](graph-coloring.md) - simulated-annealing cell-type
  assignment.
- [`physicell.md`](physicell.md) - exporting tissues to PhysiCell and
  reading ABM snapshots back in.
- [`convergence.md`](convergence.md) - ADF, Mann-Kendall, rolling CV,
  `find_convergence_time`, `MultiMetricConvergence`
  (`tissue_simulator.convergence`; module-level docstrings until a
  dedicated page lands).
- [`power-analysis.md`](power-analysis.md) - Cohen's d, required
  replicates, power curves
  (`tissue_simulator.power_analysis`; module-level docstrings until a
  dedicated page lands).
- [Changelog](../changelog.md) - release notes, including
  the v0.1.2 `seed=` introduction.
