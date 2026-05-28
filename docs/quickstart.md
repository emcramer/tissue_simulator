# Quick Start

## What is `tissue_simulator`?

`tissue_simulator` is a Python package for generating 3D simulated biological
tissue sections with random sphere packing, slicing them into 2D histological
sections, and analyzing the resulting cell-cell spatial networks. It also
ships a replicate generator that matches target spatial statistics, a graph
coloring routine that assigns cell types via simulated annealing, convergence
and power-analysis utilities for ABM trajectories, a bridge to PhysiCell,
and an MCP server so an LLM can drive the whole pipeline. For the full
project pitch see [Home](index.md).

## Install

```bash
# Eventual PyPI install (the package is not yet published; use "from source" for now).
pip install tissue_simulator

# From source, editable install (recommended today):
pip install -e .

# With the optional MCP server (v0.1.4+ moved this out of requirements.txt):
pip install -e ".[mcp]"

# For contributors (pytest, pytest-cov, build, twine):
pip install -e ".[dev]"
# Or use requirements-dev.txt, which additionally includes mcp:
pip install -r requirements-dev.txt
```

Verify the install:

```bash
python verify_installation.py
```

## Five-minute tour

The snippet below generates a small 3-cell-type tissue, takes a 2D slice,
builds a spatial network, prints a couple of stats, and writes CSVs.
Copy-paste runnable; runs in well under 30 seconds.

```python
from tissue_simulator import (
    TissueSection,
    TissueSlicer,
    SpatialNetworkAnalyzer,
)

# step 1: generate a small 3D tissue with three cell types (seeded for reproducibility)
tissue = TissueSection(
    height=200,
    width=200,
    thickness=50,
    cell_radii={
        "tumor": (8.0, 12.0),
        "immune": (5.0, 8.0),
        "stroma": (6.0, 10.0),
    },
    seed=20260513,
)
n_cells = tissue.generate_cells(max_attempts=1500, min_spacing=0.5)
stats = tissue.get_cell_statistics()
print(f"Generated {n_cells} cells (packing fraction = {stats['packing_fraction']:.3f})")
by_type = {str(k): v for k, v in stats["cell_types"].items()}
print(f"  by type: {by_type}")

# step 2: extract a horizontal 2D slice at z = 25 um
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(z_position=25.0)
print(f"Slice captured {len(slice_cells)} cell cross-sections at z=25 um")

# step 3: build a radius-mode spatial network on the 3D tissue
analyzer = SpatialNetworkAnalyzer()
analyzer.build_network_from_tissue(tissue, mode="radius", radius=25.0)

# step 4: print a couple of network stats
g = analyzer.compute_global_statistics()
print(f"Network: {g.total_nodes} nodes, {g.total_edges} edges")
print(f"  avg degree = {g.avg_degree:.2f}, clustering = {g.avg_clustering:.3f}")

# step 5: save CSV outputs (3D tissue, 2D slice, and 3 network-statistics files)
tissue.export_to_csv("tissue.csv")
slicer.export_slice_csv("slice.csv")
analyzer.export_statistics_csv("network_stats")
print("Wrote: tissue.csv, slice.csv, network_stats_{global,cell_types,interactions}.csv")
```

Expected output (line-for-line shape; exact counts depend on packing run):

```text
Generated <N> cells (packing fraction = 0.xxx)
  by type: {'tumor': ..., 'immune': ..., 'stroma': ...}
Slice captured <M> cell cross-sections at z=25 um
Network: <N> nodes, <E> edges
  avg degree = ..., clustering = ...
Wrote: tissue.csv, slice.csv, network_stats_{global,cell_types,interactions}.csv
```

`export_statistics_csv("network_stats")` writes three sibling files:
`network_stats_global.csv`, `network_stats_cell_types.csv`, and
`network_stats_interactions.csv`.

## Headline workflows

- **Generate replicates matching target spatial statistics** -> see
  [`api/replicate-generation.md`](api/replicate-generation.md) and the
  end-to-end demo in
  [examples/result1_replicate_demonstration.py](https://github.com/emcramer/tissue_simulator/blob/main/examples/result1_replicate_demonstration.py).
- **Assign cell types via simulated annealing on a graph** -> see
  [`api/graph-coloring.md`](api/graph-coloring.md). The `GraphColorizer` and
  evaluation helpers (`js_divergence`, `evaluate_graph_coloring`) are the
  entry points.
- **Bridge to a PhysiCell ABM** -> use `PhysiCellExporter` /
  `export_to_physicell` to seed PhysiCell from a `TissueSection`, then
  `PhysiCellReader` / `read_physicell_output` and `stats_to_target_statistics`
  to feed snapshots back into the `ReplicateGenerator`. The reference
  implementations live in
  [tissue_simulator/physicell_export.py](https://github.com/emcramer/tissue_simulator/blob/main/tissue_simulator/physicell_export.py)
  and
  [tissue_simulator/physicell_reader.py](https://github.com/emcramer/tissue_simulator/blob/main/tissue_simulator/physicell_reader.py).
- **Run convergence and power analysis on ABM trajectories** ->
  `tissue_simulator.convergence` (`adf_test`, `mann_kendall_test`,
  `rolling_cv`, `find_convergence_time`, `MultiMetricConvergence`) and
  `tissue_simulator.power_analysis` (`cohens_d`, `required_replicates`,
  `power_curve`, `compare_initialization_variance`). No dedicated guide
  yet; the modules are docstring-driven.
- **Drive the package from an LLM (MCP)** -> see
  [`guides/mcp.md`](guides/mcp.md) for a 5-minute walk-through and
  [`api/mcp.md`](api/mcp.md) for the full tool reference.
- **Interactive GUI** -> `python -m tissue_simulator.gui` or the installed
  console script `tissue-simulator`.

## Reproducibility

Since v0.1.2 every random draw in `TissueSection`, `SpherePacker`, and
`ReplicateGenerator` is routed through an explicit seed. Passing `seed=N`
to `TissueSection(...)` (or `ReplicateGenerator(..., seed=N)`) yields
bit-identical tissues across Python processes, independent of
`PYTHONHASHSEED`. `ReplicateGenerator` derives a deterministic per-replicate
seed from `(seed, replicate_id)` without touching the global NumPy RNG, so
mixing replicate generation with other randomness in the same process is
safe.

Note on divergence semantics: as of v0.1.2, `ReplicateStatistics.divergence_score`
is `nan` (not `0`) when both the target and the measured value are zero for
a given cell-type pair; the aggregate uses `np.nanmean`, so empty-signal
pairs no longer masquerade as perfect matches. See
[Changelog](changelog.md) for the v0.1.2 entry.

## Troubleshooting

### macOS: GUI doesn't appear

```bash
export MPLBACKEND=TkAgg
python -m tissue_simulator.gui
```

### Linux: missing Qt platform plugin

```bash
sudo apt-get install python3-pyqt5
```

### Windows: DLL load failed

Install the Microsoft Visual C++ Redistributable.

### `ImportError: No module named statsmodels`

`statsmodels` is required by the `tissue_simulator.convergence` and
`tissue_simulator.power_analysis` modules (since v0.1.1). Install with:

```bash
pip install statsmodels
```

It is pinned in `requirements.txt`, so this usually only bites users who
installed without the requirements file.

## Next steps

- [`api/core.md`](api/core.md) - `TissueSection`, `Cell`, `SpherePacker` reference
- [`api/slicing.md`](api/slicing.md) - 2D slicing and serial sections
- [`api/spatial-analysis.md`](api/spatial-analysis.md) - network construction and statistics
- [`api/graph-coloring.md`](api/graph-coloring.md) - simulated-annealing cell-type assignment
- [`api/replicate-generation.md`](api/replicate-generation.md) - target-matched replicate batches
- [`api/mcp.md`](api/mcp.md) - MCP tool reference; pair with [`guides/mcp.md`](guides/mcp.md)
- [`guides/complete-workflow.md`](guides/complete-workflow.md) - end-to-end tutorial
- [Changelog](changelog.md) - release history and migration notes
- [examples/](https://github.com/emcramer/tissue_simulator/tree/main/examples) - runnable scripts for every major feature
