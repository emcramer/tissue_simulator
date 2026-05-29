# Tissue Simulator

**Generate 3D simulated biological tissue sections, slice them, analyze cell-cell spatial networks, and feed the results into agent-based models — all from one Python package.**

`tissue_simulator` packs the loop from synthetic histology to ABM initialization
into a small set of composable building blocks, with an MCP server so an LLM
can drive the whole pipeline conversationally.

## What it does

- **3D sphere packing** — Random Sequential Addition placement of single- or
  multi-cell-type spheres in a configurable tissue volume, with collision
  detection and optional boundary cells.
- **2D slicing** — Extract planar histological sections at any angle, with
  per-cell intersection geometry and CSV export.
- **Network-based spatial analysis** — Build NetworkX graphs in `contact` or
  `radius` mode and compute global, per-cell-type, and pairwise interaction
  statistics.
- **Graph coloring for cell-type assignment** — Use simulated annealing to
  assign cell types to a network that match target node counts, edge counts,
  and neighbor distributions.
- **Replicate generation** — Iteratively tune parameters to produce batches of
  tissues that match a target spatial-statistics profile (e.g. from a real
  reference tissue).
- **ABM bridge + MCP server** — Export tissues into PhysiCell, read snapshots
  back, run convergence and power-analysis on trajectories, and drive any of
  the above from an LLM via the Model Context Protocol.

## Install

```bash
# From source, editable (recommended today; not yet on PyPI):
pip install -e .

# With the optional MCP server:
pip install -e ".[mcp]"

# For contributors (pytest, build, twine):
pip install -e ".[dev]"
```

Verify with `python verify_installation.py`.

## Five-minute tour

```python
from tissue_simulator import (
    TissueSection,
    TissueSlicer,
    SpatialNetworkAnalyzer,
)

# 1. Generate a small 3D tissue with three cell types (seeded for reproducibility)
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

# 2. Take a horizontal 2D slice at z = 25 um
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(z_position=25.0)

# 3. Build a radius-mode spatial network
analyzer = SpatialNetworkAnalyzer()
analyzer.build_network_from_tissue(tissue, mode="radius", radius=25.0)
g = analyzer.compute_global_statistics()
print(f"{n_cells} cells, network: {g.total_nodes} nodes / {g.total_edges} edges")
```

A copy-paste-runnable extension (with CSV export and slice counts) lives in
the [Quickstart](quickstart.md).

## Where to next

- **[📊 Slide deck](slides/index.md)** — a 13-slide visual tour of the package
  end-to-end (code + matplotlib output side by side).
- [Quickstart](quickstart.md) — install, the five-minute tour, headline
  workflows, reproducibility notes, and troubleshooting.
- [Complete Workflow guide](guides/complete-workflow.md) — end-to-end tutorial
  from sphere packing through evaluated cell-type assignment.
- [API Guide](api/core.md) — narrative reference for each module: core,
  slicing, spatial analysis, graph coloring, replicate generation,
  convergence, power analysis, the PhysiCell bridge, and the MCP server.
- [API Reference](reference/tissue.md) — auto-generated mkdocstrings pages for
  every public module.
- [Changelog](changelog.md) — release history and migration notes.

## Community & contributions

`tissue_simulator` is developed in the open on
[GitHub](https://github.com/emcramer/tissue_simulator). Bug reports, feature
requests, and pull requests are welcome on the
[issue tracker](https://github.com/emcramer/tissue_simulator/issues); design
notes, deeper background, and discussion live in the
[project wiki](https://github.com/emcramer/tissue_simulator/wiki).
