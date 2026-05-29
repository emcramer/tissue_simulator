# Graph Coloring Integration Guide

## Overview

The tissue simulator now includes powerful graph-based cell type assignment capabilities using simulated annealing. This integration allows you to:

1. Generate 3D tissue structures
2. Extract 2D slices
3. Build network graphs from cell spatial relationships
4. Assign cell types based on target network statistics
5. Visualize and evaluate results

## When to use what

`GraphColorizer` and `ReplicateGenerator` solve different problems. Pick the
right entry point first, then tune.

- **`GraphColorizer` — label-only assignment on a fixed graph.** Given a
  graph (positions already fixed) and a target adjacency structure
  (`node_counts`, `edge_counts`, `neighbor_dist`), it assigns cell-type
  labels to existing nodes via simulated annealing to match the target.
  It does not move or repack cells. Use it when you already have a packed
  tissue (or any spatial network) and only need to decide which node gets
  which cell-type label.

- **`ReplicateGenerator` — unstructured replicate generation by repacking.**
  Given a `TargetStatistics` (contact statistics + proportions + optional
  density), it generates new tissues by repacking positions and biasing
  cell-type sampling toward the target proportions. It does NOT run
  simulated annealing and does NOT internally compose `GraphColorizer`
  today. For strongly structured multi-type targets (e.g. a tumor disc
  with a fibroblast ring and a CD8 annulus), its single-stage
  `rng.choice + radius adjust` strategy will not converge — there is no
  optimizer that pushes labels into a particular spatial arrangement.

- **Canonical path for structured multi-type targets: two-stage workflow.**
  Use `TissueWorkflow` to compose the two stages — `generate_cells` for
  positions, then `assign_cell_types` (which wraps `GraphColorizer`) for
  labels. The convenience entry points `TissueNetworkWorkflow.run_complete_workflow`
  and `quick_workflow` both forward a `seed=` to the colorize step.

Note: folding `GraphColorizer` into `ReplicateGenerator` so that the
replicate generator itself can converge on structured targets is
intentionally out of scope here; the present answer is to document and
use the two-stage `TissueWorkflow` path.

## Quick Start

```python
from tissue_simulator import (
    TissueSection, 
    SpherePacker,
    quick_workflow
)

# Create tissue
tissue = TissueSection(height=200, width=200, thickness=50)
packer = SpherePacker(tissue)
packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1000)

# Run complete workflow
workflow = quick_workflow(
    tissue=tissue,
    cell_types=['cancer', 'immune', 'stroma'],
    target_stats_file="target_statistics.csv",
    network_radius=50.0,
    output_dir="results"
)
```

## Complete Workflow

### 1. Create a 3D Tissue

```python
from tissue_simulator import TissueSection, SpherePacker

# Define tissue dimensions
tissue = TissueSection(
    height=200,      # micrometers
    width=200,       # micrometers
    thickness=50     # micrometers
)

# Define cell types (initial assignment, will be reassigned)
cell_types = {
    'placeholder': [8, 12]  # min and max radius in micrometers
}

# Generate cells
packer = SpherePacker(tissue)
num_cells = packer.pack_cells(
    cell_types=cell_types,
    max_attempts=1000,
    min_spacing=0.5,
    allow_boundary_cells=True
)
```

### 2. Set Up Target Statistics

You can either create target statistics programmatically or load from CSV.

#### Option A: Create Programmatically

```python
target_stats = {
    'node_counts': {
        'cancer': 40,
        'immune': 30,
        'stroma': 30
    },
    'edge_counts': {
        'cancer-cancer': 45,
        'cancer-immune': 35,
        'cancer-stroma': 25,
        'immune-immune': 20,
        'immune-stroma': 15,
        'stroma-stroma': 30
    },
    'neighbor_dist': {
        'cancer': {'cancer': 2.5, 'immune': 1.8, 'stroma': 1.2},
        'immune': {'cancer': 2.0, 'immune': 1.3, 'stroma': 1.0},
        'stroma': {'cancer': 1.5, 'immune': 1.0, 'stroma': 2.0}
    }
}
```

#### Option B: Load from CSV

Create a CSV file with columns like:
- `nodes_cancer`, `nodes_immune`, `nodes_stroma`
- `edges_cancer-cancer`, `edges_cancer-immune`, etc.

```python
from tissue_simulator import load_graph_stats_csv

target_stats = load_graph_stats_csv(
    "target_statistics.csv",
    color_names=['cancer', 'immune', 'stroma']
)
```

### 3. Run the Workflow

```python
from tissue_simulator import TissueNetworkWorkflow

# Create workflow manager
workflow = TissueNetworkWorkflow()

# Run complete workflow
evaluation = workflow.run_complete_workflow(
    tissue=tissue,
    z_position=tissue.thickness / 2,  # Middle slice
    network_radius=50.0,              # Connection radius in μm
    target_stats_dict=target_stats,
    cell_types=['cancer', 'immune', 'stroma'],
    annealing_params={
        'initial_temp': 1000.0,
        'final_temp': 0.01,
        'cooling_rate': 0.998,
        'max_iterations': 15000
    },
    export_dir="results",
    visualize=True
)
```

### 4. Analyze Results

```python
# Get final statistics
final_stats = workflow.get_statistics()

# Compare with target
differences = workflow.compare_statistics(verbose=True)

# Comprehensive evaluation
evaluation = workflow.evaluate(print_report=True)

# Visualization
workflow.visualize_slice(save_path="slice.png")
workflow.visualize_network(save_path="network.png")
```

### 5. Export Results

```python
# Export all results
workflow.export_all(base_dir="results", prefix="tissue")

# Or export individually
workflow.export_slice_csv("slice.csv")
workflow.export_network("network.graphml")
workflow.export_statistics_csv("statistics.csv")
```

## Step-by-Step Manual Workflow

For more control, you can run each step manually:

```python
from tissue_simulator import (
    TissueSlicer,
    SpatialNetworkAnalyzer,
    GraphColorizer,
    evaluate_graph_coloring
)

# 1. Create slice
slicer = TissueSlicer(tissue)
slice_cells = slicer.slice_plane(z_position=tissue.thickness / 2)

# 2. Build network
analyzer = SpatialNetworkAnalyzer()
graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=50.0)

# 3. Assign cell types
colorizer = GraphColorizer(
    target_graph=graph,
    colors=['cancer', 'immune', 'stroma'],
    target_statistics=target_stats
)

cell_type_assignment = colorizer.colorize(
    initial_temp=1000.0,
    final_temp=0.01,
    cooling_rate=0.998,
    max_iterations=15000
)

# 4. Apply to slice cells
for i, slice_cell in enumerate(slicer.slice_cells):
    slice_cell.cell_type = cell_type_assignment[i]

# 5. Evaluate
from tissue_simulator import calculate_graph_statistics
final_stats = calculate_graph_statistics(
    graph, 
    cell_type_assignment, 
    ['cancer', 'immune', 'stroma']
)
```

## Simulated Annealing Parameters

The cell type assignment uses simulated annealing optimization. Key parameters:

- **initial_temp**: Starting temperature (default: 100.0)
  - Higher = more exploration
  - Typical range: 100-1000
  
- **final_temp**: Stopping temperature (default: 0.1)
  - Lower = more refinement
  - Typical range: 0.01-1.0
  
- **cooling_rate**: Temperature decrease rate (default: 0.995)
  - Closer to 1 = slower cooling
  - Typical range: 0.99-0.999
  
- **max_iterations**: Maximum iterations (default: 100000)
  - More = better results but slower
  - Typical range: 5000-50000

## Reproducibility (seed)

`GraphColorizer.__init__` accepts an optional `seed: Optional[int] = None`.
When provided, the colorizer routes every stochastic decision (the initial
coloring shuffle, the per-step pair sampling, and the metropolis
acceptance draw) through an instance-bound `random.Random(seed)`, making
`colorize(...)` bit-reproducible across Python processes — independent
of `PYTHONHASHSEED`. When `seed=None` (default), behavior is unchanged
and the unseeded stdlib `random` module is used. This matches the v0.1.2
reproducibility plumbing already documented for `TissueSection` and
`ReplicateGenerator` in [`../quickstart.md`](../quickstart.md).

The seed is also threaded through the workflow entry points:
`TissueNetworkWorkflow.assign_cell_types(..., seed=N)`,
`TissueNetworkWorkflow.run_complete_workflow(..., seed=N)`, and the
top-level `quick_workflow(..., seed=N)` all forward it to the underlying
`GraphColorizer`.

```python
from tissue_simulator import GraphColorizer

# Two colorizers with the same seed and identical inputs produce
# bit-identical assignments.
c1 = GraphColorizer(
    target_graph=graph,
    colors=['cancer', 'immune', 'stroma'],
    target_statistics=target_stats,
    seed=42,
)
a1 = c1.colorize(initial_temp=500.0, final_temp=0.01,
                 cooling_rate=0.998, max_iterations=5000, verbose=False)

c2 = GraphColorizer(
    target_graph=graph,
    colors=['cancer', 'immune', 'stroma'],
    target_statistics=target_stats,
    seed=42,
)
a2 = c2.colorize(initial_temp=500.0, final_temp=0.01,
                 cooling_rate=0.998, max_iterations=5000, verbose=False)

assert a1 == a2  # identical dicts
```

### Tuning Tips

- For quick testing: `initial_temp=100, max_iterations=5000`
- For better results: `initial_temp=1000, max_iterations=20000, cooling_rate=0.998`
- For high accuracy: `initial_temp=2000, max_iterations=50000, cooling_rate=0.999`

## Network Building Modes

Two modes for defining cell connections:

### Contact Mode
Cells are connected if they're touching:
```python
graph = analyzer.build_network_from_slice(slicer, mode="contact")
```

### Radius Mode
Cells are connected if within a distance threshold:
```python
graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=50.0)
```

**Choosing radius:**
- Small radius (20-40 μm): Local interactions
- Medium radius (40-60 μm): Neighborhood structure
- Large radius (60-100 μm): Broader patterns

## Evaluation Metrics

The evaluation provides multiple metrics:

- **Mean Absolute Error**: Average absolute difference in edge counts
- **Root Mean Squared Error**: Square root of average squared differences
- **Cosine Similarity**: Directional similarity (0-1, higher is better)
- **Cosine Distance**: 1 - cosine similarity (0-2, lower is better)
- **Jensen-Shannon Divergence**: Statistical divergence (0-1, lower is better)
- **Percent Differences**: Per-edge-type differences

```python
evaluation = workflow.evaluate(print_report=True)

# Access specific metrics
mae = evaluation['mean_absolute_error']
js_div = evaluation['js_divergence']
avg_diff = evaluation['avg_percent_difference']
```

## Output Files

The workflow generates several output files:

1. **tissue_slice.csv**: 2D slice cell data
   - Columns: x_2d, y_2d, intersection_radius, cell_type, etc.

2. **tissue_network.graphml**: Network graph (NetworkX compatible)
   - Nodes: cells with positions and types
   - Edges: spatial connections

3. **tissue_statistics.csv**: Network statistics
   - Node counts per cell type
   - Edge counts between cell types

4. **tissue_slice.png**: Visualization of 2D slice
5. **tissue_network.png**: Visualization of network graph

## Advanced Usage

### Custom Target Statistics from Existing Tissue

```python
# Extract statistics from a reference tissue
from tissue_simulator import calculate_graph_statistics

reference_graph = analyzer.build_network_from_slice(reference_slicer)
reference_coloring = {i: cell.cell_type for i, cell in enumerate(reference_slicer.slice_cells)}
target_stats = calculate_graph_statistics(reference_graph, reference_coloring, cell_types)
```

### Batch Processing Multiple Slices

```python
from tissue_simulator import create_standard_slices

# Create multiple slices
slicers = create_standard_slices(tissue, num_slices=5)

# Process each slice
for i, slicer in enumerate(slicers):
    workflow = TissueNetworkWorkflow()
    workflow.slicer = slicer
    workflow.build_network(mode="radius", radius=50.0)
    # ... continue workflow
```

### Custom Color Palettes for Visualization

```python
color_palette = {
    'cancer': '#d62728',   # Red
    'immune': '#1f77b4',   # Blue
    'stroma': '#2ca02c'    # Green
}

workflow.visualize_network(color_palette=color_palette)
```

When `color_palette` is **not** provided, `visualize_colored_graph` and
`visualize_graph_comparison` build a default palette by sorting cell types
alphabetically before mapping them to the `tab10` colors, so the same set of
types always yields the same colors regardless of input order or the Python
interpreter's `PYTHONHASHSEED`. When you do pass a `color_palette` dict, it
is honored as-is — color ordering is then under your control. Added in
v0.1.12.

## Troubleshooting

### Issue: Poor convergence (high JS divergence)

**Solutions:**
- Increase `max_iterations`
- Increase `initial_temp`
- Decrease `cooling_rate` (e.g., 0.998 → 0.999)
- Check that target statistics are achievable given network structure

### Issue: Slow performance

**Solutions:**
- Reduce `max_iterations`
- Use smaller network (smaller radius or fewer cells)
- Reduce tissue size

### Issue: Memory errors with large networks

**Solutions:**
- Process in batches
- Reduce network density (smaller radius)
- Use smaller tissue sections

## References

For more information, see:
- `examples/complete_workflow_example.py`: Complete working example
- `tissue_simulator/tissue_workflow.py`: Workflow implementation
- `tissue_simulator/graph_coloring.py`: Graph coloring algorithm
- `tissue_simulator/evaluation.py`: Evaluation metrics

## Citation

If you use this workflow in your research, please cite:

```
[Add appropriate citation information]
```
