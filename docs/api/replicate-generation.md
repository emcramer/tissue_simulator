# Replicate Generation Documentation

## Overview

The replicate generation module allows you to generate multiple tissue samples that match specific spatial interaction patterns. This is useful for:

- **Statistical analysis**: Generate multiple samples for robust statistical comparisons
- **Validation studies**: Create control tissues with known properties
- **Hypothesis testing**: Generate tissues matching theoretical interaction patterns
- **Simulation studies**: Create ensembles of tissues for computational experiments

> **When to reach for `GraphColorizer` instead.** `ReplicateGenerator` matches
> contact statistics and proportions through repacking; it does NOT perform
> simulated-annealing label assignment, and it will not converge on strongly
> structured multi-type targets (e.g. a tumor disc with a fibroblast ring or
> a CD8 annulus). For those, use the two-stage `TissueWorkflow` path —
> `generate_cells` for positions, then `assign_cell_types` for labels via
> `GraphColorizer`. See the "When to use what" section of
> [`graph-coloring.md`](graph-coloring.md#when-to-use-what).
>
> **Geometric vs. colored replicates.** `ReplicateGenerator` produces
> *geometric* replicates — each sample has a different cell packing. To instead
> hold one geometry fixed and draw several diverse *type labelings* of it, use
> `TissueNetworkWorkflow.generate_colored_replicates(n)`, documented under
> [Generating colored replicates](graph-coloring.md#generating-colored-replicates).

## Key Features

- **Target-based generation**: Generate tissues matching specified spatial statistics
- **Multiple sources**: Load targets from CSV files or existing tissues
- **Iterative optimization**: Automatically tunes parameters to match targets
- **Flexible constraints**: Control cell types, proportions, and packing
- **Batch processing**: Generate multiple replicates efficiently
- **Comprehensive export**: Save tissues and statistics in multiple formats
- **MCP integration**: Fully accessible via LLM coding assistant

## Generation methods

`ReplicateGenerator` supports two strategies via the `method` argument:

### `method="radius_tuning"` (default)

The original approach: each replicate repacks a tissue and a heuristic nudges
per-type radii to steer cell-type **proportions**. Because cell types are
assigned at packing time and the only lever on spatial **interaction** patterns
is an indirect radius proxy, this converges slowly and inconsistently for
interaction targets. Kept as the default for backward compatibility and for
cases where you specifically want to tune geometry/density.

Two knobs improve it:

- `radius_optimizer="differential_evolution"` replaces the sqrt-ratio heuristic
  with a gradient-free SciPy optimizer over per-type radius multipliers against
  a fixed-seed (deterministic) proportion objective. Gradient methods are
  deliberately *not* offered — the radius→cell-count map is integer-valued and
  stochastic, so finite-difference gradients are mostly zero. DE is more robust
  but slower.
- `patience=N` on `generate_replicates` / `generate_single_replicate` stops the
  tuning loop early once the best divergence plateaus.

### `method="graph_coloring"` (recommended for interaction targets)

Matching interaction statistics is fundamentally a **labeling** problem on a
fixed neighbor graph, not a geometry problem. This mode packs geometry **once**
per replicate (fresh seed → geometric diversity), builds the neighbor graph,
then assigns cell types with the simulated-annealing
[`GraphColorizer`](graph-coloring.md) to match the target interaction
statistics. Cell-type proportions are locked exactly (the SA swap-moves
preserve node counts), and convergence is far more consistent.

```python
gen = ReplicateGenerator(
    target_stats=target,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'cancer': (8, 12), 'immune': (5, 8), 'fibroblast': (6, 10)},
    network_mode="radius", network_radius=30.0, seed=42,
    method="graph_coloring",
    n_restarts=3,                       # keep best of 3 SA runs per replicate
    coloring_params={'max_iterations': 8000, 'patience': 2000},
)
replicates = gen.generate_replicates(num_replicates=10, parallel=True)
```

Supporting features:

- **`n_restarts`**: run several independent SA colorings per replicate and keep
  the lowest-cost one (hardens against bad local minima).
- **Adaptive stopping**: pass `patience` inside `coloring_params` to stop SA once
  the cost plateaus (uses the cost trajectory; see `colorize(return_history=True)`
  and `convergence.find_convergence_time`).
- **`parallel=True`**: replicates are independent and deterministically seeded,
  so `generate_replicates(parallel=True)` runs them across processes with
  identical results to the serial path.

### Measuring consistency

`consistency_report` quantifies run-to-run variability and compares methods
(via the `power_analysis` module), so you can *prove* an approach is more
consistent:

```python
report = gen.consistency_report({
    "radius_tuning": radius_replicates,
    "graph_coloring": colored_replicates,
})
# report["per_method"][name] -> {n, mean, std, cv};  report["pairwise"] -> Cohen's d, required N
```

Note: the coefficient of variation is `std/|mean|` and inflates when the mean
sits near zero (a method matching the target almost perfectly), so read `cv`
alongside the absolute `std`/`mean`.

### Related algorithm families

Matching target spatial statistics is well-studied. For the **labeling**
subproblem (used here) the relevant family is Potts/Markov-random-field models
optimized by simulated annealing — which is what `GraphColorizer` implements,
optionally hardened by parallel tempering / population annealing. For the
**geometry** subproblem (if you need the point pattern itself to match spatial
summary functions) the relevant families are Gibbs/Markov point processes
(Strauss, area-interaction, multitype Gibbs via MCMC birth-death-move), SA
reconstruction to pair-correlation `g(r)` / Ripley's K / nearest-neighbor
distributions, and cluster processes (Matérn / Thomas / log-Gaussian Cox) for
clustered arrangements.

## Installation

The replicate generator requires pandas and NetworkX:

```bash
pip install pandas networkx
```

## Quick Start

### Generate from Existing Tissue

```python
from tissue_simulator import (
    TissueSection,
    load_target_statistics_from_tissue,
    ReplicateGenerator
)

# Create reference tissue
reference = TissueSection(400, 400, 100, 
                         cell_radii={'type_a': (8, 12), 'type_b': (5, 8)})
reference.generate_cells(max_attempts=1000)

# Extract spatial statistics
target_stats = load_target_statistics_from_tissue(
    reference, 
    network_mode="contact"
)

# Setup generator
generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'type_a': (8, 12), 'type_b': (5, 8)},
    network_mode="contact"
)

# Generate replicates
replicates = generator.generate_replicates(num_replicates=10)
```

### Generate from CSV File

```python
from tissue_simulator import (
    load_target_statistics_from_csv,
    ReplicateGenerator
)

# Load target statistics from CSV
target_stats = load_target_statistics_from_csv("statistics.csv")

# Setup and generate
generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'type_a': (8, 12), 'type_b': (5, 8)},
    network_mode="contact"
)

replicates = generator.generate_replicates(num_replicates=10)
```

### Generate from a Coordinate CSV

`load_target_statistics_from_csv` (above) expects a *precomputed interaction
table* and leaves `cell_type_proportions` and `target_density` unset. When you
instead have **raw cell coordinates** — a measured sample, or a layout exported
from another tool to the `export_to_csv` schema (`x, y, z, radius, cell_type,
is_boundary`) — use `load_target_statistics_from_coordinates`, which returns a
*fully populated* `TargetStatistics` (interactions **plus** proportions and
density). Added in v0.1.7.

```python
from tissue_simulator import (
    load_target_statistics_from_coordinates,
    ReplicateGenerator,
)

# Full statistics straight from a coordinate CSV.
target_stats = load_target_statistics_from_coordinates(
    "cells.csv",
    network_mode="radius",
    network_radius=20.0,
)

generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'type_a': (8, 12), 'type_b': (5, 8)},
    network_mode="radius",
    network_radius=20.0,
)
replicates = generator.generate_replicates(num_replicates=10)
```

It is exactly equivalent to
`load_target_statistics_from_tissue(load_tissue_from_csv(path), ...)`; reach for
it when your source is coordinates rather than an interaction table. See
[`core.md`](core.md) for `load_tissue_from_csv` and `TissueSection.from_cells`.

## CSV Format for Target Statistics

The CSV file should contain interaction statistics with these columns:

```csv
type_a,type_b,num_interactions,normalized_interactions,avg_distance,median_distance
cancer,cancer,45,0.12,0.0,0.0
cancer,immune,38,0.15,0.0,0.0
cancer,fibroblast,42,0.14,0.0,0.0
immune,immune,28,0.18,0.0,0.0
immune,fibroblast,35,0.16,0.0,0.0
fibroblast,fibroblast,30,0.11,0.0,0.0
```

**Required columns:**
- `type_a`: First cell type
- `type_b`: Second cell type  
- `normalized_interactions`: Normalized interaction frequency (0-1)

**Optional columns:**
- `num_interactions`: Raw interaction count
- `avg_distance`: Average interaction distance
- `median_distance`: Median interaction distance

> **Two CSV shapes.** This interaction-table format is what
> `load_target_statistics_from_csv` reads. A *coordinate* CSV
> (`x, y, z, radius, cell_type, is_boundary`, written by
> `TissueSection.export_to_csv`) is a different shape — load it with
> `load_target_statistics_from_coordinates` instead.

## Core Classes

### TargetStatistics

Defines the spatial patterns to match:

```python
from tissue_simulator import TargetStatistics, InteractionStatistics

target = TargetStatistics(
    interaction_stats=[
        InteractionStatistics(
            type_a='cancer',
            type_b='immune',
            num_interactions=40,
            normalized_interactions=0.15,
            avg_distance=0.0,
            median_distance=0.0
        )
    ],
    cell_type_proportions={'cancer': 0.6, 'immune': 0.4},
    target_cell_count=200,
    target_density=0.45
)
```

**Attributes:**
- `interaction_stats`: List of InteractionStatistics objects
- `cell_type_proportions`: Dict of target proportions (optional)
- `target_cell_count`: Target total cell count (optional)
- `target_density`: Target packing fraction (optional)

### ReplicateGenerator

Main class for generating tissue replicates:

```python
generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(height, width, thickness),
    base_cell_radii={'type_a': (min_r, max_r), ...},
    network_mode="contact",  # or "radius"
    network_radius=None,     # required if mode="radius"
    seed=None                # for reproducibility
)
```

**Key Methods:**

#### generate_single_replicate()
```python
tissue, stats = generator.generate_single_replicate(
    replicate_id=0,
    max_attempts=1000,
    min_spacing=0.5,
    allow_boundary=True,
    max_iterations=5,
    tolerance=0.15
)
```

Generates a single replicate with iterative parameter adjustment.

**Parameters:**
- `replicate_id`: Unique identifier
- `max_attempts`: Max cell packing attempts
- `min_spacing`: Minimum cell spacing (μm)
- `allow_boundary`: Allow cells beyond bounds
- `max_iterations`: Max optimization iterations
- `tolerance`: Acceptable divergence threshold

**Returns:**
- `tissue`: TissueSection object
- `stats`: ReplicateStatistics object

#### generate_replicates()
```python
replicates = generator.generate_replicates(
    num_replicates=10,
    max_attempts=1000,
    min_spacing=0.5,
    allow_boundary=True,
    max_iterations=5,
    tolerance=0.15
)
```

Generates multiple replicates.

**Returns:**
- List of (TissueSection, ReplicateStatistics) tuples

### ReplicateStatistics

Contains statistics for a generated replicate:

```python
stats = ReplicateStatistics(
    replicate_id=0,
    num_cells=195,
    cell_type_counts={'cancer': 120, 'immune': 75},
    packing_fraction=0.438,
    interaction_stats=[...],
    divergence_score=0.083
)
```

**Attributes:**
- `replicate_id`: Unique identifier
- `num_cells`: Total cell count
- `cell_type_counts`: Dict of counts per type
- `packing_fraction`: Volume fraction
- `interaction_stats`: Measured interactions
- `divergence_score`: Divergence from target (lower = better)

## Export Functions

### Export Statistics

```python
generator.export_replicate_statistics(replicates, "output/stats")
```

Creates two CSV files:
- `output/stats_summary.csv`: Overall replicate statistics
- `output/stats_interactions.csv`: Detailed interaction data

### Export Tissues

```python
generator.export_replicate_tissues(replicates, "output/tissues")
```

Creates individual CSV files:
- `output/tissues/replicate_000_tissue.csv`
- `output/tissues/replicate_001_tissue.csv`
- etc.

## Algorithm Details

### Optimization Process

The generator uses an iterative approach:

1. **Initial Generation**: Create tissue with base parameters
2. **Analysis**: Measure spatial interactions
3. **Divergence Calculation**: Compare to target statistics
4. **Parameter Adjustment**: Tune cell radii and proportions
5. **Iteration**: Repeat until tolerance met or max iterations reached

### Divergence Metric

Divergence is calculated as the mean relative difference in normalized interactions:

```
divergence = mean(|measured - target| / target)
```

Lower values indicate better matches to target statistics.

### Parameter Adjustment

The generator adjusts:
- **Cell radii**: Scaled to change cell type proportions
- **Sampling**: Biased toward under-represented types

## Advanced Usage

### Custom Optimization

```python
# More aggressive optimization
replicates = generator.generate_replicates(
    num_replicates=5,
    max_attempts=2000,        # More packing attempts
    max_iterations=10,        # More optimization steps
    tolerance=0.05            # Stricter tolerance
)
```

### Reproducibility

```python
# Use seed for reproducible results
generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii=cell_radii,
    network_mode="contact",
    seed=42  # Reproducible
)
```

### Network Modes

**Contact Mode** (default):
```python
network_mode="contact"
```
- Connects cells that are touching
- Best for direct cell-cell interactions
- Faster computation

**Radius Mode**:
```python
network_mode="radius"
network_radius=30.0  # micrometers
```
- Connects cells within distance
- Captures proximity effects
- Useful for paracrine signaling

## MCP Integration

All functionality is available through the MCP server for LLM assistants.

### Workflow

1. **Load statistics**:
```
load_target_statistics(csv_filepath="stats.csv")
```
or
```
load_target_statistics(use_current_tissue=True)
```

2. **Setup generator**:
```
setup_replicate_generator(
    height=400,
    width=400,
    thickness=100,
    cell_radii={"cancer": [8, 12], "immune": [5, 8]}
)
```

3. **Generate replicates**:
```
generate_replicates(
    num_replicates=5,
    tolerance=0.15
)
```

4. **Export results**:
```
export_replicate_statistics(base_filename="output")
export_replicate_tissues(output_dir="tissues")
```

5. **Get summary**:
```
get_replicate_summary()
```

## Examples

### Example 1: Match Existing Tissue

See `examples/replicate_generation_example.py`

### Example 2: Load from CSV

See `examples/replicate_from_csv_example.py`

### Example 3: Batch Analysis

```python
# Generate multiple batches with different parameters
for tolerance in [0.05, 0.10, 0.15, 0.20]:
    generator = ReplicateGenerator(
        target_stats=target_stats,
        tissue_dimensions=(400, 400, 100),
        base_cell_radii=cell_radii,
        network_mode="contact"
    )
    
    replicates = generator.generate_replicates(
        num_replicates=10,
        tolerance=tolerance
    )
    
    generator.export_replicate_statistics(
        replicates, 
        f"output/tolerance_{int(tolerance*100)}"
    )
```

## Performance Considerations

### Generation Time

- **Single replicate**: 10-60 seconds
- **10 replicates**: 2-10 minutes
- **Factors**: Tissue size, cell count, tolerance

### Optimization Tips

1. **Start with relaxed tolerance** (0.15-0.20)
2. **Increase max_iterations** if divergence is high
3. **Use appropriate network mode** (contact is faster)
4. **Consider tissue size** (smaller = faster)

### Memory Usage

- Each tissue: ~1-10 MB
- 100 replicates: ~100-1000 MB
- Export regularly for large batches

## Interpretation Guide

### Divergence Scores

- **< 0.05**: Excellent match
- **0.05 - 0.10**: Good match
- **0.10 - 0.20**: Acceptable match
- **> 0.20**: Poor match (increase iterations)

### Cell Count Variation

Expect ±10-20% variation in cell counts across replicates due to:
- Random packing dynamics
- Iterative optimization
- Stochastic cell placement

### Interaction Patterns

Compare normalized interactions:
- **Within 5%**: Very similar patterns
- **Within 10%**: Similar patterns
- **Within 20%**: Broadly similar
- **> 20%**: Different patterns

## Troubleshooting

### Issue: High Divergence Scores

**Solutions:**
- Increase `max_iterations` (e.g., 10-15)
- Increase `tolerance` (e.g., 0.20)
- Check that cell types in target match generator config
- Verify target statistics are achievable

### Issue: Few Cells Generated

**Solutions:**
- Increase `max_attempts` (e.g., 2000-3000)
- Decrease `min_spacing` (e.g., 0.3-0.4)
- Reduce cell radii ranges
- Allow boundary cells

### Issue: Long Generation Time

**Solutions:**
- Reduce `max_iterations` (e.g., 3-5)
- Use contact mode instead of radius
- Reduce tissue dimensions
- Generate fewer replicates per batch

### Issue: Inconsistent Cell Type Proportions

**Solutions:**
- Specify `cell_type_proportions` in TargetStatistics
- Increase `max_iterations` for better tuning
- Check that base cell radii are reasonable
- Verify target proportions are achievable

## Best Practices

1. **Start Simple**: Test with small tissues first
2. **Validate Targets**: Ensure target statistics are realistic
3. **Monitor Progress**: Check divergence scores during generation
4. **Export Regularly**: Save results incrementally
5. **Use Seeds**: Enable reproducibility when needed
6. **Batch Processing**: Generate replicates in manageable batches
7. **Quality Control**: Review divergence scores before analysis

## See Also

- [Spatial Analysis Documentation](spatial-analysis.md)
- [Tissue Simulator Core API](core.md)
- [MCP Integration Guide](mcp.md)
- Example scripts in `examples/` directory
