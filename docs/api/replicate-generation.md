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

## Neighbour-graph rule (`network_mode`)

`network_mode` is `"contact"`, `"radius"` or `"mechanical"` (edge iff
`d <= interaction_factor * (r_i + r_j)`, default factor 1.5 after PhysiCell's
mechanics interaction distance). `ReplicateGenerator(...)` still defaults to
`"contact"`; **`ReplicateGenerator.from_coordinates` now defaults to
`network_mode="mechanical"`** (`network_radius` is ignored unless the mode is
`"radius"`). Pass `network_mode="radius", network_radius=20.0` to reproduce the
earlier behaviour. `interaction_factor` is accepted by `ReplicateGenerator`,
`from_coordinates`, `load_target_statistics_from_tissue` and
`load_target_statistics_from_coordinates`.

`interaction_factor="auto"` learns the factor from the source: the first
minimum of its size-normalised pair correlation `g(s)`, `s = d / (r_i + r_j)`
(bins of 0.05 in `s`, Gaussian smoothing with sigma 0.1, first minimum after the
first peak and significantly below it: at least 50 cells, a peak of at least 1.2,
a minimum at least 15% below the peak, and a peak-to-minimum difference of more
than 4 Poisson standard errors). If there is none (a Poisson source), the factor
is 1.5 and a warning is raised. A dense random-sequential-addition source can
have a genuine shell and then returns a value. `"auto"` measures the x-y
projection, so it is meant for 2-D or thin-slab sources.
The factor is learned once, on the source. `from_coordinates` and the
`load_target_statistics_from_*` functions record it in
`TargetStatistics.network_rule`, and `ReplicateGenerator(..., interaction_factor="auto")`
reads it from there, so every replicate graph uses the source's number and is
never refit. If the target statistics carry no mechanical factor (for example
CSV tables), the generator warns and uses 1.5. `gen.interaction_factor` is
always a float.

Targets and replicates must be measured with the same rule.
`TargetStatistics.network_rule` records the rule
(`{"mode", "radius", "interaction_factor", "interaction_factor_source"}`, `None`
for CSV tables); the source is `"fixed"`, `"auto"` or `"auto_fallback"`.
`ReplicateGenerator` warns when the mode, radius or numeric factor differs from
its own `network_rule`; the source label is not compared.

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
- **Step budget**: the default schedule cools from 100 to 0.1 at a rate of
  0.995 per step, which ends after about 1,400 swaps whatever
  `max_iterations` is. To use a larger budget with the legacy strategy, pass a
  matching `cooling_rate`, for example `(0.1 / 100) ** (1 / max_iterations)`.
  With `strategy="adaptive"` and a density model the generator does this
  itself: the budget is the larger of 20,000 and ten swaps per cell, unless
  `coloring_params` sets `cooling_rate` or `max_iterations`.
- **Adaptive stopping**: pass `patience` inside `coloring_params` to stop SA once
  the cost plateaus (uses the cost trajectory; see `colorize(return_history=True)`
  and `convergence.find_convergence_time`).
- **`parallel=True`**: replicates are independent and deterministically seeded,
  so `generate_replicates(parallel=True)` runs them across processes with
  identical results to the serial path.

### Density-aware scaffolds

The uniform scaffold spreads cells evenly, so replicates of a real region
lose its dense tumor nests, sparse stroma and immune margins. When the source
region's coordinates are available, fit a
[`DensityModel`](core.md#densitymodel) and pass it to the generator: each
replicate is then packed on its own density layout, and the annealer also
matches the layout's expected composition in 40 µm bins.

```python
from tissue_simulator import ReplicateGenerator

# Fits target statistics and the density model from one coordinate CSV.
gen = ReplicateGenerator.from_coordinates(
    "region.csv", network_mode="radius", network_radius=20.0,
    seed=42, layout="resample",
)
replicates = gen.generate_replicates(num_replicates=30, parallel=True)
tissue, stats = replicates[0]
print(stats.composition_error, stats.packing_report["bin_correlation"], stats.layout_flags)
```

Or build the pieces yourself with
`DensityModel.from_tissue(region)` / `fit_density_model_from_coordinates(path)`
and `ReplicateGenerator(..., method="graph_coloring", density_model=model)`.

- **`layout="resample"`** (default) draws a new arrangement of compartments
  per replicate, with the region's compartment areas, within-compartment
  density distribution, composition and patch structure. Replicates differ
  in where the nests are, so the ensemble keeps initial-condition variance.
- **`layout="copy"`** reuses the region's own maps; replicates differ only
  below the smoothing bandwidth. Use it for regions flagged `"trend"` or
  `"patch_length_at_upper_bound"`, where resampling assumes a stationarity
  the region does not have.
- **`composition_weight`** scales the composition term. Under the legacy
  strategy it is multiplied by the squared mean degree of each replicate graph
  (default 4.0, unchanged). Under `strategy="adaptive"` the term is calibrated
  per replicate so the value is its size relative to the edge-count term for a
  shuffled labeling (default 1.0; a multi-seed sweep on twelve synthetic
  scenarios showed 4.0 over-constrains strongly clustered samples under the
  mechanical graph). Larger values trade pair-fraction accuracy for
  composition accuracy.
- Regions that are no more heterogeneous than a uniform packing
  (`model.homogeneous`) get uniform layouts.
- For a 2D source, use a thin slab (thickness about 1 µm) so replicate graphs
  stay planar.

#### Adaptive strategy (opt-in)

`strategy="legacy"` is the default and reproduces earlier releases exactly;
switching the default is a separate reviewed change. `strategy="adaptive"`
fits the model with data-derived bandwidths, per-type bandwidths and a
planar/radial organization trend, uses multi-scale composition targets and a
size-compatibility term, and records fidelity diagnostics.

```python
gen = ReplicateGenerator.from_coordinates(
    "region.csv", network_mode="radius", network_radius=20.0,
    seed=42, strategy="adaptive",
)
tissue, stats = gen.generate_single_replicate(0)
print(stats.layout_organization, stats.fidelity["size_ks_by_type"])
```

Keyword-only generator options: `strategy`, `composition_scales` (list of bin
sides in µm; default `"auto"` under adaptive), `composition_weight`,
`size_weight` (default 1.0 adaptive, 0 legacy), `diagnostics` (default on for
adaptive), `max_proposals` (organization re-draws, default 20).
`DensityModel.fit` accepts `strategy`, `bandwidth_range` (`(lo, hi)`, `"auto"`
or None), `per_type_bandwidth` and `organization`. It also takes `voids` (`"auto"`/`"none"`; default `"auto"` under
adaptive): lumens and holes in the source are inferred and re-placed in each
layout, reported as `layout_voids` on the replicate statistics. `units` (`"auto"`/`"none"`;
default `"auto"` under adaptive) detects compact units (nests, follicles, glomeruli) by persistent
homology and re-places them as germ-grain units with fitted radial profiles; they appear as
`layout_units` and take precedence over voids inside them and over a radial trend centered in a
unit. Criteria and limits are in
[the design notes](../notes/density-aware-packing.md).

Design, ablations and known limits are in
[the design notes](../notes/density-aware-packing.md).

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
    seed=None,               # for reproducibility
    method="radius_tuning",  # or "graph_coloring"
    density_model=None,      # DensityModel: density-aware scaffold (graph_coloring only)
    layout="resample",       # or "copy"; used with density_model
    composition_weight=4.0,  # composition term weight (times squared mean degree)
    composition_bin=40.0,    # composition bin side in µm
    packing_params=None,     # extra InhomogeneousPacker arguments
)
```

`packing_params` is passed to `InhomogeneousPacker`. Useful keys for adaptive
models: `refine_shell` (None, the default, turns first-shell refinement on for
adaptive layouts; False turns it off), `refine_sweeps` (maximum sweeps, default
20), `refine_cap` (largest distance in µm a cell may move during refinement,
default one median radius) and `radius_assignment` (`"deck"` or
`"per_candidate"`). Refinement runs only for thin slabs (thickness at most twice
the median radius) of at least 50 cells, and only when the model carries a
first-shell profile; see
[First-shell fidelity](../notes/density-aware-packing.md#first-shell-fidelity).

`ReplicateGenerator.from_coordinates(path, ...)` builds a density-aware
generator directly from a coordinate CSV; see
[Density-aware scaffolds](#density-aware-scaffolds).

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
- `packing_report`: `PackingReport.to_dict()` of the replicate's scaffold
  (density-aware replicates only)
- `composition_error`: Fraction of cells whose type would have to move between
  composition bins to match the layout (density-aware replicates only)
- `layout_flags`: Layout mode (for example `"mode:resample"`) followed by the
  density model's flags (density-aware replicates only)
- `requested_cell_type_counts`, `achieved_cell_type_counts`: per-type cell
  quota from the layout and final counts (density-aware replicates)
- `layout_organization`: the layout's organization dict (`model`, `geometry`,
  `direction`/`center`, `proposals_tried`, `accepted`, `fallback`)
- `layout_voids`: placed voids of the layout (`centers`, `radii`, `n_requested`,
  `n_placed`, `anchored`); None when the model has none
- `layout_units`: placed units of the layout (`centers`, `outer_radii`,
  `inner_radii`, `kinds`, `n_requested`, `n_placed`, `shortfall`, `anchored`);
  None when the model has none
- `fidelity["anneal_budget"]`: the maximum number of swaps the annealer was
  allowed (adaptive strategy). A user-supplied `cooling_rate` or `patience` in
  `coloring_params` can stop the run earlier, so it is not the number of swaps
  run.
- `fidelity`: diagnostics against the source (`size_nll`, `size_ks_by_type`,
  `nn_distance_quantiles`, `mixing_index`, `organization_rmse`,
  `interface_fraction`, `n_components`, `n_holes`, `persistence_distance`
  per type: H0 Wasserstein-1 of source vs replicate KDE maps); None unless
  `diagnostics` is on. Evidence, not calibrated intervals.
- `fidelity["first_shell"]` (adaptive replicates with diagnostics on):
  `{"factor": 1.5, "replicate": summary, "source": summary, "ratios": {...}, "same_window": bool}`.
  `same_window` says whether the replicate window has the source's size
  (None for models saved before this field existed); the ratios are only
  comparable, and the warning is only raised, when it is true.
  A summary holds `mean_degree` (mechanical graph at the factor),
  `overlap_pairs_per_cell` (pairs closer than `r_i + r_j`), `median_radius` and
  `area_fraction`; `ratios` holds replicate over source for each. A warning is
  raised when an adaptive replicate's `mean_degree` ratio is below 0.9. The
  source's value assumes the replicate window has the source's cell density.
- `separation`: nearest-neighbor `clearance_quantiles`,
  `normalized_distance_quantiles` and `n_nearest_neighbour`

`packing_report` (`PackingReport.to_dict()`) also carries `bin_shortfall`,
`quota_floor`, `quota_floor_source` (`"layout"` or `"legacy"`),
`clearance_quantiles`, `normalized_distance_quantiles` and
`dense_bin_fraction_short`, `radius_assignment` (`"deck"` for adaptive
layouts, `"per_candidate"` otherwise), `refinement` (None, or `sweeps`,
`energy_before`, `energy_after`, `n_moved`, `mean_displacement`,
`max_displacement`, `accept_rate`, `s_target`, `seconds`) and `first_shell`.

`layout_organization` also records `failed` (`"coverage"` or `"composition"`)
and `best_coverage` when no proposal met the criteria and the best-coverage
fallback was used. Replicate windows much smaller than the source region
cannot span the source's trend range, so the fallback is likely there. With a
trend in the model, use a window at least as large as the source region.

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

   For a density-aware scaffold add `density_layout` (`"resample"` or
   `"copy"`, needs source coordinates loaded with
   `load_target_statistics_from_coordinates`) and optionally `strategy`
   (`"legacy"` default or `"adaptive"`), `diagnostics`, `composition_weight`
   and `size_weight`. The result echoes `strategy` and `diagnostics`; adaptive
   fits the source model with `strategy="adaptive"`.

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
