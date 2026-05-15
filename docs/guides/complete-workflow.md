# Complete Workflow Guide: Graph-Based Cell Type Assignment

This guide demonstrates the complete workflow for generating simulated tissue sections with cell types assigned based on network statistics using simulated annealing.

## Table of Contents
1. [Quick Start](#quick-start)
2. [Complete Workflow](#complete-workflow)
3. [Understanding Target Statistics](#understanding-target-statistics)
4. [Tuning Simulated Annealing](#tuning-simulated-annealing)
5. [Network Building Modes](#network-building-modes)
6. [Evaluation Metrics](#evaluation-metrics)
7. [Examples](#examples)
8. [Troubleshooting](#troubleshooting)

---

## Quick Start

The fastest way to run the complete workflow:

```python
from tissue_simulator import TissueSection, SpherePacker, quick_workflow

# 1. Create tissue
tissue = TissueSection(height=300, width=300, thickness=80)
packer = SpherePacker(tissue)
packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1500)

# 2. Run complete workflow
workflow = quick_workflow(
    tissue=tissue,
    cell_types=['cancer', 'immune', 'stroma'],
    target_stats_file="target_statistics.csv",
    network_radius=50.0,
    output_dir="results"
)
```

---

## Complete Workflow

### Step 1: Generate 3D Tissue

First, create a 3D tissue structure using random sphere packing:

```python
from tissue_simulator import TissueSection, SpherePacker

# Define tissue dimensions (micrometers)
tissue = TissueSection(
    height=300,      # Y-axis
    width=300,       # X-axis
    thickness=80     # Z-axis
)

# Create packer and generate cells
packer = SpherePacker(tissue)
num_cells = packer.pack_cells(
    cell_types={'placeholder': [8, 12]},  # Initial assignment (will be reassigned)
    max_attempts=1500,
    min_spacing=0.5,
    allow_boundary_cells=True
)

print(f"Generated {num_cells} cells")
```

**Note**: The initial cell type assignment is just a placeholder. The graph coloring algorithm will reassign cell types based on network statistics.

### Step 2: Extract 2D Slice

Extract a 2D cross-section from the 3D tissue:

```python
from tissue_simulator import TissueSlicer

slicer = TissueSlicer(tissue)

# Option A: Horizontal slice at specified z-position
slice_cells = slicer.slice_plane(z_position=tissue.thickness / 2)

# Option B: Angled slice
slice_cells = slicer.slice_plane(angle_x=15.0, angle_y=30.0)

print(f"Slice contains {len(slice_cells)} cells")
```

### Step 3: Build Network Graph

Create a spatial relationship graph from the 2D slice:

```python
from tissue_simulator import SpatialNetworkAnalyzer

analyzer = SpatialNetworkAnalyzer()

# Option A: Contact mode (cells must touch)
graph = analyzer.build_network_from_slice(slicer, mode="contact")

# Option B: Radius mode (cells within distance threshold)
graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=50.0)

print(f"Network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
```

### Step 4: Define Target Statistics

Specify the desired spatial organization:

```python
target_stats = {
    'node_counts': {
        'cancer': 40,
        'immune': 35,
        'stroma': 25
    },
    'edge_counts': {
        'cancer-cancer': 60,    # Cancer cells cluster together
        'cancer-immune': 45,    # Moderate cancer-immune interaction
        'cancer-stroma': 30,    # Some cancer-stroma interaction
        'immune-immune': 25,    # Immune cells moderately clustered
        'immune-stroma': 20,    # Immune-stroma interaction
        'stroma-stroma': 15     # Stroma less clustered
    },
    'neighbor_dist': {
        'cancer': {'cancer': 3.0, 'immune': 2.2, 'stroma': 1.5},
        'immune': {'cancer': 2.5, 'immune': 1.8, 'stroma': 1.3},
        'stroma': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.2}
    }
}
```

### Step 5: Assign Cell Types

Use simulated annealing to assign cell types:

```python
from tissue_simulator import GraphColorizer

colorizer = GraphColorizer(
    target_graph=graph,
    colors=['cancer', 'immune', 'stroma'],
    target_statistics=target_stats
)

cell_type_assignment = colorizer.colorize(
    initial_temp=1000.0,
    final_temp=0.01,
    cooling_rate=0.998,
    max_iterations=10000,
    verbose=True
)

# Apply to slice cells
for i, slice_cell in enumerate(slicer.slice_cells):
    slice_cell.cell_type = cell_type_assignment[i]
```

### Step 6: Visualize Results

Create visualizations:

```python
from tissue_simulator import visualize_colored_graph

# Visualize network
visualize_colored_graph(
    graph,
    cell_type_assignment,
    layout="spring",
    title="Tissue Network with Assigned Cell Types",
    save_path="network.png"
)

# Visualize 2D slice
slicer.visualize_slice_2d()
```

### Step 7: Evaluate Against Target

Compare results with target statistics:

```python
from tissue_simulator import (
    calculate_graph_statistics,
    evaluate_graph_coloring,
    print_evaluation_report
)

# Calculate final statistics
final_stats = calculate_graph_statistics(
    graph, 
    cell_type_assignment, 
    ['cancer', 'immune', 'stroma']
)

# Convert target stats to comparison format
target_stats_formatted = {}
for color in ['cancer', 'immune', 'stroma']:
    target_stats_formatted[f'nodes_{color}'] = target_stats['node_counts'][color]
for key, value in target_stats['edge_counts'].items():
    target_stats_formatted[f'edges_{key}'] = value

# Comprehensive evaluation
evaluation = evaluate_graph_coloring(target_stats_formatted, final_stats)
print_evaluation_report(target_stats_formatted, final_stats, evaluation)
```

### Step 8: Export Results

Save all results:

```python
from tissue_simulator import export_colored_graph_statistics

# Export slice data
slicer.export_slice_csv("slice_data.csv", include_3d=True)

# Export network
analyzer.export_network("network.graphml", format="graphml")

# Export statistics
export_colored_graph_statistics(
    graph,
    cell_type_assignment,
    ['cancer', 'immune', 'stroma'],
    "statistics.csv"
)
```

---

## Understanding Target Statistics

Target statistics define the desired spatial organization with three components:

### 1. Node Counts
How many cells of each type:
```python
'node_counts': {
    'cancer': 40,   # 40 cancer cells
    'immune': 35,   # 35 immune cells
    'stroma': 25    # 25 stromal cells
}
```

### 2. Edge Counts
How many connections between cell type pairs:
```python
'edge_counts': {
    'cancer-cancer': 60,    # 60 connections between cancer cells
    'cancer-immune': 45,    # 45 connections between cancer and immune
    'cancer-stroma': 30,    # 30 connections between cancer and stroma
    'immune-immune': 25,    # etc.
    'immune-stroma': 20,
    'stroma-stroma': 15
}
```

**High edge counts** indicate clustering or strong spatial association.
**Low edge counts** indicate dispersion or weak spatial association.

### 3. Neighbor Distribution
Average number of neighbors of each type for each cell type:
```python
'neighbor_dist': {
    'cancer': {
        'cancer': 3.0,   # Average cancer cell has 3.0 cancer neighbors
        'immune': 2.2,   # Average cancer cell has 2.2 immune neighbors
        'stroma': 1.5    # Average cancer cell has 1.5 stromal neighbors
    },
    # ... similar for other cell types
}
```

This captures the local neighborhood composition around each cell type.

### Creating Target Statistics from Data

#### From CSV File:
```python
from tissue_simulator import load_graph_stats_csv

target_stats = load_graph_stats_csv(
    "reference_statistics.csv",
    cell_types=['cancer', 'immune', 'stroma']
)
```

#### From Existing Tissue:
```python
from tissue_simulator import calculate_graph_statistics

# Assuming you have a reference graph with known cell types
reference_coloring = {i: cell.cell_type for i, cell in enumerate(reference_slice.cells)}
target_stats_raw = calculate_graph_statistics(
    reference_graph, 
    reference_coloring, 
    ['cancer', 'immune', 'stroma']
)

# Convert to GraphColorizer format
target_stats = {
    'node_counts': {},
    'edge_counts': {},
    'neighbor_dist': {}
}
# ... (conversion logic)
```

---

## Tuning Simulated Annealing

The simulated annealing algorithm has four key parameters:

### 1. Initial Temperature (`initial_temp`)
- **What it does**: Controls how much random exploration happens at the start
- **Higher values**: More exploration, can escape local minima
- **Lower values**: Less exploration, faster convergence
- **Typical range**: 100 - 2000
- **Default**: 100.0

### 2. Final Temperature (`final_temp`)
- **What it does**: Determines when to stop optimizing
- **Lower values**: More refinement, better final quality
- **Higher values**: Stops earlier, faster runtime
- **Typical range**: 0.001 - 1.0
- **Default**: 0.1

### 3. Cooling Rate (`cooling_rate`)
- **What it does**: How quickly temperature decreases
- **Closer to 1.0**: Slower cooling, more iterations at each temperature
- **Closer to 0.9**: Faster cooling, fewer iterations
- **Typical range**: 0.99 - 0.999
- **Default**: 0.995

### 4. Maximum Iterations (`max_iterations`)
- **What it does**: Hard limit on number of optimization steps
- **Higher values**: Better results but slower
- **Lower values**: Faster but may not converge
- **Typical range**: 5,000 - 50,000
- **Default**: 100,000

### Preset Configurations

#### Fast Testing (5-10 seconds)
```python
colorizer.colorize(
    initial_temp=100.0,
    final_temp=0.1,
    cooling_rate=0.995,
    max_iterations=3000
)
```

#### Balanced (30-60 seconds)
```python
colorizer.colorize(
    initial_temp=500.0,
    final_temp=0.01,
    cooling_rate=0.998,
    max_iterations=10000
)
```

#### High Quality (2-5 minutes)
```python
colorizer.colorize(
    initial_temp=2000.0,
    final_temp=0.001,
    cooling_rate=0.999,
    max_iterations=20000
)
```

### When to Adjust Parameters

**Poor convergence (high JS divergence > 0.3)**:
- Increase `initial_temp` to 1000-2000
- Decrease `cooling_rate` to 0.998-0.999
- Increase `max_iterations` to 15000-30000

**Slow performance**:
- Decrease `max_iterations` to 3000-5000
- Increase `cooling_rate` to 0.99-0.995
- Increase `final_temp` to 0.1-1.0

**Already good results but want refinement**:
- Keep `initial_temp` the same
- Decrease `final_temp` to 0.001-0.01
- Increase `max_iterations` by 50%

---

## Network Building Modes

Two modes for defining spatial connections:

### Contact Mode
Cells are connected if they're touching (intersecting):

```python
graph = analyzer.build_network_from_slice(slicer, mode="contact")
```

**Pros**:
- Captures direct physical contact
- No parameter tuning needed
- Biologically meaningful for contact-dependent interactions

**Cons**:
- Sparse networks (fewer edges)
- Sensitive to cell spacing
- May miss nearby but non-touching cells

**Best for**:
- Adherent cells (epithelial tissues)
- Contact-dependent signaling
- Tight cell packing

### Radius Mode
Cells are connected if their centers are within a distance threshold:

```python
graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=50.0)
```

**Pros**:
- Captures broader spatial relationships
- Adjustable sensitivity via radius parameter
- More robust to spacing variations

**Cons**:
- Requires radius parameter tuning
- May include spurious long-range connections

**Best for**:
- Paracrine signaling
- Diffusible factor interactions
- Sparse tissues

### Choosing the Radius

**Small radius (20-40 μm)**:
- Immediate neighbors only
- Local interactions
- Similar to contact mode but more forgiving

**Medium radius (40-60 μm)**:
- Neighborhood structure
- 2-3 cell diameters
- **Recommended starting point**

**Large radius (60-100 μm)**:
- Broader spatial patterns
- Long-range organization
- Risk of over-connection

**Rule of thumb**: Set radius to 2-3× the average cell diameter.

---

## Evaluation Metrics

The evaluation provides multiple metrics to assess quality:

### 1. Jensen-Shannon Divergence (JS Divergence)
- **Range**: 0 to 1
- **Lower is better**: 0 = perfect match
- **Interpretation**:
  - < 0.1: Excellent match
  - 0.1 - 0.2: Good match
  - 0.2 - 0.3: Fair match
  - > 0.3: Poor match

### 2. Cosine Similarity
- **Range**: 0 to 1
- **Higher is better**: 1 = perfect match
- **Interpretation**:
  - > 0.95: Excellent
  - 0.90 - 0.95: Good
  - 0.85 - 0.90: Fair
  - < 0.85: Poor

### 3. Mean Absolute Error (MAE)
- **Range**: 0 to ∞
- **Lower is better**: 0 = perfect match
- **Interpretation**: Average absolute difference in edge counts

### 4. Root Mean Squared Error (RMSE)
- **Range**: 0 to ∞
- **Lower is better**: 0 = perfect match
- **Interpretation**: Penalizes larger errors more than MAE

### 5. Percent Differences
- **Per edge type**: Percent difference for each cell type pair
- **Average**: Mean percent difference across all edge types
- **Interpretation**:
  - < 10%: Excellent
  - 10-20%: Good
  - 20-30%: Fair
  - > 30%: Poor

### Accessing Metrics

```python
evaluation = evaluate_graph_coloring(target_stats, final_stats)

# Access individual metrics
js_div = evaluation['js_divergence']
cos_sim = evaluation['cosine_similarity']
mae = evaluation['mean_absolute_error']
rmse = evaluation['root_mean_squared_error']
avg_diff = evaluation['avg_percent_difference']
max_diff = evaluation['max_percent_difference']

# Print comprehensive report
print_evaluation_report(target_stats, final_stats, evaluation)
```

---

## Examples

### Example 1: Complete Workflow with TissueNetworkWorkflow

```python
from tissue_simulator import TissueSection, SpherePacker, TissueNetworkWorkflow

# Create tissue
tissue = TissueSection(height=300, width=300, thickness=80)
packer = SpherePacker(tissue)
packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1500)

# Define target statistics
target_stats = {
    'node_counts': {'cancer': 40, 'immune': 35, 'stroma': 25},
    'edge_counts': {
        'cancer-cancer': 60, 'cancer-immune': 45, 'cancer-stroma': 30,
        'immune-immune': 25, 'immune-stroma': 20, 'stroma-stroma': 15
    },
    'neighbor_dist': {
        'cancer': {'cancer': 3.0, 'immune': 2.2, 'stroma': 1.5},
        'immune': {'cancer': 2.5, 'immune': 1.8, 'stroma': 1.3},
        'stroma': {'cancer': 2.0, 'immune': 1.5, 'stroma': 1.2}
    }
}

# Run workflow
workflow = TissueNetworkWorkflow()
evaluation = workflow.run_complete_workflow(
    tissue=tissue,
    z_position=40,
    network_radius=50.0,
    target_stats_dict=target_stats,
    cell_types=['cancer', 'immune', 'stroma'],
    annealing_params={
        'initial_temp': 1000.0,
        'final_temp': 0.01,
        'cooling_rate': 0.998,
        'max_iterations': 10000
    },
    export_dir="results",
    visualize=True
)
```

### Example 2: Batch Processing Multiple Slices

```python
from tissue_simulator import TissueSection, SpherePacker, TissueNetworkWorkflow

# Create tissue
tissue = TissueSection(height=300, width=300, thickness=100)
packer = SpherePacker(tissue)
packer.pack_cells(cell_types={'placeholder': [8, 12]}, max_attempts=1500)

# Process multiple slices
z_positions = [20, 40, 60, 80]
cell_types = ['cancer', 'immune', 'stroma']

for i, z_pos in enumerate(z_positions):
    print(f"Processing slice {i+1}/{len(z_positions)} at z={z_pos}")
    
    workflow = TissueNetworkWorkflow()
    workflow.set_tissue(tissue)
    workflow.create_slice(z_position=z_pos)
    workflow.build_network(mode="radius", radius=50.0)
    workflow.load_target_statistics(statistics=target_stats, cell_types=cell_types)
    workflow.assign_cell_types(max_iterations=5000, verbose=False)
    
    # Export
    workflow.export_all(base_dir=f"results/slice_{i+1}", prefix=f"slice_{i+1}")
    
    # Evaluate
    evaluation = workflow.evaluate(print_report=False)
    print(f"  JS Divergence: {evaluation['js_divergence']:.4f}")
```

### Example 3: Comparing Network Modes

```python
# Same tissue, different network modes
modes = [
    ('contact', None),
    ('radius', 30.0),
    ('radius', 50.0),
    ('radius', 70.0)
]

for mode_type, radius in modes:
    print(f"\nTesting mode: {mode_type}" + (f" (radius={radius})" if radius else ""))
    
    workflow = TissueNetworkWorkflow()
    workflow.set_tissue(tissue)
    workflow.create_slice(z_position=40)
    
    if mode_type == 'contact':
        workflow.build_network(mode="contact")
    else:
        workflow.build_network(mode="radius", radius=radius)
    
    workflow.load_target_statistics(statistics=target_stats, cell_types=cell_types)
    workflow.assign_cell_types(max_iterations=5000, verbose=False)
    
    evaluation = workflow.evaluate(print_report=False)
    print(f"Network edges: {workflow.graph.number_of_edges()}")
    print(f"JS Divergence: {evaluation['js_divergence']:.4f}")
```

---

## Troubleshooting

### Issue: High JS Divergence (> 0.3)

**Possible causes**:
1. Target statistics are not achievable with the given network structure
2. Insufficient optimization (too few iterations)
3. Parameters cooling too quickly

**Solutions**:
- Increase `max_iterations` to 20000-50000
- Increase `initial_temp` to 1000-2000
- Decrease `cooling_rate` to 0.998-0.999
- Check that node counts in target match graph size
- Try different network radius

### Issue: Slow Performance

**Possible causes**:
1. Too many iterations
2. Large network (many nodes/edges)
3. Slow cooling rate

**Solutions**:
- Reduce `max_iterations` to 3000-5000
- Use smaller network (smaller radius, fewer cells)
- Increase `cooling_rate` to 0.99-0.995
- Reduce tissue size

### Issue: Poor Match for Specific Cell Types

**Possible causes**:
1. Target statistics are imbalanced
2. Network structure doesn't support desired patterns

**Solutions**:
- Check that node counts are reasonable given network size
- Verify edge counts are consistent with network connectivity
- Try different network radius to get more/fewer edges
- Adjust target statistics to be more achievable

### Issue: Converges Quickly but Poor Quality

**Possible causes**:
1. Cooling too fast
2. Initial temperature too low
3. Getting stuck in local minimum

**Solutions**:
- Decrease `cooling_rate` to 0.998-0.999
- Increase `initial_temp` to 1000-2000
- Increase `max_iterations`
- Run multiple times and take best result

### Issue: Memory Errors with Large Networks

**Possible causes**:
1. Network too large
2. Statistics calculation too expensive

**Solutions**:
- Process in smaller batches
- Reduce network density (smaller radius)
- Use smaller tissue sections
- Reduce number of nodes

---

## Advanced Topics

### Custom Color Palettes
```python
color_palette = {
    'cancer': '#d62728',   # Red
    'immune': '#1f77b4',   # Blue  
    'stroma': '#2ca02c'    # Green
}

visualize_colored_graph(graph, cell_type_assignment, color_palette=color_palette)
```

### Extracting Statistics from Reference Tissue
```python
from tissue_simulator import calculate_graph_statistics

# Build network from reference tissue
ref_analyzer = SpatialNetworkAnalyzer()
ref_graph = ref_analyzer.build_network_from_slice(ref_slicer, mode="radius", radius=50.0)

# Get existing cell types
ref_coloring = {i: cell.cell_type for i, cell in enumerate(ref_slicer.slice_cells)}

# Calculate statistics
ref_stats = calculate_graph_statistics(ref_graph, ref_coloring, cell_types)

# Use as target for new tissue
# (Convert format as needed)
```

### Parameter Sweep to Find Optimal Settings
```python
param_grid = {
    'initial_temp': [100, 500, 1000, 2000],
    'cooling_rate': [0.995, 0.998, 0.999]
}

best_params = None
best_js_div = float('inf')

for init_temp in param_grid['initial_temp']:
    for cool_rate in param_grid['cooling_rate']:
        # Run workflow
        workflow = TissueNetworkWorkflow()
        # ... setup ...
        workflow.assign_cell_types(
            initial_temp=init_temp,
            cooling_rate=cool_rate,
            max_iterations=5000,
            verbose=False
        )
        
        evaluation = workflow.evaluate(print_report=False)
        js_div = evaluation['js_divergence']
        
        if js_div < best_js_div:
            best_js_div = js_div
            best_params = {'initial_temp': init_temp, 'cooling_rate': cool_rate}

print(f"Best parameters: {best_params}, JS Divergence: {best_js_div:.4f}")
```

---

## References

See also:
- `examples/complete_graph_coloring_workflow.py`: Comprehensive examples
- `tissue_simulator/tissue_workflow.py`: Workflow implementation
- `tissue_simulator/graph_coloring.py`: Simulated annealing algorithm
- `tissue_simulator/evaluation.py`: Evaluation metrics
- `docs/api/graph-coloring.md`: Additional documentation

---

**For questions or issues, please check the documentation or file an issue on the repository.**
