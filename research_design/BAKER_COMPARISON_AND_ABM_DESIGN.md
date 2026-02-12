# Comprehensive Analysis: Baker et al. vs. Tissue Simulator for ABM Spatial Power Analysis

## Executive Summary

This document provides:
1. Detailed analysis of Baker et al.'s approach to spatial power analysis
2. Comparison with the `tissue_simulator` package capabilities
3. Novel experimental design framework for agent-based models (ABMs) with spatial metrics
4. Recommendations for rigorous ABM studies that account for temporal dynamics

---

## Part I: Baker et al. Paper Structure & Methodology

### 1. Problem Definition

**Core Issue Identified:**
- Spatial profiling technologies (CODEX, osmFISH, HDST, Visium) lack statistical frameworks for experimental design
- Unknown: How many samples, FOVs, or cells needed to detect spatial features at given confidence?
- Existing power analysis tools designed for bulk or single-cell (dissociated) data don't account for spatial structure

**Key Parameters Affecting Spatial Experiments:**
1. **Sample size** (number of cells/FOVs)
2. **Cellular composition** (cell type proportions)
3. **FOV size** (area captured per field)
4. **Number of FOVs** (spatial coverage)
5. **Spatial distribution** (how cell types are organized)
6. **Spatial resolution** (single-cell vs. binned data)

### 2. Solution: In Silico Tissue (IST) Framework

**Core Innovation:**
Generate synthetic tissues that recapitulate spatial statistics from real data, enabling power analysis without requiring massive datasets.

**Two-Stage Process:**

**Stage 1: Blank Tissue Scaffold**
- Random circle packing algorithm
- Creates planar graph where nodes = cells, edges = adjacencies
- Similar to your `SpherePacker` but 2D focused

**Stage 2: Cell-Type Label Assignment**

Two methods:
1. **Heuristic assignment** (faster)
   - Grid-based neighborhood sampling
   - Sample cell types from multinomial based on proportions
   - Sample neighbors based on adjacency matrix H
   - Iterative swapping to match targets

2. **Optimization-based assignment**
   - Formulated as inverse optimization problem
   - Minimize ||H(B) - H(B*))||² 
   - Where H(B) is achieved adjacency matrix given assignment B
   - Subject to constraints on cell type proportions
   - Relaxed to continuous optimization via augmented Lagrangian
   - Note: Related to NP-complete graph coloring problems

**Key Input Parameters:**
- **p**: Vector of cell type proportions (K cell types)
- **H**: K×K matrix of pairwise adjacency probabilities
  - H_ij = probability that cell type i is adjacent to cell type j
  
**Regional/Macrostructure Handling:**
- Segment tissue into morphological zones
- Estimate p and H separately for each zone
- Generate tiles for each zone
- Stitch tiles together to create full tissue

### 3. Experiments Demonstrating the Problem

**Experiment 1: Cell Type Detection - Spatial vs. Non-Spatial Sampling**

*Setup:*
- Generated small ISTs (2,186 cells, ~500×500 μm)
- Three spatial configurations:
  1. Random rare cell (3% abundance, randomly distributed)
  2. Self-preference structure (one cell type clusters)
  3. Unstructured null model (random placement)

*Models Developed:*
- **Beta-binomial model**: Cells needed for cell-type detection (accounts for spatial overdispersion)
- **Gamma-Poisson model**: FOVs needed for cell-type detection
- **Binomial model**: Non-spatial single-cell sampling (baseline)

*Key Findings:*
- Spatial structure dramatically impacts sampling requirements
- Breast cancer (unstructured): 80% probability of detecting T cell (3% abundance) requires:
  - 1 FOV at 5% tissue size (~500 cells) with spatial sampling
  - ~100 cells with non-spatial sampling
- Mouse cortex (highly structured, non-repetitive): Detecting L6 pyramidal neurons (9% abundance):
  - 2 FOVs at 5% (~650 cells total)
  - 17 cells with non-spatial sampling
- Mouse spleen (repetitive structures): Detecting megakaryocytes (0.1% abundance):
  - 1 FOV at 5% (~4,300 cells) - 80% probability
  - OR 4 FOVs at 0.5% each (~1,700 cells total)
  - ~1,270 cells with non-spatial sampling

*Critical Insight:*
Power analysis based only on overall cell frequencies vastly underestimates FOV requirements for spatial experiments.

**Experiment 2: Impact of FOV Size on Cell-Cell Adjacency Detection**

*Setup:*
- Generated ISTs with specific enriched adjacencies
- Varied FOV size from 0.5% to 10% of tissue
- Applied permutation test to detect significant adjacencies

*Key Findings:*
- Smaller FOVs less impacted by spatial overdispersion
- FOV size at which inflection occurs reflects the length scale of spatial organization
- For mouse spleen CD4+/CD8+ T cell adjacency:
  - Requires >7.5% tissue size (~123×123 μm, ~5,600 cells) at 80% probability
  - Sharp inflection point indicates macroscale organization
  - TMAs of insufficient size may never capture feature of interest

**Experiment 3: Comparing Tissues/Cohorts via Adjacency Enrichment**

*Setup:*
- Created two mouse spleen tissues:
  1. Original tissue
  2. Modified tissue with 37% reduction in CD4+/CD8+ adjacency (but same cell type proportions)
- Drew FOVs of varying sizes (5%, 7.5%, 10%)
- Calculated Adjacency Enrichment Statistic (AES) for each FOV

*AES Definition:*
```
Σ = 2·f_A·f_B·|E|  (expected edges between type A and B)
AES = (N_AB / Σ) - 1
```
- AES = 0: no enrichment
- AES > 0: enrichment
- AES < 0: depletion

*Key Findings:*
- With 5% FOV: Cannot distinguish tissues (Z-test P = 0.41)
- With 7.5% FOV: Can distinguish (Z-test P = 0.018)
- With 10% FOV: Clearly distinguishable (Z-test P << 0.01)
- Required FOVs to distinguish at P = 0.05: ~1,000, ~100, ~50 respectively

**Experiment 4: Spatial Resolution Impact**

*Setup:*
- Spatially binned CODEX spleen and HDST breast cancer to Visium-like resolution (55 μm spots)
- Applied same power analysis models
- Assumed 10% threshold for cell type deconvolution

*Key Findings:*
- Lower resolution increases sampling requirements:
  - Spleen: 7 FOVs (1% each) vs. 3 FOVs with single-cell resolution
  - Breast cancer: 5 FOVs (1% each) vs. 2 FOVs with single-cell resolution
- Sampling from multiple tissue samples showed no benefit over multiple FOVs from one sample
  - Suggests narrow spatial heterogeneity between samples in their datasets

**Experiment 5: Unknown Features (Cohort-Level Analysis)**

*Setup:*
- Assembled 3 real mouse spleen tissues
- Estimated parameters from 1 tissue
- Generated 20 ISTs with shuffled macrostructures (to avoid spurious patterns)
- Called significant adjacencies (P < 0.01, permutation test)

*Key Findings:*
- 729 possible pairwise adjacencies
- 44 significant in all 20 ISTs
- 50 significant in all 3 real tissues
- 37 of these 50 (84%) overlap with the 44 from ISTs
- Other 13 were boundary artifacts

*Predictive Model:*
To detect an adjacency at 80% probability that occurs in 5/20 ISTs, need ≥6 tissue samples.

### 4. Experimental Design Demonstrated in Paper

The paper's structure follows a classic applied methods paper:

```
1. Introduction → Problem Definition
   ├─ Spatial profiling landscape
   ├─ Existing power analysis limitations
   └─ Unique challenges in spatial experiments

2. Results
   ├─ IST Framework Introduction
   │  ├─ Tissue scaffold generation
   │  ├─ Label assignment (heuristic + optimization)
   │  └─ Regional/macrostructure handling
   │
   ├─ Cell-Type Detection Power Analysis
   │  ├─ Statistical models (beta-binomial, gamma-Poisson)
   │  ├─ Spatial vs. non-spatial sampling
   │  └─ Impact of tissue structure
   │
   ├─ Cell-Cell Adjacency Detection
   │  ├─ Permutation testing framework
   │  ├─ FOV size impact
   │  └─ Length scale of organization
   │
   ├─ Differential Adjacency Testing
   │  ├─ AES statistic development
   │  ├─ Z-test framework
   │  └─ Required sample sizes
   │
   ├─ IST Validation Across Modalities
   │  ├─ HDST breast cancer (unstructured)
   │  ├─ osmFISH mouse cortex (highly structured)
   │  └─ CODEX mouse spleen (repetitive structures)
   │
   └─ Special Considerations
      ├─ Spatial resolution effects
      └─ Unknown features (exploratory analysis)

3. Discussion
   ├─ Framework applicability
   ├─ Limitations
   └─ Future directions

4. Methods
   └─ Detailed algorithms and statistical models
```

---

## Part II: Comparison with Tissue Simulator Package

### Similarities

| Feature | Baker et al. | Tissue Simulator |
|---------|-------------|------------------|
| **Core Algorithm** | 2D random circle packing | 3D random sphere packing |
| **Graph Representation** | Planar graph of adjacencies | Network graph (contact or radius-based) |
| **Label Assignment** | Heuristic + optimization | Simulated annealing optimization |
| **Input Parameters** | Proportions (p) + adjacency matrix (H) | Cell type counts + edge counts + neighbor distributions |
| **Regional Handling** | Tile-based with zone-specific parameters | Not yet implemented |
| **Validation** | Compare IST to real tissue | Evaluation metrics (JS divergence, cosine similarity) |
| **Export** | CSV output | CSV output (3D tissue, 2D slices, network statistics) |

### Key Differences

| Aspect | Baker et al. | Tissue Simulator |
|--------|-------------|------------------|
| **Dimensionality** | 2D focus (with some 3D) | Native 3D with 2D slicing |
| **Time Dynamics** | ❌ Static snapshots only | ✅ Designed for ABM temporal analysis |
| **Primary Use Case** | Power analysis for experimental design | ABM output analysis + synthetic tissue generation |
| **Statistical Framework** | Beta-binomial, gamma-Poisson models | Spatial network statistics |
| **Optimization Method** | Augmented Lagrangian (continuous relaxation) | Simulated annealing (discrete) |
| **Target Statistics** | Adjacency matrix H (probabilities) | Edge counts + neighbor distributions (absolute) |
| **Macrostructure** | Explicit regional annotations | Can be added via workflow |
| **Replicate Generation** | ✅ Core feature (20+ ISTs) | ✅ Recently added feature |
| **ABM Integration** | ❌ Not designed for ABM | ✅ Explicit design goal |

### Tissue Simulator Advantages for ABM Analysis

1. **3D Native Structure**
   - ABMs often simulate 3D space
   - Can slice at any angle for analysis
   - True volumetric analysis

2. **Temporal Dynamics Ready**
   - Can generate tissues at multiple timepoints
   - Track spatial statistics over time
   - Analyze emergent spatial patterns

3. **Network Flexibility**
   - Contact-based or radius-based networks
   - Configurable distance thresholds
   - Multiple network export formats

4. **Comprehensive Spatial Metrics**
   - Global: degree, density, clustering, path lengths
   - Per cell type: degree distributions, centrality
   - Pairwise: interaction strengths, distances
   - Can be extended for temporal tracking

5. **Graph Coloring Integration**
   - Target-based cell type assignment
   - Evaluation framework included
   - Iterative optimization for matching real tissues

### Tissue Simulator Gaps (Opportunities for Development)

1. **Statistical Power Analysis Models**
   - ❌ No beta-binomial / gamma-Poisson models
   - ❌ No predictive sampling requirement calculations
   - ❌ No formal power curves

2. **Regional/Macrostructure Handling**
   - ❌ No built-in segmentation tools
   - ❌ No zone-specific parameter estimation
   - ❌ No automatic tile stitching

3. **Permutation Testing Framework**
   - ❌ No built-in permutation test for adjacencies
   - ❌ No significance calling infrastructure

4. **Cohort-Level Analysis**
   - ❌ No tools for comparing tissue sets
   - ❌ No AES-based statistical testing

5. **Temporal Analysis Framework**
   - ❌ No time series handling
   - ❌ No convergence diagnostics
   - ❌ No transient vs. steady-state separation

---

## Part III: Novel Experimental Design for ABM with Spatial Metrics

### The Core Challenge: Temporal Dynamics in ABM

**Traditional ABM Practice (Problematic):**
```
1. Run simulation for T timesteps
2. Discard first N timesteps (e.g., first 48 hours)
3. Analyze remaining timesteps as "converged/settled"
4. Assumption: Dynamics have reached steady-state
```

**Why This Is Insufficient:**
- Arbitrary choice of N
- No formal convergence criteria
- Loses information about transient dynamics
- Spatial patterns may have different convergence rates
- May miss important emergent behaviors

### New Framework: Spatial-Temporal Power Analysis for ABM

#### 1. Define Spatial Features of Interest

**Categories of Spatial Features:**

**A. Cell-Type Based:**
- Proportions over time
- Spatial distributions (clustered vs. dispersed)
- Migration patterns

**B. Adjacency-Based:**
- Cell-cell contact frequencies
- Neighborhood compositions
- Mixing indices

**C. Higher-Order Structures:**
- Multicellular patterns (e.g., tumor-immune interfaces)
- Spatial domains (regions of similar composition)
- Gradient structures

**D. Network-Based:**
- Degree distributions
- Clustering coefficients
- Community structures
- Network evolution rates

#### 2. Establish Convergence Criteria

**Multi-Metric Convergence Framework:**

Instead of arbitrary burn-in, test convergence of each spatial metric separately:

```python
class SpatialConvergenceAnalyzer:
    """
    Determines when spatial metrics have converged in ABM simulation.
    """
    
    def __init__(self, metrics_of_interest):
        """
        Parameters:
        -----------
        metrics_of_interest : dict
            Dict of {metric_name: metric_function}
            e.g., {'avg_degree': compute_avg_degree,
                   'clustering': compute_clustering,
                   'cd4_cd8_adjacency': compute_cd4_cd8_adjacency}
        """
        self.metrics = metrics_of_interest
        self.time_series = {metric: [] for metric in metrics_of_interest}
        
    def add_timepoint(self, tissue, t):
        """Record all metrics for tissue at time t."""
        for metric_name, metric_func in self.metrics.items():
            value = metric_func(tissue)
            self.time_series[metric_name].append((t, value))
    
    def test_convergence(self, window=100, tolerance=0.05):
        """
        Test if metrics have converged using sliding window.
        
        Criteria:
        - Coefficient of variation < tolerance over window
        - Augmented Dickey-Fuller test (stationarity)
        - Mann-Kendall test (no trend)
        
        Returns:
        --------
        dict : {metric_name: (is_converged, convergence_time)}
        """
        results = {}
        for metric_name, series in self.time_series.items():
            if len(series) < window:
                results[metric_name] = (False, None)
                continue
                
            # Extract recent window
            recent = [val for (t, val) in series[-window:]]
            
            # Test 1: Coefficient of variation
            mean = np.mean(recent)
            std = np.std(recent)
            cv = std / mean if mean != 0 else float('inf')
            
            # Test 2: Stationarity (ADF test)
            from statsmodels.tsa.stattools import adfuller
            adf_result = adfuller(recent)
            is_stationary = adf_result[1] < 0.05  # p-value
            
            # Test 3: No trend (Mann-Kendall)
            from scipy.stats import kendalltau
            times = list(range(len(recent)))
            tau, p_value = kendalltau(times, recent)
            no_trend = p_value > 0.05
            
            # All criteria must pass
            is_converged = (cv < tolerance) and is_stationary and no_trend
            conv_time = series[-window][0] if is_converged else None
            
            results[metric_name] = (is_converged, conv_time)
        
        return results
    
    def get_burn_in_period(self):
        """
        Returns the latest convergence time across all metrics.
        This is the recommended burn-in period.
        """
        convergence = self.test_convergence()
        conv_times = [t for (converged, t) in convergence.values() 
                      if converged and t is not None]
        
        if not conv_times:
            return None  # Not converged yet
        
        return max(conv_times)
```

**Example Usage:**
```python
# Define spatial metrics
def compute_cd4_cd8_adjacency(tissue):
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_tissue(tissue, mode="contact")
    stats = analyzer.compute_pairwise_statistics()
    return stats['cd4']['cd8']['edge_count']

metrics = {
    'avg_degree': lambda t: compute_avg_degree(t),
    'clustering': lambda t: compute_clustering(t),
    'cd4_cd8_adjacency': compute_cd4_cd8_adjacency
}

# Run simulation
convergence = SpatialConvergenceAnalyzer(metrics)

for t in range(max_timesteps):
    tissue = abm.step()  # Run one timestep of ABM
    convergence.add_timepoint(tissue, t)
    
    if t % 10 == 0:  # Check every 10 timesteps
        results = convergence.test_convergence()
        if all(converged for (converged, _) in results.values()):
            print(f"All metrics converged by timestep {t}")
            burn_in = convergence.get_burn_in_period()
            print(f"Recommended burn-in: {burn_in} timesteps")
            break
```

#### 3. Control Initial Conditions

**Problem with Random Initialization:**
- High variance between simulation runs
- Unknown impact on convergence time
- Difficult to compare experimental conditions

**Solution: Controlled Initialization Using Tissue Simulator**

**Strategy A: Match Reference Tissue**
```python
# Generate initial condition matching real tissue spatial statistics
from tissue_simulator import (
    load_target_statistics_from_csv,
    ReplicateGenerator
)

# Load real tissue statistics
target_stats = load_target_statistics_from_csv("patient_biopsy_stats.csv")

# Generate matching initial conditions
generator = ReplicateGenerator(
    target_stats=target_stats,
    tissue_dimensions=(400, 400, 100),
    base_cell_radii={'cancer': (8, 12), 'immune': (5, 8), 'stroma': (10, 15)},
    network_mode="contact"
)

# Generate 10 initial conditions for 10 simulation replicates
initial_conditions = generator.generate_replicates(num_replicates=10)

# Initialize ABM with each initial condition
for i, ic in enumerate(initial_conditions):
    abm = AgentBasedModel()
    abm.initialize_from_tissue(ic)
    abm.run(max_time=1000)
    results[i] = abm.get_final_state()
```

**Strategy B: Systematic Initial Condition Design**
```python
# Design initial conditions spanning a parameter space

initial_conditions = []

# Vary spatial organization
for clustering_level in [0.0, 0.25, 0.5, 0.75, 1.0]:
    # clustering_level controls self-preference in adjacency matrix
    
    # Define target statistics
    H_target = construct_adjacency_matrix(clustering_level)
    target_stats = {
        'node_counts': {'cancer': 100, 'immune': 50},
        'neighbor_dist': H_target
    }
    
    # Generate 5 replicates at this clustering level
    generator = ReplicateGenerator(target_stats=target_stats, ...)
    replicates = generator.generate_replicates(num_replicates=5)
    
    initial_conditions.extend([
        (clustering_level, rep) for rep in replicates
    ])

# Run factorial experiment
for (clustering, initial_tissue) in initial_conditions:
    for treatment in ['control', 'drug_A', 'drug_B']:
        abm = AgentBasedModel()
        abm.initialize_from_tissue(initial_tissue)
        abm.apply_treatment(treatment)
        abm.run(max_time=1000)
        
        # Record results
        results[(clustering, treatment)] = abm.get_final_state()
```

#### 4. Temporal Sampling Strategy

**Problem:** At what timesteps should we measure spatial metrics?

**Solution A: Adaptive Sampling Based on Change Rate**
```python
class AdaptiveSpatialSampler:
    """
    Sample spatial metrics more frequently when they're changing rapidly,
    less frequently when stable.
    """
    
    def __init__(self, initial_interval=1, max_interval=100):
        self.interval = initial_interval
        self.max_interval = max_interval
        self.last_values = {}
        self.last_sample_time = 0
    
    def should_sample(self, t, tissue, metrics):
        """
        Decide whether to sample at time t.
        """
        if t - self.last_sample_time < self.interval:
            return False
        
        # Sample and compare to last values
        current_values = {name: func(tissue) 
                          for name, func in metrics.items()}
        
        if not self.last_values:
            self.last_values = current_values
            self.last_sample_time = t
            return True
        
        # Compute relative change
        max_rel_change = 0
        for name, current in current_values.items():
            last = self.last_values[name]
            rel_change = abs(current - last) / (abs(last) + 1e-10)
            max_rel_change = max(max_rel_change, rel_change)
        
        # Adjust interval based on change rate
        if max_rel_change > 0.1:  # Rapid change
            self.interval = max(1, self.interval // 2)
        elif max_rel_change < 0.01:  # Slow change
            self.interval = min(self.max_interval, self.interval * 2)
        
        self.last_values = current_values
        self.last_sample_time = t
        return True

# Usage in simulation
sampler = AdaptiveSpatialSampler()
spatial_time_series = []

for t in range(max_timesteps):
    tissue = abm.step()
    
    if sampler.should_sample(t, tissue, metrics):
        # Record all spatial metrics
        spatial_time_series.append({
            'time': t,
            **{name: func(tissue) for name, func in metrics.items()}
        })
```

**Solution B: Phase-Based Sampling**
```python
class PhaseDependentSampler:
    """
    Different sampling strategies for different phases of simulation.
    """
    
    def __init__(self):
        self.phases = {
            'initialization': {'duration': 100, 'interval': 1},
            'transient': {'duration': 400, 'interval': 10},
            'equilibration': {'duration': 500, 'interval': 50},
            'steady_state': {'duration': float('inf'), 'interval': 100}
        }
        self.current_phase = 'initialization'
        self.phase_start = 0
    
    def update_phase(self, t):
        """Advance to next phase if duration exceeded."""
        time_in_phase = t - self.phase_start
        phase_info = self.phases[self.current_phase]
        
        if time_in_phase > phase_info['duration']:
            # Advance phase
            phase_order = list(self.phases.keys())
            current_idx = phase_order.index(self.current_phase)
            if current_idx < len(phase_order) - 1:
                self.current_phase = phase_order[current_idx + 1]
                self.phase_start = t
    
    def should_sample(self, t):
        """Sample according to current phase interval."""
        self.update_phase(t)
        interval = self.phases[self.current_phase]['interval']
        return (t - self.phase_start) % interval == 0
```

#### 5. Power Analysis for Temporal Spatial Features

**Key Question:** How many simulation replicates are needed to detect a significant difference in spatial organization between conditions?

**Framework:**

```python
class TemporalSpatialPowerAnalysis:
    """
    Power analysis for ABM studies with spatial metrics evolving over time.
    """
    
    def __init__(self, tissue_generator, abm_class, metrics):
        self.tissue_generator = tissue_generator
        self.abm_class = abm_class
        self.metrics = metrics
    
    def run_pilot_study(self, n_replicates=10, conditions=['control', 'treatment']):
        """
        Run pilot simulations to estimate effect sizes and variance.
        
        Returns:
        --------
        dict : {metric_name: {'effect_size': float, 
                              'pooled_std': float,
                              'time_to_diverge': int}}
        """
        results = {condition: [] for condition in conditions}
        
        for condition in conditions:
            for rep in range(n_replicates):
                # Generate initial condition
                initial_tissue = self.tissue_generator.generate_replicates(1)[0]
                
                # Run ABM
                abm = self.abm_class()
                abm.initialize_from_tissue(initial_tissue)
                abm.set_condition(condition)
                
                # Run with temporal sampling
                time_series = self.run_with_sampling(abm, max_time=1000)
                results[condition].append(time_series)
        
        # Analyze to determine effect sizes
        return self.compute_effect_sizes(results)
    
    def compute_effect_sizes(self, results):
        """
        For each metric, compute Cohen's d at each timepoint.
        """
        effect_sizes = {}
        
        for metric_name in self.metrics.keys():
            # Extract time series for this metric
            control_series = [ts[metric_name] for ts in results['control']]
            treatment_series = [ts[metric_name] for ts in results['treatment']]
            
            # Compute Cohen's d at each timepoint
            timepoints = control_series[0].keys()  # Assuming same sampling
            cohens_d = {}
            
            for t in timepoints:
                control_vals = [series[t] for series in control_series]
                treatment_vals = [series[t] for series in treatment_series]
                
                mean_diff = np.mean(treatment_vals) - np.mean(control_vals)
                pooled_std = np.sqrt(
                    (np.var(control_vals) + np.var(treatment_vals)) / 2
                )
                
                cohens_d[t] = mean_diff / (pooled_std + 1e-10)
            
            # Find when effect becomes detectable (|d| > 0.5)
            time_to_diverge = None
            for t, d in sorted(cohens_d.items()):
                if abs(d) > 0.5:
                    time_to_diverge = t
                    break
            
            # Find maximum effect size (at steady state)
            max_d = max(cohens_d.values(), key=abs)
            
            effect_sizes[metric_name] = {
                'cohens_d_timeline': cohens_d,
                'max_effect_size': max_d,
                'time_to_diverge': time_to_diverge,
                'pooled_std': pooled_std
            }
        
        return effect_sizes
    
    def calculate_required_replicates(self, effect_sizes, alpha=0.05, power=0.8):
        """
        Calculate required sample size for each metric.
        
        Uses standard power analysis:
        n = 2 * (Z_α/2 + Z_β)² * σ² / d²
        
        Where:
        - Z_α/2 = critical value for two-tailed test at significance α
        - Z_β = critical value for power (1-β)
        - σ² = pooled variance
        - d = effect size (difference in means)
        """
        from scipy.stats import norm
        
        z_alpha = norm.ppf(1 - alpha/2)
        z_beta = norm.ppf(power)
        
        required_n = {}
        
        for metric_name, stats in effect_sizes.items():
            d = stats['max_effect_size']
            sigma = stats['pooled_std']
            
            # Standard formula
            n_per_group = 2 * ((z_alpha + z_beta) ** 2) * (sigma ** 2) / (d ** 2)
            n_per_group = int(np.ceil(n_per_group))
            
            required_n[metric_name] = {
                'n_per_group': n_per_group,
                'total_n': 2 * n_per_group,
                'time_to_diverge': stats['time_to_diverge']
            }
        
        return required_n
    
    def run_full_study(self, n_replicates_per_condition):
        """
        Run the full properly-powered study.
        """
        # Similar to pilot but with required sample size
        pass
```

**Example Usage:**
```python
# Set up power analysis
from tissue_simulator import ReplicateGenerator

generator = ReplicateGenerator(...)
power_analyzer = TemporalSpatialPowerAnalysis(
    tissue_generator=generator,
    abm_class=MyABM,
    metrics={
        'avg_degree': compute_avg_degree,
        'cd4_cd8_adjacency': compute_cd4_cd8_adjacency,
        'clustering': compute_clustering
    }
)

# Run pilot study
print("Running pilot study with 10 replicates per condition...")
effect_sizes = power_analyzer.run_pilot_study(n_replicates=10)

# Calculate required sample sizes
required_n = power_analyzer.calculate_required_replicates(
    effect_sizes,
    alpha=0.05,
    power=0.8
)

# Report results
for metric, req in required_n.items():
    print(f"\nMetric: {metric}")
    print(f"  Required replicates per group: {req['n_per_group']}")
    print(f"  Time until divergence: {req['time_to_diverge']} timesteps")
    print(f"  Recommendation: Run {req['n_per_group']} sims per condition")
    print(f"                  for {req['time_to_diverge'] + 500} timesteps")
```

---

## Part IV: Recommended Experimental Design Pipeline for ABM Studies

### Complete Workflow

```
┌─────────────────────────────────────────────────────────────┐
│                    PHASE 1: PILOT STUDY                      │
└─────────────────────────────────────────────────────────────┘

Step 1: Define Research Question
├─ What spatial feature(s) distinguish conditions?
├─ What is the biological mechanism of interest?
└─ What are the relevant timescales?

Step 2: Identify Spatial Metrics
├─ Cell type proportions
├─ Cell-cell adjacencies
├─ Clustering indices
├─ Network properties
└─ Custom domain-specific metrics

Step 3: Generate Realistic Initial Conditions
├─ Option A: Match existing tissue data
│  └─ Use tissue_simulator.ReplicateGenerator with target statistics
├─ Option B: Systematic parameter exploration
│  └─ Generate tissues spanning range of spatial organizations
└─ Generate 3-5 initial conditions per parameter setting

Step 4: Run Small-Scale Simulations (5-10 replicates/condition)
├─ Implement adaptive or phase-based temporal sampling
├─ Monitor convergence of spatial metrics
└─ Record full spatial time series

Step 5: Analyze Pilot Results
├─ Determine burn-in period (convergence analysis)
├─ Calculate effect sizes (Cohen's d) for each metric
├─ Identify time when conditions diverge
└─ Estimate variance components

┌─────────────────────────────────────────────────────────────┐
│              PHASE 2: POWER CALCULATION                      │
└─────────────────────────────────────────────────────────────┘

Step 6: Formal Power Analysis
├─ For each spatial metric:
│  ├─ Effect size (from pilot)
│  ├─ Desired power (typically 0.8)
│  ├─ Significance level (typically 0.05)
│  └─ Calculate required sample size
└─ Take maximum across all metrics of interest

Step 7: Design Full Study
├─ Number of replicates (from power analysis)
├─ Simulation duration (convergence time + observation period)
├─ Temporal sampling strategy
└─ Conditions/treatments to test

┌─────────────────────────────────────────────────────────────┐
│                PHASE 3: FULL STUDY                           │
└─────────────────────────────────────────────────────────────┘

Step 8: Execute Full Study
├─ Generate all initial conditions
├─ Run all simulation replicates
├─ Collect spatial metrics at predetermined timepoints
└─ Export results for analysis

Step 9: Statistical Analysis
├─ Pre-registered analysis plan (avoid p-hacking)
├─ Multiple comparison correction (Bonferroni, FDR)
├─ Effect size reporting (not just p-values)
└─ Visualization of spatial-temporal patterns

Step 10: Validation & Robustness
├─ Sensitivity analysis to initial conditions
├─ Sensitivity to model parameters
├─ Comparison to experimental data (if available)
└─ Cross-validation across different spatial metrics

┌─────────────────────────────────────────────────────────────┐
│            PHASE 4: PUBLICATION & REPORTING                  │
└─────────────────────────────────────────────────────────────┘

Step 11: Comprehensive Reporting
├─ All spatial metrics tracked (not just significant ones)
├─ Full temporal evolution plots
├─ Initial condition distributions
├─ Power analysis results
├─ Effect sizes with confidence intervals
└─ Raw data and analysis code (reproducibility)
```

### Checklist for Rigorous ABM Spatial Studies

**Initial Conditions:**
- [ ] Initial conditions controlled or systematically varied
- [ ] Initial condition generation method documented
- [ ] Multiple initial conditions tested per condition
- [ ] Initial spatial statistics reported

**Convergence:**
- [ ] Convergence criteria defined a priori
- [ ] Multiple spatial metrics tested for convergence
- [ ] Burn-in period determined empirically, not arbitrarily
- [ ] Convergence diagnostics reported for all metrics

**Temporal Sampling:**
- [ ] Sampling strategy justified (adaptive, phase-based, or fixed)
- [ ] Sufficient timepoints to capture dynamics
- [ ] Temporal resolution appropriate for phenomena of interest

**Statistical Power:**
- [ ] Pilot study conducted
- [ ] Effect sizes calculated for all spatial metrics
- [ ] Power analysis performed
- [ ] Required sample size justified

**Analysis:**
- [ ] Pre-registered analysis plan
- [ ] Multiple comparison correction applied
- [ ] Effect sizes reported (not just p-values)
- [ ] Confidence intervals provided
- [ ] Sensitivity analyses conducted

**Reproducibility:**
- [ ] Random seeds documented
- [ ] Initial conditions preserved
- [ ] All parameter values reported
- [ ] Code and data publicly available
- [ ] Computational environment documented

---

## Part V: Implementation Roadmap for Tissue Simulator

To enable the experimental design framework described above, the following features should be added to `tissue_simulator`:

### 1. Temporal Extensions Module

**New Module: `tissue_simulator/temporal.py`**

```python
class TemporalTissueSequence:
    """
    Container for tissue states over time from ABM simulation.
    """
    
    def __init__(self):
        self.timepoints = []
        self.tissues = []
    
    def add_timepoint(self, t, tissue):
        """Add tissue state at time t."""
        pass
    
    def get_spatial_metrics_over_time(self, metrics):
        """Compute spatial metrics for all timepoints."""
        pass
    
    def export_temporal_csv(self, filename):
        """Export time series of spatial metrics."""
        pass
    
    def visualize_evolution(self, metric_name):
        """Plot evolution of specific spatial metric."""
        pass


class ConvergenceAnalyzer:
    """
    Implements convergence testing for spatial metrics.
    """
    
    def __init__(self, metrics, window=100, tolerance=0.05):
        pass
    
    def test_convergence(self, time_series):
        """
        Run multiple convergence tests:
        - Coefficient of variation
        - ADF test (stationarity)
        - Mann-Kendall (trend)
        """
        pass
    
    def get_burn_in_period(self, time_series):
        """Return recommended burn-in based on latest convergence."""
        pass
    
    def plot_convergence_diagnostics(self):
        """Visualize convergence for all metrics."""
        pass
```

### 2. Power Analysis Module

**New Module: `tissue_simulator/power_analysis.py`**

```python
class SpatialPowerAnalyzer:
    """
    Power analysis for spatial experiments (static and temporal).
    """
    
    def __init__(self, metrics):
        pass
    
    def estimate_effect_sizes(self, condition1_tissues, condition2_tissues):
        """
        Compute Cohen's d for spatial metrics between conditions.
        """
        pass
    
    def calculate_required_replicates(self, effect_size, alpha, power):
        """
        Standard power calculation for spatial metrics.
        """
        pass
    
    def run_pilot_study(self, generator, abm_class, conditions, n_pilot=10):
        """
        Automated pilot study to estimate parameters.
        """
        pass
    
    def generate_power_curves(self, effect_sizes):
        """
        Plot power as function of sample size for each metric.
        """
        pass


class TemporalPowerAnalyzer(SpatialPowerAnalyzer):
    """
    Extended power analysis for temporal spatial data.
    """
    
    def estimate_temporal_effect_sizes(self, condition1_sequences, 
                                       condition2_sequences):
        """
        Compute effect sizes at each timepoint.
        """
        pass
    
    def identify_divergence_time(self, effect_size_timeline, threshold=0.5):
        """
        Find when conditions begin to differ significantly.
        """
        pass
```

### 3. Statistical Testing Module

**New Module: `tissue_simulator/statistical_tests.py`**

```python
class PermutationTest:
    """
    Permutation testing for spatial features.
    (Similar to Baker et al.)
    """
    
    def __init__(self, n_permutations=1000):
        pass
    
    def test_adjacency_enrichment(self, tissue, cell_type_a, cell_type_b):
        """
        Test if A-B adjacency is enriched vs. null (shuffled labels).
        """
        pass
    
    def test_all_adjacencies(self, tissue):
        """
        Test all pairwise adjacencies, return P-values.
        """
        pass
    
    def correct_multiple_comparisons(self, p_values, method='fdr'):
        """
        Apply Bonferroni, FDR, or other corrections.
        """
        pass


class AdjacencyEnrichmentStatistic:
    """
    Implements AES from Baker et al.
    """
    
    @staticmethod
    def compute_aes(tissue, cell_type_a, cell_type_b):
        """
        AES = (N_AB / Σ) - 1
        where Σ = 2·f_A·f_B·|E|
        """
        pass
    
    @staticmethod
    def compare_tissues_ztest(tissues_1, tissues_2, cell_type_a, cell_type_b):
        """
        Z-test comparing AES distributions between tissue sets.
        """
        pass
```

### 4. Regional/Macrostructure Module

**New Module: `tissue_simulator/regions.py`**

```python
class RegionalTissue:
    """
    Tissue with multiple distinct spatial regions.
    """
    
    def __init__(self, tissue, region_annotations):
        """
        Parameters:
        -----------
        tissue : TissueSection
        region_annotations : dict
            Maps (x,y,z) ranges to region labels
        """
        pass
    
    def segment_by_regions(self):
        """Split tissue into sub-tissues by region."""
        pass
    
    def estimate_regional_parameters(self):
        """
        Estimate p and H (or equivalent) for each region separately.
        """
        pass
    
    def generate_regional_replicate(self):
        """
        Generate new tissue with same regional structure but
        new random instantiation of each region.
        """
        pass


class TileStitcher:
    """
    Stitch together regional tiles into full tissue.
    """
    
    def __init__(self, region_tiles):
        pass
    
    def stitch(self, layout):
        """
        Combine tiles according to layout specification.
        """
        pass
```

### 5. ABM Integration Module

**New Module: `tissue_simulator/abm_integration.py`**

```python
class ABMTissueInterface:
    """
    Interface between tissue_simulator and ABM frameworks.
    """
    
    @staticmethod
    def tissue_to_abm_initial_condition(tissue, abm_class):
        """
        Convert TissueSection to ABM initial state.
        """
        pass
    
    @staticmethod
    def abm_to_tissue(abm_state):
        """
        Convert ABM state to TissueSection for analysis.
        """
        pass
    
    @staticmethod
    def run_abm_with_temporal_sampling(abm, sampler, max_time):
        """
        Run ABM and collect tissue snapshots based on sampling strategy.
        
        Returns:
        --------
        TemporalTissueSequence
        """
        pass


class AdaptiveSampler:
    """
    Adaptive temporal sampling based on rate of spatial change.
    """
    
    def __init__(self, metrics, initial_interval=1, max_interval=100):
        pass
    
    def should_sample(self, t, current_tissue):
        """Decide whether to sample at time t."""
        pass


class PhaseBasedSampler:
    """
    Different sampling rates for different simulation phases.
    """
    
    def __init__(self, phases):
        """
        Parameters:
        -----------
        phases : dict
            {phase_name: {'duration': int, 'interval': int}}
        """
        pass
    
    def should_sample(self, t):
        """Sample based on current phase."""
        pass
```

### 6. Documentation & Examples

**New Tutorial: `docs/ABM_EXPERIMENTAL_DESIGN.md`**
- Complete guide to designing rigorous ABM studies
- Worked examples with code
- Interpretation of results
- Common pitfalls

**New Examples:**
- `examples/abm_power_analysis.py` - Full pipeline demo
- `examples/convergence_analysis.py` - Temporal convergence testing
- `examples/controlled_initial_conditions.py` - Systematic IC generation
- `examples/regional_tissue_generation.py` - Macrostructure handling

---

## Part VI: Key Takeaways & Recommendations

### For Rigorous ABM Studies:

1. **Don't Use Arbitrary Burn-In**
   - Implement formal convergence testing
   - Test each spatial metric separately
   - Report convergence diagnostics

2. **Control Initial Conditions**
   - Generate realistic ICs matching real tissue
   - Or systematically vary ICs to test robustness
   - Always report IC distribution

3. **Conduct Pilot Studies**
   - Estimate effect sizes before full study
   - Use pilot to determine required sample size
   - Don't rely on intuition for power

4. **Temporal Sampling Strategy**
   - Match sampling to rate of spatial change
   - Don't oversample (computational waste)
   - Don't undersample (miss dynamics)

5. **Report Comprehensively**
   - All spatial metrics, not just significant ones
   - Effect sizes with confidence intervals
   - Full temporal evolution
   - Code and data for reproducibility

### For Tissue Simulator Development:

1. **Priority Additions:**
   - Temporal tissue sequence handling
   - Convergence analysis tools
   - Power analysis framework
   - Statistical testing suite

2. **Integration with Existing Tools:**
   - Current replicate generation is excellent foundation
   - Extend to temporal case
   - Add formal statistical framework

3. **Differentiation from Baker et al.:**
   - Focus on ABM integration (they don't have this)
   - Temporal dynamics (they're static)
   - 3D native with slicing (they're mostly 2D)

### For Publications:

**Proposed Papers:**

1. **Methods Paper:** "Statistical Power Analysis for Agent-Based Models with Spatial Structure"
   - Novel temporal-spatial power framework
   - Convergence testing methodology
   - Implementation in tissue_simulator
   - Validation on 2-3 ABM systems

2. **Application Paper:** "Controlled Initial Conditions Improve Reproducibility of Spatial Agent-Based Models"
   - Demonstrate IC impact on results
   - Systematic IC design framework
   - Comparison across ABM platforms

3. **Software Paper:** "tissue_simulator: A Python Package for Spatial Analysis of Agent-Based Models"
   - Full package description
   - Comparison to Baker et al. and other tools
   - Tutorial for ABM researchers

---

## Appendix: Code Integration Examples

### A. Integrating Convergence Analysis

```python
# In tissue_simulator/temporal.py

from scipy import stats
import numpy as np
from statsmodels.tsa.stattools import adfuller

class SpatialConvergenceAnalyzer:
    """
    Test convergence of spatial metrics in ABM simulations.
    """
    
    def __init__(self, metrics, window=100, cv_threshold=0.05, 
                 adf_alpha=0.05, mk_alpha=0.05):
        """
        Parameters:
        -----------
        metrics : dict
            {metric_name: metric_function} where metric_function 
            takes a TissueSection and returns float
        window : int
            Number of recent timepoints to test for convergence
        cv_threshold : float
            Maximum coefficient of variation for convergence
        adf_alpha : float
            Significance level for ADF stationarity test
        mk_alpha : float
            Significance level for Mann-Kendall trend test
        """
        self.metrics = metrics
        self.window = window
        self.cv_threshold = cv_threshold
        self.adf_alpha = adf_alpha
        self.mk_alpha = mk_alpha
        
        # Storage for time series
        self.time_series = {name: [] for name in metrics.keys()}
    
    def add_timepoint(self, t, tissue):
        """
        Compute and store all metrics for tissue at time t.
        
        Parameters:
        -----------
        t : int
            Timestep
        tissue : TissueSection
            Current tissue state
        """
        for metric_name, metric_func in self.metrics.items():
            value = metric_func(tissue)
            self.time_series[metric_name].append((t, value))
    
    def test_convergence(self):
        """
        Test if all metrics have converged using multiple criteria.
        
        Returns:
        --------
        dict : {metric_name: {
                   'is_converged': bool,
                   'convergence_time': int or None,
                   'cv': float,
                   'adf_pvalue': float,
                   'mk_pvalue': float
               }}
        """
        results = {}
        
        for metric_name, series in self.time_series.items():
            if len(series) < self.window:
                results[metric_name] = {
                    'is_converged': False,
                    'convergence_time': None,
                    'cv': None,
                    'adf_pvalue': None,
                    'mk_pvalue': None,
                    'reason': f'Insufficient data (need {self.window} points)'
                }
                continue
            
            # Extract recent window
            recent = [val for (t, val) in series[-self.window:]]
            times = [t for (t, val) in series[-self.window:]]
            
            # Test 1: Coefficient of Variation
            mean = np.mean(recent)
            std = np.std(recent)
            cv = std / abs(mean) if mean != 0 else float('inf')
            cv_pass = cv < self.cv_threshold
            
            # Test 2: Augmented Dickey-Fuller (stationarity)
            try:
                adf_result = adfuller(recent, autolag='AIC')
                adf_pvalue = adf_result[1]
                adf_pass = adf_pvalue < self.adf_alpha
            except Exception as e:
                adf_pvalue = None
                adf_pass = False
            
            # Test 3: Mann-Kendall (no trend)
            try:
                tau, mk_pvalue = stats.kendalltau(range(len(recent)), recent)
                mk_pass = mk_pvalue > self.mk_alpha
            except Exception as e:
                mk_pvalue = None
                mk_pass = False
            
            # All criteria must pass
            is_converged = cv_pass and adf_pass and mk_pass
            
            # If converged, convergence time is start of current window
            conv_time = times[0] if is_converged else None
            
            results[metric_name] = {
                'is_converged': is_converged,
                'convergence_time': conv_time,
                'cv': cv,
                'cv_pass': cv_pass,
                'adf_pvalue': adf_pvalue,
                'adf_pass': adf_pass,
                'mk_pvalue': mk_pvalue,
                'mk_pass': mk_pass
            }
        
        return results
    
    def get_burn_in_period(self):
        """
        Return the latest convergence time across all metrics.
        This is the recommended burn-in period.
        
        Returns:
        --------
        int or None : Latest convergence time, or None if not all converged
        """
        conv_results = self.test_convergence()
        
        conv_times = [
            result['convergence_time'] 
            for result in conv_results.values() 
            if result['is_converged'] and result['convergence_time'] is not None
        ]
        
        if len(conv_times) < len(self.metrics):
            return None  # Not all metrics converged
        
        return max(conv_times)
    
    def plot_convergence_diagnostics(self, save_path=None):
        """
        Create diagnostic plots for convergence analysis.
        
        Parameters:
        -----------
        save_path : str, optional
            Path to save figure
        """
        import matplotlib.pyplot as plt
        
        n_metrics = len(self.metrics)
        fig, axes = plt.subplots(n_metrics, 2, figsize=(12, 4*n_metrics))
        
        if n_metrics == 1:
            axes = axes.reshape(1, -1)
        
        conv_results = self.test_convergence()
        
        for idx, (metric_name, series) in enumerate(self.time_series.items()):
            times = [t for (t, val) in series]
            values = [val for (t, val) in series]
            
            # Left plot: Time series with convergence marker
            ax_left = axes[idx, 0]
            ax_left.plot(times, values, 'b-', alpha=0.7)
            ax_left.set_xlabel('Time')
            ax_left.set_ylabel(metric_name)
            ax_left.set_title(f'{metric_name} - Time Series')
            ax_left.grid(True, alpha=0.3)
            
            # Mark convergence time if available
            result = conv_results[metric_name]
            if result['is_converged']:
                conv_t = result['convergence_time']
                ax_left.axvline(conv_t, color='g', linestyle='--', 
                               label=f'Converged at t={conv_t}')
                ax_left.legend()
            
            # Right plot: Rolling statistics
            ax_right = axes[idx, 1]
            
            # Compute rolling CV
            window = self.window
            if len(values) >= window:
                rolling_cv = []
                rolling_times = []
                for i in range(window, len(values)):
                    window_vals = values[i-window:i]
                    mean = np.mean(window_vals)
                    std = np.std(window_vals)
                    cv = std / abs(mean) if mean != 0 else float('inf')
                    rolling_cv.append(cv)
                    rolling_times.append(times[i])
                
                ax_right.plot(rolling_times, rolling_cv, 'r-', alpha=0.7)
                ax_right.axhline(self.cv_threshold, color='g', linestyle='--',
                                label=f'CV threshold = {self.cv_threshold}')
                ax_right.set_xlabel('Time')
                ax_right.set_ylabel('Rolling CV')
                ax_right.set_title(f'{metric_name} - Coefficient of Variation')
                ax_right.legend()
                ax_right.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def export_time_series_csv(self, filepath):
        """
        Export all metric time series to CSV.
        
        Parameters:
        -----------
        filepath : str
            Path to output CSV file
        """
        import pandas as pd
        
        # Build DataFrame
        all_data = {'time': []}
        for metric_name in self.metrics.keys():
            all_data[metric_name] = []
        
        # Get all unique timepoints
        all_times = set()
        for series in self.time_series.values():
            all_times.update(t for (t, val) in series)
        all_times = sorted(all_times)
        
        # Fill in data
        for t in all_times:
            all_data['time'].append(t)
            for metric_name, series in self.time_series.items():
                # Find value at time t
                val = None
                for (time, value) in series:
                    if time == t:
                        val = value
                        break
                all_data[metric_name].append(val)
        
        df = pd.DataFrame(all_data)
        df.to_csv(filepath, index=False)
        print(f"Time series data exported to {filepath}")


# Example usage:
if __name__ == "__main__":
    from tissue_simulator import TissueSection, SpatialNetworkAnalyzer
    
    # Define spatial metrics
    def compute_avg_degree(tissue):
        analyzer = SpatialNetworkAnalyzer()
        analyzer.build_network_from_tissue(tissue, mode="contact")
        stats = analyzer.compute_global_statistics()
        return stats.avg_degree
    
    def compute_clustering(tissue):
        analyzer = SpatialNetworkAnalyzer()
        analyzer.build_network_from_tissue(tissue, mode="contact")
        stats = analyzer.compute_global_statistics()
        return stats.avg_clustering
    
    metrics = {
        'avg_degree': compute_avg_degree,
        'clustering': compute_clustering
    }
    
    # Set up convergence analyzer
    analyzer = SpatialConvergenceAnalyzer(metrics, window=50)
    
    # Simulate ABM (placeholder)
    for t in range(500):
        # In real usage, this would be: tissue = abm.step()
        tissue = TissueSection(400, 400, 100, cell_radii=(8, 12))
        tissue.generate_cells(max_attempts=1000)
        
        analyzer.add_timepoint(t, tissue)
        
        # Check convergence every 10 steps
        if t % 10 == 0 and t >= 50:
            results = analyzer.test_convergence()
            all_converged = all(r['is_converged'] for r in results.values())
            
            if all_converged:
                burn_in = analyzer.get_burn_in_period()
                print(f"\nAll metrics converged by timestep {t}")
                print(f"Recommended burn-in period: {burn_in} timesteps\n")
                break
            else:
                not_converged = [name for name, r in results.items() 
                               if not r['is_converged']]
                print(f"t={t}: Not converged: {not_converged}")
    
    # Export and visualize
    analyzer.export_time_series_csv("convergence_time_series.csv")
    analyzer.plot_convergence_diagnostics("convergence_diagnostics.png")
```

This comprehensive document provides:
1. Deep analysis of Baker et al.'s methodology
2. Comparison with your package
3. Novel framework for ABM experimental design
4. Concrete implementation examples
5. Roadmap for package development

Would you like me to:
1. Implement any of these modules?
2. Create specific example scripts?
3. Design experiments for validation?
4. Draft a methods paper outline?
