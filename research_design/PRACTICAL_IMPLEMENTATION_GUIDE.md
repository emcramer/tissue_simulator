# Practical Implementation Guide: From Baker et al. to Tissue Simulator

## Quick Reference: Baker et al. Paper Structure

### Paper Organization

The Baker et al. paper follows a classical applied methods structure that we should emulate:

**1. Introduction (Problem Setup)**
- Current state of spatial omics technologies
- Gap: lack of power analysis frameworks
- Why existing tools don't work (designed for bulk/single-cell, not spatial)
- Preview of solution (IST framework)

**2. Results (Core Contribution)**

**Section 2.1: IST Framework**
- How it works (circle packing → label assignment)
- Two algorithms (heuristic vs. optimization)
- Handling complex tissues (regional parameters)

**Section 2.2: Cell-Type Detection**
- Problem: How many cells/FOVs needed to detect rare cell type?
- Three statistical models developed
- Validation across 3 tissue types with different structures

**Section 2.3: Cell-Cell Adjacency Detection**
- Problem: Sampling requirements for detecting enriched adjacencies
- Permutation test framework
- FOV size impact analysis

**Section 2.4: Comparing Tissues/Cohorts**
- Problem: Distinguishing tissues by spatial organization
- AES statistic development
- Z-test framework for differential testing

**Section 2.5: Resolution Effects**
- Problem: How does spatial resolution impact sampling?
- Visium-like binning analysis
- Increased requirements at lower resolution

**Section 2.6: Unknown Features**
- Problem: Exploratory studies without pre-specified features
- Cohort-level analysis (3 real + 20 ISTs)
- Predicting required cohort size

**3. Discussion**
- Framework applicability
- Comparison to existing approaches
- Limitations and future work

**4. Methods**
- Detailed algorithms
- Statistical model derivations
- Parameter estimation procedures

---

## Experiments That Demonstrated the Problem

### Experiment Design Pattern

All Baker et al. experiments follow this pattern:

```
1. Generate synthetic tissues (ISTs) with known ground truth
2. Sample from these tissues in various ways
3. Measure: What fraction of samples detect the feature?
4. Compare: Does this match theoretical predictions?
5. Validate: Does this hold in real data?
```

### Key Experiments

**Experiment 1: The Spatial Structure Problem**

*Hypothesis:* Spatial structure impacts sampling requirements

*Design:*
- Generate 3 tissue types:
  - Unstructured (breast cancer)
  - Highly structured (brain cortex)
  - Repetitive structure (spleen)
- Sample with varying FOV sizes and counts
- Measure probability of detecting rare cell type

*Key Finding:*
- Same cell type abundance (e.g., 3%) requires:
  - Unstructured: 1 FOV (5% area)
  - Structured: 2 FOVs (5% area each)
  - Repetitive: 1 FOV (5% area) but different considerations

*Implication:*
Power analysis that ignores spatial structure underestimates requirements by orders of magnitude.

**Experiment 2: The FOV Size Problem**

*Hypothesis:* FOV size affects which spatial features can be detected

*Design:*
- Generate mouse spleen ISTs
- Vary FOV size (0.5%, 1%, 5%, 7.5%, 10% of tissue)
- Measure probability of detecting CD4+/CD8+ T cell adjacency

*Key Finding:*
- Sharp inflection point at ~7.5% tissue size
- Below this: very low probability
- Above this: high probability
- Inflection reflects length scale of spatial organization

*Implication:*
TMAs must be sized appropriately for the spatial feature of interest, or they will never capture it.

**Experiment 3: The Cohort Comparison Problem**

*Hypothesis:* FOV size affects ability to distinguish cohorts

*Design:*
- Create 2 mouse spleens differing in one adjacency (37% reduction in CD4+/CD8+)
- Draw many FOVs of different sizes
- Calculate AES for each FOV
- Test if distributions differ (Z-test)

*Key Finding:*
- 5% FOV: Cannot distinguish (P = 0.41)
- 7.5% FOV: Can distinguish (P = 0.018)
- 10% FOV: Clearly different (P << 0.01)

*Implication:*
Required sample size depends critically on FOV size—can reduce from 1000s to 10s of samples.

---

## How to Structure Our Papers

### Paper 1: Methods Paper on ABM Spatial Power Analysis

**Title:** "Statistical Power Analysis for Agent-Based Models with Evolving Spatial Structure"

**Abstract Structure:**
```
[Problem] Agent-based models (ABMs) are increasingly used to study spatial 
tissue dynamics, but rigorous experimental design principles are lacking. 
Current practice relies on arbitrary burn-in periods and sample sizes.

[Gap] Spatial structure evolves over time in ABMs, requiring frameworks that 
account for both spatial organization and temporal dynamics.

[Solution] We introduce a statistical framework for power analysis of ABM 
studies with spatial metrics, including convergence testing, controlled 
initial condition generation, and temporal effect size estimation.

[Validation] We demonstrate our framework on [2-3 ABM systems], showing that:
(1) Convergence times vary dramatically across spatial metrics
(2) Controlled initial conditions reduce required replicates by 40%
(3) Temporal sampling strategy affects power significantly

[Impact] Our framework and open-source implementation (tissue_simulator) 
enable rigorous, reproducible ABM studies.
```

**Results Section Structure:**

**Result 1: Spatial Metrics Converge at Different Rates**

*Setup:*
- Run [ABM system 1] with standard parameters
- Track 5 spatial metrics over time
- Apply convergence testing framework

*Finding:*
- Cell-type proportions converge quickly (t=50)
- Network clustering converges slowly (t=300)
- Adjacency patterns slowest (t=500)

*Implication:*
Cannot use single burn-in period—must test each metric.

*Figure:*
- Panel A: Time series of all metrics
- Panel B: Rolling CV for each metric
- Panel C: Convergence diagnostic statistics
- Panel D: Recommended burn-in for each metric

**Result 2: Controlled Initial Conditions Reduce Variance**

*Setup:*
- Generate ABM initial conditions three ways:
  1. Random placement (current standard)
  2. Match reference tissue
  3. Systematic parameter exploration
- Run 20 replicates each
- Measure variance in spatial metrics at t=final

*Finding:*
- Random: High variance (CV = 0.45)
- Matched: Medium variance (CV = 0.22)
- Systematic: Low variance (CV = 0.18)

*Implication:*
Controlled ICs can reduce required replicates by ~40% for equivalent power.

*Figure:*
- Panel A: Example initial conditions from each method
- Panel B: Final spatial organization from each method
- Panel C: Distribution of spatial metrics across replicates
- Panel D: Power curves for each IC method

**Result 3: Temporal Effect Sizes Guide Study Design**

*Setup:*
- Run pilot study with [ABM system 2], two conditions (control vs treatment)
- N=10 replicates per condition
- Compute Cohen's d for spatial metrics at each timepoint

*Finding:*
- Conditions diverge at different times for different metrics:
  - Cell proportions: t=100 (d=0.6)
  - Adjacency patterns: t=300 (d=0.8)
  - Network clustering: t=500 (d=1.2)

*Implication:*
Can optimize study duration to capture effect of interest, avoiding unnecessarily long simulations.

*Figure:*
- Panel A: Time series of spatial metrics for both conditions
- Panel B: Cohen's d over time for each metric
- Panel C: Required sample size vs. observation time
- Panel D: Power curves for different observation windows

**Result 4: Framework Generalizes Across ABM Platforms**

*Setup:*
- Apply framework to 3 different ABM systems:
  1. Lattice-based (cellular automaton)
  2. Agent-based (off-lattice)
  3. Hybrid (continuum + discrete)

*Finding:*
- Framework applies to all
- Specific convergence times differ
- Power analysis principles identical

*Implication:*
Framework is platform-agnostic, can be adopted widely.

*Figure:*
- Panel A: Example tissues from each ABM type
- Panel B: Convergence diagnostics for each
- Panel C: Power curves for spatial metrics
- Panel D: Required sample sizes comparison

**Result 5: Open-Source Implementation**

*Description:*
- tissue_simulator package with temporal extensions
- Full workflow: IC generation → convergence testing → power analysis
- Integration with common ABM frameworks
- Documentation and tutorials

*Figure/Table:*
- Table 1: Feature comparison with existing tools
- Table 2: Performance benchmarks
- Supp Figure: Software architecture

---

## Implementing the Experiments

### Experiment Template for Our Papers

```python
"""
Template for running experiments following Baker et al. design pattern.
"""

from tissue_simulator import TissueSection, ReplicateGenerator
from tissue_simulator.temporal import SpatialConvergenceAnalyzer
from tissue_simulator.power_analysis import TemporalPowerAnalyzer
import numpy as np
import matplotlib.pyplot as plt

class ExperimentTemplate:
    """
    Base class for experiments demonstrating ABM spatial power analysis.
    """
    
    def __init__(self, name):
        self.name = name
        self.results = {}
    
    def generate_ground_truth(self):
        """
        Create synthetic data with known properties.
        Override in subclasses.
        """
        raise NotImplementedError
    
    def run_sampling_experiment(self):
        """
        Sample from ground truth in various ways.
        Override in subclasses.
        """
        raise NotImplementedError
    
    def analyze_results(self):
        """
        Compute statistics and test hypotheses.
        Override in subclasses.
        """
        raise NotImplementedError
    
    def visualize_results(self, save_dir):
        """
        Create figures for publication.
        Override in subclasses.
        """
        raise NotImplementedError
    
    def run_full_experiment(self, save_dir):
        """
        Execute complete experiment pipeline.
        """
        print(f"Running experiment: {self.name}")
        
        print("Step 1: Generating ground truth...")
        self.generate_ground_truth()
        
        print("Step 2: Running sampling experiments...")
        self.run_sampling_experiment()
        
        print("Step 3: Analyzing results...")
        self.analyze_results()
        
        print("Step 4: Creating visualizations...")
        self.visualize_results(save_dir)
        
        print(f"Experiment complete! Results saved to {save_dir}")
        
        return self.results


class ConvergenceRateExperiment(ExperimentTemplate):
    """
    Experiment 1: Demonstrate that spatial metrics converge at different rates.
    
    This addresses the problem of arbitrary burn-in periods in ABM studies.
    """
    
    def __init__(self):
        super().__init__("Differential Convergence Rates")
        self.metrics = None
        self.abm = None
        self.convergence_analyzer = None
    
    def generate_ground_truth(self):
        """
        Set up ABM and metrics to track.
        """
        # Import ABM (placeholder - replace with actual ABM)
        # from my_abm import TissueABM
        # self.abm = TissueABM(...)
        
        # Define spatial metrics to track
        from tissue_simulator import SpatialNetworkAnalyzer
        
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
        
        def compute_cell_type_prop(tissue):
            # Compute proportion of cell type A
            type_a_count = sum(1 for cell in tissue.cells 
                              if cell.cell_type == 'type_a')
            return type_a_count / len(tissue.cells)
        
        def compute_adjacency(tissue):
            analyzer = SpatialNetworkAnalyzer()
            analyzer.build_network_from_tissue(tissue, mode="contact")
            stats = analyzer.compute_pairwise_statistics()
            if 'type_a' in stats and 'type_b' in stats['type_a']:
                return stats['type_a']['type_b']['edge_count']
            return 0
        
        self.metrics = {
            'cell_type_proportion': compute_cell_type_prop,
            'avg_degree': compute_avg_degree,
            'clustering': compute_clustering,
            'type_a_type_b_adjacency': compute_adjacency
        }
    
    def run_sampling_experiment(self):
        """
        Run ABM and track convergence of each metric.
        """
        from tissue_simulator.temporal import SpatialConvergenceAnalyzer
        
        self.convergence_analyzer = SpatialConvergenceAnalyzer(
            metrics=self.metrics,
            window=100,
            cv_threshold=0.05
        )
        
        max_timesteps = 1000
        for t in range(max_timesteps):
            # Run one timestep of ABM
            # tissue = self.abm.step()
            
            # For now, generate random tissue as placeholder
            tissue = TissueSection(
                height=400, width=400, thickness=100,
                cell_radii={'type_a': (8, 12), 'type_b': (8, 12)}
            )
            tissue.generate_cells(max_attempts=1000)
            
            # Record spatial metrics
            self.convergence_analyzer.add_timepoint(t, tissue)
            
            # Check convergence every 50 steps
            if t % 50 == 0 and t >= 100:
                results = self.convergence_analyzer.test_convergence()
                all_converged = all(r['is_converged'] 
                                  for r in results.values())
                
                if all_converged:
                    print(f"All metrics converged by timestep {t}")
                    break
        
        self.results['convergence_results'] = \
            self.convergence_analyzer.test_convergence()
        self.results['burn_in_period'] = \
            self.convergence_analyzer.get_burn_in_period()
    
    def analyze_results(self):
        """
        Compute statistics on convergence times.
        """
        conv_results = self.results['convergence_results']
        
        # Extract convergence times
        conv_times = {
            name: result['convergence_time']
            for name, result in conv_results.items()
            if result['is_converged']
        }
        
        # Summary statistics
        self.results['summary'] = {
            'fastest_convergence': min(conv_times.items(), 
                                      key=lambda x: x[1]),
            'slowest_convergence': max(conv_times.items(), 
                                      key=lambda x: x[1]),
            'mean_convergence_time': np.mean(list(conv_times.values())),
            'range_convergence_times': (min(conv_times.values()), 
                                       max(conv_times.values()))
        }
        
        print("\n=== Convergence Analysis Summary ===")
        print(f"Burn-in period: {self.results['burn_in_period']} timesteps")
        print(f"Fastest: {self.results['summary']['fastest_convergence']}")
        print(f"Slowest: {self.results['summary']['slowest_convergence']}")
        print(f"Range: {self.results['summary']['range_convergence_times']}")
    
    def visualize_results(self, save_dir):
        """
        Create publication-quality figures.
        """
        import os
        os.makedirs(save_dir, exist_ok=True)
        
        # Main convergence diagnostics figure
        fig = self.convergence_analyzer.plot_convergence_diagnostics(
            save_path=os.path.join(save_dir, "convergence_diagnostics.png")
        )
        plt.close(fig)
        
        # Export time series data
        self.convergence_analyzer.export_time_series_csv(
            os.path.join(save_dir, "convergence_time_series.csv")
        )
        
        # Create summary table
        import pandas as pd
        conv_results = self.results['convergence_results']
        
        summary_data = []
        for metric_name, result in conv_results.items():
            summary_data.append({
                'Metric': metric_name,
                'Converged': result['is_converged'],
                'Convergence Time': result['convergence_time'],
                'Final CV': result['cv'],
                'ADF p-value': result['adf_pvalue'],
                'MK p-value': result['mk_pvalue']
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv(os.path.join(save_dir, "convergence_summary.csv"), 
                 index=False)
        
        print(f"\nResults saved to {save_dir}")


class InitialConditionExperiment(ExperimentTemplate):
    """
    Experiment 2: Demonstrate impact of controlled initial conditions.
    
    This addresses the problem of high variance in ABM studies.
    """
    
    def __init__(self):
        super().__init__("Controlled Initial Conditions")
        self.ic_methods = ['random', 'matched', 'systematic']
        self.n_replicates = 20
    
    def generate_ground_truth(self):
        """
        Create reference tissue to match.
        """
        # Generate a reference tissue with known spatial properties
        self.reference_tissue = TissueSection(
            height=400, width=400, thickness=100,
            cell_radii={'cancer': (8, 12), 'immune': (5, 8)}
        )
        self.reference_tissue.generate_cells(max_attempts=1000)
        
        # Extract spatial statistics
        from tissue_simulator import load_target_statistics_from_tissue
        self.target_stats = load_target_statistics_from_tissue(
            self.reference_tissue,
            network_mode="contact"
        )
    
    def run_sampling_experiment(self):
        """
        Generate initial conditions using three methods and run ABM.
        """
        from tissue_simulator import ReplicateGenerator
        
        self.results['final_states'] = {method: [] 
                                       for method in self.ic_methods}
        
        for method in self.ic_methods:
            print(f"\nGenerating initial conditions: {method}")
            
            if method == 'random':
                # Standard random placement
                initial_conditions = []
                for i in range(self.n_replicates):
                    tissue = TissueSection(
                        height=400, width=400, thickness=100,
                        cell_radii={'cancer': (8, 12), 'immune': (5, 8)}
                    )
                    tissue.generate_cells(max_attempts=1000)
                    initial_conditions.append(tissue)
            
            elif method == 'matched':
                # Match reference tissue spatial statistics
                generator = ReplicateGenerator(
                    target_stats=self.target_stats,
                    tissue_dimensions=(400, 400, 100),
                    base_cell_radii={'cancer': (8, 12), 'immune': (5, 8)},
                    network_mode="contact"
                )
                initial_conditions = generator.generate_replicates(
                    num_replicates=self.n_replicates
                )
            
            elif method == 'systematic':
                # Systematic parameter exploration (simplified for demo)
                initial_conditions = []
                clustering_levels = np.linspace(0.2, 0.8, self.n_replicates)
                for clustering in clustering_levels:
                    # Modify target stats based on clustering level
                    # (Implementation would adjust H matrix here)
                    generator = ReplicateGenerator(
                        target_stats=self.target_stats,
                        tissue_dimensions=(400, 400, 100),
                        base_cell_radii={'cancer': (8, 12), 'immune': (5, 8)},
                        network_mode="contact"
                    )
                    ic = generator.generate_replicates(num_replicates=1)[0]
                    initial_conditions.append(ic)
            
            # Run ABM from each initial condition
            for ic in initial_conditions:
                # Placeholder: run ABM simulation
                # final_state = self.run_abm(ic, max_time=1000)
                
                # For demo, just use the IC as "final state"
                final_state = ic
                
                self.results['final_states'][method].append(final_state)
    
    def analyze_results(self):
        """
        Compute variance in spatial metrics across replicates.
        """
        from tissue_simulator import SpatialNetworkAnalyzer
        
        def compute_metric(tissue):
            analyzer = SpatialNetworkAnalyzer()
            analyzer.build_network_from_tissue(tissue, mode="contact")
            stats = analyzer.compute_global_statistics()
            return stats.avg_clustering
        
        self.results['variance_analysis'] = {}
        
        for method in self.ic_methods:
            final_states = self.results['final_states'][method]
            
            # Compute metric for each replicate
            values = [compute_metric(tissue) for tissue in final_states]
            
            # Compute statistics
            mean = np.mean(values)
            std = np.std(values)
            cv = std / mean if mean != 0 else float('inf')
            
            self.results['variance_analysis'][method] = {
                'mean': mean,
                'std': std,
                'cv': cv,
                'values': values
            }
        
        # Print summary
        print("\n=== Variance Analysis ===")
        for method, stats in self.results['variance_analysis'].items():
            print(f"{method:12s}: CV = {stats['cv']:.3f}, "
                  f"mean = {stats['mean']:.3f} ± {stats['std']:.3f}")
    
    def visualize_results(self, save_dir):
        """
        Create comparison figures.
        """
        import os
        os.makedirs(save_dir, exist_ok=True)
        
        # Create violin plot comparing distributions
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        data = [self.results['variance_analysis'][method]['values']
                for method in self.ic_methods]
        
        parts = ax.violinplot(data, showmeans=True, showmedians=True)
        ax.set_xticks(range(1, len(self.ic_methods) + 1))
        ax.set_xticklabels(self.ic_methods)
        ax.set_ylabel('Clustering Coefficient')
        ax.set_title('Distribution of Spatial Metrics Across Replicates')
        ax.grid(True, alpha=0.3)
        
        # Add CV annotations
        for i, method in enumerate(self.ic_methods):
            cv = self.results['variance_analysis'][method]['cv']
            ax.text(i + 1, ax.get_ylim()[1] * 0.95,
                   f'CV = {cv:.3f}',
                   ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, "variance_comparison.png"), 
                   dpi=300)
        plt.close()
        
        print(f"\nResults saved to {save_dir}")


# Main execution script
if __name__ == "__main__":
    import sys
    
    # Run Experiment 1: Convergence Rates
    print("=" * 60)
    print("EXPERIMENT 1: Differential Convergence Rates")
    print("=" * 60)
    
    exp1 = ConvergenceRateExperiment()
    results1 = exp1.run_full_experiment("results/experiment1")
    
    # Run Experiment 2: Initial Conditions
    print("\n" + "=" * 60)
    print("EXPERIMENT 2: Controlled Initial Conditions")
    print("=" * 60)
    
    exp2 = InitialConditionExperiment()
    results2 = exp2.run_full_experiment("results/experiment2")
    
    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 60)
```

---

## Creating Publication-Quality Figures

### Figure Design Guidelines from Baker et al.

**Multi-Panel Layout:**
- Use 2-4 panels per main figure
- Each panel tells one part of the story
- Consistent color schemes across panels
- Clear panel labels (a, b, c, d)

**Typical Figure Structure:**

```
Figure 1: Problem Demonstration
├─ Panel A: Schematic of problem
├─ Panel B: Example data showing issue
├─ Panel C: Quantification across conditions
└─ Panel D: Summary statistics

Figure 2: Solution Framework
├─ Panel A: Method overview schematic
├─ Panel B: Step-by-step example
├─ Panel C: Validation on synthetic data
└─ Panel D: Comparison to baseline

Figure 3: Application to Real Data
├─ Panel A: Real tissue image
├─ Panel B: IST generation result
├─ Panel C: Power curve
└─ Panel D: Sampling recommendations

Figure 4: Generalization
├─ Panel A: Application to system 1
├─ Panel B: Application to system 2
├─ Panel C: Application to system 3
└─ Panel D: Unified framework summary
```

### Example Figure Generation Code

```python
def create_publication_figure_1(results, save_path):
    """
    Create Figure 1: Differential convergence rates across spatial metrics.
    
    4-panel figure:
    A) Time series of all metrics
    B) Rolling CV showing convergence
    C) Convergence time comparison
    D) Burn-in recommendations
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Panel A: Time series
    ax_a = fig.add_subplot(gs[0, 0])
    time_series = results['time_series']
    for metric_name, series in time_series.items():
        times = [t for (t, val) in series]
        values = [val for (t, val) in series]
        ax_a.plot(times, values, label=metric_name, linewidth=2, alpha=0.7)
    
    ax_a.set_xlabel('Time (simulation steps)', fontsize=12)
    ax_a.set_ylabel('Metric Value', fontsize=12)
    ax_a.set_title('A. Temporal Evolution of Spatial Metrics', 
                   fontsize=14, fontweight='bold')
    ax_a.legend(frameon=True, fontsize=10)
    ax_a.grid(True, alpha=0.3)
    ax_a.text(-0.15, 1.05, 'A', transform=ax_a.transAxes,
             fontsize=16, fontweight='bold', va='top')
    
    # Panel B: Rolling CV
    ax_b = fig.add_subplot(gs[0, 1])
    window = 100
    for metric_name, series in time_series.items():
        values = [val for (t, val) in series]
        times = [t for (t, val) in series]
        
        if len(values) >= window:
            rolling_cv = []
            rolling_times = []
            for i in range(window, len(values)):
                window_vals = values[i-window:i]
                mean = np.mean(window_vals)
                std = np.std(window_vals)
                cv = std / abs(mean) if mean != 0 else 0
                rolling_cv.append(cv)
                rolling_times.append(times[i])
            
            ax_b.plot(rolling_times, rolling_cv, label=metric_name, 
                     linewidth=2, alpha=0.7)
    
    ax_b.axhline(0.05, color='red', linestyle='--', linewidth=2,
                label='Convergence threshold')
    ax_b.set_xlabel('Time (simulation steps)', fontsize=12)
    ax_b.set_ylabel('Coefficient of Variation', fontsize=12)
    ax_b.set_title('B. Rolling Variability Over Time',
                  fontsize=14, fontweight='bold')
    ax_b.legend(frameon=True, fontsize=10)
    ax_b.grid(True, alpha=0.3)
    ax_b.text(-0.15, 1.05, 'B', transform=ax_b.transAxes,
             fontsize=16, fontweight='bold', va='top')
    
    # Panel C: Convergence time comparison
    ax_c = fig.add_subplot(gs[1, 0])
    conv_results = results['convergence_results']
    metric_names = list(conv_results.keys())
    conv_times = [conv_results[name]['convergence_time'] 
                 for name in metric_names]
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(metric_names)))
    bars = ax_c.barh(metric_names, conv_times, color=colors, alpha=0.8)
    
    ax_c.set_xlabel('Convergence Time (steps)', fontsize=12)
    ax_c.set_title('C. Metric-Specific Convergence Times',
                  fontsize=14, fontweight='bold')
    ax_c.grid(True, alpha=0.3, axis='x')
    ax_c.text(-0.15, 1.05, 'C', transform=ax_c.transAxes,
             fontsize=16, fontweight='bold', va='top')
    
    # Panel D: Burn-in recommendations
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis('off')
    
    # Create text summary
    burn_in = results['burn_in_period']
    summary_text = f"""
    BURN-IN RECOMMENDATIONS:
    
    Recommended burn-in period: {burn_in} steps
    
    This represents the time required for ALL
    spatial metrics to reach steady-state behavior.
    
    Metric-specific convergence times:
    """
    
    for name, time in sorted(zip(metric_names, conv_times), 
                            key=lambda x: x[1]):
        summary_text += f"\n  • {name}: {time} steps"
    
    summary_text += f"""
    
    IMPLICATIONS:
    
    • Arbitrary burn-in periods (e.g., 48 hours)
      may be insufficient for slower metrics
      
    • Different analyses may require different
      observation windows
      
    • Convergence testing should be standard
      practice in ABM studies
    """
    
    ax_d.text(0.05, 0.95, summary_text, transform=ax_d.transAxes,
             fontsize=10, verticalalignment='top', family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_d.text(-0.05, 1.05, 'D', transform=ax_d.transAxes,
             fontsize=16, fontweight='bold', va='top')
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Figure saved to {save_path}")
```

---

## Next Steps

### Immediate Implementation Priorities

1. **Create convergence analysis module** (highest priority)
   - Implements all statistical tests
   - Integrates with existing tissue_simulator
   - Well-documented with examples

2. **Implement experiment templates**
   - Ready-to-run scripts following Baker et al. pattern
   - Clear documentation of each experiment
   - Easy to adapt for different ABM systems

3. **Develop figure generation utilities**
   - Consistent styling across all figures
   - Publication-quality by default
   - Easy customization

4. **Write comprehensive documentation**
   - Tutorial for ABM researchers
   - Mathematical derivations in appendix
   - Code examples for each use case

5. **Validation studies**
   - Run on 2-3 real ABM systems
   - Compare to literature results
   - Demonstrate improvements

### Writing Timeline

**Week 1-2: Implementation**
- Core modules (convergence, power analysis)
- Example experiments
- Figure generation

**Week 3-4: Validation**
- Run experiments on real ABMs
- Generate all figures
- Statistical analysis

**Week 5-6: Writing**
- Draft Introduction + Methods
- Draft Results
- Create all figures

**Week 7-8: Revision**
- Internal review
- Revise based on feedback
- Prepare supplementary materials

**Week 9: Submission**
- Format for target journal
- Final checks
- Submit

---

This document provides a practical roadmap for:
1. Understanding Baker et al.'s experimental design
2. Structuring our own papers
3. Implementing the experiments
4. Creating publication-quality figures
5. Timeline for completion

The key insight is to follow their proven pattern:
**Problem → Solution → Validation → Generalization**

Each experiment should clearly demonstrate a specific problem,
show how our framework solves it, and validate on real data.
