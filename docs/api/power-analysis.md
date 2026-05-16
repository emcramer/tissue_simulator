# Power Analysis

## Overview

`tissue_simulator.power_analysis` provides effect-size, coefficient-of-variation,
and statistical-power tools for planning the number of ABM replicates needed to
detect a given effect. The intended use case is comparing endpoint distributions
produced by different initialization strategies (for example random placement
versus `tissue_simulator`-driven placement) on the same downstream ABM, in
support of the paper roadmap's variance-reduction claim: lower inter-replicate
variance means fewer replicates are required to reach a target power. The
sample-size and power-curve functions delegate to `statsmodels.stats.power.TTestIndPower`
and therefore require the optional `statsmodels` dependency
(`pip install statsmodels`); `cohens_d` and `coefficient_of_variation` are
pure NumPy.

For detecting burn-in on the trajectories that produce these endpoints see
[`convergence.md`](convergence.md). For the broader project tour see
[`../quickstart.md`](../quickstart.md).

## Public API

The public names below match `tissue_simulator/__init__.py` and can be imported
either from the top-level package or from `tissue_simulator.power_analysis`.

### `cohens_d(group_a, group_b) -> float`

Pooled-standard-deviation Cohen's d effect size on two 1D samples. Uses
`ddof=1`. Returns `float`. Positive when `group_a` has the larger mean;
returns `inf` (or `-inf`) when the pooled standard deviation is zero but the
means differ; returns `0.0` when both groups are constant and equal. Raises
`ValueError` if either group has fewer than two observations.

Use when summarizing the magnitude of the mean shift between two replicate
sets in standard-deviation units.

### `coefficient_of_variation(values) -> float`

Sample coefficient of variation, `std(ddof=1) / |mean|`. Returns `float`;
returns `inf` when the sample mean is exactly zero. Raises `ValueError` for
fewer than two observations.

Use as a unit-free dispersion summary for a single initialization method's
endpoints.

### `required_replicates(effect_size, alpha=0.05, power=0.8, alternative="two-sided") -> int`

Solve for the smallest per-group sample size of an independent two-sample
t-test that achieves at least the requested power for the given Cohen's d
effect size. Returns `int` (rounded up). Raises `ValueError` if
`effect_size == 0` and `ImportError` if `statsmodels` is not installed.

Use to translate an observed (or hypothesized) effect size into a replicate
budget per group.

### `power_curve(effect_sizes, n_range, alpha=0.05, alternative="two-sided") -> numpy.ndarray`

Evaluate achieved power on a 2D grid: rows over `effect_sizes`, columns over
per-group sample sizes `n_range`. Returns a `numpy.ndarray` of shape
`(len(effect_sizes), len(n_range))` with entries in `[0, 1]`. Requires
`statsmodels`.

Use to sweep the (d, N) plane when picking a replicate count by hand.

### `compare_initialization_variance(endpoints_by_method, alpha=0.05, power=0.8) -> dict`

For each named initialization method, compute `n`, `mean`, `std`, and `cv`.
For each unordered pair of methods, compute Cohen's d and the per-group N
required to detect that effect. Returns a `dict` of the form:

```text
{
    "per_method": {
        "<name>": {"n": int, "mean": float, "std": float, "cv": float},
        ...
    },
    "pairwise": [
        {
            "method_a": str,
            "method_b": str,
            "cohens_d": float,
            "required_n_per_group": int | None,
            "alpha": float,
            "power": float,
        },
        ...
    ],
}
```

`required_n_per_group` is `None` when the effect size is zero or
non-finite, or when `statsmodels` is unavailable.

Use as the single entry point for the variance-reduction comparison.

### `summarize_power_analysis(comparison) -> str`

Format a multi-line, fixed-width report from the output of
`compare_initialization_variance`. Returns `str`. Use for logging or pasting
into a notebook.

## Worked example

The snippet below feeds three synthetic endpoint arrays (one per
"initialization method") to `compare_initialization_variance` and prints the
resulting summary. It is reproducible via `numpy.random.default_rng(seed=...)`.

```python
import numpy as np

from tissue_simulator.power_analysis import (
    compare_initialization_variance,
    summarize_power_analysis,
)

# step 1: synthesize endpoint values for three "initialization methods".
# Each array is one scalar endpoint per ABM replicate.
rng = np.random.default_rng(seed=20260515)

# step 2: random_init has the largest spread (highest CV).
random_init = rng.normal(loc=10.0, scale=2.0, size=30)

# step 3: tissue_sim_init has the same mean but a much tighter spread.
tissue_sim_init = rng.normal(loc=10.0, scale=0.5, size=30)

# step 4: grid_init has a small mean shift and a moderate spread.
grid_init = rng.normal(loc=11.0, scale=1.0, size=30)

# step 5: run the comparison and print the human-readable report.
comparison = compare_initialization_variance(
    {
        "random_init": random_init,
        "tissue_sim_init": tissue_sim_init,
        "grid_init": grid_init,
    },
    alpha=0.05,
    power=0.8,
)

print(summarize_power_analysis(comparison))
```

Expected output (line-for-line shape; exact numbers depend on the seed):

```text
============================================================
INITIALIZATION VARIANCE / POWER ANALYSIS REPORT
============================================================

--- Per-Method Summary ---
method                   n         mean          std         cv
random_init             30      <mean>       <std>      <cv>
tissue_sim_init         30      <mean>       <std>      <cv>
grid_init               30      <mean>       <std>      <cv>

--- Pairwise Effect Sizes ---
method_a             method_b               cohens_d  req_n/grp
random_init          tissue_sim_init           <d>       <N>
random_init          grid_init                 <d>       <N>
tissue_sim_init      grid_init                 <d>       <N>

(alpha=0.050, power=0.80, two-sided t-test)
============================================================
```

`tissue_sim_init` is constructed to have a much smaller CV than `random_init`;
the test
`tests/test_power_analysis.py::test_compare_initialization_variance_structure`
asserts the analogous inequality for a closely related comparison.

## Interpretation cheatsheet

- Cohen's d magnitudes follow the conventional labels:
  - `|d| ~ 0.2`: small effect.
  - `|d| ~ 0.5`: medium effect.
  - `|d| ~ 0.8`: large effect.
- `required_replicates(d, alpha=0.05, power=0.8)` is the per-group N for a
  two-sample independent t-test (equal group sizes, `ratio=1.0`). The total
  number of replicates across both groups is therefore `2 * N`.
- The coefficient of variation is unit-free; halving the CV of an
  initialization method roughly halves the per-group N needed to detect a
  fixed mean shift, because Cohen's d scales as mean-difference divided by
  pooled standard deviation.
- `compare_initialization_variance` returns `required_n_per_group=None` when
  the two methods have identical means (effect size zero), since sample size
  is undefined in that case.

## References

- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*
  (2nd ed.). Lawrence Erlbaum Associates. The conventional small/medium/large
  effect-size thresholds (0.2 / 0.5 / 0.8) come from this work.
- `statsmodels` reference for the underlying power solver:
  <https://www.statsmodels.org/stable/generated/statsmodels.stats.power.TTestIndPower.html>.

See also [`convergence.md`](convergence.md) for burn-in detection on the
trajectories that produce these endpoints, and
[`../quickstart.md`](../quickstart.md) for the project-wide tour.
