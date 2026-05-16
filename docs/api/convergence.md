# Convergence Diagnostics

## Overview

`tissue_simulator.convergence` provides stationarity and convergence-detection
utilities for agent-based-model (ABM) trajectories. The intended use case is
identifying when a simulated tissue ABM has left its initial transient and
reached steady state, so that subsequent samples can be treated as draws from
the long-run distribution. This module targets the burn-in-elimination claim
in the paper roadmap: by combining the Augmented Dickey-Fuller (ADF) test, the
Mann-Kendall trend test, and a rolling coefficient-of-variation criterion, it
gives a per-trajectory burn-in length that can be discarded before downstream
analysis. The ADF code path requires the optional `statsmodels` dependency
(`pip install statsmodels`); the Mann-Kendall and rolling-CV paths depend only
on NumPy and SciPy.

For replicate-count planning on top of these converged samples see
[`power-analysis.md`](power-analysis.md). For the broader project tour see
[`../quickstart.md`](../quickstart.md).

## Public API

The public names below match `tissue_simulator/__init__.py` and can be imported
either from the top-level package or from `tissue_simulator.convergence`.

### `adf_test(series, regression="c") -> dict`

Run the Augmented Dickey-Fuller test for unit-root non-stationarity. `series`
is a 1D numeric sequence; `regression` is passed through to
`statsmodels.tsa.stattools.adfuller` (`"c"`, `"ct"`, `"ctt"`, or `"n"`).

Returns a `dict` with keys:

- `statistic` (float) - the ADF test statistic.
- `pvalue` (float) - two-sided p-value.
- `is_stationary` (bool) - True when `pvalue < 0.05` (null of unit root
  rejected).
- `critical_values` (dict[str, float]) - mapping of significance level (e.g.
  `"1%"`, `"5%"`, `"10%"`) to its critical value.
- `lags_used` (int) - number of lags selected by AIC.

Use when you need a formal stationarity test on a single window. Raises
`ImportError` if `statsmodels` is not installed.

### `mann_kendall_test(series, alpha=0.05) -> dict`

Two-sided Mann-Kendall trend test with the standard tie-adjusted variance and
continuity-corrected z-score.

Returns a `dict` with keys:

- `trend` (str) - one of `"increasing"`, `"decreasing"`, or `"no trend"`.
- `S` (float) - Mann-Kendall S statistic.
- `z` (float) - normal-approximation z-score after continuity correction.
- `pvalue` (float) - two-sided p-value from the standard normal.
- `is_significant` (bool) - True when `pvalue < alpha`.

Use as a non-parametric alternative to ADF when the trajectory may have a
monotonic trend but is not necessarily Gaussian. Pure NumPy/SciPy; no
`statsmodels` dependency.

### `rolling_cv(series, window) -> numpy.ndarray`

Rolling coefficient of variation, `std / |mean|`, computed with `ddof=0`.
Returns a NumPy array the same length as `series`; the first `window - 1`
entries are `NaN`, and any window whose mean is zero or non-finite is `NaN`.

Use as a cheap, distribution-free "is the signal still wandering" check.
Raises `ValueError` if `window < 2`.

### `find_convergence_time(series, window=20, cv_threshold=0.05, require_stationary=True) -> int | None`

Return the first index `t` at which both criteria hold:

1. The rolling CV over `series[t - window + 1 : t + 1]` is below
   `cv_threshold`.
2. If `require_stationary=True`, the ADF test on that same window rejects the
   unit-root null. Degenerate (constant) windows automatically pass.

Returns `None` if no such index exists. With `require_stationary=False` the
ADF step is skipped, which removes the `statsmodels` dependency. This is the
primary entry point for burn-in detection.

### `MultiMetricConvergence(metrics, window=20, cv_threshold=0.05, require_stationary=True)`

Aggregate `find_convergence_time` across several named metrics. `metrics` is a
dict mapping metric name to a 1D series. After construction:

- `mmc.convergence_times` (dict[str, int | None]) - per-metric convergence
  index, or `None` for metrics that never converged.
- `mmc.global_convergence_time` (int | None) - the maximum over all per-metric
  times, or `None` if any metric returned `None`.
- `mmc.summary()` (str) - formatted multi-line report.

Use when an ABM tracks several scalar outputs (e.g., cell counts per type,
mean nearest-neighbor distance) and you want a single burn-in length that
satisfies all of them.

## Worked example

The snippet below builds a synthetic trajectory consisting of a noisy ramp
transient followed by a stationary AR(1) signal, then asks
`find_convergence_time` where the burn-in ends. It is reproducible via
`numpy.random.default_rng(seed=...)`.

```python
import numpy as np

from tissue_simulator.convergence import find_convergence_time

# step 1: build a synthetic trajectory: noisy transient ramp then AR(1) steady state.
rng = np.random.default_rng(seed=20260515)

# step 2: 80 transient samples that decay linearly from 20 -> 5 with heavy noise.
transient = np.linspace(20.0, 5.0, 80) + rng.normal(0.0, 5.0, size=80)

# step 3: 200 stationary AR(1) samples centered at 10 (small residual sigma).
phi = 0.4
eps = rng.normal(0.0, 0.05, size=200)
steady = np.zeros(200)
steady[0] = 10.0
for i in range(1, 200):
    steady[i] = 10.0 + phi * (steady[i - 1] - 10.0) + eps[i]

trajectory = np.concatenate([transient, steady])

# step 4: detect the first index where rolling CV is small AND ADF says stationary.
t_conv = find_convergence_time(
    trajectory,
    window=20,
    cv_threshold=0.05,
    require_stationary=True,
)

# step 5: report the result.
if t_conv is None:
    print("Trajectory never converged within the recorded horizon.")
else:
    print(f"Converged at t = {t_conv} (out of {trajectory.size} timesteps)")
    print(f"Value at convergence: {trajectory[t_conv]:.4f}")
    print(f"Burn-in length to discard: {t_conv + 1} samples")
```

Expected output (line-for-line shape; exact `t` value depends on the seed):

```text
Converged at t = <t> (out of 280 timesteps)
Value at convergence: <value near 10.0>
Burn-in length to discard: <t + 1> samples
```

By construction the steady-state window starts at index 80, so the reported
convergence index will satisfy `t >= 99` (the earliest `window=20` slice that
lies entirely inside the steady-state segment). The corresponding test in
`tests/test_convergence.py::test_find_convergence_time_transient_then_stable`
asserts the analogous lower bound for a closely related trajectory.

## References

- Hipel, K. W. and McLeod, A. I. (1994). *Time Series Modelling of Water
  Resources and Environmental Systems*. Elsevier. Chapter 23 covers the
  Mann-Kendall S statistic and its tie-adjusted variance, which is what
  `mann_kendall_test` implements.
- Dickey, D. A. and Fuller, W. A. (1979). "Distribution of the Estimators for
  Autoregressive Time Series with a Unit Root." *Journal of the American
  Statistical Association*, 74(366a), 427-431. The ADF test used by
  `adf_test` is the augmented version of this statistic.
- `statsmodels` reference for the underlying ADF implementation:
  <https://www.statsmodels.org/stable/generated/statsmodels.tsa.stattools.adfuller.html>.

See also [`power-analysis.md`](power-analysis.md) for replicate-count planning
on top of the post-burn-in samples, and [`../quickstart.md`](../quickstart.md)
for the project-wide tour.
