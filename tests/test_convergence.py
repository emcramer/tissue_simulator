"""
Tests for the convergence diagnostics module.
"""

from typing import Optional

import numpy as np
import pytest

from tissue_simulator.convergence import (
    MultiMetricConvergence,
    adf_test,
    find_convergence_time,
    mann_kendall_test,
    rolling_cv,
)


RNG = np.random.default_rng(seed=42)


def _ar1(n: int, phi: float = 0.5, sigma: float = 1.0,
         rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """Generate a stationary AR(1) series."""
    if rng is None:
        rng = np.random.default_rng(seed=42)
    eps = rng.normal(0.0, sigma, size=n)
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + eps[i]
    return x


def _random_walk(n: int, sigma: float = 1.0,
                 rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """Generate a non-stationary random walk."""
    if rng is None:
        rng = np.random.default_rng(seed=42)
    return np.cumsum(rng.normal(0.0, sigma, size=n))


def test_adf_lazy_import_works():
    """The lazy statsmodels import inside adf_test should succeed."""
    series = _ar1(200, phi=0.3, rng=np.random.default_rng(seed=42))
    result = adf_test(series)
    assert "statistic" in result
    assert "pvalue" in result
    assert "is_stationary" in result
    assert "critical_values" in result
    assert "lags_used" in result
    assert isinstance(result["critical_values"], dict)


def test_adf_flags_stationary_ar1():
    """A stationary AR(1) process should be flagged stationary."""
    series = _ar1(500, phi=0.3, rng=np.random.default_rng(seed=42))
    result = adf_test(series)
    assert result["is_stationary"] is True
    assert result["pvalue"] < 0.05


def test_adf_flags_random_walk_nonstationary():
    """A random walk should not be flagged stationary."""
    series = _random_walk(500, rng=np.random.default_rng(seed=42))
    result = adf_test(series)
    assert result["is_stationary"] is False
    assert result["pvalue"] >= 0.05


def test_mann_kendall_detects_increasing_trend():
    """Linear trend plus noise should be flagged as increasing."""
    rng = np.random.default_rng(seed=42)
    n = 200
    series = np.linspace(0.0, 10.0, n) + rng.normal(0.0, 0.5, size=n)
    result = mann_kendall_test(series)
    assert result["trend"] == "increasing"
    assert result["is_significant"] is True
    assert result["S"] > 0
    assert result["z"] > 0
    assert result["pvalue"] < 0.05


def test_mann_kendall_no_trend_on_white_noise():
    """White noise should generally show no significant trend."""
    rng = np.random.default_rng(seed=42)
    series = rng.normal(0.0, 1.0, size=300)
    result = mann_kendall_test(series)
    assert result["trend"] == "no trend"
    assert result["is_significant"] is False


def test_mann_kendall_handles_short_series():
    """Series shorter than 3 points should return a no-trend stub."""
    result = mann_kendall_test([1.0, 2.0])
    assert result["trend"] == "no trend"
    assert result["is_significant"] is False


def test_rolling_cv_shape_and_nan_prefix():
    """rolling_cv must match the input length and pad with NaNs."""
    rng = np.random.default_rng(seed=42)
    series = rng.normal(5.0, 0.1, size=50)
    window = 10
    cv = rolling_cv(series, window)
    assert cv.shape == series.shape
    assert np.all(np.isnan(cv[: window - 1]))
    assert np.all(np.isfinite(cv[window - 1 :]))
    # CV should be small relative to mean for tight series.
    assert np.nanmax(cv) < 0.1


def test_rolling_cv_invalid_window():
    """Window smaller than 2 must raise."""
    with pytest.raises(ValueError):
        rolling_cv([1.0, 2.0, 3.0], window=1)


def test_find_convergence_time_transient_then_stable():
    """A noisy transient followed by a tight steady state should converge."""
    rng = np.random.default_rng(seed=42)
    transient = rng.normal(0.0, 5.0, size=80) + np.linspace(20.0, 5.0, 80)
    steady = 10.0 + rng.normal(0.0, 0.02, size=200)
    series = np.concatenate([transient, steady])

    t_conv = find_convergence_time(
        series, window=20, cv_threshold=0.05, require_stationary=True
    )
    assert t_conv is not None
    # Must occur after the transient ends.
    assert t_conv >= 80


def test_find_convergence_time_never_converges():
    """A deterministic linear ramp violates both CV and stationarity."""
    series = np.linspace(1.0, 50.0, 300)
    t_conv = find_convergence_time(
        series, window=20, cv_threshold=0.001, require_stationary=True
    )
    assert t_conv is None


def test_find_convergence_time_without_stationarity():
    """Disabling ADF should still find a low-CV window."""
    rng = np.random.default_rng(seed=42)
    series = 10.0 + rng.normal(0.0, 0.02, size=100)
    t_conv = find_convergence_time(
        series, window=20, cv_threshold=0.05, require_stationary=False
    )
    assert t_conv is not None
    assert t_conv >= 19


def test_multi_metric_convergence_aggregates():
    """MultiMetricConvergence should report per-metric times and global max."""
    rng = np.random.default_rng(seed=42)
    steady_a = 5.0 + rng.normal(0.0, 0.01, size=120)
    steady_b = np.concatenate([
        rng.normal(0.0, 2.0, size=40) + np.linspace(10.0, 5.0, 40),
        5.0 + rng.normal(0.0, 0.01, size=120),
    ])
    mmc = MultiMetricConvergence(
        {"metric_a": steady_a, "metric_b": steady_b},
        window=20,
        cv_threshold=0.05,
        require_stationary=True,
    )
    assert set(mmc.convergence_times) == {"metric_a", "metric_b"}
    assert mmc.convergence_times["metric_a"] is not None
    assert mmc.convergence_times["metric_b"] is not None
    assert mmc.global_convergence_time == max(
        mmc.convergence_times["metric_a"],
        mmc.convergence_times["metric_b"],
    )


def test_multi_metric_global_none_if_any_fails():
    """If any metric never converges, the global time must be None."""
    rng = np.random.default_rng(seed=42)
    steady = 5.0 + rng.normal(0.0, 0.001, size=120)
    # A deterministic linear ramp is neither low-CV (over a 20-window the
    # std/|mean| stays above 0.01) nor stationary, so it never converges.
    ramp = np.linspace(1.0, 50.0, 120)
    mmc = MultiMetricConvergence(
        {"good": steady, "bad": ramp},
        window=20,
        cv_threshold=0.01,
        require_stationary=True,
    )
    assert mmc.convergence_times["good"] is not None
    assert mmc.convergence_times["bad"] is None
    assert mmc.global_convergence_time is None


def test_multi_metric_summary_is_string():
    """summary() should return a non-empty string mentioning each metric."""
    rng = np.random.default_rng(seed=42)
    steady = 5.0 + rng.normal(0.0, 0.01, size=120)
    mmc = MultiMetricConvergence(
        {"metric_a": steady}, window=20, cv_threshold=0.05
    )
    report = mmc.summary()
    assert isinstance(report, str)
    assert "metric_a" in report
    assert "Global convergence time" in report
