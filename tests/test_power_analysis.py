"""
Tests for tissue_simulator.power_analysis.
"""

import numpy as np
import pytest

from tissue_simulator.power_analysis import (
    cohens_d,
    coefficient_of_variation,
    required_replicates,
    power_curve,
    compare_initialization_variance,
    summarize_power_analysis,
)


SEED = 42


def test_cohens_d_known_value():
    rng = np.random.default_rng(seed=SEED)
    n = 5000
    a = rng.normal(loc=0.0, scale=1.0, size=n)
    b = rng.normal(loc=0.5, scale=1.0, size=n)

    d = cohens_d(b, a)
    assert d == pytest.approx(0.5, abs=0.05)

    d_signed = cohens_d(a, b)
    assert d_signed == pytest.approx(-d, abs=1e-9)


def test_cohens_d_zero_when_identical():
    rng = np.random.default_rng(seed=SEED)
    a = rng.normal(0, 1, size=100)
    d = cohens_d(a, a)
    assert d == pytest.approx(0.0, abs=1e-12)


def test_cohens_d_requires_two_samples():
    with pytest.raises(ValueError):
        cohens_d([1.0], [2.0, 3.0])


def test_cv_constant_array_is_zero():
    values = np.full(10, 3.14)
    assert coefficient_of_variation(values) == pytest.approx(0.0, abs=1e-12)


def test_cv_zero_mean_returns_inf():
    values = np.array([-1.0, 1.0, -1.0, 1.0])
    assert np.isinf(coefficient_of_variation(values))


def test_cv_known_value():
    rng = np.random.default_rng(seed=SEED)
    x = rng.normal(loc=10.0, scale=2.0, size=10000)
    cv = coefficient_of_variation(x)
    assert cv == pytest.approx(0.2, abs=0.02)


def test_required_replicates_monotonic_in_shrinking_effect():
    effect_sizes = [1.0, 0.8, 0.5, 0.3, 0.2]
    ns = [required_replicates(es) for es in effect_sizes]
    assert all(ns[i] <= ns[i + 1] for i in range(len(ns) - 1))
    assert all(n > 0 for n in ns)


def test_required_replicates_known_benchmark():
    # Cohen's d = 0.5, alpha = 0.05, power = 0.8 -> ~64 per group.
    n = required_replicates(0.5, alpha=0.05, power=0.8)
    assert 60 <= n <= 70


def test_required_replicates_rejects_zero_effect():
    with pytest.raises(ValueError):
        required_replicates(0.0)


def test_power_curve_shape_and_bounds():
    effect_sizes = np.array([0.2, 0.5, 0.8])
    n_range = np.array([10, 25, 50, 100])
    grid = power_curve(effect_sizes, n_range)

    assert grid.shape == (3, 4)
    assert np.all(grid >= 0.0)
    assert np.all(grid <= 1.0)

    # Power should be non-decreasing in N (for fixed effect size).
    for i in range(grid.shape[0]):
        assert np.all(np.diff(grid[i, :]) >= -1e-9)

    # Power should be non-decreasing in effect size (for fixed N).
    for j in range(grid.shape[1]):
        assert np.all(np.diff(grid[:, j]) >= -1e-9)


def test_compare_initialization_variance_structure():
    rng = np.random.default_rng(seed=SEED)
    endpoints = {
        "random_init": rng.normal(loc=10.0, scale=2.0, size=30),
        "tissue_sim_init": rng.normal(loc=10.0, scale=0.5, size=30),
        "grid_init": rng.normal(loc=11.0, scale=1.0, size=30),
    }

    result = compare_initialization_variance(endpoints)

    assert "per_method" in result
    assert "pairwise" in result
    assert set(result["per_method"].keys()) == set(endpoints.keys())

    for name, stats in result["per_method"].items():
        for key in ("n", "mean", "std", "cv"):
            assert key in stats
        assert stats["n"] == 30
        assert stats["std"] >= 0.0

    # Three methods -> C(3, 2) = 3 unordered pairs.
    assert len(result["pairwise"]) == 3
    for rec in result["pairwise"]:
        for key in (
            "method_a",
            "method_b",
            "cohens_d",
            "required_n_per_group",
            "alpha",
            "power",
        ):
            assert key in rec
        assert rec["method_a"] != rec["method_b"]
        assert isinstance(rec["cohens_d"], float)

    # tissue_sim_init should have a lower CV than random_init.
    assert (
        result["per_method"]["tissue_sim_init"]["cv"]
        < result["per_method"]["random_init"]["cv"]
    )


def test_summarize_power_analysis_returns_str():
    rng = np.random.default_rng(seed=SEED)
    endpoints = {
        "a": rng.normal(0.0, 1.0, size=20),
        "b": rng.normal(0.5, 1.0, size=20),
    }
    comparison = compare_initialization_variance(endpoints)
    report = summarize_power_analysis(comparison)

    assert isinstance(report, str)
    assert "INITIALIZATION VARIANCE" in report
    assert "a" in report and "b" in report
    assert "cohens_d" in report
