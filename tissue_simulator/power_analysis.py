"""
Power analysis utilities for ABM replicate planning.

This module provides effect-size, coefficient-of-variation, and
statistical-power tools used to demonstrate that tissue_simulator
initialization reduces inter-replicate variance, lowering the number
of simulation replicates required to achieve a target power.

The functions here are intended to compare endpoint distributions
produced by different ABM initialization strategies (e.g., random
placement versus tissue_simulator-driven placement) on identical
downstream models.
"""

from typing import Dict, Mapping, Sequence, Union

import numpy as np


ArrayLike = Union[Sequence[float], np.ndarray]


def cohens_d(group_a: ArrayLike, group_b: ArrayLike) -> float:
    """
    Compute the pooled-standard-deviation Cohen's d effect size.

    Uses sample standard deviations with ``ddof=1`` (Bessel's correction)
    and the standard pooled-variance formula.

    Parameters
    ----------
    group_a : array_like
        First sample of endpoint values.
    group_b : array_like
        Second sample of endpoint values.

    Returns
    -------
    float
        Cohen's d. Positive values indicate ``group_a`` has the larger
        mean. Returns 0.0 when both groups are constant and identical;
        returns ``inf`` when the pooled standard deviation is zero but
        the means differ.
    """
    a = np.asarray(group_a, dtype=float).ravel()
    b = np.asarray(group_b, dtype=float).ravel()

    if a.size < 2 or b.size < 2:
        raise ValueError("Each group must contain at least two observations.")

    n_a, n_b = a.size, b.size
    var_a = np.var(a, ddof=1)
    var_b = np.var(b, ddof=1)

    pooled_var = ((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2)
    pooled_std = np.sqrt(pooled_var)

    mean_diff = np.mean(a) - np.mean(b)

    if pooled_std == 0:
        if mean_diff == 0:
            return 0.0
        return float(np.inf) if mean_diff > 0 else float(-np.inf)

    return float(mean_diff / pooled_std)


def coefficient_of_variation(values: ArrayLike) -> float:
    """
    Compute the coefficient of variation: std / |mean|.

    Uses sample standard deviation (``ddof=1``).

    Parameters
    ----------
    values : array_like
        Sample of endpoint values.

    Returns
    -------
    float
        Coefficient of variation. Returns ``inf`` if the sample mean is
        zero (CV is undefined in that case).
    """
    x = np.asarray(values, dtype=float).ravel()
    if x.size < 2:
        raise ValueError("At least two observations are required.")

    mean = np.mean(x)
    if mean == 0:
        return float(np.inf)

    std = np.std(x, ddof=1)
    return float(std / abs(mean))


def required_replicates(
    effect_size: float,
    alpha: float = 0.05,
    power: float = 0.8,
    alternative: str = "two-sided",
) -> int:
    """
    Solve for the per-group sample size of a two-sample t-test.

    Returns the smallest integer N per group such that an independent
    two-sample t-test achieves at least the requested power for the
    given Cohen's d effect size.

    Parameters
    ----------
    effect_size : float
        Cohen's d. The absolute value is used; sign is irrelevant for
        sample-size calculation.
    alpha : float, optional
        Type-I error rate. Defaults to 0.05.
    power : float, optional
        Desired statistical power (1 - beta). Defaults to 0.8.
    alternative : {"two-sided", "larger", "smaller"}, optional
        Alternative hypothesis as accepted by statsmodels. Defaults to
        ``"two-sided"``.

    Returns
    -------
    int
        Required sample size per group, rounded up.

    Raises
    ------
    ImportError
        If statsmodels is not installed.
    ValueError
        If ``effect_size`` is zero (sample size is undefined).
    """
    if effect_size == 0:
        raise ValueError("effect_size must be non-zero to solve for N.")

    try:
        from statsmodels.stats.power import TTestIndPower
    except ImportError as exc:
        raise ImportError(
            "required_replicates requires statsmodels. "
            "Install it with: pip install statsmodels"
        ) from exc

    analysis = TTestIndPower()
    n = analysis.solve_power(
        effect_size=abs(float(effect_size)),
        alpha=alpha,
        power=power,
        ratio=1.0,
        alternative=alternative,
    )
    return int(np.ceil(n))


def power_curve(
    effect_sizes: ArrayLike,
    n_range: ArrayLike,
    alpha: float = 0.05,
    alternative: str = "two-sided",
) -> np.ndarray:
    """
    Build a 2D grid of achieved power across effect sizes and sample sizes.

    Parameters
    ----------
    effect_sizes : array_like
        Cohen's d values to evaluate (rows of the output).
    n_range : array_like
        Per-group sample sizes to evaluate (columns of the output).
    alpha : float, optional
        Type-I error rate. Defaults to 0.05.
    alternative : {"two-sided", "larger", "smaller"}, optional
        Alternative hypothesis. Defaults to ``"two-sided"``.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(len(effect_sizes), len(n_range))`` containing
        the achieved power for each (effect size, N) pair.

    Raises
    ------
    ImportError
        If statsmodels is not installed.
    """
    try:
        from statsmodels.stats.power import TTestIndPower
    except ImportError as exc:
        raise ImportError(
            "power_curve requires statsmodels. "
            "Install it with: pip install statsmodels"
        ) from exc

    es = np.asarray(effect_sizes, dtype=float).ravel()
    ns = np.asarray(n_range, dtype=float).ravel()

    analysis = TTestIndPower()
    grid = np.empty((es.size, ns.size), dtype=float)
    for i, d in enumerate(es):
        for j, n in enumerate(ns):
            grid[i, j] = analysis.power(
                effect_size=abs(float(d)),
                nobs1=float(n),
                alpha=alpha,
                ratio=1.0,
                alternative=alternative,
            )
    return grid


def compare_initialization_variance(
    endpoints_by_method: Mapping[str, ArrayLike],
    alpha: float = 0.05,
    power: float = 0.8,
) -> Dict[str, object]:
    """
    Compare endpoint variability across initialization methods.

    For each method computes mean, sample standard deviation, and
    coefficient of variation. For each unordered pair of methods,
    computes Cohen's d on the endpoint samples and the per-group N
    needed to detect that effect at the requested alpha and power.

    Parameters
    ----------
    endpoints_by_method : dict[str, array_like]
        Mapping of method name to a 1D array of endpoint values (one
        value per replicate of that method).
    alpha : float, optional
        Type-I error rate used for the required-N calculation.
    power : float, optional
        Desired statistical power used for the required-N calculation.

    Returns
    -------
    dict
        Structured results with two top-level keys:

        ``"per_method"`` : dict[str, dict]
            Per-method ``n``, ``mean``, ``std``, ``cv``.
        ``"pairwise"`` : list[dict]
            One record per unordered method pair, each with
            ``method_a``, ``method_b``, ``cohens_d``,
            ``required_n_per_group``, ``alpha``, ``power``.
            ``required_n_per_group`` is ``None`` when the effect size
            is zero or statsmodels is unavailable.
    """
    per_method: Dict[str, Dict[str, float]] = {}
    arrays: Dict[str, np.ndarray] = {}

    for name, values in endpoints_by_method.items():
        arr = np.asarray(values, dtype=float).ravel()
        arrays[name] = arr
        per_method[name] = {
            "n": int(arr.size),
            "mean": float(np.mean(arr)) if arr.size else float("nan"),
            "std": float(np.std(arr, ddof=1)) if arr.size > 1 else float("nan"),
            "cv": (
                coefficient_of_variation(arr) if arr.size > 1 else float("nan")
            ),
        }

    names = list(arrays.keys())
    pairwise = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a_name, b_name = names[i], names[j]
            d = cohens_d(arrays[a_name], arrays[b_name])

            required_n: Union[int, None]
            if d == 0 or not np.isfinite(d):
                required_n = None
            else:
                try:
                    required_n = required_replicates(
                        effect_size=d, alpha=alpha, power=power
                    )
                except ImportError:
                    required_n = None

            pairwise.append(
                {
                    "method_a": a_name,
                    "method_b": b_name,
                    "cohens_d": float(d),
                    "required_n_per_group": required_n,
                    "alpha": float(alpha),
                    "power": float(power),
                }
            )

    return {"per_method": per_method, "pairwise": pairwise}


def summarize_power_analysis(comparison: Dict[str, object]) -> str:
    """
    Format a human-readable report from ``compare_initialization_variance``.

    Parameters
    ----------
    comparison : dict
        Output of :func:`compare_initialization_variance`.

    Returns
    -------
    str
        Multi-line report covering per-method summary statistics and
        all pairwise effect sizes and required sample sizes.
    """
    per_method = comparison.get("per_method", {})
    pairwise = comparison.get("pairwise", [])

    lines = []
    lines.append("=" * 60)
    lines.append("INITIALIZATION VARIANCE / POWER ANALYSIS REPORT")
    lines.append("=" * 60)

    lines.append("")
    lines.append("--- Per-Method Summary ---")
    header = f"{'method':20s} {'n':>5s} {'mean':>12s} {'std':>12s} {'cv':>10s}"
    lines.append(header)
    for name, stats in per_method.items():
        lines.append(
            f"{name:20s} {stats['n']:>5d} "
            f"{stats['mean']:>12.4f} {stats['std']:>12.4f} "
            f"{stats['cv']:>10.4f}"
        )

    lines.append("")
    lines.append("--- Pairwise Effect Sizes ---")
    pair_header = (
        f"{'method_a':20s} {'method_b':20s} "
        f"{'cohens_d':>10s} {'req_n/grp':>10s}"
    )
    lines.append(pair_header)
    for rec in pairwise:
        req_n = rec["required_n_per_group"]
        req_n_str = "n/a" if req_n is None else f"{req_n:d}"
        lines.append(
            f"{rec['method_a']:20s} {rec['method_b']:20s} "
            f"{rec['cohens_d']:>10.4f} {req_n_str:>10s}"
        )

    if pairwise:
        lines.append("")
        lines.append(
            f"(alpha={pairwise[0]['alpha']:.3f}, "
            f"power={pairwise[0]['power']:.2f}, two-sided t-test)"
        )

    lines.append("=" * 60)
    return "\n".join(lines)
