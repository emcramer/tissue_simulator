"""
Convergence diagnostics for agent-based-model trajectories.

This module provides stationarity testing and convergence-detection utilities
intended for evaluating how quickly an ABM trajectory reaches steady state.
It combines the Augmented Dickey-Fuller test, the Mann-Kendall trend test,
and a rolling coefficient-of-variation criterion.
"""

from typing import Dict, List, Optional, Union

import numpy as np
from scipy.stats import norm


ArrayLike = Union[np.ndarray, List[float]]


def adf_test(series: ArrayLike, regression: str = "c") -> dict:
    """
    Run the Augmented Dickey-Fuller test for unit-root non-stationarity.

    Parameters
    ----------
    series : array_like
        1D numeric time series.
    regression : str, optional
        Trend/constant specification passed through to ``adfuller``
        (``"c"``, ``"ct"``, ``"ctt"`` or ``"n"``). Default ``"c"``.

    Returns
    -------
    dict
        Dictionary with keys ``statistic``, ``pvalue``, ``is_stationary``,
        ``critical_values`` and ``lags_used``. ``is_stationary`` is True
        when ``pvalue < 0.05`` (null of unit root rejected).

    Raises
    ------
    ImportError
        If statsmodels is not installed.
    """
    try:
        from statsmodels.tsa.stattools import adfuller
    except ImportError as exc:  # pragma: no cover - exercised only when missing
        raise ImportError(
            "adf_test requires statsmodels. Install it with: "
            "pip install statsmodels"
        ) from exc

    arr = np.asarray(series, dtype=float)
    statistic, pvalue, lags_used, _nobs, crit_values, _icbest = adfuller(
        arr, regression=regression
    )

    return {
        "statistic": float(statistic),
        "pvalue": float(pvalue),
        "is_stationary": bool(pvalue < 0.05),
        "critical_values": {k: float(v) for k, v in crit_values.items()},
        "lags_used": int(lags_used),
    }


def mann_kendall_test(series: ArrayLike, alpha: float = 0.05) -> dict:
    """
    Two-sided Mann-Kendall trend test with tie correction.

    Computes the S statistic as the sum of signs of all pairwise differences,
    applies the standard tie-adjusted variance, and converts to a z-score
    using the continuity correction. The two-sided p-value is obtained from
    the standard normal distribution.

    Parameters
    ----------
    series : array_like
        1D numeric time series.
    alpha : float, optional
        Significance level. Default 0.05.

    Returns
    -------
    dict
        Dictionary with keys ``trend`` (``"increasing"``, ``"decreasing"``
        or ``"no trend"``), ``S``, ``z``, ``pvalue`` and ``is_significant``.
    """
    x = np.asarray(series, dtype=float)
    n = x.size

    if n < 3:
        return {
            "trend": "no trend",
            "S": 0.0,
            "z": 0.0,
            "pvalue": 1.0,
            "is_significant": False,
        }

    # S statistic: sum of signs of all pairs (j > i).
    diffs = x[None, :] - x[:, None]
    iu = np.triu_indices(n, k=1)
    S = float(np.sum(np.sign(diffs[iu])))

    # Tie-corrected variance. Group ties by exact value.
    _unique, counts = np.unique(x, return_counts=True)
    tie_sum = float(np.sum(counts[counts > 1] * (counts[counts > 1] - 1)
                           * (2 * counts[counts > 1] + 5)))
    var_S = (n * (n - 1) * (2 * n + 5) - tie_sum) / 18.0

    # Continuity-corrected z-statistic.
    if var_S <= 0:
        z = 0.0
    elif S > 0:
        z = (S - 1) / np.sqrt(var_S)
    elif S < 0:
        z = (S + 1) / np.sqrt(var_S)
    else:
        z = 0.0

    pvalue = float(2.0 * norm.sf(abs(z)))
    is_significant = bool(pvalue < alpha)

    if is_significant and z > 0:
        trend = "increasing"
    elif is_significant and z < 0:
        trend = "decreasing"
    else:
        trend = "no trend"

    return {
        "trend": trend,
        "S": S,
        "z": float(z),
        "pvalue": pvalue,
        "is_significant": is_significant,
    }


def rolling_cv(series: ArrayLike, window: int) -> np.ndarray:
    """
    Rolling coefficient of variation (std / |mean|).

    Parameters
    ----------
    series : array_like
        1D numeric time series.
    window : int
        Window length (must be >= 2).

    Returns
    -------
    numpy.ndarray
        Array of the same length as ``series``. The first ``window - 1``
        entries are NaN. Windows whose mean is zero produce NaN.
    """
    if window < 2:
        raise ValueError("window must be at least 2")

    x = np.asarray(series, dtype=float)
    n = x.size
    out = np.full(n, np.nan, dtype=float)

    if n < window:
        return out

    for i in range(window - 1, n):
        w = x[i - window + 1 : i + 1]
        mean = np.mean(w)
        std = np.std(w, ddof=0)
        if mean == 0 or not np.isfinite(mean):
            out[i] = np.nan
        else:
            out[i] = std / abs(mean)

    return out


def find_convergence_time(
    series: ArrayLike,
    window: int = 20,
    cv_threshold: float = 0.05,
    require_stationary: bool = True,
) -> Optional[int]:
    """
    Find the first timestep at which the trajectory has converged.

    A timestep ``t`` is declared converged when (a) the rolling CV computed
    over ``series[t - window + 1 : t + 1]`` is below ``cv_threshold`` and
    (b) if ``require_stationary`` is True, the ADF test on the same window
    rejects the unit-root null.

    Parameters
    ----------
    series : array_like
        1D numeric time series.
    window : int, optional
        Window length used for both CV and ADF. Default 20.
    cv_threshold : float, optional
        Maximum allowable rolling CV. Default 0.05.
    require_stationary : bool, optional
        If True, also require the ADF test to flag the window as stationary.
        Default True.

    Returns
    -------
    int or None
        Index of the first converged timestep, or None if never converged.
    """
    x = np.asarray(series, dtype=float)
    n = x.size

    if n < window:
        return None

    cv = rolling_cv(x, window)

    for t in range(window - 1, n):
        cv_t = cv[t]
        if not np.isfinite(cv_t) or cv_t >= cv_threshold:
            continue

        if require_stationary:
            w = x[t - window + 1 : t + 1]
            # Skip degenerate (constant) windows: ADF cannot be evaluated,
            # but a flat window with low CV is trivially converged.
            if np.ptp(w) == 0:
                return int(t)
            try:
                result = adf_test(w)
            except Exception:
                continue
            if not result["is_stationary"]:
                continue

        return int(t)

    return None


class MultiMetricConvergence:
    """
    Aggregate convergence diagnostics across multiple ABM metrics.

    Parameters
    ----------
    metrics : dict
        Mapping from metric name to a 1D series (array or list).
    window : int, optional
        Window length used for each metric. Default 20.
    cv_threshold : float, optional
        Rolling CV threshold. Default 0.05.
    require_stationary : bool, optional
        Whether each metric must also pass ADF. Default True.

    Attributes
    ----------
    convergence_times : dict
        Per-metric convergence index (or None).
    global_convergence_time : int or None
        Maximum convergence index across all metrics, or None if any
        metric never converged.
    """

    def __init__(
        self,
        metrics: Dict[str, ArrayLike],
        window: int = 20,
        cv_threshold: float = 0.05,
        require_stationary: bool = True,
    ) -> None:
        self.metrics = {k: np.asarray(v, dtype=float) for k, v in metrics.items()}
        self.window = window
        self.cv_threshold = cv_threshold
        self.require_stationary = require_stationary

        self.convergence_times: Dict[str, Optional[int]] = {
            name: find_convergence_time(
                series,
                window=window,
                cv_threshold=cv_threshold,
                require_stationary=require_stationary,
            )
            for name, series in self.metrics.items()
        }

    @property
    def global_convergence_time(self) -> Optional[int]:
        """Maximum convergence time across metrics, or None if any failed."""
        if not self.convergence_times:
            return None
        times = list(self.convergence_times.values())
        if any(t is None for t in times):
            return None
        return int(max(times))  # type: ignore[arg-type]

    def summary(self) -> str:
        """
        Build a formatted multi-metric convergence report.

        Returns
        -------
        str
            Human-readable summary suitable for logging.
        """
        lines = []
        lines.append("=" * 60)
        lines.append("MULTI-METRIC CONVERGENCE REPORT")
        lines.append("=" * 60)
        lines.append(f"window={self.window}, cv_threshold={self.cv_threshold}, "
                     f"require_stationary={self.require_stationary}")
        lines.append("-" * 60)
        for name, t in self.convergence_times.items():
            label = "never" if t is None else f"t={t}"
            lines.append(f"  {name:30s} {label}")
        lines.append("-" * 60)
        gct = self.global_convergence_time
        gct_label = "never" if gct is None else f"t={gct}"
        lines.append(f"Global convergence time: {gct_label}")
        lines.append("=" * 60)
        return "\n".join(lines)
