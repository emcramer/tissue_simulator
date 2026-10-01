"""
Spatial organization models for single-sample tissue replicates (private).

A segmented tissue section often has structure that a stationary model cannot
express: layers (a tumor | margin | stroma sequence across the window) or
center-to-periphery progressions (necrotic core, proliferating rim, ...).
This module detects such a *trend* in cell density and/or cell-type
composition along one spatial coordinate and summarizes it as profiles that
:class:`~tissue_simulator.density.DensityModel` stores in its
``organization`` dictionary.

Module layout (WP1 adds the estimation half; later work packages add the
sampling half to this same module):

* Estimation (this file, WP1): :func:`fit_organization`, the candidate
  fitters, :func:`profile_along`, :func:`coordinate_along` and
  :func:`source_trend_maps`. Pure numpy/scipy, RNG-free, deterministic.
* Sampling (WP2/WP3, to be added below the marked section): drawing a new
  direction or center, evaluating the trend on a replicate window and
  accepting proposals. Nothing in this file draws random numbers.

Method
------
Coordinates. ``planar``: ``s = x cos(theta) + y sin(theta)`` (theta in
``[0, pi)``). ``radial``: ``s = |p - c|``. ``s`` is standardized to ``z`` over
the tissue-mask pixels before fitting.

Candidate trend models (each fitted for both geometries):

* ``density``: inhomogeneous Poisson, ``log lambda(s) = a + b z + c z**2``.
  The log-likelihood gain over a homogeneous pattern is
  ``sum_i g(z_i) - n log mean_pixels exp(g)`` (normalization on the mask
  grid, thinned to at most ``MAX_INTEGRAL_PIXELS`` pixels by a fixed stride).
* ``composition``: multinomial logit, ``logit p_t = a_t + b_t z`` (type 0 is
  the reference), against constant proportions.
* ``both``: sum of the two log-likelihood gains.

Model selection by BIC. ``delta_bic = 2 * delta_loglik - k log n`` is the BIC
advantage over "no trend" with ``k`` extra parameters (density 2, composition
``T - 1`` slopes, plus 1 for a planar angle or 2 for a radial center). A trend
is accepted only when the best candidate has ``delta_bic >= BIC_THRESHOLD``
(10, "very strong" evidence on the usual Kass-Raftery scale). A radial model
must beat the best planar candidate by the same margin; otherwise the radial
model is only recorded as ``weak_candidate = "radial"`` and the best planar
candidate (or no trend) is used.

The geometry (angle or center) is optimized once on the ``both`` objective and
shared by the three variants, which keeps the search cost to one 1-D (planar:
12 starts + bounded refinement) or 2-D (radial: 5x5 starts + Nelder-Mead, at
most ``CENTER_MAX_ITER`` iterations) optimization.

Profiles. For the selected model, density and composition are re-estimated
non-parametrically along ``s``: ``n_knots = clip(n // 40, 5, 25)`` equal-count
bins (knot = median ``s`` of the bin's cells); density is cells per unit mask
area relative to the overall mean; composition is
``(count + 0.5) / (n_bin + 0.5 T)``. Components that were not selected get a
flat profile. The transition width of a type is the ``s``-distance over which
its fraction rises from 25 % to 75 % of its range across knots.

Reported ``r2`` is the share of the binned (12 equal-count bins) log-likelihood
gain that the parametric trend captures, clipped to [0, 1].
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import optimize, special

# ---------------------------------------------------------------------------
# Named criteria (see the module docstring)
# ---------------------------------------------------------------------------

BIC_THRESHOLD = 10.0           # minimum delta BIC to accept a trend / radial over planar
MIN_CELLS = 60                 # fewer cells: no trend is attempted
MAX_INTEGRAL_PIXELS = 4000     # mask pixels used for the density normalization
N_ANGLE_STARTS = 12            # planar angle grid over [0, pi)
CENTER_GRID = 5                # radial center start grid (per axis)
CENTER_MAX_ITER = 200          # Nelder-Mead iteration cap
CELLS_PER_KNOT = 40
MIN_KNOTS, MAX_KNOTS = 5, 25
PROFILE_SMOOTHING = 0.5        # additive count per type in composition profiles
BINNED_R2_BINS = 12
TRANSITION_LOW, TRANSITION_HIGH = 0.25, 0.75
MIN_TRANSITION_RANGE = 0.02    # fraction range below which a width is undefined
_PARAM_BOUND = 30.0            # |density coefficient| bound
_LOGIT_BOUND = 20.0            # |composition coefficient| bound
_SLOPE_RIDGE = 1e-4            # tiny ridge on composition slopes (separation guard)

MODELS = ("none", "planar_density", "planar_composition", "planar_both",
          "radial_density", "radial_composition", "radial_both")


# ---------------------------------------------------------------------------
# Data container and coordinates
# ---------------------------------------------------------------------------

class _Data:
    """Cells and mask pixels used by all candidate fits."""

    def __init__(self, xs, ys, type_idx, n_types, mask, step, width, height):
        self.xs = np.asarray(xs, dtype=float)
        self.ys = np.asarray(ys, dtype=float)
        self.type_idx = np.asarray(type_idx, dtype=int)
        self.n_types = int(n_types)
        self.n = int(self.xs.size)
        self.step = float(step)
        self.width = float(width)
        self.height = float(height)
        rows, cols = np.nonzero(mask)
        self.pix_x = (cols + 0.5) * self.step
        self.pix_y = (rows + 0.5) * self.step
        stride = max(1, int(math.ceil(rows.size / MAX_INTEGRAL_PIXELS)))
        self.sub_x = self.pix_x[::stride]
        self.sub_y = self.pix_y[::stride]
        self.counts = np.bincount(self.type_idx, minlength=self.n_types).astype(float)
        self.null_comp_ll = float(sum(c * math.log(c / self.n) for c in self.counts if c > 0))


def coordinate_along(organization: Dict, x, y) -> np.ndarray:
    """Trend coordinate ``s`` of points for a fitted organization dict."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if organization.get("geometry") == "planar":
        dx, dy = organization["direction"]
        return x * dx + y * dy
    cx, cy = organization["center"]
    return np.hypot(x - cx, y - cy)


def _coordinate(geometry: str, param, x, y) -> np.ndarray:
    if geometry == "planar":
        return x * math.cos(param) + y * math.sin(param)
    return np.hypot(x - param[0], y - param[1])


def _standardized(data: _Data, geometry: str, param):
    s_pix = _coordinate(geometry, param, data.sub_x, data.sub_y)
    mu, sd = float(s_pix.mean()), float(s_pix.std())
    if sd <= 1e-12:
        return None
    z_pix = (s_pix - mu) / sd
    z_cell = (_coordinate(geometry, param, data.xs, data.ys) - mu) / sd
    return z_cell, z_pix


# ---------------------------------------------------------------------------
# Candidate fitters
# ---------------------------------------------------------------------------

def _density_gain(z_cell: np.ndarray, z_pix: np.ndarray) -> Tuple[float, np.ndarray]:
    """Log-likelihood gain of a quadratic log-intensity trend over homogeneity."""
    phi_cell = np.stack([z_cell, z_cell ** 2])
    phi_pix = np.stack([z_pix, z_pix ** 2])
    n, n_pix = z_cell.size, z_pix.size
    total = phi_cell.sum(axis=1)

    def neg(beta):
        g = beta @ phi_pix
        lse = special.logsumexp(g) - math.log(n_pix)
        w = np.exp(g - special.logsumexp(g))
        return -(float(beta @ total) - n * lse), -(total - n * (phi_pix @ w))

    res = optimize.minimize(neg, np.zeros(2), jac=True, method="L-BFGS-B",
                            bounds=[(-_PARAM_BOUND, _PARAM_BOUND)] * 2)
    return max(0.0, -float(res.fun)), res.x


def _composition_gain(z_cell: np.ndarray, type_idx: np.ndarray, n_types: int,
                      null_ll: float) -> Tuple[float, np.ndarray]:
    """Log-likelihood gain of multinomial-logit slopes over constant proportions."""
    if n_types < 2:
        return 0.0, np.zeros((0, 2))
    n = z_cell.size
    onehot = np.zeros((n, n_types))
    onehot[np.arange(n), type_idx] = 1.0
    design = np.column_stack([np.ones(n), z_cell])

    def neg(flat):
        coef = flat.reshape(n_types - 1, 2)
        eta = np.zeros((n, n_types))
        eta[:, 1:] = design @ coef.T
        lse = special.logsumexp(eta, axis=1)
        ll = float((eta * onehot).sum() - lse.sum())
        p = np.exp(eta - lse[:, None])
        grad = ((onehot - p)[:, 1:].T @ design)
        ridge = _SLOPE_RIDGE * float((coef[:, 1] ** 2).sum())
        grad[:, 1] -= 2.0 * _SLOPE_RIDGE * coef[:, 1]
        return -(ll) + ridge, -grad.ravel()

    start = np.zeros((n_types - 1, 2))
    props = np.maximum(np.bincount(type_idx, minlength=n_types), 0.5) / n
    start[:, 0] = np.log(props[1:] / props[0])
    res = optimize.minimize(neg, start.ravel(), jac=True, method="L-BFGS-B",
                            bounds=[(-_LOGIT_BOUND, _LOGIT_BOUND)] * (2 * (n_types - 1)))
    coef = res.x.reshape(n_types - 1, 2)
    eta = np.zeros((n, n_types))
    eta[:, 1:] = design @ coef.T
    ll = float((eta * onehot).sum() - special.logsumexp(eta, axis=1).sum())
    return max(0.0, ll - null_ll), coef


def _objective(data: _Data, geometry: str, param, components: str) -> float:
    """Log-likelihood gain used while optimizing the geometry."""
    zs = _standardized(data, geometry, param)
    if zs is None:
        return 0.0
    z_cell, z_pix = zs
    gain = 0.0
    if components in ("density", "both"):
        gain += _density_gain(z_cell, z_pix)[0]
    if components in ("composition", "both") and data.n_types > 1:
        gain += _composition_gain(z_cell, data.type_idx, data.n_types, data.null_comp_ll)[0]
    return gain


def _search_planar(data: _Data, components: str) -> float:
    angles = np.arange(N_ANGLE_STARTS) * math.pi / N_ANGLE_STARTS
    scores = [_objective(data, "planar", a, components) for a in angles]
    best = float(angles[int(np.argmax(scores))])
    half = math.pi / N_ANGLE_STARTS
    res = optimize.minimize_scalar(lambda a: -_objective(data, "planar", a, components),
                                   bounds=(best - half, best + half), method="bounded",
                                   options={"xatol": 0.01, "maxiter": 20})
    refined = float(res.x)
    if _objective(data, "planar", refined, components) < max(scores):
        refined = best
    return refined % math.pi


def _search_radial(data: _Data, components: str) -> Tuple[float, float]:
    fr = np.linspace(0.1, 0.9, CENTER_GRID)
    starts = [(fx * data.width, fy * data.height) for fy in fr for fx in fr]
    scores = [_objective(data, "radial", c, components) for c in starts]
    best = np.array(starts[int(np.argmax(scores))])

    def clipped(c):
        return (float(np.clip(c[0], 0.0, data.width)), float(np.clip(c[1], 0.0, data.height)))

    res = optimize.minimize(lambda c: -_objective(data, "radial", clipped(c), components),
                            best, method="Nelder-Mead",
                            options={"maxiter": CENTER_MAX_ITER, "xatol": data.step,
                                     "fatol": 1e-3,
                                     "initial_simplex": best + data.step * 4.0 * np.array(
                                         [[0, 0], [1, 0], [0, 1]], dtype=float)})
    refined = clipped(res.x)
    if _objective(data, "radial", refined, components) < max(scores):
        refined = (float(best[0]), float(best[1]))
    return refined


def _binned_gain(data: _Data, z_cell, z_pix, components: str) -> float:
    """Log-likelihood gain of a 12-bin step profile (denominator of ``r2``)."""
    order = np.argsort(z_cell, kind="stable")
    groups = np.array_split(order, min(BINNED_R2_BINS, max(1, data.n // 5)))
    edges = [z_cell[g[-1]] for g in groups[:-1]]
    pix_bin = np.searchsorted(edges, z_pix, side="left")
    gain = 0.0
    n_pix = z_pix.size
    for b, g in enumerate(groups):
        nb = g.size
        if components in ("density", "both"):
            ab = int((pix_bin == b).sum())
            if ab > 0 and nb > 0:
                gain += nb * math.log((nb / data.n) / (ab / n_pix))
        if components in ("composition", "both") and data.n_types > 1 and nb > 0:
            c = np.bincount(data.type_idx[g], minlength=data.n_types).astype(float)
            for t in range(data.n_types):
                if c[t] > 0:
                    gain += c[t] * math.log((c[t] / nb) / (data.counts[t] / data.n))
    return gain


def _evaluate_candidate(data: _Data, geometry: str, param, components: str) -> Optional[Dict]:
    zs = _standardized(data, geometry, param)
    if zs is None:
        return None
    z_cell, z_pix = zs
    use_d = components in ("density", "both")
    use_c = components in ("composition", "both") and data.n_types > 1
    delta = 0.0
    k = 1 if geometry == "planar" else 2
    if use_d:
        delta += _density_gain(z_cell, z_pix)[0]
        k += 2
    if use_c:
        delta += _composition_gain(z_cell, data.type_idx, data.n_types, data.null_comp_ll)[0]
        k += data.n_types - 1
    if not (use_d or use_c):
        return None
    delta_bic = 2.0 * delta - k * math.log(data.n)
    binned = _binned_gain(data, z_cell, z_pix,
                          "both" if (use_d and use_c) else ("density" if use_d else "composition"))
    cand = {
        "model": f"{geometry}_{components}", "geometry": geometry, "components": components,
        "delta_loglik": float(delta), "n_params": int(k),
        "bic": float(-delta_bic), "delta_bic": float(delta_bic),
        "delta_bic_per_cell": float(delta_bic / data.n),
        "r2": float(np.clip(delta / binned, 0.0, 1.0)) if binned > 1e-9 else 0.0,
    }
    if geometry == "planar":
        cand["theta"] = float(param)
        cand["direction"] = [math.cos(param), math.sin(param)]
    else:
        cand["center"] = [float(param[0]), float(param[1])]
    return cand


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------

def _transition_width(knots: np.ndarray, f: np.ndarray) -> Optional[float]:
    """Distance over which ``f`` rises from 25 % to 75 % of its range (None if flat)."""
    lo, hi = float(f.min()), float(f.max())
    if hi - lo < MIN_TRANSITION_RANGE:
        return None
    g = f if f[-1] >= f[0] else lo + hi - f
    l25, l75 = lo + TRANSITION_LOW * (hi - lo), lo + TRANSITION_HIGH * (hi - lo)

    def cross(i, level):  # interpolated crossing in segment (i, i+1)
        d = g[i + 1] - g[i]
        return knots[i] if abs(d) < 1e-12 else knots[i] + (level - g[i]) / d * (knots[i + 1] - knots[i])

    j = int(np.argmax(g >= l75))
    s75 = knots[0] if j == 0 else cross(j - 1, l75)
    below = np.flatnonzero(g[:j] <= l25)
    if below.size == 0:
        return float(abs(s75 - knots[0]))
    i = int(below[-1])
    s25 = cross(i, l25) if i + 1 < g.size else knots[i]
    return float(abs(s75 - s25))


def _profiles(data: _Data, s_cell: np.ndarray, s_all_pix: np.ndarray,
              use_d: bool, use_c: bool, mask_area_pixels: int) -> Dict:
    n_knots = int(np.clip(data.n // CELLS_PER_KNOT, MIN_KNOTS, MAX_KNOTS))
    order = np.argsort(s_cell, kind="stable")
    groups = np.array_split(order, n_knots)
    knots = np.array([np.median(s_cell[g]) for g in groups])
    edges = [0.5 * (s_cell[a[-1]] + s_cell[b[0]]) for a, b in zip(groups[:-1], groups[1:])]
    pix_bin = np.searchsorted(edges, s_all_pix, side="left")
    density = np.ones(n_knots)
    comp = np.tile(data.counts / data.n, (n_knots, 1))
    mean_density = data.n / mask_area_pixels
    for b, g in enumerate(groups):
        if use_d:
            area = int((pix_bin == b).sum())
            if area > 0:
                density[b] = (g.size / area) / mean_density
        if use_c:
            c = np.bincount(data.type_idx[g], minlength=data.n_types).astype(float)
            comp[b] = (c + PROFILE_SMOOTHING) / (g.size + PROFILE_SMOOTHING * data.n_types)
    widths = [_transition_width(knots, comp[:, t]) if use_c else None
              for t in range(data.n_types)]
    return {
        "s_knots": knots.tolist(),
        "density_profile": density.tolist(),
        "composition_profile": comp.tolist(),
        "s_range": [float(s_all_pix.min()), float(s_all_pix.max())],
        "s_mid": float(0.5 * (knots[0] + knots[-1])),
        "n_knots": n_knots,
        "transition_widths": widths,
        "gradient_strength": {
            "density_log_span": float(math.log(density.max() / max(density.min(), 1e-12)))
            if use_d else 0.0,
            "composition_span": float(np.ptp(comp, axis=0).max()) if use_c else 0.0,
            "s_extent": float(knots[-1] - knots[0]),
        },
    }


def profile_along(organization: Dict, s) -> Tuple[np.ndarray, np.ndarray]:
    """Density (relative to the mean) and composition at coordinate values ``s``.

    Piecewise-linear in the stored knots, clamped outside them. Returns
    ``(density (m,), composition (m, T))``; composition rows sum to 1.
    """
    s = np.atleast_1d(np.asarray(s, dtype=float))
    knots = np.asarray(organization["s_knots"], dtype=float)
    dens = np.interp(s, knots, np.asarray(organization["density_profile"], dtype=float))
    comp_k = np.asarray(organization["composition_profile"], dtype=float)
    comp = np.column_stack([np.interp(s, knots, comp_k[:, t]) for t in range(comp_k.shape[1])])
    comp = comp / comp.sum(axis=1, keepdims=True)
    return dens, comp


def source_trend_maps(organization: Dict, shape, step: float) -> Tuple[np.ndarray, np.ndarray]:
    """Trend on the source grid: ``(relative density (ny, nx), composition (T, ny, nx))``."""
    ny, nx = shape
    yy, xx = np.meshgrid((np.arange(ny) + 0.5) * step, (np.arange(nx) + 0.5) * step,
                         indexing="ij")
    dens, comp = profile_along(organization, coordinate_along(organization, xx.ravel(), yy.ravel()))
    return dens.reshape(shape), np.moveaxis(comp.reshape(shape + (comp.shape[1],)), -1, 0)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _none(reason: str, candidates: Optional[List[Dict]] = None,
          weak: Optional[str] = None) -> Dict:
    return {"model": "none", "candidates": candidates or [], "weak_candidate": weak,
            "fallback": reason}


def _candidate_search(data: _Data, comps: str, variants) -> List[Dict]:
    """Geometry search (planar angle, radial centre) and every candidate's fit."""
    candidates: List[Dict] = []
    p_theta = _search_planar(data, comps)
    r_center = _search_radial(data, comps)
    for geometry, param in (("planar", p_theta), ("radial", r_center)):
        for variant in variants:
            cand = _evaluate_candidate(data, geometry, param, variant)
            if cand is not None:
                candidates.append(cand)
    return candidates


def trend_statistic(xs, ys, type_idx, n_types, mask, grid_step, width, height) -> float:
    """Maximum ``delta_bic`` over all candidates (the trend-null test statistic).

    Runs exactly the candidate search of :func:`fit_organization` (same
    geometry search, same variants) but estimates no profiles. Returns
    ``-inf`` when the data are too few or no candidate is defined.
    """
    data = _Data(xs, ys, type_idx, n_types, mask, grid_step, width, height)
    if data.n < MIN_CELLS:
        return float("-inf")
    comps = "both" if data.n_types > 1 else "density"
    variants = ("density", "composition", "both") if data.n_types > 1 else ("density",)
    candidates = _candidate_search(data, comps, variants)
    return max((c["delta_bic"] for c in candidates), default=float("-inf"))


def fit_organization(xs, ys, type_idx, n_types, mask, grid_step, width, height) -> Dict:
    """Select and summarize a planar or radial trend in density and composition.

    Args:
        xs, ys: Cell centers in µm, relative to the window origin.
        type_idx: Integer cell-type index per cell.
        n_types: Number of cell types.
        mask: Boolean tissue mask on the model grid (``[row, col] = [y, x]``).
        grid_step, width, height: Grid pixel size and window size in µm.

    Returns:
        A JSON-serializable dict. Keys: ``model`` (one of :data:`MODELS`),
        ``geometry``, ``components``, ``direction``/``theta`` or ``center``,
        ``candidates`` (every fitted candidate with ``delta_loglik``,
        ``n_params``, ``bic`` (BIC minus the no-trend BIC, lower is better),
        ``delta_bic`` (its negative), ``r2`` and geometry), profile keys
        (``s_knots``, ``density_profile``, ``composition_profile``,
        ``s_range``, ``s_mid``, ``n_knots``, ``transition_widths``,
        ``gradient_strength``), ``weak_candidate`` (``"radial"`` or None),
        ``fallback`` (reason string when the model is ``"none"`` or a radial
        candidate was downgraded, else None) and ``bic_threshold``.
    """
    data = _Data(xs, ys, type_idx, n_types, mask, grid_step, width, height)
    if data.n < MIN_CELLS:
        return _none("too_few_cells")
    comps = "both" if data.n_types > 1 else "density"
    variants = ("density", "composition", "both") if data.n_types > 1 else ("density",)

    candidates = _candidate_search(data, comps, variants)
    if not candidates:
        return _none("degenerate_window")

    best = max(candidates, key=lambda c: c["delta_bic"])
    weak, fallback = None, None
    if best["delta_bic"] < BIC_THRESHOLD:
        chosen, fallback = None, "no_candidate_reached_bic_threshold"
    elif best["geometry"] == "radial":
        planar = [c for c in candidates if c["geometry"] == "planar"]
        best_planar = max(planar, key=lambda c: c["delta_bic"]) if planar else None
        if best_planar is None or best["delta_bic"] - best_planar["delta_bic"] >= BIC_THRESHOLD:
            chosen = best
        else:
            weak = "radial"
            if best_planar["delta_bic"] >= BIC_THRESHOLD:
                chosen, fallback = best_planar, "radial_not_decisive_over_planar"
            else:
                chosen, fallback = None, "no_candidate_reached_bic_threshold"
    else:
        chosen = best
    if chosen is None:
        return _none(fallback, candidates, weak)
    geometry = chosen["geometry"]
    param = chosen["theta"] if geometry == "planar" else tuple(chosen["center"])
    s_cell = _coordinate(geometry, param, data.xs, data.ys)
    s_pix = _coordinate(geometry, param, data.pix_x, data.pix_y)
    use_d = chosen["components"] in ("density", "both")
    use_c = chosen["components"] in ("composition", "both") and data.n_types > 1
    out = {
        "model": chosen["model"], "geometry": geometry, "components": chosen["components"],
        "candidates": candidates, "weak_candidate": weak, "fallback": fallback,
        "bic_threshold": BIC_THRESHOLD, "delta_bic": chosen["delta_bic"],
        "delta_bic_per_cell": chosen["delta_bic_per_cell"], "r2": chosen["r2"],
        "window_centroid": [float(data.pix_x.mean()), float(data.pix_y.mean())],
    }
    if geometry == "planar":
        out["theta"], out["direction"] = chosen["theta"], chosen["direction"]
    else:
        out["center"] = chosen["center"]
    out.update(_profiles(data, s_cell, s_pix, use_d, use_c, int(data.pix_x.size)))
    return out


# ---------------------------------------------------------------------------
# Sampling (added by later work packages) goes below this line.
# ---------------------------------------------------------------------------

MIN_COVERAGE = 0.9             # fraction of the source knot range the new window must span
PROPORTION_ABS_TOL = 0.02      # accepted |p'_t - p_t| floor (absolute)
PROPORTION_REL_TOL = 0.15      # accepted |p'_t - p_t| <= max(abs, rel * p_t)
DEFAULT_MAX_PROPOSALS = 20     # direction re-draws before the best-coverage fallback
_LOG_FLOOR = 1e-12             # floor before taking logs of residual maps
CENTER_MARGIN = 0.10           # radial centers are drawn from the window shrunk by this fraction per side


def sample_direction(rng: np.random.Generator) -> float:
    """Draw a planar trend angle ``theta' ~ U(0, 2 pi)`` (one RNG draw)."""
    return float(rng.uniform(0.0, 2.0 * math.pi))


def sample_center(rng: np.random.Generator, shape, step: float) -> Tuple[float, float]:
    """Draw a radial trend centre ``c' ~ U(window shrunk by CENTER_MARGIN per side)``.

    One RNG call (two values), made before the noise fields like the planar angle.
    """
    ny, nx = shape
    w, h = nx * step, ny * step
    cx, cy = rng.uniform([CENTER_MARGIN * w, CENTER_MARGIN * h],
                         [(1 - CENTER_MARGIN) * w, (1 - CENTER_MARGIN) * h])
    return float(cx), float(cy)


def _evaluate_on(organization: Dict, s: np.ndarray, shape) -> Tuple[np.ndarray, np.ndarray, float]:
    dens, comp = profile_along(organization, s.ravel())
    log_density = np.log(np.maximum(dens, _LOG_FLOOR)).reshape(shape)
    logits = np.moveaxis(np.log(np.maximum(comp, _LOG_FLOOR)).reshape(shape + (comp.shape[1],)),
                         -1, 0)
    knots = np.asarray(organization["s_knots"], dtype=float)
    span = float(knots[-1] - knots[0])
    if span <= 0:
        return log_density, logits, 1.0
    overlap = min(float(s.max()), knots[-1]) - max(float(s.min()), knots[0])
    return log_density, logits, float(np.clip(overlap / span, 0.0, 1.0))


def evaluate_radial_trend(organization: Dict, center, shape, step: float
                          ) -> Tuple[np.ndarray, np.ndarray, float]:
    """Source radial trend around a new centre ``c'`` on a new window.

    ``s' = |p - c'|`` (no offset: radial coordinates are anchored at the centre,
    so the source ``s`` axis maps directly). Returns the same triple as
    :func:`evaluate_trend`; coverage is the fraction of the source knot range
    ``[s_min, s_max]`` represented in the window's range of ``s'``.
    """
    ny, nx = shape
    yy, xx = np.meshgrid((np.arange(ny) + 0.5) * step, (np.arange(nx) + 0.5) * step,
                         indexing="ij")
    return _evaluate_on(organization, np.hypot(xx - center[0], yy - center[1]), shape)


def evaluate_trend(organization: Dict, theta: float, shape, step: float
                   ) -> Tuple[np.ndarray, np.ndarray, float]:
    """Source trend re-expressed along a new planar direction on a new window.

    The coordinate is ``s' = (p - c') . d' + s_mid`` with ``c'`` the window
    centre and ``d' = (cos theta, sin theta)``: the source profile's mid-knot
    lands on the window centroid, so physical gradient strength and transition
    widths (in µm) are preserved.

    Returns:
        ``(log_density (ny, nx) [log of density relative to the mean],
        logits (T, ny, nx) [log composition profile], coverage)`` where
        coverage is the fraction of the source knot range ``[k_0, k_last]``
        lying within the window's range of ``s'``.
    """
    ny, nx = shape
    yy, xx = np.meshgrid((np.arange(ny) + 0.5) * step, (np.arange(nx) + 0.5) * step,
                         indexing="ij")
    cx, cy = 0.5 * nx * step, 0.5 * ny * step
    s = (xx - cx) * math.cos(theta) + (yy - cy) * math.sin(theta) + organization["s_mid"]
    dens, comp = profile_along(organization, s.ravel())
    log_density = np.log(np.maximum(dens, _LOG_FLOOR)).reshape(shape)
    logits = np.moveaxis(np.log(np.maximum(comp, _LOG_FLOOR)).reshape(shape + (comp.shape[1],)),
                         -1, 0)
    knots = np.asarray(organization["s_knots"], dtype=float)
    span = float(knots[-1] - knots[0])
    if span <= 0:
        return log_density, logits, 1.0
    overlap = min(float(s.max()), knots[-1]) - max(float(s.min()), knots[0])
    return log_density, logits, float(np.clip(overlap / span, 0.0, 1.0))


def combine_trend(log_density, trend_logits, resid_intensity, resid_composition,
                  proportions, n_target: int, step: float) -> Tuple[np.ndarray, np.ndarray]:
    """Multiply the trend onto the residual maps and normalize.

    ``intensity = exp(log_density) * residual intensity`` (scaled to
    ``n_target`` cells); ``composition = softmax(trend logits + log residual
    composition - log proportions)``. The residual composition describes the
    detrended data and equals the global proportions when it has no structure,
    so the ``- log proportions`` term avoids counting the prior twice.
    Composition is finite and sums to 1 per pixel.
    """
    intensity = np.exp(log_density) * resid_intensity
    mass = intensity.sum() * step ** 2
    if mass > 0:
        intensity = intensity * (n_target / mass)
    z = (trend_logits + np.log(np.maximum(resid_composition, _LOG_FLOOR))
         - np.log(np.maximum(np.asarray(proportions, dtype=float), _LOG_FLOOR))[:, None, None])
    z = z - z.max(axis=0, keepdims=True)
    e = np.exp(z)
    return intensity, e / e.sum(axis=0, keepdims=True)


def expected_proportions(intensity: np.ndarray, composition: np.ndarray) -> np.ndarray:
    """Expected type proportions, ``sum(intensity * composition) / sum(intensity)``."""
    w = (composition * intensity[None]).reshape(composition.shape[0], -1).sum(axis=1)
    total = w.sum()
    return w / total if total > 0 else np.full(composition.shape[0], 1.0 / composition.shape[0])


def accept_proposal(coverage: float, expected, target) -> bool:
    """Window-aware acceptance.

    Accept when ``coverage >= MIN_COVERAGE`` (0.9) and every type satisfies
    ``|p'_t - p_t| <= max(0.02, 0.15 p_t)``.
    """
    expected = np.asarray(expected, dtype=float)
    target = np.asarray(target, dtype=float)
    tol = np.maximum(PROPORTION_ABS_TOL, PROPORTION_REL_TOL * target)
    return bool(coverage >= MIN_COVERAGE and np.all(np.abs(expected - target) <= tol))


def propose_layout(organization: Dict, theta0: float, rng: np.random.Generator, shape,
                   step: float,
                   resid_intensity, resid_composition, proportions, n_target: int,
                   max_proposals: int = DEFAULT_MAX_PROPOSALS, resid_fn=None):
    """Draw directions (planar) or centres (radial) until one is accepted.

    ``theta0`` is the first direction (planar) or centre ``(cx, cy)`` (radial), drawn by the caller *before* the noise
    fields; rejected proposals re-draw only the direction/centre from ``rng`` (the
    residual maps are reused).

    ``resid_fn(theta) -> (resid_intensity, resid_composition, n_target)``, when
    given, replaces the fixed residual maps (void placement that depends on the
    proposed centre, so every proposal re-places the anchored hole).

    Returns:
        ``(intensity, composition, info)``; ``info`` has ``theta`` (the angle, or the centre for radial),
        ``proposals_tried``, ``accepted``, ``fallback`` (``None`` or
        ``"best_of_proposals"``, the best-coverage proposal then being used).
    """
    radial = organization.get("geometry") != "planar"
    best = None
    tried = 0
    for k in range(max(1, int(max_proposals))):
        if k == 0:
            theta = theta0
        else:
            theta = sample_center(rng, shape, step) if radial else sample_direction(rng)
        tried += 1
        if resid_fn is not None:
            resid_intensity, resid_composition, n_target = resid_fn(theta)
        if radial:
            log_d, logits, coverage = evaluate_radial_trend(organization, theta, shape, step)
        else:
            log_d, logits, coverage = evaluate_trend(organization, theta, shape, step)
        intensity, composition = combine_trend(log_d, logits, resid_intensity,
                                               resid_composition, proportions, n_target, step)
        ok = accept_proposal(coverage, expected_proportions(intensity, composition),
                             proportions)
        if ok:
            return intensity, composition, {"theta": theta, "proposals_tried": tried,
                                            "accepted": True, "fallback": None}
        if best is None or coverage > best[0]:
            best = (coverage, theta, intensity, composition)
    return best[2], best[3], {"theta": best[1], "proposals_tried": tried,
                              "accepted": False, "fallback": "best_of_proposals"}
