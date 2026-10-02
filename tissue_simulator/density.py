"""
Density-aware layouts for tissue replicates.

A :class:`DensityModel` summarizes where a source tissue region is dense or
sparse and which cell types dominate where. It is fitted once from the
region's cells and then samples one :class:`Layout` per replicate, which
:class:`~tissue_simulator.packing.InhomogeneousPacker` fills with cells.

Fitting has four steps:

1. Kernel intensity maps per cell type (bandwidth by likelihood
   cross-validation).
2. Compartments: k-means on each pixel's per-type intensities (for example
   dense tumor nests, sparse stroma, a mixed margin), so splits fall on
   blurred boundaries. Each compartment's density and composition come from
   the cells inside it, profiled by distance to the compartment edge
   (0-10, 10-20, 20-40 and >40 µm) so margins such as an immune cuff around
   nests are kept.
3. Patch structure: the correlation length and smoothness of two latent
   Gaussian fields whose thresholded maps reproduce how often two points at a
   given distance fall in the same compartment (calibrated by simulation).
4. Radius marks and the hard core, from the observed cells.

Two layout modes are supported:

* ``"resample"`` draws a new arrangement of compartments with the region's
  compartment areas, density distribution within each compartment,
  composition and patch structure. Replicates differ in where the nests are,
  not only in cell positions.
* ``"copy"`` reuses the region's own density and composition maps, so
  replicates differ only at scales below the smoothing bandwidth.

Resampling assumes the region is statistically stationary. Regions with a
trend across the window, or patches as large as the window, are flagged in
:attr:`DensityModel.flags` and are better served by ``"copy"``.

All maps are two-dimensional (x, y); cell z positions are left to the packer.
"""

import math
import warnings
from dataclasses import dataclass, field, fields
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import fft as sp_fft
from scipy import ndimage, special
from scipy.spatial import cKDTree
from scipy.stats import rankdata

from . import _organization
from . import _units
from . import _grains
from . import _voids

_QUANTILE_LEVELS = np.linspace(0.0, 1.0, 201)
_N_RADIUS_BINS = 5
_TREND_R2_LIMIT = 0.3
_SMOOTHNESS_GRID = (0.5, 1.0, 2.0, 4.0)
_N_LENGTHS = 10
_CALIBRATION_REPEATS = 2
# Distance bands (µm) from a pixel to the edge of its compartment, used to
# profile density and composition near boundaries (margins, invasive fronts).
_BAND_EDGES_UM = (10.0, 20.0, 40.0)
_MIN_BAND_CELLS = 10
# Adaptive detrending regularization: per-cell inverse-trend weights are clipped
# to [1/cap, cap] (then renormalized to mean 1) so rare-type cells in layers
# where the trend says they are rare do not dominate residual maps, and
# per-(compartment, band) compositions are shrunk toward the global
# proportions with this many pseudo-cells (Dirichlet prior).
_DETREND_WEIGHT_CAP = 5.0
_COMPOSITION_PRIOR_CELLS = 10

# Adaptive-strategy criteria (all deterministic, RNG-free).
# Per-type bandwidths: a type with fewer than _TYPE_BW_MIN_CELLS cells, or whose
# cross-validated likelihood gain over the pooled bandwidth is below
# _TYPE_BW_MIN_GAIN nats, keeps the pooled bandwidth ("pooled_fallback");
# otherwise log h_t = (n_t log h_cv + m log h_pool) / (n_t + m) with
# m = _TYPE_BW_SHRINK_M. Evidence is "cv" when n_t / (n_t + m) >= _TYPE_BW_CV_WEIGHT
# and "shrunk" otherwise.
_TYPE_BW_MIN_CELLS = 10
_TYPE_BW_MIN_GAIN = 2.0
_TYPE_BW_SHRINK_M = 30.0
_TYPE_BW_CV_WEIGHT = 0.9
# Size model: Gaussian of log radius per (type, density bin), shrunk toward the
# type's pooled Gaussian with m = _SIZE_SHRINK_M pseudo-cells; types with fewer
# than _SIZE_MIN_TYPE_CELLS cells use the pooled-over-types Gaussian.
_SIZE_SHRINK_M = 20.0
_SIZE_MIN_TYPE_CELLS = 10
_MIN_LOG_SIGMA = 0.02
# Learned band edges: quantiles of the cells' distance to their compartment edge.
_BAND_EDGE_QUANTILES = (0.25, 0.5, 0.75)
_DETREND_FLOOR = 1e-3

UNIT_OWNS_LUMEN = True
"""Precedence: a void whose centroid lies inside a unit's outer radius is not placed
by the voids path (the unit's lumen is); diagnostics keep it."""

UNIT_PLACEMENT_TRIES = 500
"""Center draws per grain in :func:`_grains.place_grains` (cheap; the source's
smallest center separation can be a tight constraint for several large units)."""

UNIT_OWNS_RADIAL = True
"""Precedence: a radial trend centered inside a unit's outer radius is skipped for
sampling (layout.organization = {"model": "none", "fallback": "units"}); planar
trends coexist with units (grains are applied after the trend)."""


# ---------------------------------------------------------------------------
# Small numerical helpers
# ---------------------------------------------------------------------------

def _as_generator(rng) -> np.random.Generator:
    """Accept a Generator, SeedSequence, int or None."""
    return rng if isinstance(rng, np.random.Generator) else np.random.default_rng(rng)


def _grid_shape(width: float, height: float, step: float) -> Tuple[int, int]:
    return (max(1, int(math.ceil(height / step))), max(1, int(math.ceil(width / step))))


def _pixel_indices(xs, ys, step: float, shape: Tuple[int, int]):
    ny, nx = shape
    ix = np.clip(np.floor(np.asarray(xs) / step).astype(int), 0, nx - 1)
    iy = np.clip(np.floor(np.asarray(ys) / step).astype(int), 0, ny - 1)
    return iy, ix


def _smooth(image: np.ndarray, sigma: float) -> np.ndarray:
    return ndimage.gaussian_filter(image, sigma=sigma, mode='constant', truncate=4.0)


def _largest_remainder(weights, total: int) -> np.ndarray:
    """Integer allocation of ``total`` proportional to ``weights``."""
    w = np.clip(np.asarray(weights, dtype=float), 0.0, None)
    if w.sum() <= 0:
        w = np.ones_like(w)
    raw = w / w.sum() * total
    counts = np.floor(raw).astype(int)
    remainder = int(total - counts.sum())
    if remainder > 0:
        order = np.argsort(-(raw - counts), kind='stable')
        counts[order[:remainder]] += 1
    return counts


def _normal_scores(values: np.ndarray) -> np.ndarray:
    """Rank-based transform to standard-normal scores (ties share a score)."""
    ranks = rankdata(values, method='average')
    return special.ndtri((ranks - 0.5) / values.size)


def _standardize(values: np.ndarray) -> np.ndarray:
    sd = values.std()
    return (values - values.mean()) / sd if sd > 0 else np.zeros_like(values)


def _helmert_basis(n: int) -> np.ndarray:
    """Orthonormal (n-1) x n basis used for isometric log-ratio coordinates."""
    basis = np.zeros((n - 1, n))
    for k in range(1, n):
        basis[k - 1, :k] = 1.0 / k
        basis[k - 1, k] = -1.0
        basis[k - 1] *= math.sqrt(k / (k + 1.0))
    return basis


# ---------------------------------------------------------------------------
# Intensity maps and compartments
# ---------------------------------------------------------------------------

def _intensity_grids(iy, ix, type_idx, n_types, shape, step, bandwidth, mask, weights=None):
    """Edge-corrected Gaussian kernel intensity per type, in cells per µm².

    ``bandwidth`` is a scalar or one value per type; ``weights`` (per cell,
    default 1) are used to detrend counts before smoothing.
    """
    if np.ndim(bandwidth) > 0:  # one bandwidth per type
        sigmas = np.asarray(bandwidth, dtype=float) / step
        edges = {}
        grids = np.zeros((n_types,) + tuple(shape))
        for t in range(n_types):
            sigma = float(sigmas[t])
            if sigma not in edges:
                edges[sigma] = np.where(
                    mask, np.maximum(_smooth(mask.astype(float), sigma), 1e-6), 1.0)
            counts = np.zeros(shape)
            sel = type_idx == t
            np.add.at(counts, (iy[sel], ix[sel]), 1.0 if weights is None else weights[sel])
            grids[t] = np.where(mask, _smooth(counts, sigma) / edges[sigma], 0.0) / step ** 2
        return grids
    sigma = bandwidth / step
    edge = np.where(mask, np.maximum(_smooth(mask.astype(float), sigma), 1e-6), 1.0)
    grids = np.zeros((n_types,) + tuple(shape))
    for t in range(n_types):
        counts = np.zeros(shape)
        sel = type_idx == t
        np.add.at(counts, (iy[sel], ix[sel]), 1.0 if weights is None else weights[sel])
        grids[t] = np.where(mask, _smooth(counts, sigma) / edge, 0.0) / step ** 2
    return grids


def _bandwidth_scores(iy, ix, shape, step, mask, candidates):
    """Leave-one-out likelihood cross-validation score of each candidate bandwidth."""
    counts = np.zeros(shape)
    np.add.at(counts, (iy, ix), 1.0)
    floor = 1e-3 * iy.size / (mask.sum() * step ** 2)
    scores = []
    for h in candidates:
        sigma = h / step
        r = int(math.ceil(4.0 * sigma))
        delta = np.zeros((2 * r + 1, 2 * r + 1))
        delta[r, r] = 1.0
        self_weight = _smooth(delta, sigma)[r, r]
        edge = np.maximum(_smooth(mask.astype(float), sigma), 1e-6)
        smoothed = _smooth(counts, sigma)
        leave_one_out = (smoothed[iy, ix] - self_weight) / edge[iy, ix] / step ** 2
        integral = (smoothed / edge)[mask].sum()
        scores.append(float(np.log(np.maximum(leave_one_out, floor)).sum() - integral))
    return scores


def _select_bandwidth(iy, ix, shape, step, mask, bandwidth_range, n_candidates=10):
    """Likelihood cross-validation bandwidth for the total intensity.

    Returns ``(bandwidth, at_bound)`` where ``at_bound`` is ``"lower"`` or
    ``"upper"`` when the best candidate is an end of ``bandwidth_range``
    (the range may be too narrow) and None otherwise.
    """
    lo, hi = bandwidth_range
    candidates = np.geomspace(lo, hi, n_candidates)
    scores = _bandwidth_scores(iy, ix, shape, step, mask, candidates)
    best = int(np.argmax(scores))
    at_bound = "lower" if best == 0 else ("upper" if best == len(scores) - 1 else None)
    return float(candidates[best]), at_bound


def _auto_bandwidth_range(xs, ys, radii, mask, step) -> Tuple[float, float]:
    """Data-derived bandwidth search range for the adaptive strategy.

    With ``d_nn`` the median nearest-neighbor distance, ``r_med`` the median
    radius and ``L`` the shorter side of the tissue-mask bounding box:
    ``lo = max(step, d_nn, 2 r_med)``, ``hi = max(min(L / 3, 20 d_nn), 3 lo)``.
    """
    points = np.column_stack([xs, ys])
    d_nn = float(np.median(cKDTree(points).query(points, k=2)[0][:, 1]))
    r_med = float(np.median(radii))
    rows, cols = np.nonzero(mask)
    side = min((np.ptp(rows) + 1), (np.ptp(cols) + 1)) * step
    lo = max(float(step), d_nn, 2.0 * r_med)
    hi = max(min(side / 3.0, 20.0 * d_nn), 3.0 * lo)
    return float(lo), float(hi)


def _type_bandwidths(iy, ix, type_idx, n_types, shape, step, mask, bandwidth_range,
                     h_pool, n_candidates=10):
    """Regularized per-type bandwidths and the evidence behind each.

    Returns ``(array (T,), evidence list)``; see the ``_TYPE_BW_*`` constants.
    """
    lo, hi = bandwidth_range
    grid = list(np.geomspace(lo, hi, n_candidates))
    candidates = np.array(grid + [h_pool])
    out = np.full(n_types, float(h_pool))
    evidence = []
    for t in range(n_types):
        sel = type_idx == t
        n_t = int(sel.sum())
        if n_t < _TYPE_BW_MIN_CELLS:
            evidence.append("pooled_fallback")
            continue
        scores = np.array(_bandwidth_scores(iy[sel], ix[sel], shape, step, mask, candidates))
        best = int(np.argmax(scores[:-1]))
        if scores[best] - scores[-1] < _TYPE_BW_MIN_GAIN:
            evidence.append("pooled_fallback")
            continue
        w = n_t / (n_t + _TYPE_BW_SHRINK_M)
        out[t] = float(np.clip(np.exp(w * math.log(candidates[best])
                                      + (1.0 - w) * math.log(h_pool)), lo, hi))
        evidence.append("cv" if w >= _TYPE_BW_CV_WEIGHT else "shrunk")
    return out, evidence


def _size_model(radii, type_idx, n_types, bins, n_bins):
    """Gaussian of log radius per (type, density bin), shrunk toward the type's pooled fit.

    ``bins`` is each cell's density bin (the same bins as the radius marks).
    Shrinkage adds ``_SIZE_SHRINK_M`` pseudo-cells of the type-level mean and
    variance; types with fewer than ``_SIZE_MIN_TYPE_CELLS`` cells use the
    all-cell fit (per bin, shrunk to the all-cell pooled fit). Returns
    ``(mu (T, B), sigma (T, B), fallback type indices)``. RNG-free.
    """
    logr = np.log(np.maximum(np.asarray(radii, dtype=float), 1e-6))
    m = _SIZE_SHRINK_M

    def fit(x, b):
        mu0 = float(x.mean())
        var0 = float(x.var())
        mu = np.full(n_bins, mu0)
        var = np.full(n_bins, var0)
        for k in range(n_bins):
            xb = x[b == k]
            if xb.size:
                mu[k] = (xb.sum() + m * mu0) / (xb.size + m)
                var[k] = (xb.size * float(xb.var()) + m * var0) / (xb.size + m)
        return mu, np.sqrt(np.maximum(var, _MIN_LOG_SIGMA ** 2))

    all_mu, all_sigma = fit(logr, bins)
    mu_out = np.zeros((n_types, n_bins))
    sigma_out = np.zeros((n_types, n_bins))
    fallbacks = []
    for t in range(n_types):
        sel = type_idx == t
        if sel.sum() < _SIZE_MIN_TYPE_CELLS:
            mu_out[t], sigma_out[t] = all_mu, all_sigma
            fallbacks.append(t)
        else:
            mu_out[t], sigma_out[t] = fit(logr[sel], bins[sel])
    return mu_out, sigma_out, fallbacks


def _layout_features(lam_types: np.ndarray, mask: np.ndarray):
    """Density and leading-composition score of every in-mask pixel.

    These are the two latent axes of resampled layouts. Returns
    ``(density, composition_score, loading)`` where ``loading`` is the
    leading principal direction of the composition fractions.
    """
    n_types = lam_types.shape[0]
    density = lam_types.sum(axis=0)[mask]
    if n_types < 2:
        return density, np.zeros_like(density), np.zeros(0)
    c = 1e-3 * float(density.mean()) / n_types + 1e-300
    centered = ((lam_types[:, mask] + c) / (density + n_types * c)).T
    centered = centered - centered.mean(axis=0)
    loading = np.linalg.svd(centered, full_matrices=False)[2][0]
    if loading[np.argmax(np.abs(loading))] < 0:
        loading = -loading
    return density, centered @ loading, loading


def _kmeans(points: np.ndarray, k: int, rng: np.random.Generator,
            n_init: int = 4, max_iter: int = 100) -> np.ndarray:
    """Seeded k-means++ / Lloyd labels, numbered by ascending first coordinate."""
    n = points.shape[0]
    k = max(1, min(k, n))
    best = None
    for _ in range(n_init):
        centers = points[[rng.integers(n)]]
        for _ in range(1, k):
            d2 = ((points[:, None, :] - centers[None]) ** 2).sum(-1).min(axis=1)
            total = d2.sum()
            idx = rng.integers(n) if total <= 0 else rng.choice(n, p=d2 / total)
            centers = np.vstack([centers, points[idx]])
        for _ in range(max_iter):
            labels = ((points[:, None, :] - centers[None]) ** 2).sum(-1).argmin(axis=1)
            updated = np.array([points[labels == j].mean(axis=0) if np.any(labels == j)
                                else centers[j] for j in range(k)])
            if np.allclose(updated, centers):
                break
            centers = updated
        d2 = ((points[:, None, :] - centers[None]) ** 2).sum(-1)
        labels = d2.argmin(axis=1)
        inertia = d2[np.arange(n), labels].sum()
        if best is None or inertia < best[0]:
            best = (inertia, centers, labels)
    _, centers, labels = best

    used = np.flatnonzero(np.bincount(labels, minlength=k) > 0)
    centers = centers[used]
    labels = np.searchsorted(used, labels)
    order = np.lexsort((centers[:, 1], centers[:, 0]))
    remap = np.empty(order.size, dtype=int)
    remap[order] = np.arange(order.size)
    return remap[labels]


def _assign_compartments(points: np.ndarray, centers: np.ndarray,
                         fractions: np.ndarray) -> np.ndarray:
    """Nearest-center labels with offsets so label counts hit ``fractions`` exactly."""
    n, k = points.shape[0], centers.shape[0]
    cost = ((points[:, None, :] - centers[None]) ** 2).sum(-1)
    targets = _largest_remainder(fractions, n)
    offsets = np.zeros(k)
    rate = 4.0
    for _ in range(200):
        labels = (cost - offsets).argmin(axis=1)
        diff = targets - np.bincount(labels, minlength=k)
        if not diff.any():
            return labels
        offsets += rate * diff / n
        rate *= 0.98

    adjusted = cost - offsets
    labels = adjusted.argmin(axis=1)
    counts = np.bincount(labels, minlength=k)
    while True:
        over = np.flatnonzero(counts > targets)
        if not over.size:
            return labels
        i, j = over[0], np.flatnonzero(counts < targets)[0]
        members = np.flatnonzero(labels == i)
        regret = adjusted[members, j] - adjusted[members, i]
        n_move = min(counts[i] - targets[i], targets[j] - counts[j])
        moved = members[np.argsort(regret, kind='stable')[:n_move]]
        labels[moved] = j
        counts[i] -= n_move
        counts[j] += n_move


def _band_index(labels: np.ndarray, step: float, edges=_BAND_EDGES_UM) -> np.ndarray:
    """Distance band of each pixel from the edge of its own compartment (-1: no tissue).

    Band ``b`` holds pixels whose distance to the nearest pixel of another
    compartment (or of no tissue) is at most ``edges[b]`` (default
    ``_BAND_EDGES_UM``); the last band is the interior. Window edges do not
    count as compartment edges.
    """
    band = np.full(labels.shape, -1, dtype=int)
    for c in np.unique(labels[labels >= 0]):
        inside = labels == c
        distance = ndimage.distance_transform_edt(inside) * step
        band[inside] = np.searchsorted(edges, distance[inside], side='left')
    return band


def _edge_distance_map(labels: np.ndarray, step: float) -> np.ndarray:
    """Distance (µm) of each tissue pixel to the edge of its compartment (0 elsewhere)."""
    out = np.zeros(labels.shape)
    for c in np.unique(labels[labels >= 0]):
        inside = labels == c
        out[inside] = (ndimage.distance_transform_edt(inside) * step)[inside]
    return out


def _learned_band_edges(distances: np.ndarray, step: float, patch_length: float) -> Tuple[float, ...]:
    """Band edges from quantiles 0.25/0.5/0.75 of cell-to-compartment-edge distance.

    Clipped to ``[step, max(patch_length, step)]``, rounded to the grid and
    deduplicated.
    """
    q = np.quantile(distances, _BAND_EDGE_QUANTILES)
    q = np.clip(q, step, max(patch_length, step))
    q = np.maximum(np.round(q / step) * step, step)
    return tuple(float(v) for v in np.unique(q))


# ---------------------------------------------------------------------------
# Latent fields and patch calibration
# ---------------------------------------------------------------------------

def _spectrum(w2: np.ndarray, length: float, bandwidth: float, nu: float = 1.0) -> np.ndarray:
    """2D Matérn spectrum (smoothness ``nu``) seen through the Gaussian smoothing kernel."""
    return np.exp(-(bandwidth ** 2) * w2) / (length ** -2 + w2) ** (nu + 1.0)


def _noise_shape(shape, pad: int) -> Tuple[int, int]:
    return (sp_fft.next_fast_len(shape[0] + 2 * pad), sp_fft.next_fast_len(shape[1] + 2 * pad))


def _padding(shape, step: float, length: float, bandwidth: float) -> int:
    return min(int(math.ceil(2.0 * max(length, bandwidth) / step)), 4 * max(shape))


def _simulate_field(shape, step: float, length: float, bandwidth: float, nu: float,
                    noise: np.ndarray, pad: int) -> np.ndarray:
    """Stationary Gaussian random field with the model spectrum, via padded FFT.

    ``noise`` is standard-normal white noise of shape ``_noise_shape(shape, pad)``;
    reusing it across parameter values gives common random numbers.
    """
    ny, nx = shape
    big_y, big_x = noise.shape
    wy = 2 * np.pi * sp_fft.fftfreq(big_y, d=step)
    wx = 2 * np.pi * sp_fft.fftfreq(big_x, d=step)
    amplitude = np.sqrt(_spectrum(wy[:, None] ** 2 + wx[None, :] ** 2, length, bandwidth, nu))
    amplitude[0, 0] = 0.0
    simulated = sp_fft.ifft2(sp_fft.fft2(noise) * amplitude).real
    return simulated[pad:pad + ny, pad:pad + nx]


def _sample_labels(shape, step, length, nu, bandwidth, rho, centers, fractions,
                   noise_u, noise_w, pad):
    """Compartment labels (flat) and density-axis scores of one simulated layout."""
    u = _standardize(_simulate_field(shape, step, length, bandwidth, nu, noise_u, pad))
    w = _standardize(_simulate_field(shape, step, length, bandwidth, nu, noise_w, pad))
    v = rho * u + math.sqrt(1.0 - rho ** 2) * w
    u_scores = _normal_scores(u.ravel())
    points = np.column_stack([u_scores, _normal_scores(v.ravel())])
    if np.ptp(centers[:, 1]) <= 1e-12:
        points[:, 1] = centers[0, 1]
    return _assign_compartments(points, centers, fractions), u_scores


def _annulus_index(shape) -> np.ndarray:
    """Rounded lag distance, in pixels, of each entry of an FFT-ordered lag map."""
    ly = sp_fft.fftfreq(shape[0], d=1.0 / shape[0])
    lx = sp_fft.fftfreq(shape[1], d=1.0 / shape[1])
    return np.rint(np.sqrt(ly[:, None] ** 2 + lx[None, :] ** 2)).astype(int)


def _same_label_profile(labels: np.ndarray, n_labels: int, n_lags: int) -> np.ndarray:
    """Probability that two pixels ``r`` pixels apart share a label, for r < n_lags.

    Pixels labeled -1 (no tissue) are ignored; zero-padding avoids wrap-around.
    """
    ny, nx = labels.shape
    shape = (sp_fft.next_fast_len(ny + n_lags), sp_fft.next_fast_len(nx + n_lags))

    def autocorrelation(image):
        padded = np.zeros(shape)
        padded[:ny, :nx] = image
        return sp_fft.ifft2(np.abs(sp_fft.fft2(padded)) ** 2).real

    same = sum(autocorrelation(labels == k) for k in range(n_labels))
    pairs = autocorrelation(labels >= 0)
    index = _annulus_index(shape)
    sel = index < n_lags
    return (np.bincount(index[sel], weights=same[sel], minlength=n_lags)
            / np.maximum(np.bincount(index[sel], weights=pairs[sel], minlength=n_lags), 1e-12))


def _calibrate_patches(region_labels, n_labels, centers, fractions, rho, step, bandwidth,
                       lower, upper, rng):
    """Latent correlation length and smoothness matching the region's patch structure.

    Candidates are scored by the squared difference between the region's
    same-compartment probability profile (lags up to ``upper``) and the mean
    profile of a few simulated layouts that share white noise across
    candidates. Returns ``(length, nu, at_lower_bound, at_upper_bound)``.
    """
    shape = region_labels.shape
    n_lags = int(upper // step) + 1
    target = _same_label_profile(region_labels, n_labels, n_lags)[1:]
    pad = _padding(shape, step, upper, bandwidth)
    noises = [(rng.standard_normal(_noise_shape(shape, pad)),
               rng.standard_normal(_noise_shape(shape, pad)))
              for _ in range(_CALIBRATION_REPEATS)]
    lengths = np.geomspace(lower, upper, _N_LENGTHS)
    best = None
    for nu in _SMOOTHNESS_GRID:
        for i, length in enumerate(lengths):
            profiles = [
                _same_label_profile(
                    _sample_labels(shape, step, length, nu, bandwidth, rho, centers,
                                   fractions, noise_u, noise_w, pad)[0].reshape(shape),
                    n_labels, n_lags)[1:]
                for noise_u, noise_w in noises
            ]
            score = float(((np.mean(profiles, axis=0) - target) ** 2).sum())
            if best is None or score < best[0]:
                best = (score, float(length), float(nu), i)
    _, length, nu, i = best
    return length, nu, i == 0, i == lengths.size - 1


def _heterogeneity_statistics(iy, ix, type_idx, n_types, shape, step, bandwidth, mask,
                              comp_bandwidth=None):
    """Variance of log density and density-weighted composition departure.

    The density statistic always uses the pooled ``bandwidth``. When
    ``comp_bandwidth`` (one value per type, adaptive strategy) is given, the
    composition statistic is computed from per-type maps at those bandwidths
    so that compositional blobs smaller than the pooled bandwidth are not
    blurred away; None keeps the single-map (legacy) computation exactly.
    """
    lam = _intensity_grids(iy, ix, type_idx, n_types, shape, step, bandwidth, mask)
    total = lam.sum(axis=0)[mask]
    mean = total.mean()
    t_density = float(np.var(np.log(total + 0.01 * mean))) if mean > 0 else 0.0
    if n_types < 2 or total.sum() <= 0:
        return t_density, 0.0
    if comp_bandwidth is not None:
        lam = _intensity_grids(iy, ix, type_idx, n_types, shape, step, comp_bandwidth, mask)
        total = lam.sum(axis=0)[mask]
        if total.sum() <= 0:
            return t_density, 0.0
    global_p = np.bincount(type_idx, minlength=n_types) / type_idx.size
    comp = lam[:, mask] / np.maximum(total, 1e-12)
    safe = np.where(global_p > 0, global_p, 1.0)
    chi = (((comp - global_p[:, None]) ** 2) / safe[:, None]).sum(axis=0)
    return t_density, float((chi * total).sum() / total.sum())


# ---------------------------------------------------------------------------
# Public data structures
# ---------------------------------------------------------------------------

@dataclass
class RadiusMarks:
    """Empirical cell radii grouped by local-density quantile bin.

    Attributes:
        density_edges: Interior bin edges in cells per µm² (ascending).
        radii_by_bin: Radii (µm) observed in each bin; ``len(density_edges) + 1`` bins.
    """
    density_edges: np.ndarray
    radii_by_bin: List[np.ndarray]

    def draw(self, rng: np.random.Generator, local_density: float) -> float:
        """Draw a radius observed at a similar local density."""
        b = int(np.searchsorted(self.density_edges, local_density, side='right'))
        radii = self.radii_by_bin[b]
        return float(radii[rng.integers(radii.size)])

    @property
    def median_radius(self) -> float:
        return float(np.median(np.concatenate(self.radii_by_bin)))


@dataclass
class Layout:
    """Target density and composition maps for one replicate.

    Arrays are indexed ``[row, col] = [y, x]`` on a square grid of side
    ``grid_step`` µm starting at the origin.

    Attributes:
        width, height: Window size in µm (x and y extents).
        grid_step: Pixel size in µm.
        intensity: Expected cells per µm², shape ``(ny, nx)``.
        composition: Expected type fractions, shape ``(n_types, ny, nx)``.
        compartment: Compartment index per pixel (-1 outside tissue).
        cell_types: Type names, ordered as in ``composition``.
        marks: Density-conditioned radius distribution.
        kappa: Hard-core factor; centers must be at least ``kappa * (r_i + r_j)`` apart.
        n_target: Number of cells to place.
        target_overlap_fraction: Fraction of source cells whose nearest
            neighbor is closer than the hard core; relaxation stops at this level.
        mode: ``"resample"``, ``"copy"`` or ``"uniform"``.
        flags: Diagnostics inherited from the density model.
        bandwidth: Smoothing bandwidth of the maps in µm (sets packer bin size).
        strategy: ``"legacy"`` or ``"adaptive"`` (from the model).
        organization: For organized layouts, ``{model, geometry, direction,
            theta, proposals_tried, accepted, fallback}``; empty otherwise.
        voids: Placed voids ``{centers, radii, n_requested, n_placed, anchored}``
            (µm); empty when the model has none. ``"voids_unsatisfied"`` is in
            ``flags`` when fewer holes than requested could be placed.
        quota_scale: Adaptive quota bin size in µm, ``max(bandwidth / 2,
            max(2 d_nn, 2 grid_step, sqrt(4 / mean_density)))`` snapped to the
            grid; None for legacy layouts.
        units: Placed germ-grain units (``place_grains`` dict: ``centers``,
            ``outer_radii``, ``inner_radii``, ``kinds``, ``n_requested``,
            ``n_placed``, ``shortfall``, ``anchored``); empty when the model
            has no units. ``"units_unsatisfied"`` is in ``flags`` on shortfall.
    """
    width: float
    height: float
    grid_step: float
    intensity: np.ndarray
    composition: np.ndarray
    compartment: np.ndarray
    cell_types: Tuple[str, ...]
    marks: RadiusMarks
    kappa: float
    n_target: int
    target_overlap_fraction: float
    mode: str
    flags: Tuple[str, ...] = ()
    bandwidth: Optional[float] = None
    strategy: str = "legacy"
    organization: Dict = field(default_factory=dict)
    quota_scale: Optional[float] = None
    voids: Dict = field(default_factory=dict)
    units: Dict = field(default_factory=dict)

    def pixel(self, x: float, y: float) -> Tuple[int, int]:
        ny, nx = self.intensity.shape
        return (min(max(int(y // self.grid_step), 0), ny - 1),
                min(max(int(x // self.grid_step), 0), nx - 1))

    def intensity_at(self, x: float, y: float) -> float:
        return float(self.intensity[self.pixel(x, y)])

    def composition_at(self, x: float, y: float) -> np.ndarray:
        iy, ix = self.pixel(x, y)
        return self.composition[:, iy, ix]


# ---------------------------------------------------------------------------
# Density model
# ---------------------------------------------------------------------------

def _fit_maps(xs, ys, type_idx, n_types, iy, ix, shape, grid_step, mask, bandwidth,
              type_bw, bw_evidence, proportions, n_compartments, rng, org=None):
    """Per-type maps, compartments, quantiles, bands and projection for one pass.

    Shared by the legacy and adaptive strategies. When ``org`` is a selected
    trend (``org["model"] != "none"``) cells are weighted by the inverse trend
    (detrending) and the trend is multiplied back into the stored maps;
    otherwise the maps are the raw (stationary) estimates. Draws RNG only in
    the k-means step, as the previously inlined code did.

    Returns a dict of everything :meth:`DensityModel.fit` needs from the pass.
    """
    weights = None
    trend_factor = None
    if org is not None and org["model"] != "none":
        dens_map, comp_map = _organization.source_trend_maps(org, shape, grid_step)
        trend_factor = np.maximum(dens_map[None] * comp_map / proportions[:, None, None],
                                  _DETREND_FLOOR)
        dens_c, comp_c = _organization.profile_along(
            org, _organization.coordinate_along(org, xs, ys))
        weights = 1.0 / np.maximum(
            dens_c * comp_c[np.arange(xs.size), type_idx] / proportions[type_idx],
            _DETREND_FLOOR)
        weights = np.clip(weights, 1.0 / _DETREND_WEIGHT_CAP, _DETREND_WEIGHT_CAP)
        weights = weights * (weights.size / weights.sum())
    lam_types = _intensity_grids(iy, ix, type_idx, n_types, shape, grid_step,
                                 bandwidth if bw_evidence is None else type_bw, mask,
                                 weights)
    lam_total = lam_types.sum(axis=0)
    
    # Compartments: k-means on per-type intensities, which are linear in
    # cell counts, so splits fall midway between compartments on their
    # blurred boundaries. (Composition fractions or log scales would push
    # splits into the sparser side; normal scores would flatten them.)
    density, composition_score, loading = _layout_features(lam_types, mask)
    log_density = np.log(density + 0.01 * density.mean() + 1e-300)
    pixel_labels = _kmeans(lam_types[:, mask].T / max(float(density.std()), 1e-300),
                           n_compartments, rng)
    k = int(pixel_labels.max()) + 1
    rank = np.empty(k, dtype=int)
    rank[np.argsort([density[pixel_labels == c].mean() for c in range(k)], kind='stable')] = np.arange(k)
    pixel_labels = rank[pixel_labels]
    region_compartments = np.full(shape, -1, dtype=int)
    region_compartments[mask] = pixel_labels
    fractions = np.bincount(pixel_labels, minlength=k) / pixel_labels.size
    
    # Resampling works on normal scores; express the partition there.
    u = _normal_scores(log_density)
    v = (_normal_scores(composition_score) if np.ptp(composition_score) > 1e-12
         else np.zeros_like(u))
    centers = np.array([[u[pixel_labels == c].mean(), v[pixel_labels == c].mean()]
                        for c in range(k)])
    rho = (float(np.clip(np.corrcoef(u, v)[0, 1], -0.99, 0.99))
           if np.ptp(u) > 1e-12 and np.ptp(v) > 1e-12 else 0.0)
    
    # Each compartment's density level and composition come from the cells
    # inside it, which undoes the kernel's blurring across compartment
    # boundaries; the smoothed map only supplies within-compartment variation.
    cell_compartment = region_compartments[iy, ix]
    density_quantiles = np.zeros((k, _QUANTILE_LEVELS.size))
    composition = np.zeros((k, n_types))
    for c in range(k):
        sel = region_compartments == c
        in_c = cell_compartment == c
        smoothed = lam_total[sel]
        if in_c.any() and smoothed.mean() > 0:
            if weights is None:
                factor = (in_c.sum() / (sel.sum() * grid_step ** 2)) / smoothed.mean()
                composition[c] = np.bincount(type_idx[in_c], minlength=n_types) / in_c.sum()
            else:
                wc = weights[in_c]
                factor = (wc.sum() / (sel.sum() * grid_step ** 2)) / smoothed.mean()
                composition[c] = (np.bincount(type_idx[in_c], weights=wc, minlength=n_types)
                                  + _COMPOSITION_PRIOR_CELLS * proportions
                                  ) / (wc.sum() + _COMPOSITION_PRIOR_CELLS)
        else:
            factor = 0.0
            composition[c] = proportions
        lam_types[:, sel] *= factor
        density_quantiles[c] = np.quantile(smoothed * factor, _QUANTILE_LEVELS)
    lam_total = lam_types.sum(axis=0)
    
    # Boundary profiles: density and composition by distance to the
    # compartment edge, again counted from cells (falls back to the
    # compartment's values where a band holds too few cells).
    def band_stats(edges):
        band = _band_index(region_compartments, grid_step, edges)
        cell_band = band[iy, ix]
        n_bands = len(edges) + 1
        band_quantiles = np.repeat(density_quantiles[:, None, :], n_bands, axis=1)
        band_composition = np.repeat(composition[:, None, :], n_bands, axis=1)
        for c in range(k):
            for b in range(n_bands):
                sel = (region_compartments == c) & (band == b)
                in_cb = (cell_compartment == c) & (cell_band == b)
                smoothed = lam_total[sel]
                if in_cb.sum() < _MIN_BAND_CELLS or smoothed.mean() <= 0:
                    continue
                if weights is None:
                    band_density = in_cb.sum() / (sel.sum() * grid_step ** 2)
                    band_comp = np.bincount(type_idx[in_cb], minlength=n_types) / in_cb.sum()
                else:
                    wb = weights[in_cb]
                    band_density = wb.sum() / (sel.sum() * grid_step ** 2)
                    band_comp = (np.bincount(type_idx[in_cb], weights=wb, minlength=n_types)
                                 + _COMPOSITION_PRIOR_CELLS * proportions
                                 ) / (wb.sum() + _COMPOSITION_PRIOR_CELLS)
                band_quantiles[c, b] = np.quantile(
                    smoothed * (band_density / smoothed.mean()), _QUANTILE_LEVELS)
                band_composition[c, b] = band_comp
        return band_quantiles, band_composition
    
    band_quantiles, band_composition = band_stats(_BAND_EDGES_UM)
    
    adjacency = np.zeros((k, k))
    for a, b in ((region_compartments[:, :-1], region_compartments[:, 1:]),
                 (region_compartments[:-1, :], region_compartments[1:, :])):
        both = (a >= 0) & (b >= 0)
        np.add.at(adjacency, (a[both], b[both]), 1.0)
    adjacency = adjacency + adjacency.T
    adjacency /= max(adjacency.sum(), 1.0)
    
    # Put the trend back so stored maps and radius marks describe the real
    # (not detrended) density; compartments and bands above stay residual.
    region_lam = lam_types
    if trend_factor is not None:
        region_lam = lam_types * trend_factor
        region_lam *= type_idx.size / (region_lam.sum() * grid_step ** 2)

    return dict(lam_types=lam_types, lam_total=lam_total, log_density=log_density,
            region_compartments=region_compartments, fractions=fractions, centers=centers,
            rho=rho, k=k, density_quantiles=density_quantiles, composition=composition,
            band_quantiles=band_quantiles, band_composition=band_composition,
            adjacency=adjacency, region_lam=region_lam, trend_factor=trend_factor,
            weights=weights, loading=loading, band_stats=band_stats)


_MAP_KEYS = ("lam_types", "lam_total", "log_density", "region_compartments", "fractions",
             "centers", "rho", "k", "density_quantiles", "composition", "band_quantiles",
             "band_composition", "adjacency", "region_lam", "trend_factor", "loading",
             "band_stats")


@dataclass
class DensityModel:
    """Fitted density and compartment model of a source tissue region.

    Build with :meth:`fit` or :meth:`from_tissue`; sample replicate layouts
    with :meth:`sample_layout`. See the module docstring for the method.

    Attributes (selected):
        bandwidth: Kernel bandwidth of the intensity maps (µm).
        compartment_fractions: Area fraction of each compartment.
        compartment_composition: Cell-type fractions per compartment.
        patch_length, patch_smoothness: Calibrated latent-field correlation
            length (µm) and Matérn smoothness that set patch sizes in
            resampled layouts.
        kappa: Hard-core factor from nearest-neighbor spacing.
        homogeneous: True when the region is not significantly more
            heterogeneous than a uniform packing; layouts are then uniform.
        flags: ``"trend"``, ``"patch_length_at_upper_bound"`` (both mean the
            region may not be stationary) or ``"patch_length_at_lower_bound"``.
            Adaptive fits may add ``"bandwidth_at_lower_bound"``,
            ``"bandwidth_at_upper_bound"`` and ``"organization:<model>"``.
        strategy: ``"legacy"`` or ``"adaptive"`` (how the model was fitted).
        type_bandwidths: Per-type kernel bandwidths (µm); all equal to
            ``bandwidth`` unless fitted adaptively.
        band_edges: Distance-band edges (µm) from the compartment edge;
            ``(10, 20, 40)`` unless learned (adaptive).
        size_log_mu, size_log_sigma: Mean and standard deviation of log radius
            per (type, density bin of ``marks.density_edges``), shrunk toward
            the type-level fit.
        estimation: Estimation metadata (bandwidth range and whether the best
            candidate sat at its bound, per-type bandwidth evidence, type
            counts, size-model fallbacks).
        organization: Selected spatial trend model (see
            :mod:`tissue_simulator._organization`); ``{"model": "none"}`` when
            none was fitted or accepted.
        voids: Void (lumen/hole) summary from :func:`_voids.infer_voids` or
            :func:`_voids.summarize_components` (``n_holes``,
            ``equivalent_diameters``, ``min_center_separation``, ``tau``, ...);
            empty when voids are off. Layouts place ``n_holes`` discs.
    """
    width: float
    height: float
    grid_step: float
    bandwidth: float
    cell_types: Tuple[str, ...]
    n_cells: int
    proportions: np.ndarray
    mask: np.ndarray
    region_intensity: np.ndarray
    region_compartments: np.ndarray
    compartment_centers: np.ndarray
    compartment_fractions: np.ndarray
    compartment_density_quantiles: np.ndarray
    compartment_composition: np.ndarray
    compartment_adjacency: np.ndarray
    band_density_quantiles: np.ndarray
    band_composition: np.ndarray
    composition_loading: np.ndarray
    rho: float
    marks: RadiusMarks
    kappa: float
    target_overlap_fraction: float
    patch_length: float = 0.0
    patch_smoothness: float = 1.0
    flags: Tuple[str, ...] = ()
    homogeneous: bool = False
    heterogeneity: Dict[str, float] = field(default_factory=dict)
    strategy: str = "legacy"
    type_bandwidths: np.ndarray = field(default_factory=lambda: np.zeros(0))
    band_edges: Tuple[float, ...] = _BAND_EDGES_UM
    size_log_mu: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    size_log_sigma: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    estimation: Dict = field(default_factory=dict)
    organization: Dict = field(default_factory=lambda: {"model": "none"})
    voids: Dict = field(default_factory=dict)
    units: Dict = field(default_factory=dict)

    # -- construction -------------------------------------------------------

    @classmethod
    def fit(cls, x: Sequence[float], y: Sequence[float], radius: Sequence[float],
            cell_type: Sequence[str],
            bounds: Optional[Tuple[float, float, float, float]] = None,
            mask: Optional[np.ndarray] = None,
            bandwidth: Optional[float] = None,
            bandwidth_range=None,
            grid_step: float = 5.0,
            n_compartments: int = 3,
            patch_prior: Optional[Dict[str, float]] = None,
            n_null: int = 19,
            alpha: float = 0.05,
            seed=None,
            *,
            strategy: str = "legacy",
            per_type_bandwidth: Optional[bool] = None,
            organization: Optional[bool] = None,
            voids: Optional[str] = None,
            units: Optional[str] = None) -> 'DensityModel':
        """Fit a density model to a 2D point pattern of cells.

        Args:
            x, y: Cell centers in µm.
            radius: Cell radii in µm (e.g. ``sqrt(area / pi)`` from segmentation).
            cell_type: Cell type label per cell.
            bounds: ``(xmin, ymin, xmax, ymax)`` of the observation window;
                defaults to the bounding box of the cells.
            mask: Optional boolean tissue mask on the model grid
                (``ceil(height / grid_step)`` by ``ceil(width / grid_step)``).
                False marks areas with no tissue (e.g. tears or artifacts);
                they are excluded from every estimate.
            bandwidth: Kernel standard deviation in µm; chosen by likelihood
                cross-validation within ``bandwidth_range`` when None.
            bandwidth_range: ``(lo, hi)`` search range in µm, ``"auto"`` for a
                data-derived range, or None (default): ``(15, 80)`` under the
                legacy strategy, auto-derived under ``"adaptive"``.
            grid_step: Pixel size of the maps in µm.
            n_compartments: Number of compartments (k-means clusters on log
                density and the leading composition axis).
            patch_prior: Optional ``{"length": µm, "smoothness": nu}`` to use
                instead of calibrating patch structure (for example pooled
                over a cohort with :func:`pooled_patch_prior`).
            n_null: Uniform packings used to test whether the region is
                heterogeneous at all. 0 skips the test.
            alpha: Significance level of that test.
            seed: Seed for k-means, calibration and the heterogeneity test.
            strategy: ``"legacy"`` (default; results are unchanged from earlier
                versions) or ``"adaptive"``, which turns on per-type
                bandwidths, trend (organization) selection with detrending,
                and learned boundary bands. The size model and estimation
                metadata are filled in both strategies.
            per_type_bandwidth: Regularized per-type bandwidths; defaults to
                ``strategy == "adaptive"``.
            organization: Fit a planar/radial trend model and detrend before
                compartments; defaults to ``strategy == "adaptive"``. The
                result is stored in :attr:`organization` and ignored by layouts
                of homogeneous models.
            voids: ``"auto"`` or ``"none"``; defaults to ``"auto"`` under
                ``"adaptive"`` and ``"none"`` under ``"legacy"``. With ``"auto"``
                and no ``mask``, empty regions (lumens, holes) are inferred from
                the cells and become the model mask; with a ``mask`` its
                interior holes are summarized instead. Layouts then place the
                same number of holes (sizes resampled from the source) at new
                positions.
            units: ``"auto"`` or ``"none"``; defaults to ``"auto"`` under
                ``"adaptive"`` and ``"none"`` under ``"legacy"``. Detects compact
                units (blobs, rings with a core or lumen) by persistent homology
                of the per-type intensity maps, with thresholds from stationary
                null draws (drawn after all other fit draws); stored in
                :attr:`units`.
        """
        if units is None:
            units = "auto" if strategy == "adaptive" else "none"
        if units not in ("auto", "none"):
            raise ValueError(f"units must be None, 'auto' or 'none', got {units!r}.")
        if voids is None:
            voids = "auto" if strategy == "adaptive" else "none"
        if voids not in ("auto", "none"):
            raise ValueError(f"voids must be None, 'auto' or 'none', got {voids!r}.")
        if strategy not in ("legacy", "adaptive"):
            raise ValueError(f"strategy must be 'legacy' or 'adaptive', got {strategy!r}.")
        adaptive = strategy == "adaptive"
        use_per_type = adaptive if per_type_bandwidth is None else bool(per_type_bandwidth)
        use_org = adaptive if organization is None else bool(organization)
        xs = np.asarray(x, dtype=float)
        ys = np.asarray(y, dtype=float)
        radii = np.asarray(radius, dtype=float)
        labels = np.asarray(cell_type).astype(str)
        if not (xs.size == ys.size == radii.size == labels.size):
            raise ValueError("x, y, radius and cell_type must have the same length.")
        if xs.size < 10:
            raise ValueError("At least 10 cells are needed to fit a density model.")

        if bounds is None:
            bounds = (xs.min(), ys.min(), xs.max(), ys.max())
        x0, y0, x1, y1 = (float(b) for b in bounds)
        width, height = x1 - x0, y1 - y0
        if width <= 0 or height <= 0:
            raise ValueError(f"Invalid bounds {bounds!r}.")
        xs, ys = xs - x0, ys - y0

        shape = _grid_shape(width, height, grid_step)
        mask_given = mask is not None
        if mask is None:
            mask = np.ones(shape, dtype=bool)
        else:
            mask = np.asarray(mask, dtype=bool)
            if mask.shape != shape:
                raise ValueError(f"mask must have shape {shape}, got {mask.shape}.")
            if not mask.any():
                raise ValueError("mask excludes the whole window.")

        inside = (xs >= 0) & (xs <= width) & (ys >= 0) & (ys <= height)
        iy, ix = _pixel_indices(xs, ys, grid_step, shape)
        keep = inside & mask[iy, ix]
        xs, ys, radii, labels, iy, ix = (a[keep] for a in (xs, ys, radii, labels, iy, ix))
        if xs.size < 10:
            raise ValueError("Fewer than 10 cells fall inside the window and mask.")
        void_info = {}
        if voids == "auto":
            # Before anything else uses the mask: the inferred tissue mask is
            # the model mask. Pixels holding a cell stay tissue.
            if not mask_given:
                mask, void_info = _voids.infer_voids(xs, ys, radii, None, grid_step,
                                                     width, height)
                mask = mask.copy()
                # Only interior holes are re-placed in layouts, so edge-touching
                # void components stay tissue in the mask (they remain in
                # ``void_info`` diagnostics); otherwise n_target is biased.
                vlab, nvl = ndimage.label(~mask)
                for k in range(1, nvl + 1):
                    vy, vx = np.nonzero(vlab == k)
                    if (vy.min() == 0 or vx.min() == 0 or vy.max() == shape[0] - 1
                            or vx.max() == shape[1] - 1):
                        mask[vy, vx] = True
                mask[iy, ix] = True
            else:
                void_info = _voids.summarize_components(~mask, grid_step)
                diam = void_info["equivalent_diameters"]
                # No inference tau for a user mask: separate holes by one median diameter.
                void_info.update(tau=float(np.median(diam)) / 2.0 if diam else 0.0,
                                 inferred=False)

        cell_types = tuple(sorted(set(labels.tolist())))
        type_idx = np.searchsorted(cell_types, labels)
        n_types = len(cell_types)
        rng = _as_generator(seed)

        if bandwidth_range is None:
            range_source = "auto" if adaptive else "legacy_default"
            bandwidth_range = (_auto_bandwidth_range(xs, ys, radii, mask, grid_step)
                               if adaptive else (15.0, 80.0))
        elif isinstance(bandwidth_range, str):
            if bandwidth_range != "auto":
                raise ValueError(f"bandwidth_range must be None, 'auto' or (lo, hi), "
                                 f"got {bandwidth_range!r}.")
            range_source = "auto"
            bandwidth_range = _auto_bandwidth_range(xs, ys, radii, mask, grid_step)
        else:
            range_source = "user"
        bandwidth_range = (float(bandwidth_range[0]), float(bandwidth_range[1]))
        at_bound = None
        if bandwidth is None:
            bandwidth, at_bound = _select_bandwidth(iy, ix, shape, grid_step, mask,
                                                    bandwidth_range)
            bandwidth_source = "cv"
        else:
            bandwidth_source = "user"
        bandwidth = float(bandwidth)
        type_bw = np.full(n_types, bandwidth)
        bw_evidence = None
        if use_per_type and bandwidth_source == "cv":
            type_bw, bw_evidence = _type_bandwidths(iy, ix, type_idx, n_types, shape,
                                                    grid_step, mask, bandwidth_range, bandwidth)
        proportions = np.bincount(type_idx, minlength=n_types) / type_idx.size

        org = {"model": "none"}
        if use_org:
            org = _organization.fit_organization(xs, ys, type_idx, n_types, mask,
                                                 grid_step, width, height)
        maps = _fit_maps(xs, ys, type_idx, n_types, iy, ix, shape, grid_step, mask,
                         bandwidth, type_bw, bw_evidence, proportions, n_compartments,
                         rng, org=org)
        (lam_types, lam_total, log_density, region_compartments, fractions, centers, rho, k,
         density_quantiles, composition, band_quantiles, band_composition, adjacency,
         region_lam, trend_factor, loading, band_stats) = (
            maps[n] for n in ("lam_types", "lam_total", "log_density", "region_compartments",
                              "fractions", "centers", "rho", "k", "density_quantiles",
                              "composition", "band_quantiles", "band_composition",
                              "adjacency", "region_lam", "trend_factor", "loading",
                              "band_stats"))
        def build(maps, org, het=None, quiet=False):
            (lam_types, lam_total, log_density, region_compartments, fractions, centers, rho,
             k, density_quantiles, composition, band_quantiles, band_composition, adjacency,
             region_lam, trend_factor, loading, band_stats) = (
                maps[n] for n in _MAP_KEYS)
            # Radius marks and hard core from the observed cells.
            lam_at_cells = region_lam.sum(axis=0)[iy, ix] if trend_factor is not None \
                else lam_total[iy, ix]
            edges = np.quantile(lam_at_cells, np.linspace(0, 1, _N_RADIUS_BINS + 1)[1:-1])
            bins = np.searchsorted(edges, lam_at_cells, side='right')
            radii_by_bin = [radii[bins == b] if np.any(bins == b) else radii.copy()
                            for b in range(_N_RADIUS_BINS)]
            size_mu, size_sigma, size_fallback = _size_model(radii, type_idx, n_types, bins,
                                                             _N_RADIUS_BINS)
            points = np.column_stack([xs, ys])
            distances, neighbors = cKDTree(points).query(points, k=2)
            ratio = distances[:, 1] / np.maximum(radii + radii[neighbors[:, 1]], 1e-12)
            kappa = float(np.clip(np.quantile(ratio, 0.02), 0.5, 1.5))

            model = cls(
                width=width, height=height, grid_step=float(grid_step), bandwidth=bandwidth,
                cell_types=cell_types, n_cells=int(xs.size), proportions=proportions,
                mask=mask, region_intensity=region_lam, region_compartments=region_compartments,
                compartment_centers=centers, compartment_fractions=fractions,
                compartment_density_quantiles=density_quantiles,
                compartment_composition=composition, compartment_adjacency=adjacency,
                band_density_quantiles=band_quantiles, band_composition=band_composition,
                composition_loading=loading, rho=rho,
                marks=RadiusMarks(density_edges=edges, radii_by_bin=radii_by_bin),
                kappa=kappa, target_overlap_fraction=float(np.mean(ratio < kappa)),
                strategy=strategy, type_bandwidths=type_bw, size_log_mu=size_mu,
                size_log_sigma=size_sigma, organization=org,
            )
            model.voids = void_info
            model.estimation = {
                "strategy": strategy,
                "detrend_weight_cap": _DETREND_WEIGHT_CAP,
                "composition_prior_cells": _COMPOSITION_PRIOR_CELLS,
                "bandwidth": bandwidth,
                "bandwidth_source": bandwidth_source,
                "bandwidth_range": [bandwidth_range[0], bandwidth_range[1]],
                "bandwidth_range_source": range_source,
                "bandwidth_at_bound": at_bound is not None,
                "bandwidth_bound": at_bound,
                "per_type_bandwidth": bw_evidence is not None,
                "type_bandwidths": {t: float(h) for t, h in zip(cell_types, type_bw)},
                "type_bandwidth_evidence": (None if bw_evidence is None
                                            else dict(zip(cell_types, bw_evidence))),
                "type_counts": {t: int(c) for t, c in
                                zip(cell_types, np.bincount(type_idx, minlength=n_types))},
                "size_fallbacks": [cell_types[t] for t in size_fallback],
                "learned_band_edges": False,
            }
            if void_info:
                model.estimation["voids"] = {k: v for k, v in void_info.items()
                                             if k != "all_components"}
            if adaptive:
                model.estimation["d_nn"] = float(np.median(
                    cKDTree(np.column_stack([xs, ys])).query(np.column_stack([xs, ys]), k=2)[0][:, 1]))
            if het is not None:
                model.heterogeneity, model.homogeneous = het
            elif n_null > 0:
                model._test_heterogeneity(xs, ys, type_idx, n_null, alpha, rng)
            if model.homogeneous and org["model"] != "none":
                # Uniform layouts ignore any trend; keep the candidates for inspection.
                org = dict(org, model="none", fallback="homogeneous", selected_model=org["model"])
                model.organization = org
            if model.homogeneous and adaptive and at_bound is not None:
                model.flags = model.flags + (f"bandwidth_at_{at_bound}_bound",)
            if not model.homogeneous:
                model._fit_structure(log_density, patch_prior, rng)
                if adaptive and model.n_compartments > 1:
                    edges = _learned_band_edges(_edge_distance_map(region_compartments,
                                                                   grid_step)[iy, ix],
                                                grid_step, model.patch_length)
                    model.band_edges = edges
                    model.band_density_quantiles, model.band_composition = band_stats(edges)
                    model.estimation["learned_band_edges"] = True
                if org["model"] != "none":
                    model.flags = model.flags + (f"organization:{org['model']}",)
                if adaptive and at_bound is not None:
                    model.flags = model.flags + (f"bandwidth_at_{at_bound}_bound",)
                if not model.stationary and not quiet:
                    warnings.warn(
                        f"Density model flags {model.flags}: the region may not be stationary; "
                        "resampled layouts may not resemble it. Consider layout='copy'.",
                        stacklevel=2)
            return model

        def finish(model):
            if adaptive and units == "auto":
                model._detect_units(rng, xs, ys, type_idx, n_null)
                if model.units.get("units"):
                    model.units["profiles"] = _grains.fit_unit_profiles(
                        xs[inside], ys[inside], type_idx[inside], n_types, radii[inside],
                        model.units, width, height)
                    if "d_nn" not in model.units and "d_nn" in model.estimation:
                        model.units["d_nn"] = model.estimation["d_nn"]
            return model

        model = build(maps, dict(org, null=None) if use_org else org)
        if (use_org and n_null > 0 and model.organization["model"] != "none"
                and not model.homogeneous):
            # Decisive trend null: can a stationary patch process (the final
            # model's residual generator, n_null draws) reproduce the trend
            # statistic? Accept iff observed >= BIC_THRESHOLD and p <= alpha;
            # otherwise refit without detrending.
            null = model._trend_null(rng, n_null, alpha, xs, ys, type_idx, org["delta_bic"])
            if null.get("degenerate"):
                model.flags = model.flags + ("trend_null_degenerate",)
            if null["accepted"]:
                model.organization["null"] = null
                return finish(model)
            maps0 = _fit_maps(xs, ys, type_idx, n_types, iy, ix, shape, grid_step, mask,
                              bandwidth, type_bw, bw_evidence, proportions, n_compartments,
                              rng, org=None)
            org0 = dict(org, model="none", fallback="stationary_null",
                        selected_model=org["model"], null=null)
            model = build(maps0, org0, het=(model.heterogeneity, model.homogeneous))
        return finish(model)

    @classmethod
    def from_tissue(cls, tissue, **kwargs) -> 'DensityModel':
        """Fit to a :class:`~tissue_simulator.tissue.TissueSection` (x, y only)."""
        if not tissue.cells:
            raise ValueError("Tissue has no cells.")
        centers = np.array([c.center for c in tissue.cells], dtype=float)
        kwargs.setdefault("bounds", (0.0, 0.0, float(tissue.width), float(tissue.height)))
        return cls.fit(centers[:, 0], centers[:, 1],
                       [c.radius for c in tissue.cells],
                       [c.cell_type for c in tissue.cells], **kwargs)

    def _fit_structure(self, log_density, patch_prior, rng) -> None:
        """Calibrate patch structure and flag signs of non-stationarity."""
        flags: List[str] = []
        upper = min(self.width, self.height) / 3.0
        lower = min(self.grid_step, upper / 4.0)
        if patch_prior:
            self.patch_length = float(patch_prior["length"])
            self.patch_smoothness = float(patch_prior.get("smoothness", 1.0))
        elif self.n_compartments < 2:
            self.patch_length = lower
        else:
            length, nu, at_lower, at_upper = _calibrate_patches(
                self.region_compartments, self.n_compartments, self.compartment_centers,
                self.compartment_fractions, self.rho, self.grid_step, self.bandwidth,
                lower, upper, rng)
            self.patch_length, self.patch_smoothness = length, nu
            if at_upper:
                flags.append("patch_length_at_upper_bound")
            if at_lower:
                flags.append("patch_length_at_lower_bound")

        # Trend: a plane explains much of log density and changes it > 1.5-fold.
        rows, cols = np.nonzero(self.mask)
        design = np.column_stack([np.ones(rows.size), cols, rows]).astype(float)
        coef = np.linalg.lstsq(design, log_density, rcond=None)[0]
        resid = log_density - design @ coef
        total_ss = float(((log_density - log_density.mean()) ** 2).sum())
        span = abs(coef[1]) * np.ptp(cols) + abs(coef[2]) * np.ptp(rows)
        if (total_ss > 0 and span > math.log(1.5)
                and 1.0 - float((resid ** 2).sum()) / total_ss > _TREND_R2_LIMIT):
            flags.append("trend")
        self.flags = tuple(flags)

    # -- diagnostics ----------------------------------------------------------

    @property
    def density(self) -> float:
        """Mean cells per µm² over the tissue mask."""
        return self.n_cells / (self.width * self.height * float(self.mask.mean()))

    @property
    def stationary(self) -> bool:
        """False when a trend or window-sized patches were detected."""
        return not any(f in ("trend", "patch_length_at_upper_bound") for f in self.flags)

    @property
    def n_compartments(self) -> int:
        return int(self.compartment_centers.shape[0])

    def _trend_null(self, rng, n_null, alpha, xs, ys, type_idx, observed) -> dict:
        """Monte Carlo test of the selected trend against stationary layouts.

        Draws ``n_null`` layouts from this (stationary, trend-free) model with
        :meth:`_residual_maps`, samples ``n_cells`` cells from each as an
        inhomogeneous Poisson sample (pixel ~ intensity, uniform jitter inside
        the pixel, type ~ composition at the pixel) and computes
        :func:`_organization.trend_statistic` on it. ``p = (1 + #{null >=
        observed}) / (n_null + 1)``; accepted iff ``observed >=
        BIC_THRESHOLD`` and ``p <= alpha``. ``n_null == 0`` accepts on the
        BIC threshold alone (``p_value`` None).
        """
        import time
        t0 = time.perf_counter()
        out = {"n_null": int(n_null), "p_value": None, "observed_delta_bic": float(observed),
               "null_delta_bic_quantiles": None, "alpha": float(alpha)}
        if n_null <= 0:
            out.update(ambiguous=None, accepted=True, seconds=0.0)
            return out
        shape = self.mask.shape
        n_types = len(self.cell_types)
        gen = self
        stats = []
        for _ in range(n_null):
            intensity, compartment, band, _n = gen._residual_maps(rng, self.width, self.height)
            w = np.where(self.mask, np.maximum(intensity, 0.0), 0.0).ravel()
            if w.sum() <= 0:
                continue
            pix = rng.choice(w.size, size=self.n_cells, p=w / w.sum())
            piy, pix_x = np.unravel_index(pix, shape)
            cx = np.clip((pix_x + rng.random(pix.size)) * self.grid_step, 0, self.width)
            cy = np.clip((piy + rng.random(pix.size)) * self.grid_step, 0, self.height)
            comp = self._composition_map(compartment, band)[piy, pix_x]
            cum = np.cumsum(comp, axis=1)
            cum /= np.maximum(cum[:, -1:], 1e-300)
            t = np.minimum((rng.random(pix.size)[:, None] > cum).sum(axis=1), n_types - 1)
            stats.append(_organization.trend_statistic(cx, cy, t, n_types, self.mask,
                                                       self.grid_step, self.width, self.height))
        stats = np.array(stats, dtype=float)
        stats = stats[np.isfinite(stats)]
        if stats.size == 0:
            # Every null statistic was non-finite: no evidence either way, so
            # behave like n_null == 0 (accept on the BIC threshold alone).
            out.update(n_null=0, ambiguous=None, degenerate=True,
                       accepted=bool(observed >= _organization.BIC_THRESHOLD),
                       seconds=float(time.perf_counter() - t0))
            return out
        p = float((1 + np.sum(stats >= observed)) / (stats.size + 1))
        q = (np.quantile(stats, [0.5, 0.95]).tolist() + [float(stats.max())]) if stats.size \
            else [None, None, None]
        out.update(n_null=int(stats.size), p_value=p,
                   null_delta_bic_quantiles={"p50": q[0], "p95": q[1], "max": q[2]},
                   ambiguous=bool(p > alpha),
                   accepted=bool(observed >= _organization.BIC_THRESHOLD and p <= alpha),
                   seconds=float(time.perf_counter() - t0))
        return out

    def _detect_units(self, rng, xs, ys, type_idx, n_null) -> None:
        """Detect units (:func:`_units.detect_units`) and store them in :attr:`units`.

        The observed per-type maps are plain kernel intensities (inverse-trend weighted when a planar trend was accepted; a radial trend is the unit itself and is not removed) at the
        per-type bandwidths times each of ``_units.UNIT_BANDWIDTH_SCALES`` over the full grid (void pixels stay in the
        domain so lumens form basins), plus the empty-space distance map.
        Null draws (``min(n_null, UNIT_NULL_DRAWS)``) are packed uniform (CSR)
        layouts with the global proportions and the source cell count, voids
        stripped, whether or not the model is homogeneous: the fitted patch
        model already contains the structures under test and is NOT the null
        for units. A second set of draws from the residual patch model is
        advisory only (``units["ambiguous_with_patches"]``). ``rng`` is used
        after every other fit draw. Each draw samples ``n_cells`` cells
        by packing the layout with :class:`~.packing.InhomogeneousPacker` (the replicate
        pipeline, so the hard core matches the source's empty-space statistics,
        which a Poisson sample would not), draws types from the layout
        composition, and maps the cells with the same estimators as the data. ``n_null == 0`` disables unit detection (nothing to calibrate against).
        """
        import dataclasses
        shape = self.mask.shape
        n_types = len(self.cell_types)
        step = self.grid_step
        full = np.ones(shape, dtype=bool)
        bw = np.asarray(self.type_bandwidths, dtype=float)
        if bw.size != n_types:
            bw = np.full(n_types, self.bandwidth)

        def maps_of(cx, cy, t, weights=None):
            iy, ix = _pixel_indices(cx, cy, step, shape)
            return (np.stack([_intensity_grids(iy, ix, t, n_types, shape, step, bw * m, full,
                                               weights)
                              for m in _units.UNIT_BANDWIDTH_SCALES]),
                    _units.empty_space_map(cx, cy, shape, step))

        weights = None
        org = self.organization
        if org.get("model", "none") != "none" and org.get("geometry") == "planar":
            # Units are residual structure: with an accepted trend, weight each
            # cell by the inverse trend (same clipping as _fit_maps) so the
            # observed maps are detrended, like the stationary null layouts.
            dens_c, comp_c = _organization.profile_along(
                org, _organization.coordinate_along(org, xs, ys))
            weights = 1.0 / np.maximum(
                dens_c * comp_c[np.arange(xs.size), type_idx] / self.proportions[type_idx],
                _DETREND_FLOOR)
            weights = np.clip(weights, 1.0 / _DETREND_WEIGHT_CAP, _DETREND_WEIGHT_CAP)
            weights = weights * (weights.size / weights.sum())
        if int(n_null) <= 0:
            self.units = {}
            self.estimation["units"] = {"skipped": "n_null=0"}
            return
        lam, empty = maps_of(xs, ys, type_idx, weights)
        draws, advisory = [], []
        n_draws = min(int(n_null), _units.UNIT_NULL_DRAWS)
        if n_draws > 0:
            from .packing import InhomogeneousPacker
            gen = dataclasses.replace(self, voids={})

            def draw(patches):
                if patches:
                    lay = gen._resampled_layout(rng, self.width, self.height)
                    bins = None
                else:  # CSR: uniform packing, global mix, source cell count
                    lay = gen._uniform_layout(self.width, self.height, self.n_cells)
                    bins = max(self.width, self.height)
                cells = InhomogeneousPacker((self.height, self.width, 0.0), lay, seed=rng,
                                            bin_size=bins).pack()
                if len(cells) < 10:
                    return None
                centers = np.array([c.center[:2] for c in cells])
                cx, cy = centers[:, 0], centers[:, 1]
                piy, pix_x = _pixel_indices(cx, cy, step, shape)
                comp = lay.composition[:, piy, pix_x].T
                cum = np.cumsum(comp, axis=1)
                cum /= np.maximum(cum[:, -1:], 1e-300)
                t = np.minimum((rng.random(len(cells))[:, None] > cum).sum(axis=1), n_types - 1)
                return maps_of(cx, cy, t)

            draws = [d for d in (draw(False) for _ in range(n_draws)) if d is not None]
            if not self.homogeneous:  # advisory patch-model null (never drops units)
                advisory = [d for d in (draw(True) for _ in range(n_draws)) if d is not None]
        self.units = _units.detect_units(lam, self.cell_types, empty, step, self.mask, draws,
                                         d_nn=self.estimation.get("d_nn"),
                                         advisory_draws=advisory)
        if self.units:
            self.estimation["units"] = {k: self.units[k] for k in (
                "n_units", "kinds", "thresholds", "min_center_separation", "radii", "null_draws")}
            self.estimation["units"]["barcode_top"] = self.units.pop("barcode_top")

    def _test_heterogeneity(self, xs, ys, type_idx, n_null, alpha, rng) -> None:
        """Monte Carlo test against uniform packings with permuted labels."""
        from .packing import InhomogeneousPacker

        shape = self.mask.shape
        n_types = len(self.cell_types)
        iy, ix = _pixel_indices(xs, ys, self.grid_step, shape)
        comp_bw = (np.asarray(self.type_bandwidths, dtype=float)
                   if self.strategy == "adaptive" and np.size(self.type_bandwidths) == n_types
                   else None)
        observed = _heterogeneity_statistics(iy, ix, type_idx, n_types, shape,
                                             self.grid_step, self.bandwidth, self.mask,
                                             comp_bw)
        null_layout = self._uniform_layout(self.width, self.height, self.n_cells,
                                           mask=self.mask)
        null = []
        for _ in range(n_null):
            # One bin: plain random sequential addition, like the uniform scaffold.
            cells = InhomogeneousPacker((self.height, self.width, 0.0), null_layout,
                                        seed=rng, bin_size=max(self.width, self.height)).pack()
            if not cells:
                continue
            centers = np.array([c.center for c in cells])
            labels = rng.permutation(type_idx)
            if labels.size != len(cells):
                labels = rng.choice(type_idx, size=len(cells), replace=True)
            niy, nix = _pixel_indices(centers[:, 0], centers[:, 1], self.grid_step, shape)
            null.append(_heterogeneity_statistics(niy, nix, labels, n_types, shape,
                                                  self.grid_step, self.bandwidth, self.mask,
                                                  comp_bw))
        null = np.array(null).reshape(-1, 2)
        p_density = (1 + np.sum(null[:, 0] >= observed[0])) / (null.shape[0] + 1)
        p_composition = (1 + np.sum(null[:, 1] >= observed[1])) / (null.shape[0] + 1)
        self.heterogeneity = {
            "density_statistic": observed[0], "composition_statistic": observed[1],
            "p_density": float(p_density), "p_composition": float(p_composition),
            "n_null": int(null.shape[0]),
        }
        if comp_bw is not None:
            self.heterogeneity["composition_bandwidths"] = {
                t: float(h) for t, h in zip(self.cell_types, comp_bw)}
        self.homogeneous = bool(p_density > alpha and p_composition > alpha)

    # -- layouts ------------------------------------------------------------

    def sample_layout(self, rng=None, width: Optional[float] = None,
                      height: Optional[float] = None,
                      layout: str = "resample", *, max_proposals: int = 20) -> Layout:
        """Sample the target maps for one replicate.

        Args:
            rng: Generator, SeedSequence, int seed or None.
            width, height: Replicate window in µm; defaults to the region's.
                ``"copy"`` requires the region's window size.
            layout: ``"resample"`` (new compartment arrangement) or ``"copy"``.
            max_proposals: Direction proposals for organized (adaptive,
                trend-model) resampling before the best-coverage fallback.

        Homogeneous regions (see :attr:`homogeneous`) always yield a uniform
        layout, which reproduces the classic uniform packing.
        """
        if layout not in ("resample", "copy"):
            raise ValueError(f"layout must be 'resample' or 'copy', got {layout!r}.")
        rng = _as_generator(rng)
        width = self.width if width is None else float(width)
        height = self.height if height is None else float(height)
        if self._units_active() and (self.homogeneous or layout == "resample"):
            return self._units_layout(rng, width, height, layout, max_proposals)

        if self.homogeneous:
            if self._has_voids():
                shape = _grid_shape(width, height, self.grid_step)
                tissue, placed = _voids.place_voids(rng, shape, self.grid_step, self.voids,
                                                    np.ones(shape, dtype=bool))
                n = int(round(self.density * width * height * float(tissue.mean())))
                return self._with_voids(self._uniform_layout(width, height, n, mask=tissue),
                                        placed)
            return self._uniform_layout(width, height, int(round(self.density * width * height)))
        if layout == "copy":
            if _grid_shape(width, height, self.grid_step) != self.mask.shape:
                raise ValueError("layout='copy' requires the region's window size.")
            return self._copy_layout(width, height)
        # Radial geometry is not organized yet (falls back to plain resampling).
        if (self.strategy == "adaptive" and self.organization.get("model", "none") != "none"
                and self.organization.get("geometry") in ("planar", "radial")):
            return self._organized_layout(rng, width, height, max_proposals)
        return self._resampled_layout(rng, width, height)

    # -- units (germ-grain) -------------------------------------------------

    def _units_active(self) -> bool:
        return (self.strategy == "adaptive" and bool(self.units.get("units"))
                and bool(self.units.get("profiles")))

    def _unit_free_voids(self) -> Dict:
        """Copy of :attr:`voids` without holes owned by a unit (precedence rule
        :data:`UNIT_OWNS_LUMEN`): interior holes whose centroid lies inside a
        unit's outer radius are dropped from the placement pool. Diagnostics
        (``all_components``) are kept; ``self.voids`` is not modified."""
        if not self._has_voids():
            return self.voids
        us = self.units["units"]

        def owned(c):
            return any(math.hypot(c["centroid"][0] - u["center"][0],
                                  c["centroid"][1] - u["center"][1]) <= u["outer_radius"]
                       for u in us)

        comps = self.voids.get("all_components", [])
        keep = [c for c in comps if not c["edge_touching"] and not owned(c)]
        out = dict(self.voids)
        out["n_holes"] = len(keep)
        if keep:
            out["equivalent_diameters"] = [c["equivalent_diameter"] for c in keep]
        out["dropped_for_units"] = sum(1 for c in comps if not c["edge_touching"]) - len(keep)
        return out

    def _radial_center_in_unit(self) -> bool:
        """Precedence rule :data:`UNIT_OWNS_RADIAL`: a radial trend centered inside
        a unit is replaced by the unit."""
        org = self.organization
        c = org.get("center")
        if org.get("model", "none") == "none" or org.get("geometry") != "radial" or c is None:
            return False
        return any(math.hypot(c[0] - u["center"][0], c[1] - u["center"][1]) <= u["outer_radius"]
                   for u in self.units["units"])

    def _units_layout(self, rng, width, height, layout, max_proposals) -> Layout:
        """Residual layout (all existing draws) with germ-grain units on top.

        Grains are placed after every other draw (noise fields; theta/center for
        organized layouts), rasterized with the fitted profiles, and multiply the
        residual intensity (lumen pixels -> 0, compartment -1); inside grains the
        composition is the grain composition. ``n_target`` follows the tissue
        fraction (voids and lumens excluded) and the intensity is renormalized to it.
        """
        import dataclasses
        skip = self._radial_center_in_unit()
        view = dataclasses.replace(
            self, units={}, voids=self._unit_free_voids(),
            organization=({"model": "none"} if skip else self.organization))
        out = view.sample_layout(rng, width, height, layout, max_proposals=max_proposals)
        if skip:
            out.organization = {"model": "none", "fallback": "units"}
        shape = out.intensity.shape
        step = self.grid_step
        placed = _grains.place_grains(rng, shape, step, self.units,
                                      self.units.get("min_center_separation"),
                                      max_tries=UNIT_PLACEMENT_TRIES)
        factor, comp, gmask, vmask = _grains.rasterize_grains(
            placed, self.units["profiles"], shape, step, self.proportions)
        intensity = np.where(vmask, 0.0, out.intensity * factor)
        composition = np.where(gmask[None], comp, out.composition)
        composition = np.maximum(composition, 0.0)
        composition = composition / np.maximum(composition.sum(axis=0, keepdims=True), 1e-300)
        out.compartment = np.where(vmask, -1, out.compartment)
        tissue = (out.compartment >= 0)
        n_target = int(round(self.density * width * height * float(tissue.mean())))
        out.intensity = self._normalize(intensity, n_target)
        out.composition = composition
        out.n_target = n_target
        out.units = placed
        if placed["shortfall"] > 0:
            out.flags = tuple(out.flags) + ("units_unsatisfied",)
        return out

    def _layout(self, width, height, intensity, composition, compartment, n_target,
                mode, marks=None) -> Layout:
        return Layout(width=width, height=height, grid_step=self.grid_step,
                      intensity=intensity, composition=composition,
                      compartment=compartment, cell_types=self.cell_types,
                      marks=self.marks if marks is None else marks, kappa=self.kappa,
                      n_target=int(n_target),
                      target_overlap_fraction=self.target_overlap_fraction,
                      mode=mode, flags=self.flags, bandwidth=self.bandwidth,
                      strategy=self.strategy, quota_scale=self._quota_scale())

    def _has_voids(self) -> bool:
        return bool(self.voids) and int(self.voids.get("n_holes", 0)) > 0

    def _with_voids(self, layout: Layout, placed: dict) -> Layout:
        layout.voids = placed
        if placed.get("n_placed", 0) < placed.get("n_requested", 0):
            layout.flags = tuple(layout.flags) + ("voids_unsatisfied",)
        return layout

    def _composition_map(self, compartment, band) -> np.ndarray:
        """Per-pixel composition ``(ny, nx, n_types)``; void pixels (compartment
        -1) take the global proportions instead of wrapping to the last compartment."""
        if not (compartment < 0).any():
            return self.band_composition[compartment, band]
        c, b = np.maximum(compartment, 0), np.maximum(band, 0)
        return np.where((compartment < 0)[..., None], self.proportions,
                        self.band_composition[c, b])

    def _radial_anchor_hole(self) -> bool:
        """True when a source hole is centered on the fitted radial center
        (centroid within one hole radius): the largest placed hole is then
        anchored at each proposed center (a lumen at the center of rings)."""
        c = self.organization.get("center")
        if c is None:
            return False
        for comp in self.voids.get("all_components", []):
            if comp["edge_touching"]:
                continue
            r = comp["equivalent_diameter"] / 2.0
            if math.hypot(comp["centroid"][0] - c[0], comp["centroid"][1] - c[1]) <= r:
                return True
        return False

    def _quota_scale(self) -> Optional[float]:
        if self.strategy != "adaptive":
            return None
        d_nn = self.estimation.get("d_nn")
        if d_nn is None or not self.bandwidth or self.density <= 0:
            return None
        floor = max(2.0 * d_nn, 2.0 * self.grid_step, math.sqrt(4.0 / self.density))
        size = max(self.bandwidth / 2.0, floor)
        return float(max(round(size / self.grid_step), 1) * self.grid_step)

    def _uniform_layout(self, width, height, n_target, mask=None) -> Layout:
        shape = _grid_shape(width, height, self.grid_step)
        tissue = np.ones(shape, dtype=bool) if mask is None else mask
        area = width * height * float(tissue.mean())
        intensity = np.where(tissue, n_target / area, 0.0)
        composition = np.broadcast_to(self.proportions[:, None, None],
                                      (len(self.cell_types),) + shape).copy()
        all_radii = RadiusMarks(density_edges=np.zeros(0),
                                radii_by_bin=[np.concatenate(self.marks.radii_by_bin)])
        return self._layout(width, height, intensity, composition,
                            np.where(tissue, 0, -1), n_target, "uniform", marks=all_radii)

    def _copy_layout(self, width, height) -> Layout:
        total = self.region_intensity.sum(axis=0)
        composition = np.where(total > 0, self.region_intensity / np.maximum(total, 1e-300),
                               self.proportions[:, None, None])
        out = self._layout(width, height, total.copy(), composition,
                           self.region_compartments.copy(), self.n_cells, "copy")
        if self._has_voids():
            comps = self.voids.get("all_components", [])
            interior = [c for c in comps if not c["edge_touching"]]
            out.voids = {"centers": [list(c["centroid"]) for c in interior],
                         "radii": [c["equivalent_diameter"] / 2.0 for c in interior],
                         "n_requested": len(interior), "n_placed": len(interior),
                         "anchored": False}
        return out

    def _noise_fields(self, rng, shape):
        """The two latent-field draws of a resampled layout (labels, scores)."""
        pad = _padding(shape, self.grid_step, self.patch_length, self.bandwidth)
        return _sample_labels(
            shape, self.grid_step, self.patch_length, self.patch_smoothness, self.bandwidth,
            self.rho, self.compartment_centers, self.compartment_fractions,
            rng.standard_normal(_noise_shape(shape, pad)),
            rng.standard_normal(_noise_shape(shape, pad)), pad)

    def _maps_from_fields(self, labels, u_scores, shape, hole_mask=None):
        """Compartment, band and intensity maps from latent fields.

        ``hole_mask`` (True inside placed voids) sets compartment -1 and zero
        intensity there and is recomputed into the band index, so bands see the
        void edge. Intensity is returned un-normalized (quantile units).
        """
        labels = labels.copy()
        if hole_mask is not None:
            labels[hole_mask.ravel()] = -1
        compartment = labels.reshape(shape)
        band = _band_index(compartment, self.grid_step, self.band_edges)
        flat_band = band.ravel()
        intensity = np.zeros(labels.size) if hole_mask is not None else np.empty(labels.size)
        for c in range(self.n_compartments):
            for b in range(self.band_density_quantiles.shape[1]):
                sel = np.flatnonzero((labels == c) & (flat_band == b))
                if not sel.size:
                    continue
                ranks = np.argsort(np.argsort(u_scores[sel], kind='stable'), kind='stable')
                intensity[sel] = np.interp((ranks + 0.5) / sel.size, _QUANTILE_LEVELS,
                                           self.band_density_quantiles[c, b])
        return intensity.reshape(shape), compartment, band

    def _normalize(self, intensity, n_target):
        mass = intensity.sum() * self.grid_step ** 2
        if mass > 0:
            intensity = intensity * n_target / mass
        return intensity

    def _residual_maps(self, rng, width, height, theta0=None):
        """Compartment/band maps (residual structure) from two noise fields.

        Returns ``(intensity normalized to n_target, compartment, band, n_target)``.
        With voids, ``n_holes`` discs are placed after the two noise fields
        (:func:`_voids.place_voids`) and ``n_target`` counts tissue pixels only;
        the placement is left in ``self._last_placed`` for the caller.
        """
        shape = _grid_shape(width, height, self.grid_step)
        n_target = int(round(self.density * width * height))
        labels, u_scores = self._noise_fields(rng, shape)
        hole_mask = None
        self._last_placed = {}
        if self._has_voids():
            tissue, placed = _voids.place_voids(rng, shape, self.grid_step, self.voids,
                                                np.ones(shape, dtype=bool))
            hole_mask = ~tissue
            self._last_placed = placed
            n_target = int(round(self.density * width * height * float(tissue.mean())))
        intensity, compartment, band = self._maps_from_fields(labels, u_scores, shape,
                                                              hole_mask)
        intensity = intensity.copy() if hole_mask is not None else intensity
        mass = intensity.sum() * self.grid_step ** 2
        if mass > 0:
            intensity *= n_target / mass
        return intensity, compartment, band, n_target

    def _resampled_layout(self, rng, width, height) -> Layout:
        intensity, compartment, band, n_target = self._residual_maps(rng, width, height)
        placed = self._last_placed
        composition = np.moveaxis(self._composition_map(compartment, band), -1, 0).copy()
        out = self._layout(width, height, intensity, composition, compartment,
                           n_target, "resample")
        return self._with_voids(out, placed) if placed else out

    def _organized_layout(self, rng, width, height, max_proposals) -> Layout:
        """Resampled layout with the source's fitted trend in a new direction.

        ``theta'`` is drawn first, then the two noise fields as in
        :meth:`_resampled_layout`; the trend (see :mod:`._organization`) is
        multiplied onto the residual maps. Rejected directions are re-drawn
        (noise reused) up to ``max_proposals``; on exhaustion the best-coverage
        proposal is used, flagged ``organization_unsatisfied`` and warned.
        """
        org = self.organization
        radial = org.get("geometry") == "radial"
        if radial:
            shape0 = _grid_shape(width, height, self.grid_step)
            theta0 = _organization.sample_center(rng, shape0, self.grid_step)
        else:
            theta0 = _organization.sample_direction(rng)
        if not self._has_voids():
            intensity, compartment, band, n_target = self._residual_maps(rng, width, height)
            resid_comp = np.moveaxis(self._composition_map(compartment, band), -1, 0)
            intensity, composition, info = _organization.propose_layout(
                org, theta0, rng, intensity.shape, self.grid_step, intensity, resid_comp,
                self.proportions, n_target, max_proposals)
            placed = {}
        else:
            # Holes are drawn after theta0 and the two noise fields, per
            # proposal when anchored on the proposed radial center.
            shape = _grid_shape(width, height, self.grid_step)
            labels, u_scores = self._noise_fields(rng, shape)
            anchored = radial and self._radial_anchor_hole()
            states = {}

            def resid_fn(theta):
                key = tuple(theta) if radial else float(theta)
                if not anchored and states:
                    state = next(iter(states.values()))
                    states[key] = state
                    return state[0], state[1], state[2]
                tissue, pl = _voids.place_voids(
                    rng, shape, self.grid_step, self.voids, np.ones(shape, dtype=bool),
                    anchor=(theta if anchored else None))
                n_t = int(round(self.density * width * height * float(tissue.mean())))
                inten, comp, band = self._maps_from_fields(labels, u_scores, shape, ~tissue)
                inten = self._normalize(inten, n_t)
                rc = np.moveaxis(self._composition_map(comp, band), -1, 0)
                states[key] = (inten, rc, n_t, comp, pl)
                return inten, rc, n_t

            intensity, composition, info = _organization.propose_layout(
                org, theta0, rng, shape, self.grid_step, None, None,
                self.proportions, None, max_proposals, resid_fn=resid_fn)
            th = info["theta"]
            _, _, n_target, compartment, placed = states[tuple(th) if radial else float(th)]
        out = self._layout(width, height, intensity, composition, compartment,
                           n_target, "resample")
        if placed:
            self._with_voids(out, placed)
        theta = info["theta"]
        geom = ({"center": [float(theta[0]), float(theta[1])]} if radial else
                {"direction": [math.cos(theta), math.sin(theta)], "theta": theta})
        out.organization = {
            "model": org["model"], "geometry": org["geometry"], **geom,
            "proposals_tried": info["proposals_tried"], "accepted": info["accepted"],
            "fallback": info["fallback"],
        }
        if not info["accepted"]:
            out.flags = tuple(out.flags) + ("organization_unsatisfied",)
            warnings.warn(
                f"No direction/center in {info['proposals_tried']} proposals met the window coverage "
                "and composition criteria; using the best-coverage proposal.", stacklevel=3)
        return out

    # -- serialization ------------------------------------------------------

    def to_dict(self) -> Dict:
        """JSON-serializable representation (see :meth:`from_dict`)."""
        out = {"format_version": 2}
        for f in fields(self):
            value = getattr(self, f.name)
            if isinstance(value, np.ndarray):
                value = {"__ndarray__": value.tolist(), "dtype": str(value.dtype)}
            elif isinstance(value, RadiusMarks):
                value = {"density_edges": value.density_edges.tolist(),
                         "radii_by_bin": [r.tolist() for r in value.radii_by_bin]}
            elif isinstance(value, tuple):
                value = list(value)
            out[f.name] = value
        return out

    @classmethod
    def from_dict(cls, data: Dict) -> 'DensityModel':
        """Rebuild a model written by :meth:`to_dict`."""
        kwargs = {}
        for f in fields(cls):
            if f.name not in data:  # format_version 1: keep the legacy default
                continue
            value = data[f.name]
            if isinstance(value, dict) and "__ndarray__" in value:
                value = np.array(value["__ndarray__"], dtype=value["dtype"])
            elif f.name == "marks":
                value = RadiusMarks(np.array(value["density_edges"], dtype=float),
                                    [np.array(r, dtype=float) for r in value["radii_by_bin"]])
            elif f.name in ("cell_types", "flags", "band_edges"):
                value = tuple(value)
            kwargs[f.name] = value
        return cls(**kwargs)


def pooled_patch_prior(models: Sequence[DensityModel]) -> Dict[str, float]:
    """Pooled patch structure over models whose calibration hit no bound.

    Returns ``{"length": geometric mean, "smoothness": median}``, or an empty
    dict when no model qualifies. Useful as ``patch_prior`` for small or
    sparse regions from the same cohort.
    """
    usable = [m for m in models
              if not m.homogeneous and not any(f.startswith("patch_length") for f in m.flags)]
    if not usable:
        return {}
    return {"length": float(np.exp(np.mean(np.log([m.patch_length for m in usable])))),
            "smoothness": float(np.median([m.patch_smoothness for m in usable]))}
