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

def _intensity_grids(iy, ix, type_idx, n_types, shape, step, bandwidth, mask):
    """Edge-corrected Gaussian kernel intensity per type, in cells per µm²."""
    sigma = bandwidth / step
    edge = np.where(mask, np.maximum(_smooth(mask.astype(float), sigma), 1e-6), 1.0)
    grids = np.zeros((n_types,) + tuple(shape))
    for t in range(n_types):
        counts = np.zeros(shape)
        sel = type_idx == t
        np.add.at(counts, (iy[sel], ix[sel]), 1.0)
        grids[t] = np.where(mask, _smooth(counts, sigma) / edge, 0.0) / step ** 2
    return grids


def _select_bandwidth(iy, ix, shape, step, mask, bandwidth_range, n_candidates=10):
    """Likelihood cross-validation bandwidth for the total intensity."""
    lo, hi = bandwidth_range
    candidates = np.geomspace(lo, hi, n_candidates)
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
    return float(candidates[int(np.argmax(scores))])


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


def _band_index(labels: np.ndarray, step: float) -> np.ndarray:
    """Distance band of each pixel from the edge of its own compartment (-1: no tissue).

    Band ``b`` holds pixels whose distance to the nearest pixel of another
    compartment (or of no tissue) is at most ``_BAND_EDGES_UM[b]``; the last
    band is the interior. Window edges do not count as compartment edges.
    """
    band = np.full(labels.shape, -1, dtype=int)
    for c in np.unique(labels[labels >= 0]):
        inside = labels == c
        distance = ndimage.distance_transform_edt(inside) * step
        band[inside] = np.searchsorted(_BAND_EDGES_UM, distance[inside], side='left')
    return band


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


def _heterogeneity_statistics(iy, ix, type_idx, n_types, shape, step, bandwidth, mask):
    """Variance of log density and density-weighted composition departure."""
    lam = _intensity_grids(iy, ix, type_idx, n_types, shape, step, bandwidth, mask)
    total = lam.sum(axis=0)[mask]
    mean = total.mean()
    t_density = float(np.var(np.log(total + 0.01 * mean))) if mean > 0 else 0.0
    if n_types < 2 or total.sum() <= 0:
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

    # -- construction -------------------------------------------------------

    @classmethod
    def fit(cls, x: Sequence[float], y: Sequence[float], radius: Sequence[float],
            cell_type: Sequence[str],
            bounds: Optional[Tuple[float, float, float, float]] = None,
            mask: Optional[np.ndarray] = None,
            bandwidth: Optional[float] = None,
            bandwidth_range: Tuple[float, float] = (15.0, 80.0),
            grid_step: float = 5.0,
            n_compartments: int = 3,
            patch_prior: Optional[Dict[str, float]] = None,
            n_null: int = 19,
            alpha: float = 0.05,
            seed=None) -> 'DensityModel':
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
        """
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

        cell_types = tuple(sorted(set(labels.tolist())))
        type_idx = np.searchsorted(cell_types, labels)
        n_types = len(cell_types)
        rng = _as_generator(seed)

        if bandwidth is None:
            bandwidth = _select_bandwidth(iy, ix, shape, grid_step, mask, bandwidth_range)
        bandwidth = float(bandwidth)
        lam_types = _intensity_grids(iy, ix, type_idx, n_types, shape, grid_step,
                                     bandwidth, mask)
        lam_total = lam_types.sum(axis=0)
        proportions = np.bincount(type_idx, minlength=n_types) / type_idx.size

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
                factor = (in_c.sum() / (sel.sum() * grid_step ** 2)) / smoothed.mean()
                composition[c] = np.bincount(type_idx[in_c], minlength=n_types) / in_c.sum()
            else:
                factor = 0.0
                composition[c] = proportions
            lam_types[:, sel] *= factor
            density_quantiles[c] = np.quantile(smoothed * factor, _QUANTILE_LEVELS)
        lam_total = lam_types.sum(axis=0)

        # Boundary profiles: density and composition by distance to the
        # compartment edge, again counted from cells (falls back to the
        # compartment's values where a band holds too few cells).
        band = _band_index(region_compartments, grid_step)
        cell_band = band[iy, ix]
        n_bands = len(_BAND_EDGES_UM) + 1
        band_quantiles = np.repeat(density_quantiles[:, None, :], n_bands, axis=1)
        band_composition = np.repeat(composition[:, None, :], n_bands, axis=1)
        for c in range(k):
            for b in range(n_bands):
                sel = (region_compartments == c) & (band == b)
                in_cb = (cell_compartment == c) & (cell_band == b)
                smoothed = lam_total[sel]
                if in_cb.sum() < _MIN_BAND_CELLS or smoothed.mean() <= 0:
                    continue
                band_density = in_cb.sum() / (sel.sum() * grid_step ** 2)
                band_quantiles[c, b] = np.quantile(smoothed * (band_density / smoothed.mean()),
                                                   _QUANTILE_LEVELS)
                band_composition[c, b] = np.bincount(type_idx[in_cb], minlength=n_types) / in_cb.sum()

        adjacency = np.zeros((k, k))
        for a, b in ((region_compartments[:, :-1], region_compartments[:, 1:]),
                     (region_compartments[:-1, :], region_compartments[1:, :])):
            both = (a >= 0) & (b >= 0)
            np.add.at(adjacency, (a[both], b[both]), 1.0)
        adjacency = adjacency + adjacency.T
        adjacency /= max(adjacency.sum(), 1.0)

        # Radius marks and hard core from the observed cells.
        lam_at_cells = lam_total[iy, ix]
        edges = np.quantile(lam_at_cells, np.linspace(0, 1, _N_RADIUS_BINS + 1)[1:-1])
        bins = np.searchsorted(edges, lam_at_cells, side='right')
        radii_by_bin = [radii[bins == b] if np.any(bins == b) else radii.copy()
                        for b in range(_N_RADIUS_BINS)]
        points = np.column_stack([xs, ys])
        distances, neighbors = cKDTree(points).query(points, k=2)
        ratio = distances[:, 1] / np.maximum(radii + radii[neighbors[:, 1]], 1e-12)
        kappa = float(np.clip(np.quantile(ratio, 0.02), 0.5, 1.5))

        model = cls(
            width=width, height=height, grid_step=float(grid_step), bandwidth=bandwidth,
            cell_types=cell_types, n_cells=int(xs.size), proportions=proportions,
            mask=mask, region_intensity=lam_types, region_compartments=region_compartments,
            compartment_centers=centers, compartment_fractions=fractions,
            compartment_density_quantiles=density_quantiles,
            compartment_composition=composition, compartment_adjacency=adjacency,
            band_density_quantiles=band_quantiles, band_composition=band_composition,
            composition_loading=loading, rho=rho,
            marks=RadiusMarks(density_edges=edges, radii_by_bin=radii_by_bin),
            kappa=kappa, target_overlap_fraction=float(np.mean(ratio < kappa)),
        )
        if n_null > 0:
            model._test_heterogeneity(xs, ys, type_idx, n_null, alpha, rng)
        if not model.homogeneous:
            model._fit_structure(log_density, patch_prior, rng)
            if not model.stationary:
                warnings.warn(
                    f"Density model flags {model.flags}: the region may not be stationary; "
                    "resampled layouts may not resemble it. Consider layout='copy'.",
                    stacklevel=2)
        return model

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

    def _test_heterogeneity(self, xs, ys, type_idx, n_null, alpha, rng) -> None:
        """Monte Carlo test against uniform packings with permuted labels."""
        from .packing import InhomogeneousPacker

        shape = self.mask.shape
        n_types = len(self.cell_types)
        iy, ix = _pixel_indices(xs, ys, self.grid_step, shape)
        observed = _heterogeneity_statistics(iy, ix, type_idx, n_types, shape,
                                             self.grid_step, self.bandwidth, self.mask)
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
                                                  self.grid_step, self.bandwidth, self.mask))
        null = np.array(null).reshape(-1, 2)
        p_density = (1 + np.sum(null[:, 0] >= observed[0])) / (null.shape[0] + 1)
        p_composition = (1 + np.sum(null[:, 1] >= observed[1])) / (null.shape[0] + 1)
        self.heterogeneity = {
            "density_statistic": observed[0], "composition_statistic": observed[1],
            "p_density": float(p_density), "p_composition": float(p_composition),
            "n_null": int(null.shape[0]),
        }
        self.homogeneous = bool(p_density > alpha and p_composition > alpha)

    # -- layouts ------------------------------------------------------------

    def sample_layout(self, rng=None, width: Optional[float] = None,
                      height: Optional[float] = None,
                      layout: str = "resample") -> Layout:
        """Sample the target maps for one replicate.

        Args:
            rng: Generator, SeedSequence, int seed or None.
            width, height: Replicate window in µm; defaults to the region's.
                ``"copy"`` requires the region's window size.
            layout: ``"resample"`` (new compartment arrangement) or ``"copy"``.

        Homogeneous regions (see :attr:`homogeneous`) always yield a uniform
        layout, which reproduces the classic uniform packing.
        """
        if layout not in ("resample", "copy"):
            raise ValueError(f"layout must be 'resample' or 'copy', got {layout!r}.")
        rng = _as_generator(rng)
        width = self.width if width is None else float(width)
        height = self.height if height is None else float(height)

        if self.homogeneous:
            return self._uniform_layout(width, height, int(round(self.density * width * height)))
        if layout == "copy":
            if _grid_shape(width, height, self.grid_step) != self.mask.shape:
                raise ValueError("layout='copy' requires the region's window size.")
            return self._copy_layout(width, height)
        return self._resampled_layout(rng, width, height)

    def _layout(self, width, height, intensity, composition, compartment, n_target,
                mode, marks=None) -> Layout:
        return Layout(width=width, height=height, grid_step=self.grid_step,
                      intensity=intensity, composition=composition,
                      compartment=compartment, cell_types=self.cell_types,
                      marks=self.marks if marks is None else marks, kappa=self.kappa,
                      n_target=int(n_target),
                      target_overlap_fraction=self.target_overlap_fraction,
                      mode=mode, flags=self.flags, bandwidth=self.bandwidth)

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
        return self._layout(width, height, total.copy(), composition,
                            self.region_compartments.copy(), self.n_cells, "copy")

    def _resampled_layout(self, rng, width, height) -> Layout:
        shape = _grid_shape(width, height, self.grid_step)
        n_target = int(round(self.density * width * height))
        pad = _padding(shape, self.grid_step, self.patch_length, self.bandwidth)
        labels, u_scores = _sample_labels(
            shape, self.grid_step, self.patch_length, self.patch_smoothness, self.bandwidth,
            self.rho, self.compartment_centers, self.compartment_fractions,
            rng.standard_normal(_noise_shape(shape, pad)),
            rng.standard_normal(_noise_shape(shape, pad)), pad)

        compartment = labels.reshape(shape)
        band = _band_index(compartment, self.grid_step)
        flat_band = band.ravel()
        intensity = np.empty(labels.size)
        for c in range(self.n_compartments):
            for b in range(self.band_density_quantiles.shape[1]):
                sel = np.flatnonzero((labels == c) & (flat_band == b))
                if not sel.size:
                    continue
                ranks = np.argsort(np.argsort(u_scores[sel], kind='stable'), kind='stable')
                intensity[sel] = np.interp((ranks + 0.5) / sel.size, _QUANTILE_LEVELS,
                                           self.band_density_quantiles[c, b])
        intensity = intensity.reshape(shape)
        mass = intensity.sum() * self.grid_step ** 2
        if mass > 0:
            intensity *= n_target / mass
        composition = np.moveaxis(self.band_composition[compartment, band], -1, 0).copy()
        return self._layout(width, height, intensity, composition, compartment,
                            n_target, "resample")

    # -- serialization ------------------------------------------------------

    def to_dict(self) -> Dict:
        """JSON-serializable representation (see :meth:`from_dict`)."""
        out = {"format_version": 1}
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
            value = data[f.name]
            if isinstance(value, dict) and "__ndarray__" in value:
                value = np.array(value["__ndarray__"], dtype=value["dtype"])
            elif f.name == "marks":
                value = RadiusMarks(np.array(value["density_edges"], dtype=float),
                                    [np.array(r, dtype=float) for r in value["radii_by_bin"]])
            elif f.name in ("cell_types", "flags"):
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
