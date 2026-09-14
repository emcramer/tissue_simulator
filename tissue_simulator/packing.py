"""
Sphere packing algorithm for cell placement.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, fields
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

from .density import Layout
from .tissue import Cell


def _stochastic_round(weights, total: int, rng: np.random.Generator) -> np.ndarray:
    """Integer allocation of ``total`` proportional to ``weights``.

    Floors are kept and the remaining units go to entries drawn without
    replacement in proportion to their fractional parts, so ties never
    favor low indices (which would bias quotas toward one corner).
    """
    w = np.clip(np.asarray(weights, dtype=float), 0.0, None)
    if total <= 0 or w.sum() <= 0:
        return np.zeros(w.size, dtype=int)
    raw = w / w.sum() * total
    counts = np.floor(raw).astype(int)
    remainder = int(total - counts.sum())
    if remainder > 0:
        frac = raw - counts
        chosen = rng.choice(w.size, size=remainder, replace=False, p=frac / frac.sum())
        counts[chosen] += 1
    return counts


class SpatialHashGrid:
    """
    Uniform hash grid over cell centers for short-range neighbor queries.

    Space is divided into cubic buckets of side ``cell_size``. A query returns
    every index stored within ``reach`` buckets of the query point, which is a
    superset of all stored points closer than ``reach * cell_size``.
    """

    def __init__(self, cell_size: float):
        if not cell_size > 0:
            raise ValueError(f"cell_size must be positive, got {cell_size!r}")
        self.cell_size = float(cell_size)
        self._buckets: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)

    def key(self, point) -> Tuple[int, int, int]:
        """Bucket key of a 3D point."""
        s = self.cell_size
        return (math.floor(point[0] / s), math.floor(point[1] / s),
                math.floor(point[2] / s))

    def insert(self, index: int, point) -> None:
        """Store ``index`` at ``point``."""
        self._buckets[self.key(point)].append(index)

    def remove(self, index: int, point) -> None:
        """Remove ``index`` previously inserted at ``point``."""
        self._buckets[self.key(point)].remove(index)

    def move(self, index: int, old_point, new_point) -> None:
        """Re-bucket ``index`` after its point moved."""
        old_key, new_key = self.key(old_point), self.key(new_point)
        if old_key != new_key:
            self._buckets[old_key].remove(index)
            self._buckets[new_key].append(index)

    def neighbors(self, point, reach: int = 1) -> Iterator[int]:
        """Yield indices stored within ``reach`` buckets of ``point``."""
        kx, ky, kz = self.key(point)
        buckets = self._buckets
        for i in range(kx - reach, kx + reach + 1):
            for j in range(ky - reach, ky + reach + 1):
                for k in range(kz - reach, kz + reach + 1):
                    bucket = buckets.get((i, j, k))
                    if bucket:
                        yield from bucket


class SpherePacker:
    """
    Random sphere packing algorithm for placing cells in tissue.
    """

    def __init__(self, bounds: Tuple[float, float, float],
                 cell_radii_config: Dict[str, Tuple[float, float]],
                 min_spacing: float = 0.5,
                 allow_boundary_cells: bool = True,
                 seed: Optional[int] = None):
        """
        Initialize sphere packer.

        Args:
            bounds: (height, width, thickness) of tissue
            cell_radii_config: Dict mapping cell types to (min, max) radii
            min_spacing: Minimum spacing between cell surfaces
            allow_boundary_cells: Allow cells extending beyond bounds
            seed: Optional integer seed for the instance RNG. When provided,
                the packing process is deterministic; when None the RNG is
                seeded from system entropy (previous default behavior).
        """
        self.bounds = bounds
        self.cell_radii_config = cell_radii_config
        self.min_spacing = min_spacing
        self.allow_boundary_cells = allow_boundary_cells
        self.seed = seed
        self._rng = np.random.default_rng(seed)

    def _select_cell_type_and_radius(self) -> Tuple[str, float]:
        """Randomly select a cell type and radius."""
        cell_type = self._rng.choice(list(self.cell_radii_config.keys()))
        min_radius, max_radius = self.cell_radii_config[cell_type]
        radius = self._rng.uniform(min_radius, max_radius)
        return cell_type, radius

    def _generate_random_position(self, radius: float) -> np.ndarray:
        """
        Generate a random position for a cell.

        Args:
            radius: Cell radius

        Returns:
            3D position array [x, y, z]
        """
        height, width, thickness = self.bounds

        if self.allow_boundary_cells:
            # Allow centers anywhere in the tissue (cells can extend beyond)
            x = self._rng.uniform(0, width)
            y = self._rng.uniform(0, height)
            z = self._rng.uniform(0, thickness)
        else:
            # Constrain centers so cells stay within bounds
            x = self._rng.uniform(radius, width - radius)
            y = self._rng.uniform(radius, height - radius)
            z = self._rng.uniform(radius, thickness - radius)

        return np.array([x, y, z])

    def _check_collision(self, candidate: Cell, existing_cells: List[Cell]) -> bool:
        """
        Check if candidate cell collides with existing cells.

        Args:
            candidate: Cell to test
            existing_cells: List of already placed cells

        Returns:
            True if collision detected, False otherwise
        """
        for cell in existing_cells:
            distance = np.linalg.norm(candidate.center - cell.center)
            min_distance = candidate.radius + cell.radius + self.min_spacing

            if distance < min_distance:
                return True

        return False

    def _neighbor_grid(self) -> SpatialHashGrid:
        """Hash grid whose buckets span the largest possible collision distance."""
        max_radius = max(r_max for _, r_max in self.cell_radii_config.values())
        reach = 2 * max_radius + self.min_spacing
        return SpatialHashGrid(reach if reach > 0 else max(2 * max_radius, 1.0))

    def _check_collision_grid(self, candidate: Cell, existing_cells: List[Cell],
                              grid: SpatialHashGrid) -> bool:
        """``_check_collision`` restricted to the grid neighbors of ``candidate``.

        Uses the identical distance comparison, so decisions (and therefore
        seeded packings) match the exhaustive check exactly.
        """
        for index in grid.neighbors(candidate.center):
            cell = existing_cells[index]
            distance = np.linalg.norm(candidate.center - cell.center)
            min_distance = candidate.radius + cell.radius + self.min_spacing

            if distance < min_distance:
                return True

        return False

    def _is_valid_placement(self, cell: Cell) -> bool:
        """
        Check if cell placement is valid given boundary constraints.

        Args:
            cell: Cell to validate

        Returns:
            True if valid, False otherwise
        """
        if self.allow_boundary_cells:
            # Only require center to be within bounds
            return cell.intersects_bounds(self.bounds)
        else:
            # Require entire cell to be within bounds
            return cell.is_within_bounds(self.bounds)

    def pack(self, max_attempts: int = 1000) -> List[Cell]:
        """
        Pack cells into tissue using random placement.

        Args:
            max_attempts: Maximum placement attempts before stopping

        Returns:
            List of successfully placed Cell objects
        """
        return self._pack(max_attempts)

    def pack_with_progress(self, max_attempts: int = 1000,
                          callback=None) -> List[Cell]:
        """
        Pack cells with progress callback for GUI updates.

        Args:
            max_attempts: Maximum placement attempts before stopping
            callback: Function called with (cells_placed, total_attempts)

        Returns:
            List of successfully placed Cell objects
        """
        return self._pack(max_attempts, callback)

    def _pack(self, max_attempts: int, callback=None) -> List[Cell]:
        """Random sequential addition shared by ``pack`` and ``pack_with_progress``."""
        cells = []
        grid = self._neighbor_grid()
        failed_attempts = 0
        total_attempts = 0

        while failed_attempts < max_attempts:
            total_attempts += 1

            # Select cell type and radius
            cell_type, radius = self._select_cell_type_and_radius()

            # Generate random position
            position = self._generate_random_position(radius)

            # Create candidate cell
            candidate = Cell(
                center=position,
                radius=radius,
                cell_type=cell_type
            )

            # Check if placement is valid
            if not self._is_valid_placement(candidate):
                failed_attempts += 1
                continue

            # Check for collisions
            if self._check_collision_grid(candidate, cells, grid):
                failed_attempts += 1
                continue

            # Mark as boundary cell if it extends beyond bounds
            candidate.is_boundary = not candidate.is_within_bounds(self.bounds)

            # Add cell and reset failure counter
            grid.insert(len(cells), candidate.center)
            cells.append(candidate)
            failed_attempts = 0

            # Call progress callback
            if callback and len(cells) % 10 == 0:
                callback(len(cells), total_attempts)

        return cells


@dataclass
class PackingReport:
    """Diagnostics of one :class:`InhomogeneousPacker` run.

    Attributes:
        n_target: Cells requested by the layout.
        n_placed: Cells placed (always ``n_target``).
        n_rsa: Cells placed by hard-core random sequential addition.
        n_inserted: Cells inserted into saturated bins at best clearance.
        n_relaxed: Cells moved by overlap relaxation.
        relaxation_iterations: Relaxation sweeps performed.
        max_displacement: Largest relaxation displacement in µm.
        overlap_fraction: Fraction of cells whose nearest neighbor is closer
            than the hard core ``kappa * (r_i + r_j)`` after packing.
        target_overlap_fraction: The same fraction in the source region.
        saturated_bins: Bins where addition stalled.
        bin_size: Bin side in µm.
        bin_targets: Cell quota per bin, shape ``(rows, cols)``.
        bin_achieved: Cells per bin after packing.
    """
    n_target: int
    n_placed: int
    n_rsa: int
    n_inserted: int
    n_relaxed: int
    relaxation_iterations: int
    max_displacement: float
    overlap_fraction: float
    target_overlap_fraction: float
    saturated_bins: int
    bin_size: float
    bin_targets: np.ndarray
    bin_achieved: np.ndarray

    @property
    def bin_correlation(self) -> float:
        """Pearson correlation of achieved against target cells per bin."""
        target = self.bin_targets.ravel().astype(float)
        achieved = self.bin_achieved.ravel().astype(float)
        if target.std() == 0 or achieved.std() == 0:
            return float('nan')
        return float(np.corrcoef(target, achieved)[0, 1])

    def to_dict(self) -> Dict:
        out = {f.name: getattr(self, f.name) for f in fields(self)}
        out['bin_targets'] = self.bin_targets.tolist()
        out['bin_achieved'] = self.bin_achieved.tolist()
        out['bin_correlation'] = self.bin_correlation
        return out


class InhomogeneousPacker:
    """
    Pack cells so that local density follows a :class:`~tissue_simulator.density.Layout`.

    1. The window is split into square bins; each bin's quota is proportional
       to the layout intensity it covers (largest-remainder rounding).
    2. One ticket per quota cell is shuffled and filled by random sequential
       addition inside its bin, with radii from the layout's
       density-conditioned marks and hard core ``kappa * (r_i + r_j)``. A bin
       that rejects ``max_failures`` consecutive candidates is saturated.
    3. Tickets left over in saturated bins are inserted at the best-clearance
       position of ``insertion_candidates`` draws, then overlaps in those bins
       and a one-bin halo are relaxed by soft-sphere pushes. No cell moves more
       than ``displacement_cap`` from where it was placed, and relaxation stops
       once the overlap fraction reaches the source region's.

    The z coordinate is uniform through the thickness; use a thin slab
    (thickness below one cell diameter) to mirror a 2D source region.
    """

    def __init__(self, bounds: Tuple[float, float, float], layout: Layout,
                 allow_boundary_cells: bool = True,
                 seed=None,
                 bin_size: Optional[float] = None,
                 max_failures: int = 100,
                 insertion_candidates: int = 20,
                 max_relax_iterations: int = 100,
                 displacement_cap: Optional[float] = None,
                 placeholder_type: str = "default"):
        """
        Args:
            bounds: (height, width, thickness) of the tissue; height and width
                must match the layout window.
            layout: Target maps from :meth:`DensityModel.sample_layout`.
            allow_boundary_cells: Allow cells extending beyond the x/y bounds.
            seed: Integer seed, SeedSequence or Generator for the packer RNG.
            bin_size: Quota bin side in µm; defaults to
                ``max(bandwidth / 2, 10)`` rounded to the layout grid.
            max_failures: Consecutive rejections before a bin is saturated.
            insertion_candidates: Positions tried per insertion.
            max_relax_iterations: Cap on relaxation sweeps.
            displacement_cap: Largest relaxation move in µm; defaults to half
                the median radius.
            placeholder_type: Cell type given to placed cells (labels are
                normally assigned afterwards).
        """
        height, width, _ = bounds
        if abs(width - layout.width) > 1e-6 or abs(height - layout.height) > 1e-6:
            raise ValueError(
                f"bounds {bounds!r} do not match the layout window "
                f"({layout.height}, {layout.width}).")
        self.bounds = tuple(float(b) for b in bounds)
        self.layout = layout
        self.allow_boundary_cells = allow_boundary_cells
        self.max_failures = int(max_failures)
        self.insertion_candidates = max(1, int(insertion_candidates))
        self.max_relax_iterations = int(max_relax_iterations)
        self.displacement_cap = (0.5 * layout.marks.median_radius
                                 if displacement_cap is None else float(displacement_cap))
        self.placeholder_type = placeholder_type
        if bin_size is None:
            bin_size = max(0.5 * (layout.bandwidth or 20.0), 10.0)
        self.bin_pixels = max(1, int(round(bin_size / layout.grid_step)))
        self.bin_size = self.bin_pixels * layout.grid_step
        self._rng = np.random.default_rng(seed)
        self.report: Optional[PackingReport] = None

    # -- geometry helpers -----------------------------------------------------

    def _bin_quotas(self) -> np.ndarray:
        p = self.bin_pixels
        ny, nx = self.layout.intensity.shape
        rows, cols = -(-ny // p), -(-nx // p)
        padded = np.zeros((rows * p, cols * p))
        padded[:ny, :nx] = self.layout.intensity
        mass = padded.reshape(rows, p, cols, p).sum(axis=(1, 3)) * self.layout.grid_step ** 2
        return _stochastic_round(mass.ravel(), self.layout.n_target, self._rng).reshape(rows, cols)

    def _bin_pixels(self, b: int):
        """Flat indices and cumulative intensity of the positive pixels in bin ``b``."""
        cached = self._pixel_cache.get(b)
        if cached is None:
            row, col = divmod(int(b), self._cols)
            p = self.bin_pixels
            ny, nx = self.layout.intensity.shape
            block = self.layout.intensity[row * p:min((row + 1) * p, ny),
                                          col * p:min((col + 1) * p, nx)]
            rows, cols = np.nonzero(block > 0)
            cached = ((rows + row * p) * nx + cols + col * p, np.cumsum(block[rows, cols]))
            self._pixel_cache[b] = cached
        return cached

    def _bin_of(self, x: float, y: float) -> int:
        col = min(max(int(x // self.bin_size), 0), self._cols - 1)
        row = min(max(int(y // self.bin_size), 0), self._rows - 1)
        return row * self._cols + col

    def _candidate(self, b: int):
        # Pixel in proportion to intensity (never a zero-intensity void),
        # then a uniform point inside it.
        rng = self._rng
        pixels, cumulative = self._bin_pixels(b)
        k = min(int(np.searchsorted(cumulative, rng.random() * cumulative[-1], side='right')),
                pixels.size - 1)
        iy, ix = divmod(int(pixels[k]), self.layout.intensity.shape[1])
        step = self.layout.grid_step
        x = rng.uniform(ix * step, min((ix + 1) * step, self.bounds[1]))
        y = rng.uniform(iy * step, min((iy + 1) * step, self.bounds[0]))
        local_density = float(self.layout.intensity[iy, ix])
        r = self.layout.marks.draw(rng, local_density)
        height, width, thickness = self.bounds
        if self.allow_boundary_cells or thickness < 2 * r:
            z = rng.uniform(0.0, thickness) if thickness > 0 else 0.0
        else:
            z = rng.uniform(r, thickness - r)
        valid = local_density > 0 and (
            self.allow_boundary_cells or (r <= x <= width - r and r <= y <= height - r))
        return x, y, z, r, valid

    def _fits(self, x: float, y: float, z: float, r: float) -> bool:
        xs, ys, zs, rs, kappa = self._xs, self._ys, self._zs, self._rs, self.layout.kappa
        for j in self._grid.neighbors((x, y, z)):
            limit = kappa * (r + rs[j])
            if (x - xs[j]) ** 2 + (y - ys[j]) ** 2 + (z - zs[j]) ** 2 < limit * limit:
                return False
        return True

    def _clearance(self, x: float, y: float, z: float, r: float) -> float:
        xs, ys, zs, rs, kappa = self._xs, self._ys, self._zs, self._rs, self.layout.kappa
        best = math.inf
        for j in self._grid.neighbors((x, y, z)):
            gap = (math.sqrt((x - xs[j]) ** 2 + (y - ys[j]) ** 2 + (z - zs[j]) ** 2)
                   - kappa * (r + rs[j]))
            best = min(best, gap)
        return best

    def _add(self, x: float, y: float, z: float, r: float) -> None:
        self._grid.insert(len(self._xs), (x, y, z))
        self._xs.append(x)
        self._ys.append(y)
        self._zs.append(z)
        self._rs.append(r)

    def _clip(self, x: float, y: float, z: float, r: float) -> Tuple[float, float, float]:
        height, width, thickness = self.bounds
        margin = 0.0 if self.allow_boundary_cells else r
        x = min(max(x, margin), width - margin)
        y = min(max(y, margin), height - margin)
        z_margin = r if (not self.allow_boundary_cells and thickness >= 2 * r) else 0.0
        z = min(max(z, z_margin), thickness - z_margin)
        return x, y, z

    # -- stages ---------------------------------------------------------------

    def _place_by_addition(self, b: int) -> bool:
        for _ in range(self.max_failures):
            x, y, z, r, valid = self._candidate(b)
            if valid and self._fits(x, y, z, r):
                self._add(x, y, z, r)
                return True
        return False

    def _insert_at_best_clearance(self, b: int) -> None:
        best, best_score = None, -math.inf
        for _ in range(self.insertion_candidates):
            x, y, z, r, valid = self._candidate(b)
            score = self._clearance(x, y, z, r) if valid else -math.inf
            if best is None or score > best_score:
                best, best_score = (x, y, z, r), score
        self._add(*self._clip(*best[:3], best[3]), best[3])

    def _relax(self, deficit_bins) -> Tuple[int, int, float]:
        rows, cols = self._rows, self._cols
        active = set()
        for b in deficit_bins:
            row, col = divmod(int(b), cols)
            for rr in range(max(row - 1, 0), min(row + 2, rows)):
                for cc in range(max(col - 1, 0), min(col + 2, cols)):
                    active.add(rr * cols + cc)

        xs, ys, zs, rs = self._xs, self._ys, self._zs, self._rs
        grid, kappa, rng = self._grid, self.layout.kappa, self._rng
        n = len(xs)
        movable = [i for i in range(n) if self._bin_of(xs[i], ys[i]) in active]
        movable_set = set(movable)
        origin = {i: (xs[i], ys[i], zs[i]) for i in movable}
        cap = self.displacement_cap
        target = self.layout.target_overlap_fraction
        moved = set()
        check = movable
        iterations = 0
        for iterations in range(1, self.max_relax_iterations + 1):
            violators = set()
            pushes = {}
            for i in check:
                fx = fy = fz = 0.0
                for j in grid.neighbors((xs[i], ys[i], zs[i])):
                    if j == i:
                        continue
                    dx, dy, dz = xs[i] - xs[j], ys[i] - ys[j], zs[i] - zs[j]
                    limit = kappa * (rs[i] + rs[j])
                    d2 = dx * dx + dy * dy + dz * dz
                    if d2 >= limit * limit:
                        continue
                    violators.add(i)
                    violators.add(j)
                    d = math.sqrt(d2)
                    if d < 1e-9:
                        angle = rng.uniform(0.0, 2.0 * math.pi)
                        dx, dy, dz, d = math.cos(angle), math.sin(angle), 0.0, 1.0
                    share = 0.5 if j in movable_set else 1.0
                    push = (limit - d) * share / d
                    fx += dx * push
                    fy += dy * push
                    fz += dz * push
                if fx or fy or fz:
                    pushes[i] = (fx, fy, fz)
            if not pushes or len(violators) / n <= target:
                break
            for i in sorted(pushes):
                fx, fy, fz = pushes[i]
                ox, oy, oz = origin[i]
                nx_, ny_, nz_ = xs[i] + fx, ys[i] + fy, zs[i] + fz
                dist = math.sqrt((nx_ - ox) ** 2 + (ny_ - oy) ** 2 + (nz_ - oz) ** 2)
                if dist > cap:
                    scale = cap / dist
                    nx_, ny_, nz_ = ox + (nx_ - ox) * scale, oy + (ny_ - oy) * scale, oz + (nz_ - oz) * scale
                nx_, ny_, nz_ = self._clip(nx_, ny_, nz_, rs[i])
                grid.move(i, (xs[i], ys[i], zs[i]), (nx_, ny_, nz_))
                xs[i], ys[i], zs[i] = nx_, ny_, nz_
                moved.add(i)
            check = sorted(v for v in violators if v in movable_set)

        max_displacement = max(
            (math.sqrt((xs[i] - origin[i][0]) ** 2 + (ys[i] - origin[i][1]) ** 2
                       + (zs[i] - origin[i][2]) ** 2) for i in moved), default=0.0)
        return len(moved), iterations, max_displacement

    def _overlap_fraction(self) -> float:
        xs, ys, zs, rs, kappa = self._xs, self._ys, self._zs, self._rs, self.layout.kappa
        n = len(xs)
        if n == 0:
            return 0.0
        violating = 0
        for i in range(n):
            for j in self._grid.neighbors((xs[i], ys[i], zs[i])):
                limit = kappa * (rs[i] + rs[j])
                if j != i and (xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2 + (zs[i] - zs[j]) ** 2 < limit * limit:
                    violating += 1
                    break
        return violating / n

    def pack(self) -> List[Cell]:
        """Place ``layout.n_target`` cells and store a :class:`PackingReport`."""
        layout = self.layout
        r_max = max(float(r.max()) for r in layout.marks.radii_by_bin)
        self._grid = SpatialHashGrid(max(2.0 * layout.kappa * r_max, 1e-6))
        self._xs, self._ys, self._zs, self._rs = [], [], [], []

        quotas = self._bin_quotas()
        self._rows, self._cols = quotas.shape
        self._pixel_cache = {}
        tickets = np.repeat(np.arange(quotas.size), quotas.ravel())
        self._rng.shuffle(tickets)

        saturated = np.zeros(quotas.size, dtype=bool)
        deficit = []
        for b in tickets:
            if saturated[b] or not self._place_by_addition(b):
                saturated[b] = True
                deficit.append(b)
        n_rsa = len(self._xs)
        for b in deficit:
            self._insert_at_best_clearance(b)
        n_relaxed, iterations, max_displacement = (
            self._relax(set(deficit)) if deficit else (0, 0, 0.0))

        achieved = np.zeros(quotas.size, dtype=int)
        for x, y in zip(self._xs, self._ys):
            achieved[self._bin_of(x, y)] += 1

        self.report = PackingReport(
            n_target=int(layout.n_target), n_placed=len(self._xs), n_rsa=n_rsa,
            n_inserted=len(deficit), n_relaxed=n_relaxed,
            relaxation_iterations=iterations, max_displacement=float(max_displacement),
            overlap_fraction=self._overlap_fraction(),
            target_overlap_fraction=float(layout.target_overlap_fraction),
            saturated_bins=int(saturated.sum()), bin_size=self.bin_size,
            bin_targets=quotas, bin_achieved=achieved.reshape(quotas.shape),
        )

        cells = []
        for x, y, z, r in zip(self._xs, self._ys, self._zs, self._rs):
            cell = Cell(center=(x, y, z), radius=r, cell_type=self.placeholder_type)
            cell.is_boundary = not cell.is_within_bounds(self.bounds)
            cells.append(cell)
        return cells
