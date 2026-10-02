"""Germ-grain unit model: per-kind radial profiles, grain placement, rasterization.

Units (tumor nests, follicles, glomeruli) are discs with a radial profile in the
normalized distance ``s = r / outer_radius``. Only placement uses an RNG.
"""

import math

import numpy as np

S_MAX = 1.5
"""Largest normalized distance covered by a profile (s > 1 is the surroundings)."""

PROFILE_PRIOR_CELLS = 5.0
"""Pseudo-cells of Dirichlet prior pulling knot composition toward the global one."""

DEFAULT_N_KNOTS = 8
"""Default number of equal-width radial bins on s in [0, S_MAX]."""

DENSITY_FLOOR = 0.05
"""Relative density floor for non-lumen knots with zero cells (keeps the model positive)."""

AREA_GRID_PIXELS = 150
"""Pixels along the shorter window side when estimating in-window annulus areas."""

RELAX_MARGIN = 5.0
"""Extra micrometers added to the sum of outer radii when separation is relaxed."""

KINDS = ("blob", "ring_core", "ring_lumen")
"""Unit kinds understood by this module."""


def _window_density(xs, width, height, window_density):
    if window_density is not None:
        return float(window_density)
    if width <= 0 or height <= 0:
        raise ValueError("width and height must be positive")
    return len(xs) / float(width * height)


def _annulus_areas(center, R, edges, width, height):
    """In-window area of each annulus of s (grid estimate, window-clipped)."""
    step = min(width, height) / float(AREA_GRID_PIXELS)
    gx = (np.arange(int(math.ceil(width / step))) + 0.5) * step
    gy = (np.arange(int(math.ceil(height / step))) + 0.5) * step
    inside_x = gx < width
    inside_y = gy < height
    dx = (gx[inside_x] - center[0]) ** 2
    dy = (gy[inside_y] - center[1]) ** 2
    s = np.sqrt(dx[None, :] + dy[:, None]) / R
    s = s[s < edges[-1]]
    h, _ = np.histogram(s, bins=edges)
    return h * step * step


def fit_unit_profiles(xs, ys, type_idx, n_types, radii, units, width, height,
                      n_knots=DEFAULT_N_KNOTS, window_density=None):
    """Pool radial density/composition profiles over units of the same kind.

    ``units`` is the ``detect_units`` dict (``units["units"]`` list). Returns a
    JSON-serializable dict keyed by kind with ``s_knots``, ``density_rel``,
    ``composition`` (n_knots x n_types), ``n_cells_per_knot``, ``inner_ratio``,
    ``n_units``. ``width``/``height`` (window size, required) give the window area; annulus areas
    are clipped to the window. ``window_density`` optionally overrides the
    window-mean density (default n_cells / (width * height)). Non-lumen knots with
    zero cells get ``density_rel = DENSITY_FLOOR``. ``radii`` is unused.
    """
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    ti = np.asarray(type_idx, int)
    lam0 = _window_density(xs, width, height, window_density)
    glob = np.bincount(ti, minlength=n_types).astype(float)
    glob = glob / glob.sum() if glob.sum() > 0 else np.full(n_types, 1.0 / n_types)
    edges = np.linspace(0.0, S_MAX, n_knots + 1)
    knots = 0.5 * (edges[:-1] + edges[1:])
    by_kind = {}
    for u in (units or {}).get("units", []):
        by_kind.setdefault(u["kind"], []).append(u)
    out = {}
    for kind, us in by_kind.items():
        counts = np.zeros((n_knots, n_types))
        area = np.zeros(n_knots)
        ratios = []
        for u in us:
            R = float(u["outer_radius"])
            if R <= 0:
                continue
            ratios.append(float(u.get("inner_radius", 0.0)) / R)
            d = np.hypot(xs - u["center"][0], ys - u["center"][1]) / R
            sel = d < S_MAX
            b = np.minimum((d[sel] / S_MAX * n_knots).astype(int), n_knots - 1)
            np.add.at(counts, (b, ti[sel]), 1.0)
            area += _annulus_areas(u["center"], R, edges, width, height)
        inner = float(np.mean(ratios)) if ratios else 0.0
        n_k = counts.sum(axis=1)
        dens = counts.sum(axis=1) / np.maximum(area, 1e-12) / max(lam0, 1e-12)
        if kind != "ring_lumen":
            dens = np.where(n_k == 0, np.maximum(dens, DENSITY_FLOOR), dens)
        if kind == "ring_lumen":
            dens[knots < inner] = 0.0
        comp = counts + PROFILE_PRIOR_CELLS * glob[None, :]
        comp = comp / comp.sum(axis=1, keepdims=True)
        out[kind] = {
            "s_knots": knots.tolist(),
            "density_rel": dens.tolist(),
            "composition": comp.tolist(),
            "n_cells_per_knot": n_k.tolist(),
            "inner_ratio": inner,
            "n_units": len(us),
        }
    return out


def _sep(min_separation, r_q, R, relax):
    base = r_q + R
    if min_separation is None:
        return base
    return min(min_separation, base + RELAX_MARGIN) if relax else min_separation


def place_grains(rng, shape, grid_step, units_summary, min_separation, max_tries=50,
                 anchor=None):
    """Place the source number of units per kind with hard-core centers.

    Outer radii are drawn first (per kind, with replacement; inner radius =
    outer * kind mean inner ratio), then centers, largest unit first. Discs lie
    inside the window (radius clipped to half the window). After ``max_tries`` at the requested separation the remaining grains retry with
    it relaxed to sum of outer radii + ``RELAX_MARGIN`` (``separation_relaxed``
    set True). Separation is
    ``min_separation`` (fallback: sum of the two outer radii). A unit failing
    ``max_tries`` is skipped and reported as ``shortfall`` (= n_requested -
    n_placed). ``min_separation=None`` means the sum-of-radii fallback, 0 means no
    constraint. ``anchor`` applies to the largest unit and is clipped into the
    window.
    """
    ny, nx = shape
    W, H = nx * grid_step, ny * grid_step
    by_kind = {}
    for u in (units_summary or {}).get("units", []):
        by_kind.setdefault(u["kind"], []).append(u)
    outer, inner, kinds = [], [], []
    for kind in sorted(by_kind):
        us = by_kind[kind]
        pool = np.array([float(u["outer_radius"]) for u in us])
        ratio = float(np.mean([float(u.get("inner_radius", 0.0)) / float(u["outer_radius"])
                               for u in us]))
        draw = rng.choice(pool, size=len(us), replace=True)
        for R in draw:
            R = min(float(R), W / 2.0, H / 2.0)
            outer.append(R)
            inner.append(R * ratio)
            kinds.append(kind)
    n = len(outer)
    placed = {"centers": [], "outer_radii": [], "inner_radii": [], "kinds": [],
              "n_requested": n, "n_placed": 0, "shortfall": n, "anchored": False,
              "separation_relaxed": False}
    order = sorted(range(n), key=lambda i: -outer[i])
    centers = []
    for rank, i in enumerate(order):
        R = outer[i]
        if rank == 0 and anchor is not None:
            c = (min(max(float(anchor[0]), 0.0), W), min(max(float(anchor[1]), 0.0), H))
            placed["anchored"] = True
        else:
            c = None
            for relax in (False, True):
                for _ in range(max_tries):
                    cand = (rng.uniform(R, W - R), rng.uniform(R, H - R))
                    if all(math.hypot(cand[0] - q[0], cand[1] - q[1])
                           >= _sep(min_separation, q[2], R, relax) for q in centers):
                        c = cand
                        break
                if c is not None:
                    placed["separation_relaxed"] |= relax
                    break
            if c is None:
                continue
        centers.append((c[0], c[1], R))
        placed["centers"].append([float(c[0]), float(c[1])])
        placed["outer_radii"].append(float(R))
        placed["inner_radii"].append(float(inner[i]))
        placed["kinds"].append(kinds[i])
    placed["n_placed"] = len(centers)
    placed["shortfall"] = n - len(centers)
    return placed


def rasterize_grains(placed, profiles, shape, grid_step, proportions):
    """Rasterize placed grains to (factor, composition, grain_mask, void_mask).

    Each pixel takes the grain with the smallest normalized distance s (nearest
    grain wins); pixels with s > S_MAX are outside (factor 1, composition
    ``proportions``). Lumen pixels (s < inner/outer) get factor 0 and void True.
    """
    ny, nx = shape
    props = np.asarray(proportions, float)
    gy, gx = np.mgrid[0:ny, 0:nx]
    px = (gx + 0.5) * grid_step
    py = (gy + 0.5) * grid_step
    best = np.full(shape, np.inf)
    factor = np.ones(shape)
    comp = np.tile(props[:, None, None], (1, ny, nx))
    void = np.zeros(shape, bool)
    for c, R, r_in, kind in zip(placed["centers"], placed["outer_radii"],
                                placed["inner_radii"], placed["kinds"]):
        prof = profiles.get(kind)
        if prof is None or R <= 0:
            continue
        s = np.hypot(px - c[0], py - c[1]) / R
        win = (s <= S_MAX) & (s < best)
        if not win.any():
            continue
        knots = np.asarray(prof["s_knots"])
        dens = np.asarray(prof["density_rel"])
        cmp_ = np.asarray(prof["composition"])
        sw = s[win]
        f = np.interp(sw, knots, dens)
        lumen = (sw < (r_in / R)) if kind == "ring_lumen" else np.zeros(sw.shape, bool)
        f = np.where(lumen, 0.0, f)
        factor[win] = f
        void[win] = lumen
        for t in range(comp.shape[0]):
            comp[t][win] = np.interp(sw, knots, cmp_[:, t])
        best[win] = sw
    comp = comp / np.maximum(comp.sum(axis=0, keepdims=True), 1e-12)
    return factor, comp, np.isfinite(best), void
