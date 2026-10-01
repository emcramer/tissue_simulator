"""Void (hole) inference, summary and placement on the density grid.

Grid conventions follow ``density.py``: arrays are [row=y, col=x] and pixel
centers sit at ((col + 0.5) * step, (row + 0.5) * step).
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

TAU_NN_FACTOR = 1.5
"""tau >= TAU_NN_FACTOR * median nearest-neighbour distance.

Deviation from the plan's 2.5: with 2.5 a hole of radius ~2.5 d_nn is barely
detectable (its empty-space distance peaks at ~R), so r=25 holes among
d_nn~10 points vanish. 1.5 keeps the Poisson term dominant in practice."""
TAU_RADIUS_FACTOR = 2.0
"""tau >= TAU_RADIUS_FACTOR * median radius + grid_step."""
OPENING_STRUCTURE = np.ones((3, 3), dtype=bool)
"""Structuring element for the binary opening of the raw void mask."""
MIN_HOLE_AREA_FACTOR = 1.0
"""Components with area < MIN_HOLE_AREA_FACTOR * pi * tau^2 are discarded."""
FALLBACK_SEPARATION_TAU = 2.0
"""place_voids separation fallback (in tau) when no separation is known."""


def _grid_shape(width, height, step):
    return (max(1, int(math.ceil(height / step))), max(1, int(math.ceil(width / step))))


def summarize_components(void_mask, grid_step: float) -> dict:
    """Summarize connected void components (JSON-serializable).

    Interior holes (not touching the array border) form the diameter pool;
    if there are none, edge-touching components are used instead.
    Elongation is sqrt(major/minor eigenvalue) of the pixel covariance.
    """
    void = np.asarray(void_mask, dtype=bool)
    labels, n = ndi.label(void)  # 4-connectivity
    comps: List[dict] = []
    ny, nx = void.shape
    for k in range(1, n + 1):
        iy, ix = np.nonzero(labels == k)
        area = float(iy.size) * grid_step ** 2
        cx = float((ix.mean() + 0.5) * grid_step)
        cy = float((iy.mean() + 0.5) * grid_step)
        if iy.size > 2:
            ev = np.linalg.eigvalsh(np.cov(np.vstack([ix, iy]).astype(float)))
            elong = float(math.sqrt(ev[1] / ev[0])) if ev[0] > 1e-12 else float("inf")
            if not math.isfinite(elong):
                elong = float(iy.size)
        else:
            elong = 1.0
        touch = bool(iy.min() == 0 or ix.min() == 0 or iy.max() == ny - 1 or ix.max() == nx - 1)
        comps.append({
            "area": area,
            "equivalent_diameter": float(2.0 * math.sqrt(area / math.pi)),
            "centroid": [cx, cy],
            "elongation": elong,
            "edge_touching": touch,
        })
    interior = [c for c in comps if not c["edge_touching"]]
    pool = interior if interior else comps
    cents = [c["centroid"] for c in interior]
    min_sep = None
    if len(cents) >= 2:
        pts = np.asarray(cents)
        d = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
        min_sep = float(d[np.triu_indices(len(pts), 1)].min())
    return {
        "n_holes": len(interior),
        "n_edge_touching": len(comps) - len(interior),
        "equivalent_diameters": [c["equivalent_diameter"] for c in pool],
        "areas": [c["area"] for c in pool],
        "centroids": [c["centroid"] for c in pool],
        "elongations": [c["elongation"] for c in pool],
        "edge_touching": [c["edge_touching"] for c in pool],
        "min_center_separation": min_sep,
        "all_components": comps,
    }


def infer_voids(xs, ys, radii, mask, grid_step, width, height):
    """Infer large empty regions from cell centers.

    Returns ``(tissue_mask, info)``. A pixel is void when its distance to the
    nearest cell center exceeds ``tau = max(sqrt(ln N / (pi lambda)),
    1.5 d_nn, 2 r_med + grid_step)`` (lambda = cell density, N = tissue
    pixels). The raw void core is regrown by tau, opened (3x3) and components smaller than
    ``pi tau^2`` are dropped.
    """
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    shape = _grid_shape(width, height, grid_step)
    mask = np.ones(shape, bool) if mask is None else np.asarray(mask, bool)
    n_pix = int(mask.sum())
    area = n_pix * grid_step ** 2
    n = xs.size
    r_med = float(np.median(radii)) if len(radii) else 0.0
    pts = np.column_stack([xs, ys])
    if n >= 2:
        d_nn = float(np.median(cKDTree(pts).query(pts, k=2)[0][:, 1]))
    else:
        d_nn = 0.0
    lam = n / area if area > 0 else 0.0
    t_poisson = math.sqrt(math.log(max(n_pix, 2)) / (math.pi * lam)) if lam > 0 else 0.0
    tau = max(t_poisson, TAU_NN_FACTOR * d_nn, TAU_RADIUS_FACTOR * r_med + grid_step)
    ny, nx = shape
    gy, gx = np.mgrid[0:ny, 0:nx]
    cx = (gx[mask] + 0.5) * grid_step
    cy = (gy[mask] + 0.5) * grid_step
    void = np.zeros(shape, bool)
    if n > 0 and n_pix > 0:
        dist = cKDTree(pts).query(np.column_stack([cx, cy]))[0]
        void[mask] = dist > tau
    if void.any():
        # The raw void is only the *core* of a hole (a hole of radius R has
        # empty-space distance > tau only within R - tau of its center), so
        # grow it back by tau (restricted to the tissue mask).
        dt = ndi.distance_transform_edt(~void) * grid_step
        void = (dt <= tau) & mask
    void = ndi.binary_opening(void, structure=OPENING_STRUCTURE)
    labels, nl = ndi.label(void)
    if nl:
        sizes = ndi.sum(void, labels, index=np.arange(1, nl + 1)) * grid_step ** 2
        small = np.nonzero(sizes < MIN_HOLE_AREA_FACTOR * math.pi * tau ** 2)[0] + 1
        void[np.isin(labels, small)] = False
    info = summarize_components(void, grid_step)
    info.update(tau=float(tau), d_nn=d_nn, r_med=r_med, inferred=True)
    return mask & ~void, info


def place_voids(rng, shape, grid_step, info, tissue_mask, anchor=None, max_tries=50):
    """Place ``info['n_holes']`` discs into ``tissue_mask`` (copy returned).

    Diameters are drawn first (with replacement), then centers. Discs lie
    fully inside the window (clipped if too large) and keep pairwise center
    separation >= ``min_center_separation`` (fallback 2 tau). With ``anchor``
    the largest disc is centered there. A hole failing ``max_tries`` is skipped.
    """
    ny, nx = shape
    W, H = nx * grid_step, ny * grid_step
    out = np.array(tissue_mask, bool, copy=True)
    n = int(info.get("n_holes", 0))
    pool = info.get("equivalent_diameters") or []
    placed = {"centers": [], "radii": [], "n_requested": n, "n_placed": 0,
              "anchored": False}
    if n == 0 or not pool:
        return out, placed
    diam = np.sort(np.asarray(rng.choice(np.asarray(pool, float), size=n, replace=True)))[::-1]
    sep = info.get("min_center_separation")
    if sep is None:
        sep = FALLBACK_SEPARATION_TAU * float(info.get("tau", 0.0))
    centers, radii = [], []
    gy, gx = np.mgrid[0:ny, 0:nx]
    px = (gx + 0.5) * grid_step
    py = (gy + 0.5) * grid_step
    for i, dd in enumerate(diam):
        r = min(dd / 2.0, W / 2.0, H / 2.0)
        if i == 0 and anchor is not None:
            c = (float(anchor[0]), float(anchor[1]))
            placed["anchored"] = True
        else:
            c = None
            for _ in range(max_tries):
                cand = (rng.uniform(r, W - r), rng.uniform(r, H - r))
                if all(math.hypot(cand[0] - q[0], cand[1] - q[1]) >= sep for q in centers):
                    c = cand
                    break
            if c is None:
                continue
        centers.append(c)
        radii.append(float(r))
        out[(px - c[0]) ** 2 + (py - c[1]) ** 2 <= r * r] = False
    placed["centers"] = [[float(a), float(b)] for a, b in centers]
    placed["radii"] = radii
    placed["n_placed"] = len(centers)
    return out, placed
