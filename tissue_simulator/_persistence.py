"""Persistent homology of 2-D intensity maps (pure numpy, no RNG).

Two filtrations of a scalar field ``f`` on a pixel grid are computed with
union-find in O(N log N) (sorting dominates):

* ``superlevel_h0``: components of ``{f >= t}`` (8-connectivity), swept from
  high to low.  A component is born at its peak; when two components merge
  the younger one (lower birth) dies (elder rule).  Death = merge level.
  The oldest component of each connected region is *essential*; its death is
  defined as ``min(f)`` over the support so that it has finite persistence and
  is included in barcodes and distances.
* ``superlevel_h1``: loops of ``{f >= t}``.  In the plane (Alexander duality)
  an H1 class of the superlevel set corresponds to a *bounded* connected
  component of the complementary sublevel set ``{f < t}``.  The superlevel
  set uses 8-connectivity, so the sublevel set uses 4-connectivity (the
  dual pairing that keeps the two sets topologically consistent).  Sweeping
  the sublevel set upward, a component is born at its minimum; on a merge
  the component with the higher minimum dies (elder rule) at the merge level
  (the saddle).  A virtual boundary node, which is older than everything,
  is joined to every active pixel on the array edge or next to a pixel
  outside the mask; components that reach it are unbounded, so a basin
  touching the border never yields a bar until it merges with the boundary,
  and the boundary component itself never yields one.  A dying component
  therefore gives an H1 bar with birth = saddle level, death = basin
  minimum (the superlevel loop is born at the saddle and filled in at the
  minimum when the level is lowered), persistence = saddle - minimum.

Mask support: pixels outside ``mask`` are absent for H0 and count as part of
the boundary for H1.  Zero-persistence bars (ties) are dropped, except
essential H0 bars (flagged ``essential=True``); ``min_persistence`` filters
the rest.
"""
from __future__ import annotations

import numpy as np

__all__ = ["superlevel_h0", "superlevel_h1", "barcode_summary", "persistence_distance"]


def _find(parent, a):
    root = a
    while parent[root] != root:
        root = parent[root]
    while parent[a] != root:
        parent[a], a = root, parent[a]
    return root


def _prep(f, mask):
    f = np.asarray(f, dtype=float)
    if f.ndim != 2:
        raise ValueError("f must be 2-D")
    valid = np.ones(f.shape, bool) if mask is None else np.asarray(mask, bool)
    if valid.shape != f.shape:
        raise ValueError("mask shape must match f")
    valid = valid & np.isfinite(f)
    h, w = f.shape
    fp = np.zeros((h + 2, w + 2))
    vp = np.zeros((h + 2, w + 2), bool)
    fp[1:-1, 1:-1] = f
    vp[1:-1, 1:-1] = valid
    return fp.ravel(), vp.ravel(), w + 2, (h + 2, w + 2)


def superlevel_h0(f, mask=None, connectivity=8, min_persistence=0.0):
    """H0 bars of the superlevel filtration.

    Returns a list of dicts ``birth, death, persistence, peak (row, col),
    area_at_death`` sorted by decreasing persistence.
    """
    if connectivity not in (4, 8):
        raise ValueError("connectivity must be 4 or 8")
    fl, vl, W, shape = _prep(f, mask)
    idx = np.flatnonzero(vl)
    if idx.size == 0:
        return []
    order = idx[np.argsort(-fl[idx], kind="stable")]
    offs = [-W, W, -1, 1] + ([-W - 1, -W + 1, W - 1, W + 1] if connectivity == 8 else [])
    parent = {}
    peak, size = {}, {}
    bars = []
    for p in order.tolist():
        parent[p] = p
        peak[p], size[p] = p, 1
        v = fl[p]
        for o in offs:
            q = p + o
            if q not in parent:
                continue
            a, b = _find(parent, p), _find(parent, q)
            if a == b:
                continue
            # elder = strictly higher birth; exact ties keep the component that
            # appeared first in the sweep (q's), so a flat field's essential
            # peak is the first pixel in sweep order
            if not fl[peak[a]] > fl[peak[b]]:
                a, b = b, a
            # b is younger and dies at v
            pers = fl[peak[b]] - v
            if pers > 0 and pers >= min_persistence:
                bars.append((fl[peak[b]], v, peak[b], size[b], False))
            parent[b] = a
            size[a] += size[b]
    lo = fl[idx].min()
    for r in {_find(parent, p) for p in parent}:
        bars.append((fl[peak[r]], lo, peak[r], size[r], True))
    out = []
    for birth, death, pk, area, ess in bars:
        out.append({"birth": float(birth), "death": float(death),
                    "persistence": float(birth - death),
                    "peak": (pk // W - 1, pk % W - 1), "area_at_death": int(area),
                    "essential": ess})
    out.sort(key=lambda d: -d["persistence"])
    return out


def superlevel_h1(f, mask=None, min_persistence=0.0):
    """H1 bars of the superlevel filtration via sublevel-set duality.

    Returns dicts ``birth`` (saddle), ``death`` (basin minimum),
    ``persistence``, ``basin_min (row, col)``, ``basin_area`` (pixels of the
    basin just before it merges), sorted by decreasing persistence.
    """
    fl, vl, W, shape = _prep(f, mask)
    idx = np.flatnonzero(vl)
    if idx.size == 0:
        return []
    order = idx[np.argsort(fl[idx], kind="stable")]
    offs = [-W, W, -1, 1]
    BND = -1
    parent = {BND: BND}
    mn = {BND: -np.inf}
    mpos, size = {BND: BND}, {BND: 0}
    bars = []

    def union(a, b, level):
        a, b = _find(parent, a), _find(parent, b)
        if a == b:
            return
        if mn[a] > mn[b]:  # b is the elder; a dies
            a, b = b, a
        # now a is the elder (lower min), b dies; the boundary is always elder
        if b != BND and level - mn[b] > 0 and level - mn[b] >= min_persistence:
            bars.append((level, mn[b], mpos[b], size[b]))
        parent[b] = a
        size[a] += size[b]

    for p in order.tolist():
        parent[p] = p
        mn[p], mpos[p], size[p] = fl[p], p, 1
        v = fl[p]
        for o in offs:
            q = p + o
            if not vl[q]:
                union(p, BND, v)
            elif q in parent:
                union(p, q, v)
    out = []
    for birth, death, pk, area in bars:
        out.append({"birth": float(birth), "death": float(death),
                    "persistence": float(birth - death),
                    "basin_min": (pk // W - 1, pk % W - 1), "basin_area": int(area)})
    out.sort(key=lambda d: -d["persistence"])
    return out


def barcode_summary(bars, k=5):
    """Summary dict: ``count``, ``top`` (k largest persistences), ``total``."""
    pers = sorted((float(b["persistence"]) for b in bars), reverse=True)
    return {"count": len(pers), "top": pers[:k], "total": float(sum(pers))}


def persistence_distance(bars_a, bars_b, drop_essential=False):
    """Wasserstein-1 distance between sorted persistence vectors.

    Both vectors are sorted in decreasing order and zero-padded to equal
    length; the distance is the sum of absolute differences.  Essential H0
    bars (death = min f) are included as ordinary finite bars unless
    ``drop_essential=True``; the H0 distance is otherwise dominated by the
    field range (the essential bar's persistence is max f - min f).
    """
    def vec(bars):
        return np.sort([x["persistence"] for x in bars
                        if not (drop_essential and x.get("essential", False))])[::-1]
    a, b = vec(bars_a), vec(bars_b)
    n = max(len(a), len(b))
    a = np.pad(a, (0, n - len(a)))
    b = np.pad(b, (0, n - len(b)))
    return float(np.abs(a - b).sum())
