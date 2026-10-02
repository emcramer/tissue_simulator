"""Detection of compact tissue units (nests, follicles, glomeruli) by persistent homology.

Inputs are per-type intensity maps ``lam_types`` (n_types, ny, nx), the
empty-space distance map (distance from each pixel center to the nearest cell
center, large in lumens) and a few stationary null draws of the same maps.
``lam_types`` may also be 4-D ``(n_scales, n_types, ny, nx)`` (see
``UNIT_BANDWIDTH_SCALES``); threshold keys then carry an ``@scale`` suffix.
Features of three kinds are extracted with :mod:`._persistence`:

* H0 bars of ``lam_types[t]``: compact blobs (peaks) of type ``t``;
* H1 bars of ``lam_types[t]``: rings of type ``t`` around a basin;
* H0 bars of the empty-space map: lumens.

Significance: for each (type, degree) and for the lumen map the threshold is
the ``quantile`` (``UNIT_NULL_QUANTILE``) of the maximum bar persistence over
the null draws (the null's essential H0 bar counts, so the comparison with the
data's own essential bar is like for like).  Without null draws the threshold
falls back to ``UNIT_FALLBACK_FRACTION`` of the map's own largest persistence.

With few null draws (``fit`` uses at most ``UNIT_NULL_DRAWS`` = 9) the 95th
percentile of the null max is effectively the null maximum, and about 5 % false
positives are expected per map, scale and degree.  Multiplicity control at the
fine scale: a 0.5x blob with no 1.0x footprint containing its peak must exceed
the null *maximum* and have outer radius >= ``UNIT_MIN_BLOB_RADIUS_NN`` nearest-
neighbour distances (``d_nn``); otherwise it is kept only when paired as a core
or lumen of a strong partner (escalated: see ``FINE_ONLY_BLOBS_ALLOWED``, default
no fine-only blobs).  Blobs whose peak lies on the array border are ignored.

Primary detector: H0 blobs plus footprint pairing (a lumen or another type's
blob inside a blob's/ring's filled footprint).  H1 rings contribute only at the
finer scale in practice (significant for follicles at 0.5x, none for glomeruli
or arteries); pairing is by footprint containment, not by H1 basins.

Compactness: a unit must be small.  Features whose component (at half
persistence for blobs, filled component at the loop's birth level for rings)
covers more than ``MAX_UNIT_AREA_FRACTION`` of the window are discarded.  This
removes the stroma's own "ring" around a unit, its essential H0 bar and
large-scale trends.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from ._persistence import barcode_summary, superlevel_h0, superlevel_h1

UNIT_NULL_QUANTILE = 0.95
"""Quantile of the null max-persistence distribution used as the threshold."""
UNIT_FALLBACK_FRACTION = 0.5
"""Threshold as a fraction of the map's own max persistence when no null draws exist."""
UNIT_NULL_DRAWS = 9
"""Stationary null layouts drawn by ``DensityModel.fit`` for unit thresholds."""
MAX_UNIT_AREA_FRACTION = 0.25
"""Largest allowed unit footprint as a fraction of the window area."""
UNIT_BANDWIDTH_SCALES = (1.0, 0.5)
"""Multipliers of the per-type bandwidths at which the intensity maps are
analysed. The fitted bandwidth smooths rings narrower than itself into blobs
and sparse core types into noise, so a finer map (half bandwidth) is analysed
too, with its own null thresholds; candidates from all scales are pooled."""
FINE_ONLY_BLOBS_ALLOWED = False
"""Escalated multiplicity rule: Potts-like clustering controls still produced two
0.5x-only blobs that cleared the null max and the radius bar, so fine-scale-only
blobs are never reported on their own (only as a core or lumen paired with a
strong partner). Set True to use the null-max + radius rule alone."""
UNIT_MIN_BLOB_RADIUS_NN = 2.0
"""Minimum outer radius, in nearest-neighbour distances, of a fine-scale-only blob."""
MIN_UNIT_RADIUS_PX = 1.5
"""Features smaller than this many pixels in radius are discarded."""
UNIT_WEAK_QUANTILE = 0.5
"""Null quantile of the max persistence above which a feature is *weak*
evidence. A weak feature is kept only when paired (a lumen inside a ring or
blob, or a core inside a ring) with a partner that is strong (above
``UNIT_NULL_QUANTILE``): both being chance features is rarer than either alone,
and a few-cell ring or a lumen of a few nearest-neighbour distances is
marginal on its own."""


def empty_space_map(xs, ys, shape, grid_step) -> np.ndarray:
    """Distance (µm) from each pixel center to the nearest cell center."""
    ny, nx = shape
    gy, gx = np.mgrid[0:ny, 0:nx]
    pts = np.column_stack([(gx.ravel() + 0.5) * grid_step, (gy.ravel() + 0.5) * grid_step])
    d, _ = cKDTree(np.column_stack([xs, ys])).query(pts)
    return d.reshape(shape)


def _max_persistence(bars) -> float:
    return max((b["persistence"] for b in bars), default=0.0)


def _thresholds(own_max, null_max, quantile):
    if null_max:
        return float(np.quantile(null_max, quantile))
    return UNIT_FALLBACK_FRACTION * own_max * (quantile / UNIT_NULL_QUANTILE)


def _component_at(f, level, rc, structure=np.ones((3, 3), bool)):
    lab, _ = ndi.label(f >= level, structure=structure)
    k = lab[rc]
    return lab == k if k > 0 else np.zeros(f.shape, bool)


def _blob_geometry(f, bar):
    """Radius (px) of an H0 blob: ``min(area_at_death, area at half persistence)``.

    The half-persistence component is the one containing the peak at level
    ``birth - persistence / 2``; capping with it keeps essential bars (whose
    death area is the whole window) and bars that merge late from absorbing
    their surroundings.
    """
    comp = _component_at(f, bar["birth"] - 0.5 * bar["persistence"], tuple(bar["peak"]))
    area = min(float(comp.sum()), float(bar["area_at_death"]))
    return math.sqrt(area / math.pi), comp


def _ring_geometry(f, bar):
    """Outer radius (px) and basin mask of an H1 ring.

    Outer radius = equivalent radius of the hole-filled superlevel component
    of ``f`` at the loop's birth level that surrounds the basin.  The basin is
    the bounded component of ``{f < birth}`` (4-connectivity) holding the
    basin minimum.  Returns ``(None, basin)`` when no surrounding component is found.
    """
    rc = tuple(bar["basin_min"])
    below, _ = ndi.label(f < bar["birth"])
    basin = below == below[rc]
    ring_lab, _ = ndi.label(f >= bar["birth"], structure=np.ones((3, 3), bool))
    near = ring_lab[ndi.binary_dilation(basin)]
    near = near[near > 0]
    if near.size == 0:
        return None, basin
    k = np.bincount(near).argmax()
    filled = ndi.binary_fill_holes(ring_lab == k)
    return math.sqrt(float(filled.sum()) / math.pi), basin


def _filled(mask):
    return ndi.binary_fill_holes(mask)


def detect_units(lam_types, cell_types: Sequence[str], empty_space, grid_step: float,
                 mask, null_draws: Optional[List[Tuple[np.ndarray, np.ndarray]]] = None,
                 *, quantile: float = UNIT_NULL_QUANTILE, d_nn: Optional[float] = None, scales: Optional[Sequence[float]] = None
                 ) -> dict:
    """Detect compact units; returns ``{}`` when nothing is significant.

    Args:
        lam_types: ``(n_types, ny, nx)`` or ``(n_scales, n_types, ny, nx)``
            per-type intensity maps.
        cell_types: Type names.
        empty_space: ``(ny, nx)`` distance-to-nearest-cell map (µm).
        grid_step: Pixel size (µm).
        mask: Boolean tissue mask (unused beyond API symmetry; maps are
            filtered over the full grid so that lumens stay in-domain).
        null_draws: ``[(lam_types_null, empty_space_null), ...]``, same shapes.
        quantile: Null quantile of the max persistence used as threshold.
        scales: Labels of the scale axis (default ``UNIT_BANDWIDTH_SCALES``).

    Output (JSON-friendly)::

        {"units": [{"kind": "blob"|"ring_core"|"ring_lumen", "center": [x, y],
                    "outer_radius", "inner_radius", "ring_type", "core_type",
                    "persistence", "score", "scale", "edge_touching"}, ...],
         "n_units", "kinds": {kind: count}, "thresholds": {key: float},
         "min_center_separation", "radii": [outer radii], "barcode_top": {key: top-5 persistences},
         "null_draws": int, "quantile", "window_area_px"}

    Candidate features (all compact, significant): *blobs* (H0 of type t,
    footprint = hole-filled component at half persistence, radius
    ``min(area_at_death, footprint area)``), *rings* (H1 of type t; outer
    radius from the hole-filled component at the birth level, inner radius from
    the basin) and *lumens* (H0 of the empty-space map; radius = its peak
    distance).  A blob or ring of type t whose footprint contains a lumen peak
    (pair rule: at least one of the two is above its strong threshold, both above the weak one) is ``ring_lumen`` (inner radius = lumen radius, centered on the lumen);
    else if it contains the peak of a significant blob of another type u it is
    ``ring_core`` with ``core_type`` u (centered on the core); else it is a
    ``blob`` of type t.  Consumed children are not reported separately and
    blobs inside an accepted unit's outer radius are absorbed.  Units are
    ordered by (kind is not blob, score = persistence / threshold) and
    units whose center lies inside the outer radius of an earlier-ranked unit
    (duplicates, e.g. the same unit found at two scales, or bumps on a ring) are dropped.
    """
    lam_types = np.asarray(lam_types, dtype=float)
    if lam_types.ndim == 3:
        lam_types = lam_types[None]
    n_scales, n_types, ny, nx = lam_types.shape
    scales = tuple(scales if scales is not None else
                   (UNIT_BANDWIDTH_SCALES if n_scales == len(UNIT_BANDWIDTH_SCALES) else (1.0,) * n_scales))
    window_area = float(ny * nx)
    null_draws = list(null_draws or [])
    full = np.ones((ny, nx), bool)

    def nulls(f):
        return [f(np.asarray(ln, dtype=float).reshape(lam_types.shape), e) for ln, e in null_draws]

    thr, barcodes = {}, {}
    bars0, bars1 = {}, {}
    for si in range(n_scales):
        suffix = f"@{scales[si]:g}" if n_scales > 1 else ""
        for t, name in enumerate(cell_types):
            bars0[si, t] = superlevel_h0(lam_types[si, t], full)
            bars1[si, t] = superlevel_h1(lam_types[si, t], full)
            n0 = nulls(lambda ln, e: _max_persistence(superlevel_h0(ln[si, t], full)))
            n1 = nulls(lambda ln, e: _max_persistence(superlevel_h1(ln[si, t], full)))
            for key, bars, vals in ((f"{name}:h0{suffix}", bars0[si, t], n0),
                                    (f"{name}:h1{suffix}", bars1[si, t], n1)):
                thr[key] = _thresholds(_max_persistence(bars), vals, quantile)
                thr[key + "~max"] = float(max(vals)) if vals else thr[key]
                thr[key + "~weak"] = min(thr[key], _thresholds(_max_persistence(bars), vals,
                                                               UNIT_WEAK_QUANTILE))
            barcodes[f"{name}{suffix}"] = {"h0": barcode_summary(bars0[si, t], 5)["top"],
                                          "h1": barcode_summary(bars1[si, t], 5)["top"]}
    es = np.asarray(empty_space, dtype=float)
    lum_bars = superlevel_h0(es, full)
    nl = [_max_persistence(superlevel_h0(np.asarray(e, dtype=float), full)) for _, e in null_draws]
    thr["lumen"] = _thresholds(_max_persistence(lum_bars), nl, quantile)
    thr["lumen~weak"] = min(thr["lumen"], _thresholds(_max_persistence(lum_bars), nl, UNIT_WEAK_QUANTILE))

    def center(rc):
        return [float((rc[1] + 0.5) * grid_step), float((rc[0] + 0.5) * grid_step)]

    def too_big(r_px):
        return math.pi * r_px ** 2 > MAX_UNIT_AREA_FRACTION * window_area or r_px < MIN_UNIT_RADIUS_PX

    lumens = []
    for b in lum_bars:
        peak = tuple(b["peak"])
        if b["persistence"] > thr["lumen~weak"] and not too_big(es[peak] / grid_step):
            lumens.append({"peak": peak, "radius": float(es[peak]),
                           "strong": bool(b["persistence"] > thr["lumen"]),
                           "score": b["persistence"] / max(thr["lumen"], 1e-300),
                           "persistence": float(b["persistence"])})

    features = []  # blobs and rings
    for si in range(n_scales):
        suffix = f"@{scales[si]:g}" if n_scales > 1 else ""
        for t, name in enumerate(cell_types):
            f = lam_types[si, t]
            for bar in bars0[si, t]:
                if bar["persistence"] <= thr[f"{name}:h0{suffix}~weak"]:
                    continue
                r_px, comp = _blob_geometry(f, bar)
                pr, pc = bar["peak"]
                if too_big(r_px) or min(pr, pc, ny - 1 - pr, nx - 1 - pc) == 0:
                    # peaks on the array border are kernel edge artefacts (a
                    # truncated bump); real units peak in the interior
                    continue
                features.append({"degree": 0, "type": t, "strong": bool(bar["persistence"] > thr[f"{name}:h0{suffix}"]), "scale": scales[si], "center": center(bar["peak"]),
                                 "peak": tuple(bar["peak"]), "outer": r_px * grid_step, "inner": 0.0,
                                 "foot": _filled(comp), "persistence": float(bar["persistence"]),
                                 "score": float(bar["persistence"] / max(thr[f"{name}:h0{suffix}"], 1e-300))})
            for bar in bars1[si, t]:
                if bar["persistence"] <= thr[f"{name}:h1{suffix}~weak"]:
                    break
                outer, basin = _ring_geometry(f, bar)
                if outer is None or too_big(outer):
                    continue
                rc = tuple(bar["basin_min"])
                ring_lab, _ = ndi.label(f >= bar["birth"], structure=np.ones((3, 3), bool))
                near = ring_lab[ndi.binary_dilation(basin)]
                foot = _filled(ring_lab == np.bincount(near[near > 0]).argmax())
                features.append({"degree": 1, "type": t, "strong": bool(bar["persistence"] > thr[f"{name}:h1{suffix}"]), "scale": scales[si], "center": center(rc), "peak": rc,
                                 "outer": outer * grid_step,
                                 "inner": math.sqrt(float(basin.sum()) / math.pi) * grid_step,
                                 "foot": foot, "persistence": float(bar["persistence"]),
                                 "score": float(bar["persistence"] / max(thr[f"{name}:h1{suffix}"], 1e-300))})

    coarse = [f for f in features if f["scale"] == max(scales) and f["strong"]]
    for ft in features:  # multiplicity control for fine-scale-only blobs
        if ft["degree"] == 0 and ft["scale"] != max(scales) and ft["strong"]:
            has_coarse = any(c["foot"][ft["peak"]] for c in coarse)
            key = f"{cell_types[ft['type']]}:h0@{ft['scale']:g}~max"
            if not has_coarse and (not FINE_ONLY_BLOBS_ALLOWED
                                   or ft["persistence"] <= thr.get(key, 0.0)
                                   or (d_nn is not None and ft["outer"] < UNIT_MIN_BLOB_RADIUS_NN * d_nn)):
                ft["strong"] = False
    units, consumed = [], set()
    for i, ft in enumerate(features):
        unit = {"kind": "blob", "center": ft["center"], "outer_radius": float(ft["outer"]),
                "inner_radius": float(ft["inner"]), "ring_type": cell_types[ft["type"]],
                "core_type": None, "persistence": ft["persistence"], "score": ft["score"],
                "scale": ft["scale"]}
        lum = [l for l in lumens if ft["foot"][l["peak"]] and (ft["strong"] or l["strong"])]
        core = [(j, c) for j, c in enumerate(features) if c["degree"] == 0 and j != i
                and c["type"] != ft["type"] and ft["foot"][c["peak"]] and c["outer"] < ft["outer"]
                and (ft["strong"] or c["strong"])]
        if lum:
            l = max(lum, key=lambda x: x["score"])
            unit.update(kind="ring_lumen", center=center(l["peak"]), inner_radius=l["radius"],
                        score=max(ft["score"], l["score"]))
        elif core:
            j, c = max(core, key=lambda jc: jc[1]["score"])
            unit.update(kind="ring_core", core_type=cell_types[c["type"]], center=c["center"],
                        inner_radius=float(ft["inner"] or c["outer"]))
            consumed.add(j)
        if unit["kind"] != "blob" or ft["strong"]:  # weak features need a partner
            units.append((i, unit))
    keep = [u for i, u in units if i not in consumed]

    units = sorted(keep, key=lambda u: (u["kind"] == "blob", -u["score"]))
    kept = []
    for u in units:  # drop duplicates and units inside an earlier-ranked unit
        if all(math.hypot(u["center"][0] - v["center"][0], u["center"][1] - v["center"][1])
               >= v["outer_radius"] for v in kept):
            kept.append(u)
    width, height = nx * grid_step, ny * grid_step
    for u in kept:
        x, y = u["center"]
        u["edge_touching"] = bool(min(x, y, width - x, height - y) < u["outer_radius"])
    if not kept:
        return {}
    centers = np.array([u["center"] for u in kept])
    sep = (float(cKDTree(centers).query(centers, k=2)[0][:, 1].min()) if len(kept) > 1 else None)
    kinds: Dict[str, int] = {}
    for u in kept:
        kinds[u["kind"]] = kinds.get(u["kind"], 0) + 1
    return {"units": kept, "n_units": len(kept), "kinds": kinds, "thresholds": thr,
            "min_center_separation": sep, "radii": [u["outer_radius"] for u in kept],
            "barcode_top": barcodes, "null_draws": len(null_draws),
            "quantile": float(quantile), "window_area_px": window_area}
