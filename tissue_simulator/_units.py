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

Units require a null: without null draws (``fit`` with ``n_null=0``) nothing is
detected.  Thresholds are rank-matched (see ``RANK_DEPTH``). The null is a packed uniform (CSR) layout, never the fitted patch
model (which already contains the structures under test); ``advisory_draws``
from the patch model only set ``ambiguous_with_patches``.  Gates: footprint
circularity ``UNIT_MIN_CIRCULARITY`` and size consistency ``UNIT_SIZE_FACTOR``.

Significance: for each (type, degree) and for the lumen map the threshold is
the ``quantile`` (``UNIT_NULL_QUANTILE``) of the maximum bar persistence over
the null draws (the null's essential H0 bar counts, so the comparison with the
data's own essential bar is like for like).  

With few null draws (``fit`` uses at most ``UNIT_NULL_DRAWS`` = 9) the 95th
percentile of the null max is effectively the null maximum, and about 5 % false
positives are expected per map, scale and degree.  Multiplicity control at the
fine scale: a 0.5x blob with no 1.0x footprint containing its peak must exceed
the null *maximum* and have outer radius >= ``UNIT_MIN_BLOB_RADIUS_NN`` nearest-
neighbour distances (``d_nn``); otherwise it is kept only when paired as a core
or lumen of a strong partner (escalated: see ``FINE_ONLY_BLOBS_ALLOWED``, default
no fine-only blobs).  Peaks on the array border are kept (the maps are edge-corrected and the null sees the same artefacts); units whose footprint touches the window are flagged ``edge_touching`` (censored).

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
UNIT_MIN_LUMEN_NN = 1.6
"""A paired lumen must have radius >= this many nearest-neighbour distances: the
weak (null-median) lumen threshold alone lets chance holes pair with a ring
(a germinal-centre core was reported as a lumen)."""
UNIT_CENTER_FRACTION = 0.4
"""A lumen or core must lie within this fraction of the outer radius of the
footprint centroid to be paired (a ring's hole is central; an off-centre hole
in a nest is not part of the unit)."""
ENRICH_FACTOR = 0.5
"""Proportion-aware enrichment rule: a blob's peak intensity must be at least
``1 + ENRICH_FACTOR * (1 - p_t)`` times the mean intensity of its type OUTSIDE
the blob footprint, with ``p_t`` the type's share of the cells. A rare type
(p_t -> 0) must reach 1.5x; an abundant type (p_t = 0.5) 1.25x; the bar is
lower for abundant types because a nest of a 50 % type can at most double the
local density. Referencing the outside (not the window mean) keeps the nests
themselves from raising the bar."""
UNIT_MATRIX_SHARE = 0.62
"""A type with at least this share of the cells is the *matrix* (stroma): its
bumps between units are the complement of the units, not units, and its
peak must reach ``UNIT_MATRIX_MIN_ENRICHMENT`` times the outside mean instead of
the proportion-aware bar (the matrix bumps peak at only ~1.2-1.4x). 0.62 sits between the demo stroma (0.66) and a 50/50 matrix with nests (0.58)."""
UNIT_MATRIX_MIN_ENRICHMENT = 1.5
"""Enrichment bar for matrix-type blobs (the former window-mean rule)."""
UNIT_MAX_BANDWIDTH_NN = 3.0
"""Cap, in nearest-neighbour distances, on the per-type bandwidth used for the
unit maps (applied by ``DensityModel._detect_units`` before the scales in
``UNIT_BANDWIDTH_SCALES``). The cross-validated bandwidth of an abundant type
can exceed a nest's radius and smooth the nests away."""
UNIT_MIN_CIRCULARITY = 0.6
"""Minimum footprint circularity ``area / (pi * r_max^2)``, with ``r_max`` the
largest distance from the footprint centroid to a footprint pixel. A disc scores ~1, strips,
bands and branching clusters of a Potts-like patch process score low. Data
driven (it needs no null); 0.6 keeps discs distorted by sampling noise."""
UNIT_SIZE_FACTOR = 1.7
"""Size consistency: units of one kind must have outer radii within this factor
of the kind's median, otherwise the outliers are dropped. With
``UNIT_SIZE_FILTER_MIN_UNITS`` or fewer units kept the median is not reliable
and outliers are flagged (``size_inconsistent``) instead of dropped."""
UNIT_SIZE_FILTER_MIN_UNITS = 3
"""At most this many kept units: size outliers are flagged, not filtered."""
UNIT_MIN_BLOB_RADIUS_NN = 2.0
"""Minimum outer radius, in nearest-neighbour distances, of a fine-scale-only blob."""
MIN_UNIT_RADIUS_PX = 1.5
"""Features smaller than this many pixels in radius are discarded."""
# Unpaired-blob rules (no core/lumen partner), applied after pairing:
# (1) inhomogeneous model: the blob must exceed the ``quantile`` of the maximum
#     persistence of the fitted patch-model null draws (``advisory_draws``) for
#     its key and be confirmed at the other scale; otherwise it is patch-ambiguous
#     and needs UNIT_AMBIGUOUS_MIN_BLOBS distinct confirmed same-type, same-scale
#     blobs (the patch model can make a bump or two, not a repeated set);
# (2) homogeneous model (no advisory draws; patch null == CSR null): the blob
#     must be supported at both scales, i.e. a footprint-overlapping strong H0
#     feature of the same type at the other scale. Family-wise control: ~12
#     tests (3 types x 2 scales x 2 degrees) per fit at the 95th percentile of
#     9 draws give a chance bump at one scale; a bump at both is much rarer.
UNIT_AMBIGUOUS_MIN_BLOBS = 3
# Patch-ambiguous unpaired blobs (inhomogeneous model; do not beat the patch-null
# maximum 95th percentile, i.e. rank 0 of the patch null) need at least this many supported strong blobs of
# the same type and scale (demo nests: 4; Potts J=1.5 chance clusters: fewer).
UNIT_SUPPORT_MIN_CIRCULARITY = 0.5
# Other-scale support must have footprint circularity >= this (looser than the
# 0.6 detection gate: a real nest's fine-scale footprint scored 0.60 and was lost
# by sampling noise, a chance bump scored 0.45).
# Paired units (ring_core, ring_lumen) keep the CSR-only rule.
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


RANK_DEPTH = 8
"""Number of ranked bars compared against the null (the i-th largest data bar
is compared with the null quantile of the i-th largest bar; bars beyond the depth are ignored. Comparing only against the null's largest bar (essential bar
included) is too conservative for the second and later units."""


def _top(bars, k=RANK_DEPTH):
    pers = [b["persistence"] for b in bars][:k]
    return pers + [0.0] * (k - len(pers))


def _rank_thresholds(draws, quantile):
    return np.quantile(np.array(draws, dtype=float), quantile, axis=0)


def _at(thresholds, i):
    return float(thresholds[min(i, len(thresholds) - 1)])


def _component_at(f, level, rc, structure=np.ones((3, 3), bool)):
    lab, _ = ndi.label(f >= level, structure=structure)
    k = lab[rc]
    return lab == k if k > 0 else np.zeros(f.shape, bool)


UNIT_EDGE_MIN_SCORE = 2.0
"""A window-cut (``edge_touching``) unit needs ``score = persistence / null
threshold`` of at least this: the edge-corrected kernel is noisiest at the
border, and Potts-like patch controls produced border bumps with scores
1.3-1.75, while three real nests cut by the edge score 2.3-10."""
UNIT_FOOTPRINT_LEVEL = 0.5
"""Fraction of a blob's persistence below its birth level at which its footprint is cut."""


def _blob_geometry(f, bar):
    """Radius (px) and footprint of an H0 blob: ``min(area_at_death, footprint area)``.

    The footprint is the component containing the peak at level
    ``birth - UNIT_FOOTPRINT_LEVEL * persistence``; capping with it keeps essential bars (whose
    death area is the whole window) and bars that merge late from absorbing
    their surroundings.
    """
    comp = _component_at(f, bar["birth"] - UNIT_FOOTPRINT_LEVEL * bar["persistence"],
                         tuple(bar["peak"]))
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


def _circularity(foot, peak) -> float:
    """Footprint circularity; a footprint touching the window is mirrored across
    each touched border first (a nest cut by the window is a disc, not a half disc)."""
    if foot.any():
        if foot[0].any():
            foot = np.concatenate([foot[::-1], foot], axis=0)
        if foot[-1].any():
            foot = np.concatenate([foot, foot[::-1]], axis=0)
        if foot[:, 0].any():
            foot = np.concatenate([foot[:, ::-1], foot], axis=1)
        if foot[:, -1].any():
            foot = np.concatenate([foot, foot[:, ::-1]], axis=1)
    rr, cc = np.nonzero(foot)
    if rr.size == 0:
        return 0.0
    # Distance from the footprint centroid (the peak of a noisy bump is off-centre
    # and would penalise real discs; deviation from "distance from the peak").
    rmax = float(np.hypot(rr - rr.mean(), cc - cc.mean()).max()) + 0.5
    return float(rr.size / (math.pi * rmax ** 2))


def _in_corner(foot) -> bool:
    """True when the footprint touches two perpendicular window borders.

    The edge-corrected kernel estimate at a corner rests on a quarter of the
    kernel's cells, and chance corner bumps of the (abundant) matrix type
    otherwise pass the null; corner-cut units are not reported (a known limit)."""
    return bool((foot[0].any() or foot[-1].any()) and (foot[:, 0].any() or foot[:, -1].any()))


def _filled(mask):
    return ndi.binary_fill_holes(mask)


def detect_units(lam_types, cell_types: Sequence[str], empty_space, grid_step: float,
                 mask, null_draws: Optional[List[Tuple[np.ndarray, np.ndarray]]] = None,
                 *, quantile: float = UNIT_NULL_QUANTILE, d_nn: Optional[float] = None,
                 advisory_draws: Optional[List[Tuple[np.ndarray, np.ndarray]]] = None, scales: Optional[Sequence[float]] = None
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
    if not null_draws:  # units are defined relative to a null: none, no units
        return {}
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

    thr, barcodes, rank = {}, {}, {}
    bars0, bars1 = {}, {}
    for si in range(n_scales):
        suffix = f"@{scales[si]:g}" if n_scales > 1 else ""
        for t, name in enumerate(cell_types):
            bars0[si, t] = superlevel_h0(lam_types[si, t], full)
            bars1[si, t] = superlevel_h1(lam_types[si, t], full)
            n0 = nulls(lambda ln, e: _top(superlevel_h0(ln[si, t], full)))
            n1 = nulls(lambda ln, e: _top(superlevel_h1(ln[si, t], full)))
            for key, vals in ((f"{name}:h0{suffix}", n0), (f"{name}:h1{suffix}", n1)):
                rank[key] = _rank_thresholds(vals, quantile)
                rank[key + "~weak"] = _rank_thresholds(vals, UNIT_WEAK_QUANTILE)
                thr[key] = float(rank[key][0])
                thr[key + "~max"] = float(max(v[0] for v in vals))
                thr[key + "~weak"] = float(rank[key + "~weak"][0])
            barcodes[f"{name}{suffix}"] = {"h0": barcode_summary(bars0[si, t], 5)["top"],
                                          "h1": barcode_summary(bars1[si, t], 5)["top"]}
    es = np.asarray(empty_space, dtype=float)
    lum_bars = superlevel_h0(es, full)
    nl = [_top(superlevel_h0(np.asarray(e, dtype=float), full)) for _, e in null_draws]
    rank["lumen"], rank["lumen~weak"] = (_rank_thresholds(nl, quantile),
                                         _rank_thresholds(nl, UNIT_WEAK_QUANTILE))
    thr["lumen"], thr["lumen~weak"] = float(rank["lumen"][0]), float(rank["lumen~weak"][0])

    means = lam_types.mean(axis=(2, 3))  # (n_scales, n_types)
    props = (means / np.maximum(means.sum(axis=1, keepdims=True), 1e-300))[0]

    def center(rc):
        return [float((rc[1] + 0.5) * grid_step), float((rc[0] + 0.5) * grid_step)]

    def too_big(r_px):
        return math.pi * r_px ** 2 > MAX_UNIT_AREA_FRACTION * window_area or r_px < MIN_UNIT_RADIUS_PX

    lumens = []
    for i, b in enumerate(lum_bars):
        peak = tuple(b["peak"])
        if (b["persistence"] > _at(rank["lumen~weak"], i) and not too_big(es[peak] / grid_step)
                and (d_nn is None or es[peak] >= UNIT_MIN_LUMEN_NN * d_nn)):
            lumens.append({"peak": peak, "radius": float(es[peak]),
                           "strong": bool(b["persistence"] > _at(rank["lumen"], i)),
                           "score": b["persistence"] / max(_at(rank["lumen"], i), 1e-300),
                           "persistence": float(b["persistence"])})

    features = []  # blobs and rings
    for si in range(n_scales):
        suffix = f"@{scales[si]:g}" if n_scales > 1 else ""
        for t, name in enumerate(cell_types):
            f = lam_types[si, t]
            for i, bar in enumerate(bars0[si, t][:RANK_DEPTH]):  # bars beyond the depth are not ranked
                if bar["persistence"] <= _at(rank[f"{name}:h0{suffix}~weak"], i):
                    continue
                r_px, comp = _blob_geometry(f, bar)
                pr, pc = bar["peak"]
                if too_big(r_px) or _in_corner(comp):
                    continue
                foot0 = _filled(comp)
                outside = f[~foot0]
                out_mean = float(outside.mean()) if outside.size else float(f.mean())
                need = (UNIT_MATRIX_MIN_ENRICHMENT if props[t] >= UNIT_MATRIX_SHARE
                        else 1.0 + ENRICH_FACTOR * (1.0 - props[t]))
                if (_circularity(foot0, bar["peak"]) < UNIT_MIN_CIRCULARITY
                        or f[tuple(bar["peak"])] < need * out_mean):
                    continue
                features.append({"degree": 0, "rank": i, "type": t, "strong": bool(bar["persistence"] > _at(rank[f"{name}:h0{suffix}"], i)), "scale": scales[si], "center": center(bar["peak"]),
                                 "peak": tuple(bar["peak"]), "outer": r_px * grid_step, "inner": 0.0,
                                 "foot": _filled(comp), "persistence": float(bar["persistence"]),
                                 "score": float(bar["persistence"] / max(_at(rank[f"{name}:h0{suffix}"], i), 1e-300))})
            for i, bar in enumerate(bars1[si, t][:RANK_DEPTH]):
                if bar["persistence"] <= _at(rank[f"{name}:h1{suffix}~weak"], i):
                    break
                outer, basin = _ring_geometry(f, bar)
                if outer is None or too_big(outer) or _in_corner(basin):
                    continue
                rc = tuple(bar["basin_min"])
                ring_lab, _ = ndi.label(f >= bar["birth"], structure=np.ones((3, 3), bool))
                near = ring_lab[ndi.binary_dilation(basin)]
                foot = _filled(ring_lab == np.bincount(near[near > 0]).argmax())
                if _circularity(foot, rc) < UNIT_MIN_CIRCULARITY:
                    continue
                features.append({"degree": 1, "rank": i, "type": t, "strong": bool(bar["persistence"] > _at(rank[f"{name}:h1{suffix}"], i)), "scale": scales[si], "center": center(rc), "peak": rc,
                                 "outer": outer * grid_step,
                                 "inner": math.sqrt(float(basin.sum()) / math.pi) * grid_step,
                                 "foot": foot, "persistence": float(bar["persistence"]),
                                 "score": float(bar["persistence"] / max(_at(rank[f"{name}:h1{suffix}"], i), 1e-300))})

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
                "scale": ft["scale"], "_ft": ft,
                "key": f"{cell_types[ft['type']]}:h{ft['degree']}"
                       + (f"@{ft['scale']:g}" if n_scales > 1 else "")}
        cy_, cx_ = (float(v) for v in ndi.center_of_mass(ft["foot"]))

        def central(peak):
            return math.hypot(peak[0] - cy_, peak[1] - cx_) * grid_step \
                <= UNIT_CENTER_FRACTION * ft["outer"]

        lum = [l for l in lumens if ft["foot"][l["peak"]] and central(l["peak"])
               and (ft["strong"] or l["strong"])]
        core = [(j, c) for j, c in enumerate(features) if c["degree"] == 0 and j != i
                and c["type"] != ft["type"] and ft["foot"][c["peak"]] and central(c["peak"]) and c["outer"] < ft["outer"]
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
    width, height = nx * grid_step, ny * grid_step
    for _, u in units:
        x, y = u["center"]
        u["edge_touching"] = bool(min(x, y, width - x, height - y) < u["outer_radius"])
    homogeneous = not advisory_draws
    patch_cache: dict = {}

    def patch_thr(key, deg, si, t):
        if key not in patch_cache:
            fn = superlevel_h0 if deg == 0 else superlevel_h1
            vals = [_top(fn(np.asarray(ln, float).reshape(lam_types.shape)[si, t], full))
                    for ln, _ in advisory_draws]
            patch_cache[key] = _rank_thresholds(vals, quantile)
        return patch_cache[key]

    def other_scale_support(ft):
        """A strong (rank-matched, above the CSR null) H0 bar of the same type at
        the other scale whose peak lies in the blob's footprint (no shape gates:
        the confirmation need only be a significant bump)."""
        name = cell_types[ft["type"]]
        for sj in range(n_scales):
            if scales[sj] == ft["scale"]:
                continue
            key = f"{name}:h0@{scales[sj]:g}"
            for j, bar in enumerate(bars0[sj, ft["type"]][:RANK_DEPTH]):
                if (bar["persistence"] > _at(rank[key], j) and ft["foot"][tuple(bar["peak"])]
                        and _circularity(_filled(_blob_geometry(lam_types[sj, ft["type"]], bar)[1]),
                                         bar["peak"]) >= UNIT_SUPPORT_MIN_CIRCULARITY):
                    return True
        return False

    def unpaired_ok(u):
        """Extra evidence for unpaired blobs (``UNIT_UNPAIRED_*`` rules)."""
        ft = u["_ft"]
        if homogeneous:  # patch null == CSR null: demand support at both scales
            return other_scale_support(ft)
        si = list(scales).index(ft["scale"])
        thr_p = patch_thr(u["key"], ft["degree"], si, ft["type"])
        if ft["persistence"] > _at(thr_p, 0) and other_scale_support(ft):
            return True  # beats the patch model's null (and is confirmed at both scales)
        # Patch-ambiguous: the fitted patch model reproduces a bump this strong.
        # A patch process can make one or two; accept only a repeated set of
        # supported same-key blobs (>= UNIT_AMBIGUOUS_MIN_BLOBS).
        same = [g for g in features if g["degree"] == ft["degree"] and g["type"] == ft["type"]
                and g["scale"] == ft["scale"] and g["strong"] and other_scale_support(g)]
        distinct = []
        for g in same:  # one bump seen at several ranks/footprints counts once
            if not any(h["foot"][g["peak"]] or g["foot"][h["peak"]] for h in distinct):
                distinct.append(g)
        return len(distinct) >= UNIT_AMBIGUOUS_MIN_BLOBS

    for _, u in units:
        if u["kind"] == "blob" and "_ft" in u:
            u["_ok"] = unpaired_ok(u)
    keep = [u for i, u in units if i not in consumed
            and (u["kind"] != "blob" or u.get("_ok", True))
            and (not u["edge_touching"] or u["score"] >= UNIT_EDGE_MIN_SCORE)]
    for _, u in units:
        u.pop("_ft", None)
        u.pop("_ok", None)

    units = sorted(keep, key=lambda u: (u["kind"] == "blob", -u["score"]))
    kept = []
    for u in units:  # drop duplicates and units inside an earlier-ranked unit
        if all(math.hypot(u["center"][0] - v["center"][0], u["center"][1] - v["center"][1])
               >= v["outer_radius"] for v in kept):
            kept.append(u)
    flag_only = len(kept) <= UNIT_SIZE_FILTER_MIN_UNITS
    for kind in {u["kind"] for u in kept}:  # size consistency within a kind
        radii = [u["outer_radius"] for u in kept if u["kind"] == kind]
        med = float(np.median(radii))
        ok = lambda u: med / UNIT_SIZE_FACTOR <= u["outer_radius"] <= med * UNIT_SIZE_FACTOR
        if flag_only:
            for u in kept:
                if u["kind"] == kind and not ok(u):
                    u["size_inconsistent"] = True
        else:
            kept = [u for u in kept if u["kind"] != kind or ok(u)]
    if not kept:
        return {}
    size_flags = [i for i, u in enumerate(kept) if u.get("size_inconsistent")]
    centers = np.array([u["center"] for u in kept])
    sep = (float(cKDTree(centers).query(centers, k=2)[0][:, 1].min()) if len(kept) > 1 else None)
    kinds: Dict[str, int] = {}
    for u in kept:
        kinds[u["kind"]] = kinds.get(u["kind"], 0) + 1
    ambiguous = None
    if advisory_draws:  # advisory: top unit vs the fitted patch model's null
        top = max(kept, key=lambda u: u["score"])
        name, rest = top["key"].split(":")
        deg, _, sc = rest.partition("@")
        si = list(scales).index(float(sc)) if sc else 0
        t = list(cell_types).index(name)
        fn = superlevel_h0 if deg == "h0" else superlevel_h1
        vals = [_max_persistence(fn(np.asarray(ln, float).reshape(lam_types.shape)[si, t], full))
                for ln, _ in advisory_draws]
        ambiguous = bool(top["persistence"] <= float(np.quantile(vals, quantile)))
    return {"units": kept, "size_inconsistent": size_flags, "ambiguous_with_patches": ambiguous, "n_units": len(kept), "kinds": kinds, "thresholds": thr,
            "min_center_separation": sep, "radii": [u["outer_radius"] for u in kept],
            "barcode_top": barcodes, "null_draws": len(null_draws),
            "quantile": float(quantile), "window_area_px": window_area}
