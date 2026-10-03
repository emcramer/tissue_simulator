# Density-aware scaffolds: design notes

Maintainer notes for `tissue_simulator.density`, `InhomogeneousPacker` and
density-aware replicates (`ReplicateGenerator(density_model=...)`).

## Why

Testing for the constrained-ensemble manuscript showed that the uniform
random-sequential-addition (RSA) scaffold erases tissue-scale heterogeneity.

- **Degree:** replicates of the 40 Keren et al. TNBC regions had 0.79 of the
  region's mean degree on the 20 µm graph.
- **Cell counts:** dense regions exceeded RSA saturation and came up short.
- **Larger scales:** cross-type L functions diverged beyond the graph radius.
- **ABM:** the replicate ensemble did not beat random placement as
  PhysiCell initial conditions.

The scaffold, not the labeling, was the limit: pair fractions at one radius
cannot put dense nests and sparse stroma back into an evenly packed box.

## Method

1. **Intensity maps.** Per-type Gaussian kernel intensity on a 5 µm grid,
   edge corrected, with an optional tissue mask. The bandwidth comes from
   likelihood cross-validation over 15–80 µm.
2. **Compartments.** k-means (K = 3) on per-type intensity vectors. Each
   compartment's density level and composition are counted from the cells
   inside it, which undoes kernel blurring across boundaries. Both are also
   profiled by distance to the compartment edge (0–10, 10–20, 20–40 and
   more than 40 µm), which keeps margins such as an immune cuff around nests.
3. **Heterogeneity test.** The region is compared against 19 uniform RSA
   packings with the same cell count and permuted labels. Homogeneous
   regions get uniform layouts, which reproduces the classic scaffold.
4. **Patch structure.** Two latent stationary Gaussian fields (Matérn,
   kernel-smoothed, correlated like the region's density and composition
   axes) are thresholded into compartments with exact area fractions. Their
   correlation length and smoothness are calibrated by simulation, so that
   the probability of two points at a given distance sharing a compartment
   matches the region's.
5. **Marks and hard core.** Radii are drawn from those observed at similar
   local density. The hard core is `kappa * (r_i + r_j)`, with `kappa` the
   2nd percentile of the region's nearest-neighbor `d / (r_i + r_j)`, so
   dense regions pack as tightly as the segmentation shows.
6. **Layouts.** `resample` simulates a new compartment map per replicate and
   rank-maps density within each compartment to the region's distribution.
   `copy` reuses the region's maps.
7. **Packing.** `InhomogeneousPacker` works in three stages:
   - bin quotas by stochastic rounding (largest-remainder ties biased quotas
     toward one corner),
   - per-bin RSA,
   - best-clearance insertion and capped relaxation where RSA saturates.

   Adaptive layouts add two things around those stages: each cell gets one
   radius from a source-matched deck before placement, and a final
   first-shell refinement moves cells slightly toward the source's spacing
   (thin slabs only). Both are described in
   [First-shell fidelity](#first-shell-fidelity).
8. **Labels.** The annealer matches pair fractions as before, plus the
   layout's expected composition in 40 µm bins. That term is weighted by
   `composition_weight * mean_degree**2` and warm-started from the layout.

## Decisions and rejected alternatives

| Choice | Rejected alternative | Why |
|---|---|---|
| Calibrate patches by simulation on same-compartment probability | Whittle / correlogram fit of the density field | Finite-window mean removal shrank fitted patches; resampled nests came out too small and too many |
| k-means on per-type intensities | Log density or composition-fraction features; normal scores | Boundaries moved into the sparse side (nest area 0.40–0.46 vs true 0.28); normal scores split at the median |
| Per-compartment density from cell counts | Smoothed intensity only | Kernel blur lowered nest density and mean degree |
| Likelihood CV bandwidth | Cronie–van Lieshout | Its criterion never reached 1 on the test patterns |
| Single-bin RSA null | Binned (stratified) null | Stratified quotas made the null too regular, so uniform regions tested heterogeneous |
| Weight scaled by squared mean degree, default 4 | Fixed weight | The useful weight ranged from about 16 (degree 5) to about 256 (degree 17) |
| Data-calibrated hard core plus capped relaxation | Matérn-II thinning; birth–death MCMC; PhysiCell relaxation | Lower maximum density than RSA; cost and determinism; slow and non-deterministic |
| Density and composition profiled by distance to the compartment edge | Constant composition per compartment; a hard-coded CD8 margin band | Constant composition lost margins; a generic profile keeps any interface without naming cell types |
| Candidates drawn from pixels in proportion to intensity | Uniform position within a bin | Respects gradients finer than a bin and never places cells in zero-density voids |
| 40 µm composition bins | 20 µm bins | 20 µm improved kNN composition but worsened L at 30–40 µm and the failure count |

## Validation

- **Legacy path:** `SpherePacker` output is bit-identical to the exhaustive
  check and about 16x faster at 850 cells. Legacy `graph_coloring` and
  `radius_tuning` replicates are hash-identical to the pre-change `main`.
- **Tests:** `tests/test_density.py`, `tests/test_packing.py` and
  `tests/test_density_replicates.py`. They cover homogeneity, compartment
  areas, exact resampled fractions, patch calibration, masks, trend flags,
  JSON round-trip, packer counts and bounds, incremental cost, determinism,
  parallel-equals-serial, `from_coordinates`, and MCP.

**Synthetic nest benchmark** (`examples/density_aware_replicates.py`):
20 truth samples, 3 targets × 10 replicates, 45 held-out statistics, |z|
against the truth ensemble.

| Row | mean \|z\| | \|z\| > 1 | kNN | L ≤ 20 µm | L 30–40 µm | L 60–80 µm | Degree ratio |
|---|---|---|---|---|---|---|---|
| Target itself (floor) | 0.59 | 8.3/45 | 0.55 | 0.54 | 0.57 | 0.69 | — |
| Uniform scaffold | 2.30 | 27.7/45 | 2.93 | 3.26 | 1.97 | 1.18 | 0.71 |
| Density-aware (constant composition per compartment) | 1.12 | 22.0/45 | 1.59 | 1.35 | 0.94 | 0.73 | 0.92 |
| Density-aware (boundary profiles, current) | 1.07 | 19.3/45 | 1.45 | 1.38 | 0.88 | 0.66 | 0.91 |
| Density-aware, 20 µm composition bins | 1.08 | 23.7/45 | 1.18 | 1.44 | 0.95 | 0.75 | 0.91 |

Replicate diversity was 0.89 of truth–truth, and positions shared within
3 µm were at chance.

**Tumor cells with a CD8 T cell within 10 µm**, as the replicate mean over
the target's value (6 replicates per target):

| Target (its contact) | Uniform | Density-aware, constant composition | Density-aware, boundary profiles |
|---|---|---|---|
| Nest sample 0 (0.21) | 0.27 | 0.73 | 0.76 |
| Nest sample 1 (0.19) | 0.29 | 0.93 | 0.88 |
| Nest sample 2 (0.24) | 0.23 | 0.70 | 0.75 |
| Dense lattice nests (0.38) | 0.18 | 0.64 | 0.64 |

Replicate SDs were 0.04–0.20 of the target value. Composition weights of 0,
1 and 4 changed these ratios by less than one SD.

**Composition-weight ablation** (3 replicates each; divergence / composition error):

| Weight (× mean degree²) | 0 | 0.25 | 1 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| Thinned nests (degree 5.3) | 0.094 / 0.19 | 0.089 / 0.07 | 0.080 / 0.06 | 0.068 / 0.05 | 0.079 / 0.05 | 0.119 / 0.05 |
| Dense nests (degree 17) | 0.239 / 0.12 | 0.237 / 0.09 | 0.229 / 0.08 | 0.199 / 0.05 | 0.189 / 0.04 | 0.163 / 0.04 |

## Review checkpoint

A mathematician, an ABM modeler/engineer and a pathologist reviewed the first
complete version. Each finding was checked against code and data:

| Finding | Decision |
|---|---|
| No margin ring; composition constant within a compartment | Accepted: boundary-distance profiles |
| Fallback insertion could land in zero-density voids | Accepted: pixel-weighted candidates |
| `get_cell_statistics` crashed for thickness 0 | Accepted: packing fraction is NaN |
| Null test cannot tell clustering from heterogeneity | Rejected: compartments are the intended model of a clustered layout, and short-range errors did halve |
| Patch smoothness and length not identifiable | Rejected: calibration matches the same-compartment profile, which the tests show is recovered |
| Relaxation scope too narrow in dense regions | Rejected: packer overlap (0.014) was below target and relaxation never ran; PhysiCell-scale overlap reflects local density and the missing lattice order (below) |

## Adaptive strategy (opt-in)

`strategy="adaptive"` (on `DensityModel.fit`, `ReplicateGenerator`,
`ReplicateGenerator.from_coordinates` and the MCP `setup_replicate_generator`
tool) turns on, together: data-derived bandwidth range, per-type bandwidths,
organization (trend) fitting with detrended residual maps, learned boundary
bands, multi-scale composition targets, a size-compatibility term and
fidelity diagnostics. `strategy="legacy"` is the default and draws no new
random numbers and takes no new arithmetic path; golden tests pin it. Switching
the default is a separate reviewed change. Individual pieces can be forced with
`fit(per_type_bandwidth=..., organization=..., bandwidth_range=...)`,
`composition_scales=`, `composition_bin=`, `size_weight=` and `diagnostics=`.

## Learned scales and estimation metadata

- Auto bandwidth range (`_auto_bandwidth_range`), with `d_nn` the median
  nearest-neighbor distance, `r_med` the median radius, `L` the shorter side
  of the mask bounding box: `lo = max(grid_step, d_nn, 2 r_med)`,
  `hi = max(min(L/3, 20 d_nn), 3 lo)`. Legacy keeps `(15, 80)`.
- Per-type bandwidths: a type keeps the pooled bandwidth
  (`"pooled_fallback"`) if it has fewer than `_TYPE_BW_MIN_CELLS = 10` cells or
  its cross-validated gain over the pooled bandwidth is below
  `_TYPE_BW_MIN_GAIN = 2.0` nats. Otherwise
  `log h_t = (n_t log h_cv + m log h_pool) / (n_t + m)` with
  `_TYPE_BW_SHRINK_M = 30`; evidence is `"cv"` when
  `n_t/(n_t+m) >= _TYPE_BW_CV_WEIGHT = 0.9`, else `"shrunk"`.
- Learned band edges: quantiles `_BAND_EDGE_QUANTILES = (0.25, 0.5, 0.75)` of
  the cells' distance to their compartment edge (legacy: `_BAND_EDGES_UM`).
- Size model: Gaussian in log radius per (type, density bin), shrunk toward the
  type's pooled Gaussian with `_SIZE_SHRINK_M = 20` pseudo-cells; types with
  fewer than `_SIZE_MIN_TYPE_CELLS = 10` cells use the all-type Gaussian
  (listed in `estimation["size_fallbacks"]`); `sigma >= _MIN_LOG_SIGMA = 0.02`.
  Fitted under both strategies (RNG-free); only the adaptive path uses it.
- Detrending regularization: inverse-trend weights clipped to
  `[1/_DETREND_WEIGHT_CAP, _DETREND_WEIGHT_CAP]` (5) and renormalized to mean 1;
  per-(compartment, band) compositions shrunk with
  `_COMPOSITION_PRIOR_CELLS = 10` pseudo-cells; maps floored at
  `_DETREND_FLOOR = 1e-3`.
- `DensityModel.estimation` keys: `strategy`, `detrend_weight_cap`,
  `composition_prior_cells`, `bandwidth`, `bandwidth_source`,
  `bandwidth_range`, `bandwidth_range_source`, `bandwidth_at_bound` /
  `bandwidth_bound` (also flag `bandwidth_at_{lower,upper}_bound`),
  `per_type_bandwidth`, `type_bandwidths`, `type_bandwidth_evidence`,
  `type_counts`, `size_fallbacks`, `learned_band_edges`, and (adaptive) `d_nn`.
- `Layout.quota_scale` (adaptive only): `max(bandwidth/2, 2 d_nn, 2 grid_step,
  sqrt(4/mean_density))` snapped to the grid; the packer uses it as quota bin
  side instead of the legacy `max(bandwidth/2, 10)`.

## Organization models (planar/radial)

Implemented in `tissue_simulator/_organization.py` (private; the module
docstring is the full method). Fitting is RNG-free; sampling happens in
`DensityModel.sample_layout`.

- Candidates: `density` (log-quadratic Poisson trend), `composition`
  (multinomial logit, linear), `both`; each for a planar coordinate
  `s = x cos(theta) + y sin(theta)` and a radial one `s = |p - c|`.
  Models: `none`, `planar_density|composition|both`,
  `radial_density|composition|both`.
- Selection by BIC: a candidate is accepted only if `delta_bic >=
  BIC_THRESHOLD = 10`. Radial must beat the best planar by the same margin;
  otherwise it is stored as `weak_candidate="radial"` and planar (or `none`)
  is used. No trend is attempted below `MIN_CELLS = 60`. Search: `N_ANGLE_STARTS
  = 12`, radial `CENTER_GRID = 5` per axis, `CENTER_MAX_ITER = 200`,
  `CELLS_PER_KNOT = 40`, `PROFILE_SMOOTHING = 0.5`, `MAX_INTEGRAL_PIXELS = 4000`.
- A model fitted on a region that the null test calls homogeneous is set to
  `none` with `fallback="homogeneous"`.
- Resampling: a new direction (planar) or center (radial, drawn from the window
  shrunk by `CENTER_MARGIN = 0.10` per side) is proposed up to
  `max_proposals` (`DEFAULT_MAX_PROPOSALS = 20`) times. A proposal is accepted
  if the window spans at least `MIN_COVERAGE = 0.9` of the source knot range
  and every type's expected proportion `p'_t` satisfies
  `|p'_t - p_t| <= max(PROPORTION_ABS_TOL = 0.02, PROPORTION_REL_TOL * p_t)`
  with `PROPORTION_REL_TOL = 0.15`. If none is accepted the best-coverage
  proposal is used, `Layout.organization["accepted"]` is False,
  `fallback == "best_of_proposals"`, the flag `organization_unsatisfied` is
  added and a warning is issued.
- `Layout.organization` keys: `model`, `geometry`, `direction`+`theta` or
  `center`, `proposals_tried`, `accepted`, `fallback`.

### Trend null

An adaptive fit that selects a trend (and is not homogeneous) tests it against
its own stationary residual model: `n_null` layouts (default 19) are drawn from
the final detrended model's `_residual_maps` (holes included), sampled as
inhomogeneous Poisson cell sets, and `trend_statistic` (max delta BIC over the
candidates) is recomputed on each. `p = (1 + #{null >= observed}) / (n_null + 1)`;
the trend is accepted iff `observed >= BIC_THRESHOLD` and `p <= alpha` (0.05).
With `n_null=19` the smallest attainable p is 0.05, so acceptance means no
null draw reaches the observed value. If every null statistic is non-finite
the null is reported as degenerate and the BIC threshold alone decides.
On rejection the model is rebuilt once without detrending and
`organization["fallback"] == "stationary_null"` (candidates are kept). The
result is in `organization["null"]`. A bootstrap from the fitted stationary
model was rejected: it reproduced concentric/layered trends at window scale, so
the null is deliberately the detrended one. It costs several seconds in
adaptive fits; set `n_null=0` to skip it (the BIC threshold alone then decides).

## Voids and lumens

`DensityModel.fit(..., voids="auto"|"none")` (default `"auto"` under adaptive,
`"none"` under legacy, where nothing changes and no random numbers are drawn).

- **Empty-space criterion** (`_voids.infer_voids`, run before anything uses the
  mask): a grid pixel is void when its distance to the nearest cell center
  exceeds `tau = max(sqrt(ln N / (pi lambda)), TAU_NN_FACTOR d_nn (1.5),
  TAU_RADIUS_FACTOR r_med + grid_step (2))` (largest Poisson gap, a multiple of
  the nearest-neighbour spacing, a multiple of the cell radius). The core is
  regrown by `tau`, opened with a 3x3 element (`OPENING_STRUCTURE`) and
  components smaller than `MIN_HOLE_AREA_FACTOR pi tau^2` (1) are dropped.
  Pixels holding a cell stay tissue.
- The inferred tissue mask becomes `model.mask` (edge correction, k-means pixels,
  per-tissue-area `density` and the organization fit all exclude the lumen).
  With a user `mask`, its interior holes are summarized instead (separation
  fallback: one median hole diameter) and the mask is unchanged.
- Only interior components (not touching the window edge) are removed from the
  tissue mask and re-placed in layouts; edge-touching components are reported
  but left in the mask, because they are never re-placed and removing them
  would inflate replicate cell counts. `DensityModel.voids` holds `n_holes`,
  `equivalent_diameters`, `areas`, `min_center_separation`, `tau` and
  `all_components`; `estimation["voids"]` carries the same summary without
  `all_components`.
- **Placement** (`_voids.place_voids`, in `_residual_maps` after the two noise
  fields; organized layouts draw holes after theta/center): the exact source
  interior hole count, diameters resampled with replacement, disc centers by
  rejection (<= 50 tries) keeping center separation >= the learned minimum
  (fallback `2 tau`). Inside holes compartment is -1 and intensity 0; bands are
  recomputed so they see the void edge; composition there is the global
  proportions. `n_target = round(density W H * tissue_fraction)`. Shortfalls are
  in `Layout.voids` and flagged `voids_unsatisfied`; replicates report
  `ReplicateStatistics.layout_voids`.
- **Radial anchoring**: for a radial trend whose center lies inside a source hole
  (centroid within one hole radius), the largest placed hole is put on the
  proposed center at every proposal (rings around a lumen). `combine_trend`
  multiplies, so zeros persist, and `expected_proportions` is intensity weighted,
  so void pixels do not affect acceptance.
- The packer is unchanged. Relaxation can push cells into a hole by at most
  `displacement_cap` (0.5 median radius); cells are also drawn from pixels next
  to the hole, adding up to half a pixel diagonal.

Known limits: the void edge and the compartment edge share one band axis
(band-edge conflation; `band_edge_type` would split them); holes are discs, not
ducts or elongated lumens (elongation is recorded, not used); inference needs
holes of at least ~`tau` radius.

## Persistent homology unit detection

Adaptive fits (`units="auto"`, the default; `"none"` disables) look for compact
*units*: tumor nests (blobs), follicles (a ring of one type around a core of
another) and glomeruli/arteries (a ring around an empty lumen). Detection is in
`_units.detect_units` on per-type intensity maps (`_persistence.superlevel_h0/h1`).
Rules, exactly as implemented:

- **Maps.** Edge-corrected kernel intensities at the per-type bandwidth capped
  at `UNIT_MAX_BANDWIDTH_NN = 3` nearest-neighbour distances, at the scales
  `UNIT_BANDWIDTH_SCALES = (1.0, 0.5)` of the capped value, plus the empty-space
  distance map. A planar trend is removed by inverse-trend weights; a radial
  trend is not (it is the unit itself). Nothing else is detrended.
- **Null.** `UNIT_NULL_DRAWS = 9` uniform (CSR) packings with the global type
  mix and the source cell count (no patch structure). `n_null=0` disables
  detection (`units == {}`, `estimation["units"] == {"skipped": "n_null=0"}`).
- **Significance.** Rank matched: the i-th largest bar (top `RANK_DEPTH = 8`)
  of a type/degree/scale must exceed the `UNIT_NULL_QUANTILE = 0.95` quantile
  (with 9 draws, about their maximum) of the null's i-th largest bar. Bars above
  the `UNIT_WEAK_QUANTILE = 0.5` quantile only are *weak* and kept only when
  paired with a strong partner.
- **Advisory comparison.** The same top-unit test against draws from the fitted
  patch model only sets `units["ambiguous_with_patches"]`; it never drops units.
- **Unpaired blobs** (kind `blob`; paired `ring_core`/`ring_lumen` keep the
  CSR-only rule). *Confirmation:* a strong (rank-matched, above the CSR null)
  H0 bar of the same type at the other scale, with its peak in the blob's
  footprint and footprint circularity at least
  `UNIT_SUPPORT_MIN_CIRCULARITY = 0.5` (other shape gates not applied). It is
  required in every case: with ~12 tests per fit (3 types x 2 scales x 2
  degrees) at the 95th percentile of 9 draws, a chance bump at one scale is
  expected; at both it is rare. *Homogeneous model:* the patch null equals the
  CSR null, so confirmation is the only extra rule. *Inhomogeneous model:* the
  blob must also exceed the 0.95 quantile of the patch-null maximum persistence
  (`advisory_draws`, its key); a patch-ambiguous blob is kept only when at least
  `UNIT_AMBIGUOUS_MIN_BLOBS = 3` distinct confirmed strong blobs of the same
  type and scale exist (a patch process can make a bump or two, tumor nests
  repeat). Strong Potts clustering (J=1.5) gives no units; random labels none.
- **Footprint.** Component of the map at `birth - 0.5 * persistence`
  (`UNIT_FOOTPRINT_LEVEL`). Area must be
  below `MAX_UNIT_AREA_FRACTION = 0.25` of the window and radius at least
  `MIN_UNIT_RADIUS_PX = 1.5`.
- **Circularity** `area / (pi r_max^2)` (r_max from the centroid; footprints cut
  by the window are mirrored across the border first) at least
  `UNIT_MIN_CIRCULARITY = 0.6`.
- **Enrichment** (blobs). Peak intensity at least
  `1 + ENRICH_FACTOR * (1 - p_t)` times the type's mean intensity outside the
  footprint, `ENRICH_FACTOR = 0.5`, `p_t` the type's share: 1.5x for a rare
  type, 1.25x for a 50 % type. A matrix type (share at least
  `UNIT_MATRIX_SHARE = 0.62`) needs `UNIT_MATRIX_MIN_ENRICHMENT = 1.5`, because
  its bumps between units are the complement of the units.
- **Pairing.** A lumen (empty-space peak, radius at least `UNIT_MIN_LUMEN_NN =
  1.6` nearest-neighbour distances) or another type's blob inside the footprint
  and within `UNIT_CENTER_FRACTION = 0.4` of the outer radius of its centroid
  makes `ring_lumen` / `ring_core`; otherwise `blob`.
- **No fine-only blobs.** A blob found only at 0.5x with no strong 1.0x
  footprint around its peak is reported only as a core or lumen of a strong
  partner (`FINE_ONLY_BLOBS_ALLOWED = False`).
- **Duplicates.** Units whose center lies in an earlier-ranked unit's outer
  radius are dropped.
- **Size.** Radii of one kind within `UNIT_SIZE_FACTOR = 1.7` of the kind's
  median; with `UNIT_SIZE_FILTER_MIN_UNITS = 3` or fewer units outliers are
  flagged (`size_inconsistent`, a list of unit indices) instead of dropped.
- **Edges.** Window-cut units are kept and flagged `edge_touching` when the
  outer radius reaches the border (censored: counted, to be left out of the
  radius pool). They need a score (persistence / null threshold) of at least
  `UNIT_EDGE_MIN_SCORE = 2.0`, since border bumps are the noisiest. Footprints touching two perpendicular borders (corners) are
  discarded.

Kinds: `blob`, `ring_core`, `ring_lumen`. `model.units` holds the units
(`center`, `outer_radius`, `inner_radius`, `ring_type`, `core_type`,
`persistence`, `score`, `scale`, `edge_touching`), `kinds`, `thresholds`,
`min_center_separation`, `size_inconsistent`, `ambiguous_with_patches` and,
after the profile fit, `profiles` and `d_nn`.

## Germ-grain layouts

`_grains.fit_unit_profiles` pools, per kind, the density relative to the window
mean and the type composition against the normalized distance `s = r / R_outer`
on `DEFAULT_N_KNOTS = 8` equal bins of `[0, S_MAX = 1.5]` (`PROFILE_PRIOR_CELLS = 5`
Dirichlet pseudo-cells toward the global mix; `DENSITY_FLOOR = 0.05` for empty
non-lumen knots; lumen knots have density 0). Fitting is RNG-free and runs in `fit`
after detection; profiles are stored in `model.units["profiles"]` (JSON-friendly;
`from_dict` tolerates their absence).

`sample_layout` (adaptive, units present, `layout="resample"`, homogeneous models
included) draws every existing field first (noise fields; theta/center for organized
layouts; the legacy-compatible voids) and then `_grains.place_grains` with the exact
source counts per kind, radii resampled with replacement and hard-core centers at
`model.units["min_center_separation"]` (`UNIT_PLACEMENT_TRIES = 500` draws per
grain). `rasterize_grains` yields an intensity factor, grain composition and a lumen
mask; intensity = residual x factor (lumen -> 0, compartment -1), composition = grain
composition inside grains and the residual elsewhere, `n_target =
round(density W H tissue_fraction)` and the intensity is renormalized to it.
`Layout.units` holds the placed dict; `"units_unsatisfied"` is flagged on shortfall.
`rasterize_grains` ignores the inner radius for non-lumen kinds (a core is not a
void); only `ring_lumen` grains get zero intensity inside it.

Precedence (constants `UNIT_OWNS_LUMEN`, `UNIT_OWNS_RADIAL` in `density.py`):

- Units own their lumens: a void whose centroid lies inside a unit's outer radius is
  dropped from the placement pool (a filtered copy; `model.voids` and its diagnostics
  are untouched), so holes are not placed twice.
- A radial trend whose center lies inside a unit's outer radius is skipped for
  sampling (`Layout.organization = {"model": "none", "fallback": "units"}`; the
  organization stays recorded on the model). Planar trends coexist with units: grains
  are applied after the trend.

Replicate diagnostics: `ReplicateStatistics.layout_units` and
`fidelity["persistence_distance"]` (per type, Wasserstein-1 between H0 barcodes of the
source and replicate KDE maps at the model's type bandwidths, essential bar dropped).

Known limits of units: discs only (no elongated or irregular units); fitted radii are
biased about 20 % low (the smoothed footprint is smaller than the true edge), so
placed grains are slightly small; units touching the window edge are detected but
censored (`edge_touching`) and re-placed as full discs inside the window; the null
costs `UNIT_NULL_DRAWS` packings per fit (seconds, scales with cell count); and the
source's smallest center separation can be infeasible for several large units in the
same window (for example four r = 38 blobs in 300 um), so grains may fall short and
`units_unsatisfied` is flagged.

- **Small cores.** A core of only a few cells (about five per unit in the follicle
  fixture) is detected as part of its ring but reproduced only as an enrichment:
  the per-kind profile is shrunk toward the global proportions by
  `PROFILE_PRIOR_CELLS` pseudo-cells and the annealer's composition bins are
  coarser than the core, so the replicate core is C-enriched (about 4x its
  global share) rather than C-dominant.

## Composition constraints and size compatibility

- Multi-scale composition (`ReplicateGenerator._resolve_composition_scales`).
  `composition_scales=[...]` wins; else numeric `composition_bin` is one scale;
  `"auto"` (adaptive default) uses `[max(bw/2, 2 d_nn, MIN_AUTO_SCALE_UM),
  max(patch_length, 2 first_band_edge), min(W, H)/4]`, deduplicated so
  successive scales differ by `MIN_SCALE_RATIO = 1.5`, at most
  `MAX_COMPOSITION_SCALES = 3`, `MIN_AUTO_SCALE_UM = 10`. Expected bin counts
  are the layout composition (Sinkhorn-rescaled to node counts), unshrunk.
  Passed to the annealer as `spatial_composition_scales`; each scale's weight is
  calibrated so the summed spatial term equals `composition_weight` (default
  `ADAPTIVE_COMPOSITION_WEIGHT = 1.0`; legacy keeps 4.0 on its own scale)
  times the edge-count SSE of a shuffled warm-start labeling (recorded in
  `fidelity['composition_weight_effective']` / `['composition_calibration']`).
  Unit profiles use equal-count radial bins (`PROFILE_MIN_CELLS = 10`) and a
  neighbour-bin prior (`PROFILE_PRIOR_CELLS = 2`).
- Size compatibility: the annealer's `size_compatibility` target holds
  `nll[node][color] = 0.5 ((log r - mu)/sigma)^2 + log sigma` for the cell's
  local density bin, weighted `size_weight * mean_degree` (default
  `size_weight=1.0` adaptive, 0 legacy). Radii are therefore a constraint on
  type assignment, not re-drawn.
- `cost_terms` gains `spatial_scales` and `size` entries when those targets exist.
- Requested vs achieved counts: `requested_cell_type_counts` (layout quota) and
  `achieved_cell_type_counts` are reported so rare-type shortfalls are visible.

## Diagnostics

Adaptive generators default to `diagnostics=True`; legacy to off. Per replicate
`ReplicateStatistics` gains `requested_cell_type_counts`,
`achieved_cell_type_counts`, `layout_organization`, `separation` (clearance
and normalized nearest-neighbor distance quantiles, plus count) and `fidelity`
(`fidelity_diagnostics`): `size_nll` (replicate vs source), `size_ks_by_type`,
`nn_distance_quantiles` (p5/p50/p95; `NN_QUANTILES`), `mixing_index`,
`organization_rmse`, `interface_fraction`, `n_components`, `persistence_distance`, `n_holes` (bins whose
layout-expected count exceeds `HOLE_MIN_EXPECTED = 2` but hold no cell, at the
finest scale). Entries that cannot be computed are None. Source reference
values are cached by `ReplicateGenerator._cache_source_reference(tissue)`
(called by `from_coordinates` and the MCP tool), so a generator built
directly from a model has no source KS/NLL.
`PackingReport` adds `bin_shortfall`, `quota_floor`, `quota_floor_source`,
`clearance_quantiles`, `normalized_distance_quantiles`,
`dense_bin_fraction_short` (bins with target >= 4 achieved below 0.9 of target),
`radius_assignment` (`"deck"` or `"per_candidate"`), `refinement` (None when
refinement did not run; otherwise `sweeps`, `energy_before`, `energy_after`,
`n_moved`, `mean_displacement`, `max_displacement`, `accept_rate`, `s_target`,
`seconds`) and `first_shell`. `fidelity["first_shell"]` holds the same
comparison per replicate (see [First-shell fidelity](#first-shell-fidelity)).
`layout_organization` records `failed` (`"coverage"` or `"composition"`) and
`best_coverage` when the organization fallback was used.

## Schema evolution

- `DensityModel.to_dict` writes `format_version: 2` (adds `strategy`,
  `type_bandwidths`, size model, `organization`, `estimation`).
  `from_dict` skips absent keys, so version-1 files load with legacy defaults
  (`strategy="legacy"`, `organization={"model": "none"}`) and behave as before.
- Each replicate draws layout, packing, warm start and annealing from separate
  streams spawned from `SeedSequence([seed, replicate_id])`, so adaptive
  replicates are identical serial and parallel. Legacy seeding is unchanged.
- Later additions (unreleased): `DensityModel.shell` and `Layout.shell`
  (default `{}`), and `PackingReport.radius_assignment`, `refinement` and
  `first_shell` (default None). The shell profile is measured at fit time for
  both strategies. Models saved before it existed load with an empty `shell`;
  for those, adaptive packs use the radius deck but skip refinement, because
  refinement needs the profile. Refit the model to get refinement.
- All new dataclass fields have defaults and come last; the generator stores
  only plain data, so it pickles for process pools.

## Known limits

- **Short-range structure.** On the synthetic benchmark above (legacy
  strategy), kNN composition (1.45) and L at 20 µm and below (1.38) remain
  above the sampling floor (about 0.55 truth SDs). Interfaces finer than
  10 µm are not modeled. Adaptive packs of thin slabs now match the source's
  size-normalised first-shell spacing (see
  [First-shell fidelity](#first-shell-fidelity)); that benchmark has not been
  re-run with it.
- **Where that gap comes from.** Keeping each target's true positions and
  re-annealing all labels reproduces L at 20 µm and below at the floor (0.53),
  and five times more annealing iterations changed nothing on either geometry.
  That second check was not informative: with the default cooling rate the
  schedule ended after about 1,400 swaps regardless of `max_iterations` (see
  [Annealing schedule](#annealing-schedule)). The first check still places
  the gap in the scaffold positions (local packing and density
  variation within 10–20 µm). The
  first-shell refinement is the position-refinement stage that note called
  for, but it is only applied to adaptive thin slabs.
- **Benchmark targets.** Replicates of one sample cannot be closer to the
  process than the sample itself, which already has 8.3 of 45 statistics
  beyond 1 SD. State acceptance thresholds relative to that floor.
- **Tumor–immune contact.** Tumor cells with a CD8 T cell within 10 µm
  reach 0.64–0.93 of the target's fraction (uniform scaffold: 0.18–0.29).
  Check this against the t0 contact gate before the ABM re-run.
- **Mean degree.** On the legacy strategy replicates reach about 0.9 of the
  region's mean degree on a 20 µm radius graph; the paper gate is 0.95–1.05.
  On a size-aware (mechanical, factor 1.5) graph the v0.1.18 adaptive packs
  were lower, 0.79–0.81 of the source. Adaptive thin-slab packs on this branch
  measure 0.985–1.029 at packer level on three regions; see
  [First-shell fidelity](#first-shell-fidelity) for the conditions and limits.
- **Local order.** Epithelial lattice order inside nests is not reproduced
  (g(r) is smooth).
- **Stationarity.** `resample` assumes a stationary region. Trends and
  window-sized patches are flagged; use `copy` for those regions.
- **Two dimensions.** Maps are 2D and z is uniform. For 2D sources use a
  thin slab so replicate graphs stay planar.
- **Fit cost.** A fit takes a few seconds: 19 null packings plus 80
  calibration simulations. Pass `n_null=0` or a `patch_prior` for speed.
- **Symmetric multi-nest samples.** A pattern of several similar nests can be
  selected as radial (delta BIC above `BIC_THRESHOLD = 10`) although it is not
  one center-to-periphery structure; inspect `organization["candidates"]`.
- **Auto bandwidth bound.** On near-uniform data the cross-validated bandwidth
  often sits at the upper bound of the auto range (flag
  `bandwidth_at_upper_bound`); treat the density map as nearly flat there.
- **Fidelity summaries are evidence, not calibrated intervals.** `fidelity`
  values have no null distribution; compare against replicates of other
  seeds or the source's own floor.
- **Radius definition.** Size targets reproduce whatever the source radii mean.
  State whether they are whole-cell or nuclear radii; mixing them between
  source and `cell_radii` changes packing density.
- **Radial coverage ignores masks.** The radial acceptance check uses the
  replicate window's rectangle, not a mask on the `Layout`.
- **Small windows.** When no proposal meets the coverage and proportion
  criteria the organization falls back to `best_of_proposals`
  (flag `organization_unsatisfied`); the replicate is still produced. The
  warning now says which criterion failed, and `layout.organization` records
  `failed` (`"coverage"` or `"composition"`) and `best_coverage`. Coverage is
  the fraction of the source's trend range that the window spans; a proposal
  needs at least 90%. A window much smaller than the source cannot reach that.
  On one fitted Keren field with a radial trend, the fallback fired in 1 of 300
  layouts at the source window size (798.72 µm) and the result was only
  marginally outside tolerance. It fired in 100% of layouts at a 400 µm window
  and in 74% at 600 µm, and there the fallback composition was materially wrong
  (one type about 24 percentage points over its proportion). **When the model
  has a trend, generate replicates in a window at least as large as the source
  region.** One regional measurement; it was not repeated on other datasets.
- **Unit detection power.** Persistence against a 9-draw CSR null is
  conservative: nests of an abundant type (90 % pure nests of r = 38 µm in a
  50/50 or 65/35 matrix) are found 0-3 times of 4 across seeds, and small
  window-cut units next to a corner are discarded. Nests of a rare type, follicles,
  glomeruli and arteries are found reliably. A count-based evidence stage
  would be the next step.
- **Size filter with few units.** With three units or fewer, size outliers are
  flagged, not dropped; a spurious unit can survive there.

## First-shell fidelity

### Problem

Replicates of real tissue were scored on a 20 µm radius graph, which hid a
short-range deficit. On three Keren et al. MIBI fields, v0.1.18 adaptive packs
had about 20% fewer first-shell neighbours than the source on the mechanical
graph (factor 1.5; replicate/source 0.78–0.81), about 25% fewer pairs closer
than `r_i + r_j`, a median radius 3–9% small and a disc area fraction 6–14%
low. A likely reason is that segmented tissue is often confluent, with
neighbouring discs touching or overlapping, while a random-sequential-addition
packer with a hard core leaves more space between cells.

All distances below are size-normalised: `s = d / (r_i + r_j)` for a pair at
centre distance `d`. Touching discs have `s = 1`.

### Mechanism 1: measure the source's shell at fit time

`DensityModel.fit` (both strategies) measures the source's pair profile with
`tissue_simulator/_shell.py` and stores it as `DensityModel.shell`, which is
copied to every `Layout.shell`. The measurement is deterministic.

- `edges`, `pairs_per_cell`: bins of 0.05 in `s` up to 2.2, and the mean
  number of neighbours per cell in each bin (no edge correction, so it is
  comparable to a replicate in the same window).
- `g`: an edge-corrected pair correlation in the same bins, normalised by the
  annulus area of an uncorrelated pattern with the same radii (about 1 for
  uncorrelated cells).
- `edge`: the shell edge, the first minimum after the first peak of `g`
  smoothed with a Gaussian of 0.1 in `s`. The minimum must be at `s >= 1` and
  at least 5% below the peak. It is None when there is no clear shell, as for
  Poisson-like or random-sequential-addition-like sources.
- `summary`: mean mechanical degree at factor 1.5, overlapping pairs per cell
  (`s < 1`), median radius and disc area fraction (sum of `pi r^2` over window
  area).
- `s_floor`: the hard core `kappa`.

`interaction_factor="auto"` uses `edge` as the mechanical graph's factor (see
[Mechanical neighbour graph](#mechanical-neighbour-graph)).

### Mechanism 2: a radius deck

Previously the packer drew a new radius for every candidate position, and a
small radius is accepted more often than a large one, so placed cells were
biased small (median radius 3–9% low). For adaptive layouts each cell (a
"ticket" in the bin quota) now gets one radius before placement. The pool of
source radii is sampled at stratified quantiles, so the deck reproduces the
source's size distribution, and the radii are handed out in the rank order of
a provisional density-conditioned draw, so dense areas keep their smaller
cells. Retries change only the position, never the radius.

`PackingReport.radius_assignment` is `"deck"` or `"per_candidate"`. Legacy
layouts and the fit-time null packings keep the per-candidate draw, so legacy
output and adaptive fits are unchanged. The deck matches the pooled size
distribution; per-type size distributions still depend on the annealer's size
term.

### Mechanism 3: first-shell refinement

After relaxation, an adaptive pack is refined by greedy single-cell moves.
Each sweep, every cell proposes one Gaussian step (sigma 0.25 median radius,
cut at 3 sigma). A move is kept only if it lowers

`E = sum_k (h_k - t_k)^2 / (t_k + 5)`

where `h_k` is the replicate's pair count in bin `k` of `s` (up to
`s = 2.0`) and `t_k` is the source's `pairs_per_cell` times `n / 2`. A cell
must also:

- stay in its quota bin and off zero-intensity pixels, so the layout's density
  map and bin quotas are preserved (quotas were identical with and without
  refinement in the tests below);
- end within `refine_cap` of its position before refinement;
- not create a pair below the hard core `kappa` (a pair already below `kappa`
  may only move apart). `z` never changes.

Sweeps stop after `refine_sweeps` or when a sweep lowers `E` by less than
0.1%. Cost measured on 3,000–6,000 cells: 0.3–0.9 s per pack.

Refinement runs only when all of these hold: the layout is adaptive and has a
shell profile, the pack has at least 50 cells, `refine_shell` is not False, and
the slab is thin (thickness at most twice the median radius), because the
profile is a 2-D measurement. Parameters go through `packing_params`:

| Key | Default | Meaning |
|---|---|---|
| `refine_shell` | None (on for adaptive) | False turns refinement off |
| `refine_sweeps` | 20 | Maximum sweeps |
| `refine_cap` | median radius | Largest distance in µm from the pre-refinement position |
| `radius_assignment` | `"deck"` adaptive, else `"per_candidate"` | Radius draw |

`PackingReport.refinement` records `sweeps`, `energy_before`, `energy_after`,
`n_moved`, `mean_displacement`, `max_displacement`, `accept_rate`, `s_target`
and `seconds`. The other report fields describe the refined positions.

### Reporting

`PackingReport.first_shell` and `ReplicateStatistics.fidelity["first_shell"]`
compare each replicate with the source:

```python
{"factor": 1.5, "replicate": summary, "source": summary,
 "ratios": {"mean_degree": ..., "overlap_pairs_per_cell": ...,
            "median_radius": ..., "area_fraction": ...},
 "same_window": True}
```

Each ratio is replicate over source. `same_window` is true when the replicate
window has the source's size, which the comparison assumes. An adaptive
replicate in the same window whose mean-degree ratio is below 0.9 raises a
warning. The MCP replicate summary carries the ratios as
`first_shell_ratios`.

### Results

Packer level, replicate/source, three Keren regions with 3 packs each.
"Before" is v0.1.18 and "after" is this change.

| Quantity | Before | After |
|---|---|---|
| Mean mechanical degree (factor 1.5) | 0.79–0.81 | 0.985–1.029 |
| Pairs closer than `r_i + r_j`, per cell | 0.69–0.77 | 0.971–1.079 |
| Median radius | 0.91–0.96 | 1.000 |
| Disc area fraction | 0.86–0.94 | 1.000 |
| Mean degree, 20 µm radius graph | | 0.93–1.08 |

The 20 µm graph value changed little with refinement. Bin quotas were
identical with and without refinement.

End to end through `ReplicateGenerator` (adaptive strategy, mechanical graph
at factor 1.5, `layout="resample"`, package defaults, three replicates per
region; mean, replicate over source, v0.1.18 → this change):

| Quantity | Sample 1 | Sample 3 | Sample 15 |
|---|---|---|---|
| Mean mechanical degree | 0.795 → 1.017 | 0.801 → 0.986 | 0.806 → 0.988 |
| Pairs closer than `r_i + r_j`, per cell | 0.773 → 1.050 | 0.762 → 0.983 | 0.721 → 0.978 |
| Median radius | 0.907 → 1.000 | 0.950 → 1.000 | 0.959 → 1.000 |
| Disc area fraction | 0.858 → 1.000 | 0.920 → 1.000 | 0.944 → 1.000 |
| Mean degree, 20 µm radius graph | 1.045 → 1.042 | 0.928 → 0.939 | 0.919 → 0.942 |
| Composition error | 0.360 → 0.288 | 0.385 → 0.290 | 0.079 → 0.076 |
| Divergence score | 0.313 → 0.316 | 0.196 → 0.160 | 0.785 → 0.769 |
| Seconds per replicate | 17.6 → 20.5 | 25.0 → 28.8 | 7.5 → 8.2 |

Cell counts equal the source in every replicate. One sample 1 replicate had
an overlapping-pairs ratio of 1.110, just above the 0.90–1.10 range the
requesting team proposed. The divergence score of sample 15 (95% tumour, with
type pairs of zero or one edge) varies by about 0.1 between replicates, so its
change is within noise.

These figures come from three regions of one dataset.

### Annealing schedule

Correcting the geometry exposed a second problem. The annealer's default
schedule cools from 100 to 0.1 at a rate of 0.995 per step, which reaches the
final temperature after about 1,400 swaps whatever `max_iterations` is. On a
graph of several thousand cells the labelling was therefore close to its warm
start. With v0.1.18 geometry the replicate had about 20% too few edges, and
that deficit happened to cancel the warm start's excess of contacts between
unlike types, so the divergence score looked better than the labelling
deserved. With the correct number of edges the excess showed: on sample 3
the divergence rose from 0.196 to 0.293.

The adaptive strategy now derives the cooling rate from the step budget, so
the final temperature is reached on the last step. The budget is the larger
of 20,000 and ten swaps per graph node
(`ADAPTIVE_ANNEAL_STEPS_PER_NODE`). On sample 3 (new geometry) the divergence
was 0.214 at 20,000 steps, 0.162 at 60,000 and 0.147 at 200,000. An explicit
`cooling_rate` in `coloring_params` turns this off, an explicit
`max_iterations` sets the budget, and the legacy strategy keeps the old
schedule. The step count used is in `fidelity["anneal_steps"]`.

### Limits

- Refinement applies to thin slabs only. 3-D packs get the radius deck but are
  not refined.
- The target is the source's per-cell pair rate, so it assumes the replicate
  window has the source's cell density (resample or copy layouts of the same
  window). A different density changes the expected counts per cell.
- One of nine test packs still had an overlap fraction above the source's
  (0.047 against 0.020) because relaxation reached its displacement cap.
- The deck matches the pooled size distribution, not per-type distributions.
- A packing method built for confluent tissue (collective compression, or a
  weighted Voronoi/Laguerre layout) was considered and deferred, because
  refinement met the targets on these regions. It may be needed for tissue
  that is more confluent than these three.
- There is no real-tissue validation beyond these three regions yet.

## Mechanical neighbour graph

Replicates are scored on a neighbour graph, and the contact rule (1.01 x summed
radii) is nearly empty for packings with spacing. The `"mechanical"` mode uses
`d <= 1.5 * (r_i + r_j)` (PhysiCell's default mechanics interaction distance),
so neighbourhoods scale with cell size, unlike a fixed radius. It is applied
identically to source targets and replicates (`TargetStatistics.network_rule`
is checked, with a warning on mismatch) and is the `from_coordinates` default.

With `interaction_factor="auto"` the factor is the source's shell edge, with
1.5 and a warning when there is none. `network_rule["interaction_factor_source"]`
is `"fixed"`, `"auto"` or `"auto_fallback"`, and replicates reuse the source's
number. See [First-shell fidelity](#first-shell-fidelity).
