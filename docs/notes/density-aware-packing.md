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

## Known limits

- **Short-range structure.** kNN composition (1.45) and L at 20 µm and
  below (1.38) remain above the sampling floor (about 0.55 truth SDs). Local
  packing order and interfaces finer than 10 µm are not modeled.
- **Where that gap comes from.** Keeping each target's true positions and
  re-annealing all labels reproduces L at 20 µm and below at the floor (0.53),
  and five times more annealing iterations changed nothing on either geometry.
  The gap is therefore in the scaffold positions (local packing and density
  variation within 10–20 µm), not in the annealing energy or budget. A
  position-refinement stage that matches short-range spacing statistics is the
  next step if the ABM contact check needs it.
- **Benchmark targets.** Replicates of one sample cannot be closer to the
  process than the sample itself, which already has 8.3 of 45 statistics
  beyond 1 SD. State acceptance thresholds relative to that floor.
- **Tumor–immune contact.** Tumor cells with a CD8 T cell within 10 µm
  reach 0.64–0.93 of the target's fraction (uniform scaffold: 0.18–0.29).
  Check this against the t0 contact gate before the ABM re-run.
- **Mean degree.** Replicates reach about 0.9 of the region's; the paper
  gate is 0.95–1.05.
- **Local order.** Epithelial lattice order inside nests is not reproduced
  (g(r) is smooth).
- **Stationarity.** `resample` assumes a stationary region. Trends and
  window-sized patches are flagged; use `copy` for those regions.
- **Two dimensions.** Maps are 2D and z is uniform. For 2D sources use a
  thin slab so replicate graphs stay planar.
- **Fit cost.** A fit takes a few seconds: 19 null packings plus 80
  calibration simulations. Pass `n_null=0` or a `patch_prior` for speed.
