# GPU acceleration: implementation specification

Status: proposed. No package changes are made by this document. It is a hand-off
to the team that operates the GPU server; the authors cannot run code there.

All speedups and memory figures below are **estimates extrapolated from a laptop
CPU profile**. Nothing has been run on a GPU. Section 11 lists what to measure
and send back.

## 1. Goal, scope and non-goals

Make `tissue_simulator` usable at two scales it cannot reach today:

- **Big samples**: one source region with 50,000 to 1,000,000 cells.
- **Big ensembles**: hundreds to thousands of replicates of one source.

Acceleration is **opt-in**. A plain `pip install tissue_simulator` must behave
exactly as version 0.1.18 does: same results, same random streams, no GPU
library imported.

Non-goals:

- `SpherePacker` and its golden checksums in `tests/test_packing.py` are out of
  scope and must not change.
- The `legacy` strategy and any seeded run on the NumPy backend stay
  byte-identical.
- No change to the public return types (`TissueSection`, `networkx.Graph`,
  `DensityModel`, `ReplicateStatistics`).

## 2. Measured baseline

Version 0.1.18, `strategy="adaptive"`, `network_mode="mechanical"`, one laptop
CPU core, synthetic tissue with four nests, radii 3.5 to 5.0 µm.

| Stage | 694 cells (300 µm window) | 2,702 cells (600 µm window) |
|---|---|---|
| Neighbour-graph build | 0.4 s | 6 to 12 s |
| Adaptive `DensityModel` fit | 10 s | 44 s |
| One replicate, 20,000 annealing steps | 1.6 s | 13.7 s |

Where the time goes:

- **Graph build** is quadratic. At 2,702 cells it is 12.1 s of the 13.7 s
  replicate.
- **Fit**: about 75 % is `InhomogeneousPacker.pack`, called 28 times in a row.
  Nineteen calls are the heterogeneity null (`density.py:1427`) and nine are the
  unit-detection null (`density.py:1335`). Inside `pack`, the cost is
  `SpatialHashGrid.neighbors` (`packing.py:72`).
- **Annealing** is about 0.5 s. This is misleading; see below.

Three things the profile hides:

1. The annealing budget is fixed at 20,000 steps
   (`replicate_generator.py:507`). At one million cells that touches at most
   4 % of the nodes, so a larger sample needs more steps for quality, not only
   speed.
2. `GraphColorizer.colorize` copies the whole labelling on every improving step.
3. One `Cell` object per cell and a `networkx.Graph` with a dict per node cost
   several GB and minutes at one million cells, before any arithmetic runs.

**The central finding**: the two largest costs are Python loops, not array
arithmetic. A GPU cannot accelerate them as written. They must first become
array operations, and that CPU rewrite alone is worth several orders of
magnitude on the graph. Work package P0 is therefore a prerequisite for
everything else, and it benefits users with no GPU.

To reproduce the baseline, time these three calls on a tissue built with
`SpherePacker((size, size, 1.0), {"c": (3.5, 5.0)}, min_spacing=0.3, seed=0)`:
`SpatialNetworkAnalyzer().build_network_from_tissue(t, mode="mechanical")`,
`DensityModel.from_tissue(t, strategy="adaptive", seed=0)`, and
`ReplicateGenerator(...).generate_single_replicate(replicate_id=0)` with
`coloring_params={"cooling_rate": 0.9995, "max_iterations": 20000}`.

## 3. Hotspots

"Changes results?" says whether the output or the random stream differs from
today's code. Paths are under `tissue_simulator/`.

| Hotspot | Today | Why it is slow | CPU fix (P0) | GPU formulation | Changes results? |
|---|---|---|---|---|---|
| Edge builders `_add_contact_edges`, `_add_radius_edges` and 2-D variants (`spatial_analysis.py:258`, `270`, `279`, `291`) | O(n²) | Python double loop, scalar norm, one `add_edge` per pair | `cKDTree.query_pairs`, then a vectorised per-pair filter | Sorted-grid cell list with a stencil pair kernel | No, if borderline pairs are re-checked (section 6a) |
| `compute_interaction_statistics` (`spatial_analysis.py:406`) | O(types² × edges) | Rescans all edges for each type pair | `bincount` on a pair code, grouped medians | Same | No |
| `SpatialHashGrid.neighbors` through `_fits` and `_clearance` (`packing.py:72`, `508`, `516`) | Per candidate | Generator over a tuple-keyed dict, scalar arithmetic | Array-backed cell list with the same predicate | Folded into parallel packing | No: pure predicate |
| `InhomogeneousPacker.pack`, `_place_by_addition`, `_insert_at_best_clearance` (`packing.py:641`, `543`, `551`) | Sequential | One candidate per random draw | As above | Batched propose-and-resolve (section 6b) | Yes: statistically equivalent only |
| `_relax`, `_overlap_fraction` (`packing.py:560`, `627`) | Iterations × movable cells × neighbours | Python loops | Vectorise over pair arrays | Pair list with scatter-add pushes | Summation order only |
| The 28 null packings (`density.py:1427`, `1335`) | Serial | They share one random generator | Process pool with spawned streams, opt-in | Batch dimension in the cell-list key | Yes: streams |
| Kernel maps and bandwidth search `_smooth`, `_intensity_grids`, `_bandwidth_scores`, `_type_bandwidths` (`density.py:129`, `172`, `203`, `255`) | Types × candidates × pixels | Already array code, but serial over candidates and types | Reuse the edge-correction map per bandwidth | `cupyx.scipy.ndimage.gaussian_filter`, stacked | No; about 1e-12 in double precision |
| Latent fields and patch calibration `_simulate_field`, `_sample_labels`, `_same_label_profile`, `_calibrate_patches` (`density.py:465`, `482`, `502`, `523`) | 40 candidates × 2 repeats of FFTs | Serial `scipy.fft` | `workers=-1` | `cupyx.scipy.fft`, batched over candidates | No, if noise is drawn on the host |
| Trend null `_trend_null` calling `trend_statistic` (`density.py:1278`, `_organization.py:420`; searches at `217`, `231`; optimisers at `151`, `169`) | 19 × 37 starts × 2 optimisers | SciPy optimisers driving Python closures | Draw all null inputs serially, evaluate in a pool | Batched objective over the 12 angles and 25 centres | CPU pool: no. GPU: optimiser tolerance |
| Persistence `superlevel_h0`, `superlevel_h1` (`_persistence.py:65`, `117`) | Pixels × log pixels per map | Python union-find per pixel | numba or Cython union-find, same sort and tie rule | None sensible. Keep on CPU, pool across maps | No, if the tie rule is kept |
| `empty_space_map` (`_units.py:145`) | Tree query per pixel | Fine below about one million pixels | `workers=-1` | Rasterise plus distance transform | Slightly: keep the tree as reference |
| Annealer `colorize`, `_update_statistics_incremental`, `_calculate_cost` (`graph_coloring.py:443`, `218`, `323`) | Per step | Dict copies, string keys | Integer arrays, type × type matrices | Batched chains (section 6c) | Yes: uses stdlib `random` |
| `_spatial_composition_target`, `_calibrate_scale_weights` (`replicate_generator.py:931`, `1025`) | Per cell | Python loop and one `rng.choice` per cell | Vectorise bins and the Sinkhorn step | Constrained assignment on device | Bins: no. Warm start: yes |
| `fidelity_diagnostics` (`replicate_generator.py:121`) | Per cell plus persistence | Comprehensions, per-cell lookups | Array indexing | Kernel maps on GPU, rest on CPU | No |
| `generate_replicates` (`replicate_generator.py`) | One process per replicate | Pickles whole tissue objects | Keep for small samples | Block-diagonal batch (section 6c) | Per-replicate seeds kept |

`trend_statistic` and `_organization._Data` (`_organization.py:98`) take no
random generator, so the trend null can be evaluated in a pool without changing
its stream, provided all null inputs are drawn first in the current order.

## 4. Architecture

**Backend module.** Add `tissue_simulator/_backend.py` with
`get_backend(name=None)`, returning a namespace with `xp`, `ndimage`, `fft`,
`to_device`, `to_host` and `rng(seed)`. The name resolves from the keyword
argument, then the environment variable `TISSUE_SIM_BACKEND`, then `"numpy"`.
The device comes from `TISSUE_SIM_DEVICE`. GPU libraries are imported lazily,
inside the branch that needs them.

**Public switch.** Add a keyword `backend=None` to `DensityModel.fit`,
`ReplicateGenerator` and `SpatialNetworkAnalyzer.build_network_from_tissue` /
`build_network_from_slice`. `None` means the resolved default.

**Host and device boundary.** Move data once per stage. Cell arrays (position,
radius, type) go up at the start of a stage; maps and summary numbers come down
at the end. `DensityModel` and `Layout` fields stay NumPy arrays on the host, so
`to_dict` and `from_dict` (`density.py:1821`, `1837`) are untouched.

**Arrays instead of objects.** Add two internal containers: a
structure-of-arrays for cells, and an edge list holding an `(E, 2)` integer pair
array, distances, and a compressed sparse row index. The annealer and the
interaction statistics read these arrays and never traverse NetworkX on the fast
path.

**NetworkX stays the public return type.** Build the graph from the edge arrays
in `(i, j)` order. Above roughly 200,000 nodes, build it lazily on first access.

**Keeping the default path identical.** The `backend="numpy"` branch calls
today's functions unmodified. Fast paths are new functions selected only by
backend or an explicit flag. A deterministic CPU rewrite (graph, interaction
statistics, cell-list predicate, compiled persistence) becomes the default only
after an equality test against the legacy function passes on the fixtures.

## 5. Stack and packaging

**CuPy is the only required GPU dependency.** The code is written against NumPy
and SciPy (`ndimage.gaussian_filter`, `distance_transform_edt`, `label`,
`scipy.fft`, `bincount`, `searchsorted`), and `cupyx.scipy` mirrors these
one-to-one. `cupy.RawKernel` covers the two custom kernels: the cell-list pair
search and the packing conflict resolution.

- Neighbour search uses a hand-written sorted-grid cell list. RAPIDS/cuSpatial
  is a heavy dependency and does not support the per-pair radius rule.
- Batched annealing also uses CuPy. Reserve `backend="torch"` in the accepted
  names but leave it unimplemented; DLPack interop makes it possible later.

**Packaging.** GPU support is an optional extra, never a core dependency. The
existing extras are `mcp`, `dev` and `docs`. Add:

| Extra | Contents | Purpose |
|---|---|---|
| `fast` | `numba` | Portable CPU speedups: packer cell list, persistence |
| `gpu` | `cupy-cuda12x>=13` | CUDA 12 |
| `gpu-cuda11` | `cupy-cuda11x>=13` | CUDA 11 |
| `all` | `mcp` + `fast` | Every portable extra |

```bash
pip install tissue_simulator              # unchanged, NumPy and SciPy only
pip install "tissue_simulator[fast]"      # CPU speedups
pip install "tissue_simulator[gpu]"       # CUDA 12
pip install "tissue_simulator[all,gpu]"   # everything on a CUDA 12 machine
```

`all` deliberately excludes the GPU extras. CuPy's CUDA wheels exist only for
Linux and Windows with a matching NVIDIA driver, and they are tied to a CUDA
major version. Including them would make `pip install "tissue_simulator[all]"`
fail on macOS and on any machine without CUDA.

Requesting `backend="cupy"` without the extra must raise an error that names the
extra to install. Every pure-Python fallback must keep working when `numba` is
absent. Adding extras edits `pyproject.toml`, so run
`python scripts/release.py check` afterwards, and document the install matrix in
`README.md` and `docs/quickstart.md`.

## 6. The three algorithmic changes

These are the parts that are not a mechanical port. Each carries a risk that the
acceptance tests in section 7 are designed to catch.

### 6a. Neighbour graph

The edge rule is per pair: `d <= f * (r_i + r_j)`, with `f = 1.01` for contact
mode and `1.5` for mechanical mode, or `d <= radius` for radius mode.

- **CPU**: candidate pairs from `cKDTree.query_pairs` at reach
  `f * 2 * r_max`, then the vectorised per-pair filter, then a sort by `(i, j)`.
- **GPU**: bin cells at side `f * 2 * r_max`, sort by cell index, scan the 9- or
  27-cell stencil, emit `i < j` pairs with a count pass followed by a fill pass.

Risks:

- A scalar norm and a vectorised square root of a sum can differ by one unit in
  the last place. Re-evaluate pairs within a relative `1e-9` of the limit using
  the legacy scalar expression on the host, so the edge set is identical.
- A heavy-tailed radius distribution inflates the bin side. Report stencil
  occupancy and switch to two-level binning above a threshold.
- Never form a dense pairwise block. Output is proportional to the edge count.

### 6b. Batched null packings

- **CPU**: a process pool over the null draws with `SeedSequence.spawn`. The 28
  packings share one generator today (`density.py:1403`, `1445`), so pooling
  changes the stream. It must be opt-in, for example a `null_workers` argument.
- **GPU**, per round: every pending ticket draws a candidate (`_candidate`,
  `packing.py:486`, vectorised); candidates overlapping accepted cells are
  rejected; conflicts among candidates are resolved by random priority, iterated
  to a fixed point. Checkerboard tiles of side at least `2 * kappa * 2 * r_max`
  are the simpler alternative.
- **Relaxation**: fixed Jacobi sweeps, keeping the displacement cap and the
  target-overlap stop rule.

Risk: with exact conflict resolution this is sequential addition in a permuted
order, so the result is statistically equivalent but not identical. The
saturation rule also changes, from "consecutive rejections" to "rounds without
an acceptance", which shifts behaviour near jamming. Null p-values
(`density.py:1463`) and unit thresholds can therefore move.

### 6c. Batched annealing

State on the device: integer labels; a neighbour-type histogram `H = A @ L` of
shape cells × types; the type-pair edge matrix `L.T @ A @ L`; per-scale bin ×
type counts; and the size likelihood gathered from `size_nll_matrix`
(`replicate_generator.py:56`). The change in energy for one swap is then a
handful of type-length operations on two rows of `H`, with a correction when the
two cells are neighbours.

- **Ensembles, exact**: concatenate the replicate graphs block-diagonally with a
  segment index. One swap per chain per step keeps exact Metropolis behaviour
  for each chain. Restarts are simply more chains.
- **Large samples, approximate**: several swaps per step, chosen so their closed
  neighbourhoods do not overlap.

Risk: the energy is a global sum of squared count errors (`_calculate_cost`,
`graph_coloring.py:323`), so the change from several simultaneous swaps is not
the sum of the individual changes. Bound the number of swaps per step, recompute
the exact energy each step, and roll the step back if it rises beyond the
Metropolis allowance.

Keep the sequential annealer for samples below about 20,000 cells and for every
seeded legacy run. The step count and temperature schedule must scale with
sample size; see the open questions.

## 7. Work packages and acceptance tests

Deliver in this order. Each package should merge on its own.

| Package | Scope | Acceptance |
|---|---|---|
| **P0** | A `benchmarks/` harness with stage timers and synthetic generators at 10k, 100k and 1M cells. CPU rewrites in `spatial_analysis.py`, `packing.py` (cell list), `_persistence.py`, `replicate_generator.py` (bins and Sinkhorn), `graph_coloring.py` (array state) | Edge sets identical, weights within 1e-12. Packings byte-identical for fixed seeds. Barcodes identical. Existing suite green |
| **P1** | `_backend.py`; kernel maps, bandwidth search, latent fields and patch calibration in `density.py` | Double-precision kernel maps within 1e-9 relative. Same selected bandwidth, patch length and smoothness, and compartment labels on the fixtures |
| **P2** | GPU cell-list graph, edge-list container, lazy NetworkX | Edge set identical to the P0 CPU builder at 10k and 100k cells, in 2-D and 3-D, for all three modes |
| **P3** | Parallel packing, batched nulls, batched trend objective | Two-sample Kolmogorov–Smirnov p > 0.01 over at least 200 packings for nearest-neighbour gap quantiles, overlap fraction, bin-count correlation and the heterogeneity statistics. The `homogeneous` verdict and the detected units agree with the CPU on at least 95 % of fixtures |
| **P4** | Batched annealing, block-diagonal ensembles, warm start on device | Final cost and `divergence_score` not worse than the CPU (one-sided Mann–Whitney, 100 replicates). `fidelity` values within the CPU's replicate-to-replicate spread |
| **P5** | Documentation page, the extras, a `pytest -m gpu` marker, benchmark report template | Strict documentation build passes. GPU tests run only on the server: the repository's CI has a documentation workflow and nothing else |

The fixtures for P3 and P4 are the twelve synthetic scenarios used to validate
0.1.18 (random, Potts clustering at three couplings, central domain, gradient,
rare type, artery, germinal centre, nests, glomeruli, follicles) plus the
controls in `tests/test_units.py`.

## 8. Targets

Estimates, to be confirmed on the server. "Today" extrapolates the baseline:
quadratic for the graph, linear for the fit.

| Stage | Cells | Today | After P0 (CPU) | GPU |
|---|---|---|---|---|
| Graph build | 10k | about 3 min | under 0.2 s | under 0.1 s |
| Graph build | 100k | about 5 h | about 2 s | under 0.2 s |
| Graph build | 1M | weeks | about 20 s | under 1 s |
| Fit with 28 nulls | 10k | about 3 min | about 30 s | under 10 s |
| Fit with 28 nulls | 100k | about 30 min | about 5 min | under 1 min |
| Fit with 28 nulls | 1M | about 5 h | about 1 h | under 5 min |
| One replicate, annealing scaled to size | 1M | infeasible | about 10 min | under 1 min |

Ensemble target: 1,000 replicates of 10,000 cells in under 5 minutes on one GPU.

## 9. Determinism, precision, memory and fallbacks

**Determinism.** Seed every GPU generator per stage from the existing
`SeedSequence([seed, replicate_id]).spawn(4)` (`replicate_generator.py:1123`).
Draw the white noise for latent fields on the host by default, so P1 keeps the
current stream. The guarantee is: same seed, same device and same library
versions give the same output. Across devices, and against the CPU, results are
statistically equivalent and not bit-identical. Say so in the documentation, and
record the backend, device and library versions in `DensityModel.estimation` and
in `ReplicateStatistics`. Where order affects a result, use a sort followed by a
segment reduction instead of atomic scatter.

**Precision.** Double precision is the default for geometry predicates, kernel
maps and energies. Single precision is opt-in for latent fields only. Counts
stay integer.

**Memory.**

- A 10 mm slide at the 5 µm grid is 2000 × 2000 pixels: 32 MB per
  double-precision map.
- Per-type maps at two scales over ten null draws with eight types is about
  5 GB. Chunk over draws.
- A padded complex FFT near 2400² is about 92 MB per field. The 80 calibration
  fields are about 7 GB if batched together. Chunk at 16 or fewer.
- Edge output is roughly 3 to 6 edges per cell, about 100 MB at one million
  cells including distances.
- For batched annealing, keep the total number of batched nodes at or below
  about 5 × 10⁷ per 24 GB of device memory.

**Fallbacks.** On an out-of-memory error, halve the chunk and retry; if it still
fails, fall back to the CPU fast path with a warning. A missing GPU library
raises `ImportError` only when a GPU backend was explicitly requested.

## 10. Open questions for the server team

1. Which GPU models, how much device memory, which CUDA version, and is the
   machine multi-GPU? This decides the wheel and whether ensembles shard across
   devices.
2. Is `numba` acceptable on the server and for end users, or should the compiled
   persistence be Cython?
3. Should annealing steps scale as a multiple of the cell count, and should the
   temperature be normalised by edge count? The cost is in squared-count units,
   so the current 100 to 0.1 schedule is not scale-free.
4. Is the 5 µm grid required at slide scale? Persistence over four million
   pixels, per type, scale and null draw, is the residual CPU bottleneck. Is a
   coarser grid acceptable for unit detection?
5. Once nulls are cheap, should the default `n_null` (19, `density.py:916`)
   rise, for example to 99?
6. Should a one-million-cell slide be fitted as tiles, since a single stationary
   window is a poor model at that extent?
7. What output format is wanted at scale, for example Parquet or Zarr arrays,
   in place of `TissueSection` objects and NetworkX graphs?
8. Who owns the statistical-equivalence fixtures and the sign-off thresholds in
   section 7?

## 11. Hand-off checklist

Please send back, per package:

- The timings table from section 8, filled in with measured values, plus the
  GPU model, device memory, CUDA version and CuPy version.
- The acceptance results from section 7: pass or fail, with the measured
  statistic for each test.
- Peak device memory per stage at 100k and 1M cells.
- Any deviation from this specification and the reason for it.
- Answers to the open questions in section 10.

Repository rules that apply to this work:

- Do not bump the package version as part of a feature change. Releases are
  their own commit; see `docs/notes/releasing.md`.
- Run `python scripts/release.py check` after editing `pyproject.toml`.
- Run `pytest tests/ -v` before each merge. The suite at 0.1.18 has 373 passing
  tests and 2 expected failures.
- Related design documents: `docs/design/single-sample-fidelity-spec.md` and
  `docs/notes/density-aware-packing.md`.
