# Changelog

All notable changes to `tissue_simulator` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.13] - 2026-06-04

### Added

- **`geometry: "2D" | "3D"` parameter on `PhysiCellExporter.export_slice()`**
  (`tissue_simulator/physicell_export.py`). Default `"2D"` preserves existing
  behavior byte-for-byte (slice-plane projection, supplied `z`, disk-area
  volumes). `"3D"` writes each cell's original 3D position from `center_3d`
  and sphere-volume on the cell's original radius — suited to PhysiCell 3D
  simulations where the slice was used as a spatial filter. The `z`
  argument is ignored when `geometry="3D"`.
- **Two new MCP tools** (`tissue_simulator/mcp/server.py`):
  - `export_tissue_to_physicell` — exports the current 3D tissue as a
    PhysiCell IC CSV. Always 3D.
  - `export_slice_to_physicell` — exports the current slice as a PhysiCell
    IC CSV, with the new `geometry` parameter. Closes the prior gap where
    MCP exposed no PhysiCell bridge at all.
- Four new tests in `tests/test_physicell_export.py` (backward-compat
  byte-identity for the default, 3D-semantics, validation, z-ignored).

### Changed

- **`docs/api/physicell.md`** documents the new `geometry` parameter, adds a
  worked example showing both modes, and cross-links to the new MCP tools.
- **`docs/api/mcp.md`** documents the two new tools under a new
  "PhysiCell Export" category.

## [0.1.12] - 2026-05-29

### Fixed

- **Deterministic matplotlib color assignment in `visualize_*` functions.**
  `TissueSection.visualize`, `TissueSlicer.visualize_slice_2d` /
  `visualize_slice_in_3d`, `SpatialNetworkAnalyzer.visualize_network`,
  `visualize_colored_graph`, and `visualize_graph_comparison` previously built
  their `{cell_type: color}` map by iterating an unsorted `set(...)`, whose
  order leaks CPython's per-interpreter `PYTHONHASHSEED`. Two independent
  Python processes (e.g. two CI runs) could therefore assign different colors
  to the same cell type, producing pixel-different figures despite every
  random `seed=` being pinned. The new internal helper
  `tissue_simulator._viz_utils.make_color_map` sorts cell types
  alphabetically before assignment, so the same input always yields the
  same mapping. Cross-function consistency too: "cancer" now always gets
  the same color in `visualize`, `visualize_slice_2d`, and
  `visualize_network`.

### Added

- `tissue_simulator/_viz_utils.py` (internal) — single source of truth for
  the cell-type → color map used across the package's visualizations.
- `tests/test_viz_utils.py` — covers determinism (same input in different
  orders yields the same map; same across `PYTHONHASHSEED`), palette
  passthrough, and empty/duplicate inputs.
- **Code-driven slide deck** on the docs site at
  <https://emcramer.github.io/tissue_simulator/latest/slides/>. Source: a single
  Jupyter notebook (`docs/slides/tour.ipynb`) walking through the package
  end-to-end with code AND rendered matplotlib output side by side — tissue
  generation, slicing, spatial-network analysis, replicate generation,
  simulated-annealing cell-type assignment, and seed reproducibility. Rendered
  to reveal.js via `jupyter nbconvert --to slides --execute` and deployed as
  part of the docs site (runs on every push to `main` and every release tag).
  Now that the package-level color-assignment is deterministic (above), the
  slide-deck build helper no longer needs a `PYTHONHASHSEED=0` workaround.
- `docs/slides/build_slides.py` — small build helper that wraps `nbconvert`;
  both CI and devs invoke the same entry point.
- `docs/slides/index.md` — landing page in the docs site nav linking to the
  rendered `tour.slides.html`.
- `tests/test_slides_notebook.py` — validates the notebook with `nbformat`
  (schema + slide_type metadata) so a malformed deck fails the test suite
  even before CI executes it.
- Notebook tooling (`jupyter`, `nbconvert>=7`, `ipykernel`, `matplotlib`,
  `networkx`) added to the `docs` optional-dependency group.
- CI render step in `.github/workflows/docs.yml` runs the build helper before
  `mike deploy`.

### Changed

- `mkdocs.yml` adds a `Slides` nav entry between `API Reference` and
  `Changelog`; the rendered `tour.slides.html` is listed in `not_in_nav` so
  `--strict` doesn't flag it.
- `docs/index.md` "Where to next" section adds a pointer to the slide deck.
- Per-symbol notes on the deterministic color guarantee added to
  `docs/api/core.md` (`TissueSection.visualize`), `docs/api/slicing.md`
  (both slice visualizers), `docs/api/spatial-analysis.md`
  (`visualize_network`), and `docs/api/graph-coloring.md` (the two
  `visualize_*` helpers, with the user-palette caveat).
- Wiki FAQ entry "Why do my visualizations have different colors each run?"
  pushed to the companion wiki.

## [0.1.11] - 2026-05-29

### Added

- **`type` column alias in `load_tissue_from_csv`** — the loader now accepts
  `type` as a column alias for `cell_type`, enabling direct ingestion of
  PhysiCell `cells.csv` exports (written by the PhysiCell MCP's
  `export_cells_csv`) without a bespoke schema adapter. Precedence:
  `cell_type` (canonical, preferred) → `type` (alias) → `"default"` (fallback,
  unchanged). `load_target_statistics_from_coordinates` inherits the fix
  automatically because it composes `load_tissue_from_csv`. Three new tests
  cover the alias, the precedence rule, and the unchanged fallback.

## [0.1.10] - 2026-05-28

### Added
- **MkDocs Material documentation site** at <https://emcramer.github.io/tissue_simulator/> (config: `mkdocs.yml`). Reuses the existing `docs/api/*.md` hand-written guides and adds:
  - `docs/index.md` site landing page,
  - `docs/changelog.md` (embeds `CHANGELOG.md` via include-markdown),
  - `docs/reference/*.md` (13 mkdocstrings stubs — per-symbol API reference auto-generated from the package's docstrings).
- **Versioned docs via `mike`** with a navbar version switcher; `latest` alias tracks `main`, tagged releases (`v*.*.*`) replace `latest`.
- **First GitHub Actions workflow** `.github/workflows/docs.yml` — auto-deploys on push to `main` and on release tags.
- **`docs` optional-dependency group** in `pyproject.toml`: `mkdocs-material`, `mkdocstrings[python]`, `mkdocs-include-markdown-plugin`, `mike`. Install with `pip install -e ".[docs]"`.
- **`Documentation` URL** in `pyproject.toml`'s `[project.urls]`.
- **GitHub repo wiki seeded** with Home / FAQ / Troubleshooting / Roadmap pages (community-editable, separate from the canonical Pages site).
- **`tests/test_docs_build.py`** — runs `mkdocs build --strict` in CI so broken-link / missing-page regressions fail tests.

### Changed
- `README.md` adds a Docs badge and links to the site + wiki at the top of the Documentation section.
- `docs/quickstart.md` outbound `../README.md` / `../CHANGELOG.md` / `../examples/` / `../tissue_simulator/` links rewritten to in-site pages or absolute GitHub blob URLs so the built site has no broken links.

## [0.1.9] - 2026-05-28

### Added
- **`seed: Optional[int] = None`** on `GraphColorizer.__init__`
  (`graph_coloring.py`) — when provided, the colorizer routes every
  stochastic decision (initial coloring shuffle, per-step pair sampling,
  metropolis acceptance draw) through an instance-bound
  `random.Random(seed)`, making `colorize(...)` bit-reproducible across
  Python processes independent of `PYTHONHASHSEED`. When `seed=None`
  (default) behavior is unchanged.
- **`seed=` plumbed through the workflow entry points**
  (`tissue_workflow.py`): `TissueNetworkWorkflow.assign_cell_types`,
  `TissueNetworkWorkflow.run_complete_workflow`, and the top-level
  `quick_workflow` all forward a `seed` to the underlying
  `GraphColorizer`. Matches the v0.1.2 reproducibility plumbing already
  in `TissueSection`, `SpherePacker`, and `ReplicateGenerator`.
- **`assign_cell_types`** (MCP tool, `mcp/server.py`) — slice the current
  tissue at `z_position`, build a radius-mode spatial network, load
  target statistics from a graph-coloring-format CSV (`node_counts`,
  `edge_counts`, `neighbor_dist`), and run `GraphColorizer.colorize()`
  with the given `seed`. The per-node coloring is stored on the server.
  Exposes the two-stage `generate_cells` -> `assign_cell_types`
  workflow through MCP so an LLM can drive structured multi-type
  cell-type assignment end-to-end.
- New reproducibility tests for graph coloring (same seed +
  same inputs => identical assignment dict; different seeds =>
  different assignments; `seed=None` preserves legacy behavior).

### Changed
- **`docs/api/graph-coloring.md`** documents the new `seed` parameter
  with a short reproducibility example, and adds a "When to use what"
  section clarifying the design intent: `GraphColorizer` is for
  label-only assignment on a fixed graph; `ReplicateGenerator` is for
  unstructured replicate generation by repacking and does NOT compose
  `GraphColorizer` today; for structured multi-type targets, use the
  two-stage `TissueWorkflow` (`generate_cells` -> `assign_cell_types`).
- **`docs/api/replicate-generation.md`** gains a cross-link callout
  pointing readers with structured multi-type targets at the
  "When to use what" section of `graph-coloring.md`.
- **`docs/guides/complete-workflow.md`** adds `seed=42` to the
  `quick_workflow(...)` and `TissueNetworkWorkflow.run_complete_workflow(...)`
  examples, with a brief note that the colorize step is then
  bit-reproducible.
- **`docs/api/mcp.md`** documents the new `assign_cell_types` MCP tool
  (a new "Cell Type Assignment" category and a per-tool reference entry
  matching the v0.1.8 format).

## [0.1.8] - 2026-05-27

### Added
- **`load_tissue_from_csv`** (MCP tool, `mcp/server.py`) — load a
  `TissueSection` from a coordinate CSV (`x, y, z, radius, cell_type,
  is_boundary`) and set it as the current tissue, so slicing, analysis, and
  statistics tools can run on externally sourced tissue. Dimensions are
  inferred from the coordinate bounds when omitted. The inverse of
  `export_tissue_csv`.
- **`load_target_statistics_from_coordinates`** (MCP tool, `mcp/server.py`) —
  compute full target statistics (interactions **plus**
  `cell_type_proportions` and `target_density`) straight from a coordinate
  CSV. Distinct from `load_target_statistics`'s `csv_filepath`, which reads a
  precomputed interaction table and leaves proportions/density unset.
- Together these expose the v0.1.7 external-cell ingest functions
  (`load_tissue_from_csv`, `load_target_statistics_from_coordinates`) through
  the MCP server, so an LLM can turn an externally measured or generated
  (e.g. PhysiCell) tissue into a `TargetStatistics` for `ReplicateGenerator`
  without going through the random packer.
- **`tests/test_mcp_server.py`** covering the two new MCP tools.

### Changed
- **`docs/api/mcp.md`** documents the two new tools (a new "Data Loading"
  category, per-tool reference entries, and an external-tissue workflow
  example).

## [0.1.7] - 2026-05-26

### Added
- **`TissueSection.from_cells(...)`** (classmethod, `tissue.py`) — build a
  tissue from already-positioned `Cell` objects instead of the random
  packer. Infers any omitted dimension from the cell-center bounding box
  (with a largest-diameter fallback for constant-axis 2D inputs) and derives
  per-type `(min, max)` radii when not given. The inverse of `export_to_csv`.
- **`load_tissue_from_csv(...)`** (module-level, `tissue.py`) — load a
  `TissueSection` from a coordinate CSV (`x, y, z, radius, cell_type,
  is_boundary`) written by `export_to_csv`; `radius` and `is_boundary` are
  optional. Round-trips a packed tissue's cells exactly.
- **`load_target_statistics_from_coordinates(...)`** (module-level,
  `replicate_generator.py`) — full `TargetStatistics` (interactions **plus**
  `cell_type_proportions` and `target_density`) straight from a coordinate
  CSV. Distinct from `load_target_statistics_from_csv`, which reads a
  precomputed interaction table and leaves proportions/density unset.
- Together these unblock turning an externally measured or externally
  generated (e.g. PhysiCell) tissue into a `TargetStatistics` for
  `ReplicateGenerator` without going through the random packer.
- Tests for round-trip fidelity, derived radii, 2D dimension inference, full
  target-stats validation, and the new public exports
  (`tests/test_tissue_simulator.py`, `tests/test_replicate_generator.py`).

### Changed
- **`tissue_simulator/__init__.py`** now exports `load_tissue_from_csv` and
  `load_target_statistics_from_coordinates`.
- **`docs/api/core.md`** documents `TissueSection.from_cells` and
  `load_tissue_from_csv`.
- **`docs/api/replicate-generation.md`** documents
  `load_target_statistics_from_coordinates` and clarifies the two CSV shapes
  (coordinate CSV vs. precomputed interaction table).

## [0.1.6] - 2026-05-15

### Added
- **`docs/api/convergence.md`** — full per-symbol reference for
  `adf_test`, `mann_kendall_test`, `rolling_cv`, `find_convergence_time`,
  and `MultiMetricConvergence`. Worked example verified end-to-end
  (synthetic transient + AR(1) trajectory, converges at t=99 / 280).
- **`docs/api/power-analysis.md`** — per-symbol reference for
  `cohens_d`, `coefficient_of_variation`, `required_replicates`,
  `power_curve`, `compare_initialization_variance`, and
  `summarize_power_analysis`. Worked example produces a three-method
  variance / effect-size report.
- **`docs/api/physicell.md`** — dedicated PhysiCell bridge page covering
  `PhysiCellExporter`, `export_to_physicell`, `PhysiCellReader`,
  `read_physicell_output`, and `stats_to_target_statistics`, plus the
  round-trip from tissue -> CSV -> reader -> `ReplicateGenerator`.
  Includes pyMCDS attribution.

### Changed
- **`docs/api/core.md`** trimmed and retitled from a top-level user guide
  (491 lines) to a focused per-module API reference (296 lines) for
  `TissueSection`, `Cell`, and `SpherePacker`. Installation / quickstart
  prose now lives in `docs/quickstart.md` only.
- **Decorative checkmark / cross emoji removed** from `docs/api/slicing.md`,
  `docs/api/spatial-analysis.md`, `docs/api/replicate-generation.md`, and
  `docs/guides/mcp.md`. Markdown style now matches the plain convention
  set in `CLAUDE.md`.

### Fixed
- The non-runnable return-type pseudo-code block in
  `docs/api/power-analysis.md` is now fenced as `text`, not `python`, so
  automated snippet extraction picks up only the runnable example.

### Removed
- **`PAPER_PLAN_REPORT.md` is no longer tracked on `main`.** (The v0.1.5
  attempt didn't make it into the commit; this release actually drops
  it.) The file remains as a local working-copy artifact and is listed
  in `.gitignore`.

  Historical note: the file is still present in tag history v0.1.1
  through v0.1.5. Purging it from history requires a force-rewrite
  (`git filter-repo`) plus re-pushing tags and recreating releases —
  see the project notes for that one-off operation.

### Tests
- 120 / 120 tests pass.

## [0.1.5] - 2026-05-14

### Changed
- **`docs/quickstart.md` fully rewritten** for the v0.1.x feature surface.
  Adds a runnable five-minute tour covering tissue generation, slicing,
  spatial-network construction, and CSV export; a "Headline workflows"
  section pointing at replicate generation, graph coloring, the PhysiCell
  bridge, convergence and power analysis, MCP, and the GUI; and a
  reproducibility note covering the `seed=` plumbing and the `nan`
  divergence semantics from v0.1.2. Snippet verified end-to-end on the
  real package: 495 cells, 2424 edges, runs in ~5 s.
- **`mcp_config_claude_desktop.example.json`** now uses the placeholder
  path `/ABSOLUTE/PATH/TO/...` rather than a hard-coded developer path.

### Removed
- **`PAPER_PLAN_REPORT.md` is no longer tracked.** It remains as a local
  internal-only doc (now properly listed in `.gitignore`).

## [0.1.4] - 2026-05-13

### Changed
- **Repository cleanup**: removed obsolete root-level helper scripts
  (`fix_graphml_export.py`, `quick_test.sh`, `quick_fix_mcp.sh`),
  historical "feature complete" narrative docs
  (`GRAPH_COLORING_INTEGRATION.md`, `REPLICATE_FEATURE_SUMMARY.md`,
  `docs/SLICING_SUMMARY.md`, `docs/STRUCTURE.md`,
  `docs/feature_tracking/*_COMPLETE.md`), and the now-empty
  `docs/feature_tracking/` directory. The historical signal these
  files carried is preserved in this `CHANGELOG.md`.
- **Docs reorganization**: split `docs/` into `docs/api/` (per-module
  reference), `docs/guides/` (tutorials), `docs/design/` (research /
  paper-track artifacts), and `docs/notes/` (maintainer notes such as
  the `plot_surface` regression incident). File renames are git-tracked
  so `git log --follow` still works.
- **`mcp_config_claude_desktop.json` → `mcp_config_claude_desktop.example.json`**
  to make its template nature explicit; the file contains a hard-coded
  absolute path users must adapt.
- **`configure_mcp.sh` → `scripts/configure_mcp.sh`** (plus a brief
  `scripts/README.md`).
- **`mcp` dependency moved to an optional extra**: `pip install -r requirements.txt`
  no longer pulls in `mcp`. Install with `pip install tissue_simulator[mcp]`
  or `pip install mcp` directly. This is a minor breaking change for users
  who relied on the implicit install.

### Added
- **`pyproject.toml`** with PEP 621 metadata as the authoritative source
  for package metadata; `setup.py` is now a thin shim.
- **`CHANGELOG.md`** (this file).
- **`requirements-dev.txt`** with pytest, build, twine, and the optional
  `mcp` extra for contributors.

### Fixed
- **`MANIFEST.in`** previously referenced `GUIDE.md` and `QUICKSTART.md`
  at the repo root, where they never lived; sdists silently shipped
  without those docs. The MANIFEST now uses `recursive-include docs *.md`.
- **`setup.py`** author placeholders replaced with the real maintainer.
  Project classified as `Development Status :: 4 - Beta`.
- **`.gitignore`** extended to cover `.claude/`, `.cursor/`, `.continue/`,
  `.aider*`, `.mypy_cache/`, `.ruff_cache/`, and `examples/output/`;
  `.gemini/*` glob hardened to `.gemini/`.

### Tests
- 120 / 120 tests pass.

## [0.1.3] - 2026-05-13

### Changed
- Removed dead RNG code from `TissueSection`. `SpherePacker` is now the
  single source of randomness; `seed` keyword arguments are preserved on
  the public API.
- PhysiCell `.mat` cells-matrix selection is now shape-aware: prefers
  elongated candidates with plausible signal-count axes, warns when the
  heuristic fires, and raises on close-call ambiguity.

### Tests
- 120 / 120 tests pass.

## [0.1.2] - 2026-05-13

### Added
- Explicit RNG seeding through `SpherePacker`, `TissueSection`, and
  `ReplicateGenerator`; the Result-1 demo is now bit-reproducible across
  `PYTHONHASHSEED` values.
- Real PhysiCell `.mat` parsing via pyMCDS-or-direct-scipy fallback, with
  attribution. Non-standard matrix shapes now raise instead of silently
  mis-decoding.

### Changed
- `_compute_interaction_divergence` returns `nan` (not `0`) when both
  target and measured are zero for a given pair. The aggregate uses
  `np.nanmean`, so empty-signal pairs no longer masquerade as perfect
  matches.

### Tests
- 116 / 116 tests pass.

## [0.1.1] - 2026-05-13

### Added
- **PhysiCell bridge**: `PhysiCellExporter` and `PhysiCellReader` modules
  with a `stats_to_target_statistics` adapter that lets ABM snapshots
  drive `ReplicateGenerator` directly.
- **Convergence diagnostics**: ADF, Mann-Kendall (tie-corrected), rolling
  CV, `find_convergence_time`, and `MultiMetricConvergence`.
- **Power analysis**: Cohen's d, required replicates, power curves, and
  variance comparison across initialization strategies.
- New Result-1 demo: `examples/result1_replicate_demonstration.py`.

### Changed
- `statsmodels` added as a required dependency (used by convergence and
  power-analysis modules).

### Tests
- 105 / 105 tests pass.

## [0.1.0] - 2025-10-28

### Added
- Initial public release.
- Core 3D tissue generation: `TissueSection`, `Cell`, `SpherePacker`
  (random sequential addition).
- 2D slicing: `TissueSlicer`, `SliceCell`, `create_standard_slices` for
  arbitrary-angle planar sections.
- Spatial network analysis: `SpatialNetworkAnalyzer` with contact and
  radius edge modes; global, per-cell-type, and pairwise interaction
  statistics.
- Graph-based cell-type assignment: `GraphColorizer` using simulated
  annealing to match target node, edge, and neighbor-distribution
  statistics.
- Replicate generation: `ReplicateGenerator` for batch tissues matching
  target spatial statistics.
- Unified workflow: `TissueNetworkWorkflow` and `quick_workflow`.
- Evaluation metrics: Jensen-Shannon divergence, cosine similarity,
  `evaluate_graph_coloring`.
- PyQt5 GUI: `tissue-simulator` console entry point.
- MCP server (`tissue_simulator.mcp.server`) exposing 11 tools for LLM
  integration via Model Context Protocol.
- Examples, tests, MIT license, and documentation.

[0.1.4]: https://github.com/emcramer/tissue_simulator/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/emcramer/tissue_simulator/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/emcramer/tissue_simulator/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/emcramer/tissue_simulator/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/emcramer/tissue_simulator/releases/tag/v0.1.0
