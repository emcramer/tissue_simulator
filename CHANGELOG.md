# Changelog

All notable changes to `tissue_simulator` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
