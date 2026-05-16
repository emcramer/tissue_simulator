# PhysiCell Bridge: Export and Read

## Overview

The PhysiCell bridge connects `tissue_simulator` to the
[PhysiCell](http://physicell.org) agent-based modeling framework in two
directions: **export** serializes a `TissueSection` (or `SliceCell`
list) to a PhysiCell initial-condition CSV that seeds an ABM run, and
**read** parses PhysiCell snapshots (CSV or `.mat`) back into the schema
consumed by `ReplicateGenerator`. The `.mat` reader delegates to the
optional
[PhysiCell-Tools/python-loader](https://github.com/PhysiCell-Tools/python-loader)
project (`pyMCDS` / `pcdl`, BSD-3-Clause) when importable; otherwise it
uses a direct `scipy.io.loadmat` fallback parser based on the
documented PhysiCell 1.10+ default cells-matrix layout. Neither
`pyMCDS` nor `pcdl` is a required dependency. Added in v0.1.1. Public
symbols on this page: `PhysiCellExporter`, `export_to_physicell`,
`PhysiCellReader`, `read_physicell_output`,
`stats_to_target_statistics`.

## Export

### `PhysiCellExporter(default_format="modern")`

PhysiCell historically accepts two IC CSV layouts:

- **modern** (PhysiCell 1.10+, default): header row `x,y,z,type` with
  `type` as the string name used in the PhysiCell XML config. An
  optional `volume` column may follow.
- **legacy**: headerless rows of `x,y,z,cell_type_id,volume`, where
  `cell_type_id` is the integer ID PhysiCell uses to resolve the cell
  definition in the XML. The `volume` column is mandatory.

Volumes are derived from radii: `4/3 * pi * r^3` for 3D tissues and
`pi * r^2` (disk area) for 2D slice exports.

#### `export_tissue(tissue, output_path, cell_type_mapping=None, include_volume=True, format="modern")`

Write every cell in `tissue.cells` to `output_path` and return the
written path.

- `cell_type_mapping`: dict mapping cell-type name to PhysiCell
  integer ID. Required when `format="legacy"`. Optional but validated
  for completeness when `format="modern"` (typos surface as
  `ValueError`).
- `include_volume`: append the `volume` column in modern format. In
  legacy format the volume column is mandatory; `include_volume=False`
  raises `ValueError` (tightened in v0.1.1).
- `format`: `"modern"` or `"legacy"`.

#### `export_slice(slice_cells, output_path, z=0.0, cell_type_mapping=None, include_volume=True, format="modern")`

Same arguments as `export_tissue`, but for an iterable of `SliceCell`
from `TissueSlicer`. Every row uses the supplied `z` (PhysiCell 2D
models conventionally use `z=0`) and the `volume` column carries the
2D disk area of the slice intersection.

#### `build_default_mapping(cell_types)`

Static helper that builds a stable name-to-ID mapping (alphabetical,
starting at 0) so legacy CSVs are reproducible regardless of input
ordering.

### `export_to_physicell(tissue_or_slice_cells, output_path, **kwargs)`

Module-level convenience dispatcher. A `TissueSection` is forwarded to
`export_tissue`; any other iterable is treated as `SliceCell` rows and
forwarded to `export_slice`. Accepts `cell_type_mapping`,
`include_volume`, `format`, and (for slices) `z`.

### Legacy-format example

```python
from tissue_simulator import PhysiCellExporter, TissueSection

tissue = TissueSection(
    height=200, width=200, thickness=50,
    cell_radii={"tumor": (8, 12), "immune": (5, 8)},
    seed=42,
)
tissue.generate_cells(max_attempts=1500)

# Legacy format requires an explicit name-to-id mapping.
mapping = PhysiCellExporter.build_default_mapping(["tumor", "immune"])
PhysiCellExporter().export_tissue(
    tissue, "ic_legacy.csv",
    cell_type_mapping=mapping, format="legacy",
)
```

Modern-format export is covered in the [Worked example](#worked-example)
below. The PhysiCell simulation domain must contain every cell; this
exporter does not clip or shift coordinates. Align `tissue_simulator`
bounds to PhysiCell's domain upstream of export.

## Read

### `PhysiCellReader(default_radius=8.412710547954228)`

`default_radius` is used when a snapshot has neither a `radius` nor a
`volume` column for a row (the default matches a PhysiCell
2494 cubic-micrometer cell).

#### `load_snapshot_csv(csv_path)`

Parse a single PhysiCell CSV snapshot. Accepts column aliases
(`position_x` for `x`, `total_volume` for `volume`, `cell_type` or
`type` for the type column). At minimum `x`, `y`, and a cell-type
column must be present; `z` defaults to `0` when absent. Returns a
list of cell dicts with keys `x`, `y`, `z`, `radius`, `cell_type`,
and optionally `id`, `volume`, `timestep`. When only `volume` is
present, `radius` is back-computed from the inverse sphere-volume
formula.

#### `load_snapshot_mat(mat_path, xml_path=None)`

Parse a `.mat` snapshot. Two backends are attempted in order:

1. **pyMCDS / pcdl** when importable *and* `xml_path` is supplied
   (pyMCDS needs the matching `output*.xml` to recover the
   cell-definition mapping).
2. **Direct `scipy.io.loadmat` fallback** otherwise.

The fallback parser assumes the documented PhysiCell 1.10+ default
cells-matrix column layout:

| Column | Meaning                |
| ------ | ---------------------- |
| 0      | ID                     |
| 1      | position x             |
| 2      | position y             |
| 3      | position z             |
| 4      | total_volume           |
| 5      | cell_type (integer ID) |

PhysiCell stores the matrix as `(n_signals, n_cells)`. The reader
selects the cells matrix by name first (`cells`, `cell`, or `Cells`)
and only falls back to a shape-aware heuristic when none of those
names resolve to a 2-D ndarray. The heuristic prefers elongated
candidates whose minor axis is in the plausible PhysiCell signal-count
range `[6, 250]`; if the top two candidates score within 10% of each
other the call raises `ValueError` rather than silently picking the
wrong array (shape-aware ambiguity check, v0.1.3). If axis 0 has fewer
than the six signals required by the documented layout, the parser
raises — non-standard transposed exports that earlier versions
mis-decoded are no longer silently accepted (v0.1.2).

Without an `xml_path` (or when pyMCDS is unavailable), `cell_type`
values are stringified integer IDs and a `UserWarning` is emitted.

#### `load_time_series(output_dir, pattern="snapshot_*.csv")`

Walk a directory of snapshot CSVs and return a list of
`(timestep, cells)` tuples sorted by the numeric token in each
filename (e.g. `snapshot_00010.csv` -> 10.0). Filename is the
tiebreaker when two snapshots share the same numeric token.

#### `compute_spatial_stats(cells, mode="contact", radius_threshold=None)`

Compute the `node_counts` / `edge_counts` / `neighbor_dist` triple
consumed internally by `GraphColorizer`. `mode="contact"` connects
cells whose surfaces touch (`distance <= r_i + r_j` with the same 1%
tolerance used by `SpatialNetworkAnalyzer`); `mode="radius"` connects
cells whose centers lie within `radius_threshold`.
`neighbor_dist[type_a][type_b]` is the average number of `type_b`
neighbors per `type_a` cell. Requires `scipy` (uses `KDTree`).

#### `compute_stats_time_series(time_series, mode="contact", radius_threshold=None)`

Apply `compute_spatial_stats` over the output of `load_time_series`,
returning `[(timestep, stats), ...]`.

### `read_physicell_output(path, **kwargs)`

Module-level dispatcher: a directory routes to `load_time_series`
(accepts `pattern`); a `.csv` to `load_snapshot_csv`; a `.mat` to
`load_snapshot_mat` (accepts `xml_path`). `default_radius` may also
be passed and is forwarded to the underlying `PhysiCellReader`.

## Bridging back to `ReplicateGenerator`

### `stats_to_target_statistics(stats, target_cell_count=None, target_density=None)`

The reader emits a flat dict (`node_counts`, `edge_counts`,
`neighbor_dist`) matching `GraphColorizer`'s internal schema.
`ReplicateGenerator` expects the richer
[`TargetStatistics`](replicate-generation.md) dataclass with a list of
`InteractionStatistics` plus cell-type proportions; this adapter does
the conversion. `avg_distance` and `median_distance` are not in the
reader's output and are set to `0.0` on the synthesized records —
compute them separately if needed. `target_cell_count` defaults to
`sum(node_counts.values())`.

```python
cells = read_physicell_output("output/")
stats = reader.compute_stats_time_series(cells)
target = stats_to_target_statistics(stats[-1][1])  # last timestep
gen = ReplicateGenerator(target_stats=target, ...)
```

## Worked example

Generates a small 3-cell-type tissue, exports it to the modern
PhysiCell IC CSV layout in a temp directory, reads it back, computes
contact-mode statistics, and confirms `TargetStatistics.validate()`
succeeds. Runs end-to-end in a few seconds.

```python
import os
import tempfile

from tissue_simulator import (
    PhysiCellReader,
    TissueSection,
    export_to_physicell,
    read_physicell_output,
    stats_to_target_statistics,
)

# Step 1: generate a small 3-cell-type tissue (seeded for reproducibility).
tissue = TissueSection(
    height=120, width=120, thickness=40,
    cell_radii={
        "tumor": (8.0, 12.0),
        "immune": (5.0, 8.0),
        "stroma": (6.0, 10.0),
    },
    seed=20260515,
)
n_cells = tissue.generate_cells(max_attempts=1500, min_spacing=0.5)
print(f"Generated {n_cells} cells")

# Step 2: export to PhysiCell's modern IC CSV format (header + volume col).
with tempfile.TemporaryDirectory() as tmp:
    ic_path = os.path.join(tmp, "physicell_ic.csv")
    export_to_physicell(tissue, ic_path, format="modern", include_volume=True)

    # Step 3: read it back as if it were a PhysiCell snapshot.
    reader = PhysiCellReader()
    cells = reader.load_snapshot_csv(ic_path)
    cells_via_dispatcher = read_physicell_output(ic_path)
    assert len(cells) == len(cells_via_dispatcher) == n_cells

    # Step 4: compute contact-mode spatial statistics on the snapshot.
    stats = reader.compute_spatial_stats(cells, mode="contact")
    print("node_counts:", stats["node_counts"])
    print("edge_count_total:", sum(stats["edge_counts"].values()))

    # Step 5: convert to a ReplicateGenerator-ready TargetStatistics.
    target = stats_to_target_statistics(stats)
    target.validate()  # raises if malformed; succeeds here.
    print(
        f"TargetStatistics OK: "
        f"target_cell_count={target.target_cell_count}, "
        f"n_pairs={len(target.interaction_stats)}"
    )
```

## References

- PhysiCell: http://physicell.org
- PhysiCell-Tools `python-loader` (pyMCDS / pcdl), BSD-3-Clause:
  https://github.com/PhysiCell-Tools/python-loader
- [`replicate-generation.md`](replicate-generation.md) -
  `TargetStatistics` / `ReplicateGenerator` reference.
- [`spatial-analysis.md`](spatial-analysis.md) -
  `InteractionStatistics` and contact / radius edge semantics shared
  with the reader.
- [`../../tissue_simulator/physicell_export.py`](../../tissue_simulator/physicell_export.py)
  and
  [`../../tissue_simulator/physicell_reader.py`](../../tissue_simulator/physicell_reader.py)
  for the implementations.
- [`../../CHANGELOG.md`](../../CHANGELOG.md) - v0.1.1 (bridge),
  v0.1.2 (pyMCDS attribution + non-standard-shape rejection),
  v0.1.3 (shape-aware `.mat` selection).
