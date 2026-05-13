"""
Unit tests for the PhysiCell snapshot reader.
"""

import csv
import math
import os
import tempfile
import unittest
import warnings
from typing import Dict, List

import numpy as np

from tissue_simulator.physicell_reader import (
    PhysiCellReader,
    read_physicell_output,
    stats_to_target_statistics,
)
from tissue_simulator.replicate_generator import TargetStatistics
from tissue_simulator.spatial_analysis import InteractionStatistics


def _make_physicell_cells_matrix(n_cells: int,
                                 type_ids: List[int],
                                 volume: float = 2494.0) -> np.ndarray:
    """
    Build a fake PhysiCell-style cells matrix of shape (6, n_cells).

    Rows follow the documented PhysiCell 1.10+ default layout:
    0=ID, 1=x, 2=y, 3=z, 4=total_volume, 5=cell_type.
    """
    if len(type_ids) != n_cells:
        raise ValueError("type_ids must have length n_cells")
    matrix = np.zeros((6, n_cells), dtype=float)
    matrix[0, :] = np.arange(n_cells, dtype=float)
    matrix[1, :] = np.arange(n_cells, dtype=float) * 10.0
    matrix[2, :] = np.arange(n_cells, dtype=float) * 10.0 + 1.0
    matrix[3, :] = np.arange(n_cells, dtype=float) * 10.0 + 2.0
    matrix[4, :] = float(volume)
    matrix[5, :] = np.asarray(type_ids, dtype=float)
    return matrix


def _write_cell_definitions_xml(path: str,
                                mapping: Dict[int, str]) -> None:
    """Write a tiny PhysiCell XML with <cell_definitions> entries."""
    lines = [
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
        "<PhysiCell_snapshot>",
        "  <cell_definitions>",
    ]
    for type_id, name in mapping.items():
        lines.append(
            f"    <cell_definition ID=\"{type_id}\" name=\"{name}\">"
        )
        lines.append("    </cell_definition>")
    lines.append("  </cell_definitions>")
    lines.append("</PhysiCell_snapshot>")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def _write_snapshot_csv(path: str,
                        rows: List[Dict],
                        include_volume: bool = False,
                        include_timestep: bool = False) -> None:
    """Write a small PhysiCell-like CSV snapshot to disk."""
    fieldnames = ["x", "y", "z", "radius", "cell_type"]
    if include_volume:
        fieldnames.append("volume")
    if include_timestep:
        fieldnames.append("timestep")

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


class TestSnapshotCSVLoading(unittest.TestCase):
    """Tests for load_snapshot_csv."""

    def test_loads_basic_columns(self):
        reader = PhysiCellReader()
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "snapshot.csv")
            _write_snapshot_csv(csv_path, [
                {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
                {"x": 8.0, "y": 0.0, "z": 0.0, "radius": 3.0, "cell_type": "immune"},
            ])
            cells = reader.load_snapshot_csv(csv_path)

        self.assertEqual(len(cells), 2)
        self.assertEqual(cells[0]["cell_type"], "cancer")
        self.assertAlmostEqual(cells[0]["radius"], 5.0)
        self.assertEqual(cells[1]["cell_type"], "immune")
        self.assertAlmostEqual(cells[1]["x"], 8.0)

    def test_missing_file_raises(self):
        reader = PhysiCellReader()
        with self.assertRaises(FileNotFoundError):
            reader.load_snapshot_csv("/nonexistent/path/snapshot.csv")

    def test_volume_converts_to_radius_when_missing(self):
        reader = PhysiCellReader()
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "vol.csv")
            with open(csv_path, "w", newline="") as fh:
                writer = csv.DictWriter(fh,
                                        fieldnames=["x", "y", "z",
                                                    "volume", "cell_type"])
                writer.writeheader()
                writer.writerow({"x": 0.0, "y": 0.0, "z": 0.0,
                                 "volume": 2494.0, "cell_type": "cancer"})
            cells = reader.load_snapshot_csv(csv_path)

        self.assertEqual(len(cells), 1)
        self.assertGreater(cells[0]["radius"], 8.0)
        self.assertLess(cells[0]["radius"], 9.0)

    def test_missing_required_columns_raises(self):
        reader = PhysiCellReader()
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "bad.csv")
            with open(csv_path, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=["foo", "bar"])
                writer.writeheader()
                writer.writerow({"foo": 1, "bar": 2})
            with self.assertRaises(ValueError):
                reader.load_snapshot_csv(csv_path)


class TestTimeSeriesLoading(unittest.TestCase):
    """Tests for load_time_series and read_physicell_output dispatch."""

    def test_three_snapshots_sorted_by_timestep(self):
        reader = PhysiCellReader()
        with tempfile.TemporaryDirectory() as tmp:
            # Write out of order intentionally to test sorting.
            for t in (20, 0, 10):
                fname = f"snapshot_{t:05d}.csv"
                _write_snapshot_csv(
                    os.path.join(tmp, fname),
                    [
                        {"x": float(t), "y": 0.0, "z": 0.0,
                         "radius": 5.0, "cell_type": "cancer"},
                        {"x": float(t) + 8.0, "y": 0.0, "z": 0.0,
                         "radius": 3.0, "cell_type": "immune"},
                    ],
                )

            series = reader.load_time_series(tmp, pattern="snapshot_*.csv")

        self.assertEqual(len(series), 3)
        timesteps = [t for t, _ in series]
        self.assertEqual(timesteps, [0.0, 10.0, 20.0])
        for _, cells in series:
            self.assertEqual(len(cells), 2)

    def test_read_physicell_output_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            # Single-file dispatch
            csv_path = os.path.join(tmp, "one.csv")
            _write_snapshot_csv(csv_path, [
                {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0,
                 "cell_type": "cancer"},
            ])
            single = read_physicell_output(csv_path)
            self.assertIsInstance(single, list)
            self.assertEqual(len(single), 1)
            self.assertEqual(single[0]["cell_type"], "cancer")

            # Directory dispatch
            for t in (0, 1):
                _write_snapshot_csv(
                    os.path.join(tmp, f"snapshot_{t:03d}.csv"),
                    [{"x": 0.0, "y": 0.0, "z": 0.0,
                      "radius": 5.0, "cell_type": "cancer"}],
                )
            series = read_physicell_output(tmp, pattern="snapshot_*.csv")
            self.assertEqual(len(series), 2)
            self.assertEqual(series[0][0], 0.0)
            self.assertEqual(series[1][0], 1.0)


class TestSpatialStatsSchema(unittest.TestCase):
    """Ensure the stats dict matches the TargetStatistics-compatible schema."""

    def test_stats_dict_keys(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 8.0, "y": 0.0, "z": 0.0, "radius": 3.0, "cell_type": "immune"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="contact")

        self.assertIn("node_counts", stats)
        self.assertIn("edge_counts", stats)
        self.assertIn("neighbor_dist", stats)
        self.assertIsInstance(stats["node_counts"], dict)
        self.assertIsInstance(stats["edge_counts"], dict)
        self.assertIsInstance(stats["neighbor_dist"], dict)

    def test_node_counts_correct(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 50.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 100.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "immune"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="radius",
                                             radius_threshold=10.0)
        self.assertEqual(stats["node_counts"], {"cancer": 2, "immune": 1})


class TestContactMode(unittest.TestCase):
    """Tests for contact-mode edge detection."""

    def test_two_touching_cells_produce_one_edge(self):
        reader = PhysiCellReader()
        # Two spheres of radius 5 with centers 8 micrometers apart:
        # distance (8) <= r1 + r2 (10), so they should touch.
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 8.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "immune"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="contact")

        self.assertEqual(sum(stats["edge_counts"].values()), 1)
        self.assertEqual(stats["edge_counts"].get("cancer-immune"), 1)

    def test_separated_cells_produce_no_edges(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 100.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "immune"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="contact")
        self.assertEqual(sum(stats["edge_counts"].values()), 0)

    def test_neighbor_dist_for_two_touching_cells(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 8.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "immune"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="contact")
        # Single cancer cell has 1 immune neighbor on average; vice versa.
        self.assertAlmostEqual(stats["neighbor_dist"]["cancer"]["immune"], 1.0)
        self.assertAlmostEqual(stats["neighbor_dist"]["immune"]["cancer"], 1.0)
        self.assertAlmostEqual(stats["neighbor_dist"]["cancer"]["cancer"], 0.0)


class TestRadiusMode(unittest.TestCase):
    """Tests for radius-mode edge detection."""

    def test_radius_threshold_connects_within_distance(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 15.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 100.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="radius",
                                             radius_threshold=20.0)
        self.assertEqual(stats["edge_counts"].get("cancer-cancer"), 1)

    def test_radius_mode_requires_threshold(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
        ]
        with self.assertRaises(ValueError):
            reader.compute_spatial_stats(cells, mode="radius")

    def test_same_type_neighbor_dist_counts_both_endpoints(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
            {"x": 15.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "cancer"},
        ]
        stats = reader.compute_spatial_stats(cells, mode="radius",
                                             radius_threshold=20.0)
        # Single edge between 2 cancer cells: each one has 1 cancer neighbor.
        self.assertAlmostEqual(stats["neighbor_dist"]["cancer"]["cancer"], 1.0)


class TestTimeSeriesStats(unittest.TestCase):
    """Tests for compute_stats_time_series."""

    def test_stats_applied_across_series(self):
        reader = PhysiCellReader()
        series = [
            (0.0, [
                {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0,
                 "cell_type": "cancer"},
                {"x": 8.0, "y": 0.0, "z": 0.0, "radius": 5.0,
                 "cell_type": "immune"},
            ]),
            (10.0, [
                {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0,
                 "cell_type": "cancer"},
                {"x": 100.0, "y": 0.0, "z": 0.0, "radius": 5.0,
                 "cell_type": "immune"},
            ]),
        ]
        stats_series = reader.compute_stats_time_series(series, mode="contact")
        self.assertEqual(len(stats_series), 2)
        self.assertEqual(stats_series[0][0], 0.0)
        self.assertEqual(sum(stats_series[0][1]["edge_counts"].values()), 1)
        self.assertEqual(sum(stats_series[1][1]["edge_counts"].values()), 0)


class TestMatLoading(unittest.TestCase):
    """Tests for .mat snapshot loading via the fallback parser."""

    def test_missing_mat_file_raises_file_not_found(self):
        reader = PhysiCellReader()
        with self.assertRaises(FileNotFoundError):
            reader.load_snapshot_mat("/nonexistent/output.mat")

    def test_invalid_mat_file_raises_value_error(self):
        """An unreadable .mat surfaces as ValueError (not NotImplementedError)."""
        reader = PhysiCellReader()
        with tempfile.TemporaryDirectory() as tmp:
            fake_mat = os.path.join(tmp, "fake.mat")
            with open(fake_mat, "wb") as fh:
                fh.write(b"this is not a valid MATLAB file")
            with self.assertRaises(ValueError):
                reader.load_snapshot_mat(fake_mat)

    def test_mat_loading_via_direct_parser_default_layout(self):
        """Synthesize a (6, 10) cells matrix and verify the parser output."""
        from scipy.io import savemat

        reader = PhysiCellReader()
        type_ids = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        matrix = _make_physicell_cells_matrix(10, type_ids, volume=2494.0)

        with tempfile.TemporaryDirectory() as tmp:
            mat_path = os.path.join(tmp, "out_cells.mat")
            savemat(mat_path, {"cells": matrix})
            # Suppress the "no cell-definition mapping" warning here —
            # it's covered separately in test_mat_loading_no_pymcds_no_xml_warns.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cells = reader.load_snapshot_mat(mat_path)

        self.assertEqual(len(cells), 10)
        # First cell: x=0, y=1, z=2 per the synthesizer.
        self.assertAlmostEqual(cells[0]["x"], 0.0)
        self.assertAlmostEqual(cells[0]["y"], 1.0)
        self.assertAlmostEqual(cells[0]["z"], 2.0)
        # Sixth cell (i=5): x=50, y=51, z=52
        self.assertAlmostEqual(cells[5]["x"], 50.0)
        self.assertAlmostEqual(cells[5]["y"], 51.0)
        self.assertAlmostEqual(cells[5]["z"], 52.0)
        # Radius derived from total_volume 2494: ~8.41
        expected_radius = (3.0 * 2494.0 / (4.0 * math.pi)) ** (1.0 / 3.0)
        for cell in cells:
            self.assertAlmostEqual(cell["radius"], expected_radius, places=5)
            self.assertIn(cell["cell_type"], {"0", "1"})  # stringified ints
            self.assertIn("id", cell)
            self.assertIn("volume", cell)

    def test_mat_loading_with_xml_remaps_cell_type(self):
        """When XML is supplied, cell_type IDs become human-readable names."""
        from scipy.io import savemat

        reader = PhysiCellReader()
        type_ids = [0, 1, 0, 1]
        matrix = _make_physicell_cells_matrix(4, type_ids)

        with tempfile.TemporaryDirectory() as tmp:
            mat_path = os.path.join(tmp, "out_cells.mat")
            xml_path = os.path.join(tmp, "out.xml")
            savemat(mat_path, {"cells": matrix})
            _write_cell_definitions_xml(xml_path,
                                        {0: "tumor", 1: "immune"})
            cells = reader.load_snapshot_mat(mat_path, xml_path=xml_path)

        self.assertEqual(len(cells), 4)
        types = [c["cell_type"] for c in cells]
        self.assertEqual(types, ["tumor", "immune", "tumor", "immune"])

    def test_mat_loading_rejects_non_standard_shape(self):
        """A matrix with fewer than 6 rows is rejected, not silently transposed.

        Pair B's earlier transpose-flip heuristic could mis-decode a
        (5, 1000) matrix as (1000 signals, 5 cells); the safer behavior is
        to raise so users notice the format mismatch.
        """
        from scipy.io import savemat

        reader = PhysiCellReader()
        # 4 rows × 10 cols — neither a real PhysiCell shape nor a
        # transposable one. Must reject with a clear error.
        bad = np.zeros((4, 10), dtype=float)

        with tempfile.TemporaryDirectory() as tmp:
            mat_path = os.path.join(tmp, "out_cells.mat")
            savemat(mat_path, {"cells": bad})
            with self.assertRaises(ValueError) as ctx:
                reader.load_snapshot_mat(mat_path)
        self.assertIn("axis 0", str(ctx.exception))

    def test_mat_loading_no_pymcds_no_xml_warns(self):
        """Without XML, parser still works but warns about ID stringification."""
        from scipy.io import savemat

        reader = PhysiCellReader()
        matrix = _make_physicell_cells_matrix(3, [0, 1, 0])

        with tempfile.TemporaryDirectory() as tmp:
            mat_path = os.path.join(tmp, "out_cells.mat")
            savemat(mat_path, {"cells": matrix})
            with self.assertWarns(UserWarning):
                cells = reader.load_snapshot_mat(mat_path)

        self.assertEqual(len(cells), 3)
        # Cell types are stringified integer IDs.
        self.assertEqual([c["cell_type"] for c in cells], ["0", "1", "0"])

    def test_select_cells_matrix_prefers_named_over_heuristic(self):
        """Named 'cells' wins even when a larger sibling matrix is present."""
        reader = PhysiCellReader()
        cells_mat = _make_physicell_cells_matrix(10, [0, 1] * 5)
        # A much larger near-square sibling — would dominate any size-only
        # tiebreaker, but must not be selected because 'cells' is present.
        microenv = np.ones((60, 5000), dtype=float)
        data = {"cells": cells_mat, "microenv": microenv}

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = reader._select_cells_matrix(data, "/tmp/named.mat")

        # The named 'cells' matrix is returned verbatim.
        self.assertIs(result, cells_mat)
        # No heuristic warning was emitted by _select_cells_matrix.
        heuristic_warnings = [
            w for w in caught
            if "shape-aware heuristic" in str(w.message)
        ]
        self.assertEqual(heuristic_warnings, [])

    def test_select_cells_matrix_picks_elongated_under_custom_name(self):
        """Without a standard name, the elongated candidate wins with a warning."""
        reader = PhysiCellReader()
        # Cells-like: 6 signals x 1000 cells (elongation 1000/6 ~= 166).
        cells_like = _make_physicell_cells_matrix(1000, [0] * 1000)
        # Microenv-like: nearly square (elongation 60/50 = 1.2), minor axis
        # 50 still in the plausibility window but with much lower score.
        microenv_like = np.ones((60, 50), dtype=float)

        data = {"my_cells": cells_like, "substrate_grid": microenv_like}

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = reader._select_cells_matrix(data, "/tmp/custom.mat")

        self.assertIs(result, cells_like)
        heuristic_warnings = [
            w for w in caught
            if issubclass(w.category, UserWarning)
            and "shape-aware heuristic" in str(w.message)
        ]
        self.assertEqual(len(heuristic_warnings), 1)
        self.assertIn("my_cells", str(heuristic_warnings[0].message))

    def test_select_cells_matrix_raises_on_ambiguous_heuristic(self):
        """Two custom-named candidates with nearly-equal scores raise ValueError."""
        reader = PhysiCellReader()
        # Same shape -> identical scores; both plausible-window minor axes.
        cand_a = np.ones((6, 1000), dtype=float)
        cand_b = np.ones((6, 1000), dtype=float)
        data = {"weird_a": cand_a, "weird_b": cand_b}

        with self.assertRaises(ValueError) as ctx:
            reader._select_cells_matrix(data, "/tmp/ambiguous.mat")

        msg = str(ctx.exception)
        self.assertIn("weird_a", msg)
        self.assertIn("weird_b", msg)

    def test_select_cells_matrix_no_2d_candidates_raises(self):
        """Files with no 2-D numeric ndarray of sufficient minor axis raise."""
        reader = PhysiCellReader()
        data = {
            "vector": np.arange(10, dtype=float),         # 1-D, rejected
            "tensor": np.zeros((3, 3, 3), dtype=float),   # 3-D, rejected
            "strings": np.array([["a", "b"], ["c", "d"]]),  # non-numeric
            "tiny": np.zeros((2, 100), dtype=float),       # minor axis < 6
        }

        with self.assertRaises(ValueError) as ctx:
            reader._select_cells_matrix(data, "/tmp/empty.mat")

        self.assertIn("No cells matrix", str(ctx.exception))


class TestStatsToTargetStatistics(unittest.TestCase):
    """Adapter from reader stats dict to ReplicateGenerator's TargetStatistics."""

    def _build_stats(self):
        reader = PhysiCellReader()
        cells = [
            {"x": 0.0, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "tumor"},
            {"x": 9.5, "y": 0.0, "z": 0.0, "radius": 5.0, "cell_type": "tumor"},
            {"x": 0.0, "y": 9.5, "z": 0.0, "radius": 5.0, "cell_type": "immune"},
        ]
        return reader.compute_spatial_stats(cells, mode="contact")

    def test_returns_target_statistics_instance(self):
        stats = self._build_stats()
        target = stats_to_target_statistics(stats)
        self.assertIsInstance(target, TargetStatistics)

    def test_interaction_stats_are_populated(self):
        stats = self._build_stats()
        target = stats_to_target_statistics(stats)
        # We expect at least one InteractionStatistics for each unique pair.
        pair_keys = {
            "-".join(sorted([s.type_a, s.type_b]))
            for s in target.interaction_stats
        }
        self.assertIn("immune-tumor", pair_keys)
        self.assertIn("tumor-tumor", pair_keys)
        for s in target.interaction_stats:
            self.assertIsInstance(s, InteractionStatistics)

    def test_cell_type_proportions_sum_to_one(self):
        stats = self._build_stats()
        target = stats_to_target_statistics(stats)
        self.assertAlmostEqual(sum(target.cell_type_proportions.values()), 1.0)

    def test_consumed_by_replicate_generator_validate(self):
        stats = self._build_stats()
        target = stats_to_target_statistics(stats)
        target.validate()


if __name__ == "__main__":
    unittest.main()
