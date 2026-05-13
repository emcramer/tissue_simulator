"""
Unit tests for the PhysiCell snapshot reader.
"""

import csv
import os
import tempfile
import unittest
from typing import Dict, List

from tissue_simulator.physicell_reader import (
    PhysiCellReader,
    read_physicell_output,
    stats_to_target_statistics,
)
from tissue_simulator.replicate_generator import TargetStatistics
from tissue_simulator.spatial_analysis import InteractionStatistics


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
    """Tests for graceful handling of .mat snapshots."""

    def test_missing_mat_file_raises_file_not_found(self):
        reader = PhysiCellReader()
        with self.assertRaises(FileNotFoundError):
            reader.load_snapshot_mat("/nonexistent/output.mat")

    def test_invalid_mat_file_raises_not_implemented(self):
        reader = PhysiCellReader()
        with tempfile.TemporaryDirectory() as tmp:
            fake_mat = os.path.join(tmp, "fake.mat")
            with open(fake_mat, "wb") as fh:
                fh.write(b"this is not a valid MATLAB file")
            with self.assertRaises(NotImplementedError):
                reader.load_snapshot_mat(fake_mat)


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
