"""
Tests for the PhysiCell initial-condition exporter.
"""

import csv
import math
import os
import tempfile
import unittest

import numpy as np

from tissue_simulator import Cell, TissueSection
from tissue_simulator.slicing import SliceCell
from tissue_simulator.physicell_export import (
    PhysiCellExporter,
    export_to_physicell,
)


def _make_tissue(cells):
    tissue = TissueSection(height=100, width=100, thickness=50, cell_radii=(5, 10))
    tissue.cells = list(cells)
    return tissue


def _make_slice_cell(x, y, cell_type, radius=6.0, intersection_radius=4.0):
    return SliceCell(
        center_3d=np.array([x, y, 25.0]),
        center_2d=np.array([x, y]),
        radius=radius,
        cell_type=cell_type,
        is_boundary=False,
        distance_from_plane=0.0,
        intersection_radius=intersection_radius,
    )


def _read_rows(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


class TestBuildDefaultMapping(unittest.TestCase):
    """Default name-to-id mapping helper."""

    def test_alphabetical_ids_start_at_zero(self):
        mapping = PhysiCellExporter.build_default_mapping(
            ["stroma", "cancer", "immune"]
        )
        self.assertEqual(mapping, {"cancer": 0, "immune": 1, "stroma": 2})

    def test_duplicates_collapsed(self):
        mapping = PhysiCellExporter.build_default_mapping(
            ["a", "b", "a", "c", "b"]
        )
        self.assertEqual(mapping, {"a": 0, "b": 1, "c": 2})


class TestExportTissueModern(unittest.TestCase):
    """3D exports in modern format."""

    def setUp(self):
        self.cells = [
            Cell(center=(10.0, 20.0, 30.0), radius=5.0, cell_type="cancer"),
            Cell(center=(15.0, 25.0, 35.0), radius=7.0, cell_type="immune"),
            Cell(center=(50.0, 50.0, 25.0), radius=6.0, cell_type="cancer"),
        ]
        self.tissue = _make_tissue(self.cells)

    def test_modern_default_mapping(self):
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            returned = exporter.export_tissue(self.tissue, path)
            self.assertEqual(returned, path)
            rows = _read_rows(path)

        self.assertEqual(rows[0], ["x", "y", "z", "type", "volume"])
        self.assertEqual(len(rows), len(self.cells) + 1)
        # Spot-check the first data row matches sphere volume.
        x, y, z, ct, vol = rows[1]
        self.assertEqual(ct, "cancer")
        self.assertAlmostEqual(float(x), 10.0)
        self.assertAlmostEqual(float(y), 20.0)
        self.assertAlmostEqual(float(z), 30.0)
        self.assertAlmostEqual(
            float(vol), (4.0 / 3.0) * math.pi * 5.0 ** 3, places=6
        )

    def test_modern_without_volume(self):
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            exporter.export_tissue(self.tissue, path, include_volume=False)
            rows = _read_rows(path)

        self.assertEqual(rows[0], ["x", "y", "z", "type"])
        for row in rows[1:]:
            self.assertEqual(len(row), 4)

    def test_modern_with_custom_mapping_validates_only(self):
        # Mapping is provided but format=modern still emits names. The
        # mapping is validated for completeness so typos surface early.
        mapping = {"cancer": 7, "immune": 11}
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            exporter.export_tissue(self.tissue, path, cell_type_mapping=mapping)
            rows = _read_rows(path)

        # Names, not ids, are written in modern format.
        types_written = {row[3] for row in rows[1:]}
        self.assertEqual(types_written, {"cancer", "immune"})

    def test_modern_missing_mapping_raises(self):
        mapping = {"cancer": 0}  # immune missing
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            with self.assertRaises(ValueError) as ctx:
                exporter.export_tissue(
                    self.tissue, path, cell_type_mapping=mapping
                )
        self.assertIn("immune", str(ctx.exception))


class TestExportTissueLegacy(unittest.TestCase):
    """3D exports in legacy headerless integer format."""

    def setUp(self):
        self.cells = [
            Cell(center=(10.0, 20.0, 30.0), radius=5.0, cell_type="cancer"),
            Cell(center=(15.0, 25.0, 35.0), radius=7.0, cell_type="immune"),
        ]
        self.tissue = _make_tissue(self.cells)

    def test_legacy_with_explicit_mapping(self):
        mapping = {"cancer": 1, "immune": 2}
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            exporter.export_tissue(
                self.tissue, path, cell_type_mapping=mapping, format="legacy"
            )
            rows = _read_rows(path)

        self.assertEqual(len(rows), 2)
        # No header; first row should be data.
        x, y, z, ct_id, vol = rows[0]
        self.assertEqual(int(ct_id), 1)
        self.assertAlmostEqual(float(x), 10.0)
        self.assertAlmostEqual(
            float(vol), (4.0 / 3.0) * math.pi * 5.0 ** 3, places=6
        )

        x2, y2, z2, ct_id2, vol2 = rows[1]
        self.assertEqual(int(ct_id2), 2)

    def test_legacy_requires_mapping(self):
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            with self.assertRaises(ValueError):
                exporter.export_tissue(self.tissue, path, format="legacy")

    def test_legacy_missing_type_in_mapping_raises(self):
        mapping = {"cancer": 1}  # immune missing
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            with self.assertRaises(ValueError) as ctx:
                exporter.export_tissue(
                    self.tissue, path, cell_type_mapping=mapping, format="legacy"
                )
        self.assertIn("immune", str(ctx.exception))

    def test_legacy_rejects_include_volume_false(self):
        mapping = {"cancer": 1, "immune": 2}
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ic.csv")
            with self.assertRaises(ValueError):
                exporter.export_tissue(
                    self.tissue,
                    path,
                    cell_type_mapping=mapping,
                    format="legacy",
                    include_volume=False,
                )


class TestExportSlice(unittest.TestCase):
    """2D slice exports."""

    def setUp(self):
        self.slice_cells = [
            _make_slice_cell(10.0, 20.0, "cancer", intersection_radius=3.0),
            _make_slice_cell(15.0, 25.0, "immune", intersection_radius=4.5),
            _make_slice_cell(30.0, 30.0, "cancer", intersection_radius=2.5),
        ]

    def test_slice_modern_with_volume_column(self):
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "slice.csv")
            exporter.export_slice(self.slice_cells, path, z=0.0)
            rows = _read_rows(path)

        self.assertEqual(rows[0], ["x", "y", "z", "type", "volume"])
        self.assertEqual(len(rows), len(self.slice_cells) + 1)
        # Verify 2D disk-area convention.
        _x, _y, z, _ct, vol = rows[1]
        self.assertAlmostEqual(float(z), 0.0)
        self.assertAlmostEqual(float(vol), math.pi * 3.0 ** 2, places=6)

    def test_slice_legacy(self):
        mapping = {"cancer": 0, "immune": 1}
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "slice.csv")
            exporter.export_slice(
                self.slice_cells,
                path,
                z=5.5,
                cell_type_mapping=mapping,
                format="legacy",
            )
            rows = _read_rows(path)

        self.assertEqual(len(rows), len(self.slice_cells))
        for row, sc in zip(rows, self.slice_cells):
            self.assertEqual(int(row[3]), mapping[sc.cell_type])
            self.assertAlmostEqual(float(row[2]), 5.5)

    def test_slice_no_volume_column(self):
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "slice.csv")
            exporter.export_slice(self.slice_cells, path, include_volume=False)
            rows = _read_rows(path)

        self.assertEqual(rows[0], ["x", "y", "z", "type"])
        for row in rows[1:]:
            self.assertEqual(len(row), 4)


class TestConvenienceDispatch(unittest.TestCase):
    """Top-level ``export_to_physicell`` dispatches based on input type."""

    def test_dispatch_tissue(self):
        tissue = _make_tissue(
            [Cell(center=(1.0, 2.0, 3.0), radius=4.0, cell_type="a")]
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.csv")
            export_to_physicell(tissue, path)
            rows = _read_rows(path)
        self.assertEqual(rows[0][:4], ["x", "y", "z", "type"])
        self.assertEqual(rows[1][3], "a")

    def test_dispatch_slice(self):
        slice_cells = [_make_slice_cell(1.0, 2.0, "a")]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.csv")
            export_to_physicell(slice_cells, path, z=1.25)
            rows = _read_rows(path)
        self.assertEqual(rows[0][:4], ["x", "y", "z", "type"])
        self.assertAlmostEqual(float(rows[1][2]), 1.25)

    def test_dispatch_rejects_unknown_type(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.csv")
            with self.assertRaises(TypeError):
                export_to_physicell([object()], path)


class TestRoundTrip(unittest.TestCase):
    """Read back the CSV and check counts and values."""

    def test_round_trip_modern(self):
        cells = [
            Cell(center=(1.5, 2.5, 3.5), radius=4.0, cell_type="cancer"),
            Cell(center=(11.5, 12.5, 13.5), radius=5.0, cell_type="immune"),
            Cell(center=(21.5, 22.5, 23.5), radius=6.0, cell_type="stroma"),
        ]
        tissue = _make_tissue(cells)
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "round.csv")
            exporter.export_tissue(tissue, path)

            with open(path, newline="") as f:
                reader = csv.DictReader(f)
                read_rows = list(reader)

        self.assertEqual(len(read_rows), len(cells))
        for cell, row in zip(cells, read_rows):
            self.assertAlmostEqual(float(row["x"]), cell.center[0])
            self.assertAlmostEqual(float(row["y"]), cell.center[1])
            self.assertAlmostEqual(float(row["z"]), cell.center[2])
            self.assertEqual(row["type"], cell.cell_type)
            self.assertAlmostEqual(
                float(row["volume"]),
                (4.0 / 3.0) * math.pi * cell.radius ** 3,
                places=6,
            )

    def test_round_trip_legacy(self):
        cells = [
            Cell(center=(1.0, 2.0, 3.0), radius=4.0, cell_type="cancer"),
            Cell(center=(11.0, 12.0, 13.0), radius=5.0, cell_type="immune"),
        ]
        tissue = _make_tissue(cells)
        mapping = PhysiCellExporter.build_default_mapping(
            [c.cell_type for c in cells]
        )
        exporter = PhysiCellExporter()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "round_legacy.csv")
            exporter.export_tissue(
                tissue, path, cell_type_mapping=mapping, format="legacy"
            )
            rows = _read_rows(path)

        self.assertEqual(len(rows), len(cells))
        for cell, row in zip(cells, rows):
            self.assertAlmostEqual(float(row[0]), cell.center[0])
            self.assertAlmostEqual(float(row[1]), cell.center[1])
            self.assertAlmostEqual(float(row[2]), cell.center[2])
            self.assertEqual(int(row[3]), mapping[cell.cell_type])
            self.assertAlmostEqual(
                float(row[4]),
                (4.0 / 3.0) * math.pi * cell.radius ** 3,
                places=6,
            )


if __name__ == "__main__":
    unittest.main()
