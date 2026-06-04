"""
PhysiCell initial-condition CSV exporter.

This module bridges ``tissue_simulator`` outputs (:class:`TissueSection`
instances and 2D slice cell lists from :class:`TissueSlicer`) into the
initial-condition (IC) CSV format consumed by PhysiCell.

PhysiCell historically accepted two CSV layouts for seeding initial cell
positions:

* **legacy**: headerless rows of ``x,y,z,cell_type_id,volume`` where
  ``cell_type_id`` is an integer key referencing a cell definition in the
  XML configuration. ``volume`` is the cell volume in cubic micrometers.
* **modern**: a header row of ``x,y,z,type`` (with ``type`` as the
  string name defined in the XML), optionally followed by a ``volume``
  column. These column names match the canonical PhysiCell 1.10+
  ``config/cells.csv`` loader.

For 2D PhysiCell models the ``z`` coordinate is fixed to a constant
(typically 0). PhysiCell's documentation is ambiguous about whether the
volume column for 2D models should be a surface area or a volume; this
exporter writes ``pi * r ** 2`` for 2D slices so the value can be
interpreted as either a 2D cross-sectional area or, equivalently, the
"volume" of a unit-thickness cell. Downstream users should override the
mapping if a different convention is required.

:meth:`PhysiCellExporter.export_slice` supports both 2D and 3D output
via its ``geometry`` parameter, so a slice can either be flattened to
the slice plane (for 2D PhysiCell models) or emitted at the cells'
original 3D positions (for 3D PhysiCell models where the slice acts as
a spatial filter).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Sequence, Union

from .slicing import SliceCell
from .tissue import Cell, TissueSection


_MODERN = "modern"
_LEGACY = "legacy"
_VALID_FORMATS = (_MODERN, _LEGACY)

_GEOM_2D = "2D"
_GEOM_3D = "3D"
_VALID_GEOMETRIES = (_GEOM_2D, _GEOM_3D)


@dataclass
class _Row:
    """Internal flat representation of a row before serialization."""

    x: float
    y: float
    z: float
    cell_type: str
    volume: float


class PhysiCellExporter:
    """
    Export tissue and slice data as PhysiCell initial-condition CSV files.

    The exporter supports both the legacy headerless integer-typed CSV
    format and the modern named-cell-type format. Cell volumes are derived
    from the radius using a sphere volume in 3D and a disk area in 2D.

    Parameters
    ----------
    default_format : str, optional
        Default CSV format to use when an explicit ``format`` argument is
        not provided. Must be either ``"modern"`` or ``"legacy"``.

    Notes
    -----
    PhysiCell expects all initial cells to lie inside the configured
    simulation domain. This exporter does not clip or shift coordinates;
    callers should align ``tissue_simulator`` bounds to PhysiCell's domain
    upstream of export.
    """

    def __init__(self, default_format: str = _MODERN):
        if default_format not in _VALID_FORMATS:
            raise ValueError(
                f"default_format must be one of {_VALID_FORMATS}, got "
                f"{default_format!r}"
            )
        self.default_format = default_format

    @staticmethod
    def build_default_mapping(cell_types: Sequence[str]) -> Dict[str, int]:
        """
        Build a stable name-to-id mapping for PhysiCell cell types.

        Parameters
        ----------
        cell_types : sequence of str
            Cell type names to include in the mapping. Duplicates are
            collapsed.

        Returns
        -------
        dict of str to int
            Mapping where ids are assigned alphabetically starting at 0.
            Stability under alphabetical sort makes legacy CSVs
            reproducible regardless of input ordering.
        """
        unique = sorted({str(name) for name in cell_types})
        return {name: idx for idx, name in enumerate(unique)}

    def export_tissue(
        self,
        tissue: TissueSection,
        output_path: str,
        cell_type_mapping: Optional[Dict[str, int]] = None,
        include_volume: bool = True,
        format: str = _MODERN,
    ) -> str:
        """
        Export a 3D :class:`TissueSection` to a PhysiCell IC CSV.

        Parameters
        ----------
        tissue : TissueSection
            Source tissue. Its ``cells`` attribute provides the rows to
            export.
        output_path : str
            Destination CSV path. Parent directories must already exist.
        cell_type_mapping : dict of str to int, optional
            Mapping from cell type name to PhysiCell integer id. Required
            for ``format="legacy"``; if provided for ``"modern"`` it is
            still validated so that callers can detect typos early.
        include_volume : bool, default True
            For ``"modern"`` format, append a ``volume`` column.
            For ``"legacy"`` format the volume column is mandatory and
            ``include_volume=False`` raises ``ValueError``.
        format : str, default "modern"
            Output schema. One of ``"modern"`` or ``"legacy"``.

        Returns
        -------
        str
            The ``output_path`` that was written.
        """
        fmt = self._validate_format(format)
        rows = [
            _Row(
                x=float(cell.center[0]),
                y=float(cell.center[1]),
                z=float(cell.center[2]),
                cell_type=str(cell.cell_type),
                volume=self._sphere_volume(cell.radius),
            )
            for cell in tissue.cells
        ]
        self._write(rows, output_path, cell_type_mapping, include_volume, fmt)
        return output_path

    def export_slice(
        self,
        slice_cells: Iterable[SliceCell],
        output_path: str,
        geometry: str = "2D",
        z: float = 0.0,
        cell_type_mapping: Optional[Dict[str, int]] = None,
        include_volume: bool = True,
        format: str = _MODERN,
    ) -> str:
        """
        Export a slice as a PhysiCell IC CSV.

        Parameters
        ----------
        slice_cells : iterable of SliceCell
            Cells extracted by :class:`TissueSlicer`.
        output_path : str
            Destination CSV path.
        geometry : str, default "2D"
            Output geometry mode. One of ``"2D"`` or ``"3D"``.

            * ``"2D"`` (default): the slice-plane projection. Uses each
              cell's ``center_2d`` (its position in the 2D slice frame),
              the supplied ``z`` coordinate (default ``0``), and the
              intersection-circle area (``pi * intersection_radius ** 2``)
              for the volume column. Suited to PhysiCell 2D simulations.
            * ``"3D"``: the original 3D geometry. Uses each cell's
              ``center_3d`` (its true 3D center) and sphere volume
              (``(4/3) * pi * radius ** 3``) on the cell's original 3D
              radius. Suited to PhysiCell 3D simulations where the slice
              is used only as a spatial filter. The ``z`` parameter is
              ignored in this mode.
        z : float, default 0.0
            Fixed ``z`` coordinate written for every row when
            ``geometry="2D"``. PhysiCell 2D models conventionally fix
            ``z=0``. Ignored when ``geometry="3D"``.
        cell_type_mapping : dict of str to int, optional
            Same semantics as :meth:`export_tissue`.
        include_volume : bool, default True
            Append a ``volume`` column in modern format. Required (and
            therefore must be left True) in legacy format.
        format : str, default "modern"
            Output schema.

        Returns
        -------
        str
            The ``output_path`` that was written.

        Notes
        -----
        When ``geometry="3D"`` the supplied ``z`` argument is ignored;
        the ``z`` column is populated from each cell's ``center_3d[2]``.
        """
        geom = self._validate_geometry(geometry)
        fmt = self._validate_format(format)
        if geom == _GEOM_2D:
            rows = [
                _Row(
                    x=float(sc.center_2d[0]),
                    y=float(sc.center_2d[1]),
                    z=float(z),
                    cell_type=str(sc.cell_type),
                    volume=self._disk_area(sc.intersection_radius),
                )
                for sc in slice_cells
            ]
        else:  # _GEOM_3D
            rows = [
                _Row(
                    x=float(sc.center_3d[0]),
                    y=float(sc.center_3d[1]),
                    z=float(sc.center_3d[2]),
                    cell_type=str(sc.cell_type),
                    volume=self._sphere_volume(sc.radius),
                )
                for sc in slice_cells
            ]
        self._write(rows, output_path, cell_type_mapping, include_volume, fmt)
        return output_path

    @staticmethod
    def _sphere_volume(radius: float) -> float:
        return (4.0 / 3.0) * math.pi * float(radius) ** 3

    @staticmethod
    def _disk_area(radius: float) -> float:
        return math.pi * float(radius) ** 2

    @staticmethod
    def _validate_format(format: str) -> str:
        if format not in _VALID_FORMATS:
            raise ValueError(
                f"format must be one of {_VALID_FORMATS}, got {format!r}"
            )
        return format

    @staticmethod
    def _validate_geometry(geometry: str) -> str:
        if geometry not in _VALID_GEOMETRIES:
            raise ValueError(
                f"geometry must be one of {_VALID_GEOMETRIES}, got "
                f"{geometry!r}"
            )
        return geometry

    @staticmethod
    def _validate_mapping(
        rows: Sequence[_Row],
        cell_type_mapping: Optional[Dict[str, int]],
        require: bool,
    ) -> None:
        if cell_type_mapping is None:
            if require:
                raise ValueError(
                    "cell_type_mapping is required when format='legacy'. "
                    "Use PhysiCellExporter.build_default_mapping(...) or "
                    "supply your own name-to-id dict."
                )
            return
        observed = {row.cell_type for row in rows}
        missing = sorted(observed - set(cell_type_mapping.keys()))
        if missing:
            raise ValueError(
                "cell_type_mapping is missing entries for cell types: "
                f"{missing}. Known mapping keys: "
                f"{sorted(cell_type_mapping.keys())}."
            )

    def _write(
        self,
        rows: Sequence[_Row],
        output_path: str,
        cell_type_mapping: Optional[Dict[str, int]],
        include_volume: bool,
        format: str,
    ) -> None:
        require_mapping = format == _LEGACY
        self._validate_mapping(rows, cell_type_mapping, require=require_mapping)
        if format == _LEGACY and not include_volume:
            raise ValueError(
                "include_volume=False is not supported for format='legacy'; "
                "the legacy PhysiCell IC schema requires a volume column."
            )

        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if format == _MODERN:
                self._write_modern(writer, rows, cell_type_mapping, include_volume)
            else:
                self._write_legacy(writer, rows, cell_type_mapping)

    @staticmethod
    def _write_modern(
        writer: "csv._writer",
        rows: Sequence[_Row],
        cell_type_mapping: Optional[Dict[str, int]],
        include_volume: bool,
    ) -> None:
        header = ["x", "y", "z", "type"]
        if include_volume:
            header.append("volume")
        writer.writerow(header)

        for row in rows:
            # Modern format keeps cell type names; mapping is validated only.
            record = [row.x, row.y, row.z, row.cell_type]
            if include_volume:
                record.append(row.volume)
            writer.writerow(record)

    @staticmethod
    def _write_legacy(
        writer: "csv._writer",
        rows: Sequence[_Row],
        cell_type_mapping: Dict[str, int],
    ) -> None:
        for row in rows:
            writer.writerow(
                [
                    row.x,
                    row.y,
                    row.z,
                    cell_type_mapping[row.cell_type],
                    row.volume,
                ]
            )


def export_to_physicell(
    tissue_or_slice_cells: Union[TissueSection, Iterable[SliceCell]],
    output_path: str,
    **kwargs,
) -> str:
    """
    Dispatch to the appropriate :class:`PhysiCellExporter` method.

    Parameters
    ----------
    tissue_or_slice_cells : TissueSection or iterable of SliceCell
        The data to export. If a :class:`TissueSection` is provided,
        :meth:`PhysiCellExporter.export_tissue` is invoked; otherwise the
        input is treated as a sequence of :class:`SliceCell` objects and
        passed to :meth:`PhysiCellExporter.export_slice`.
    output_path : str
        Destination CSV path.
    **kwargs
        Forwarded to the chosen export method (e.g. ``cell_type_mapping``,
        ``include_volume``, ``format``, and ``z`` for slices).

    Returns
    -------
    str
        The ``output_path`` that was written.
    """
    exporter = PhysiCellExporter()
    if isinstance(tissue_or_slice_cells, TissueSection):
        return exporter.export_tissue(tissue_or_slice_cells, output_path, **kwargs)

    # A list-of-SliceCell is the only other supported input. Materialize the
    # iterable so we can introspect it without consuming a generator twice.
    slice_cells = list(tissue_or_slice_cells)
    if slice_cells and not isinstance(slice_cells[0], SliceCell):
        raise TypeError(
            "export_to_physicell expected a TissueSection or an iterable of "
            f"SliceCell; got element of type {type(slice_cells[0]).__name__}."
        )
    return exporter.export_slice(slice_cells, output_path, **kwargs)
