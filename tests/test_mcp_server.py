"""
Tests for the MCP server tools, focusing on the two coordinate-CSV tools:
``load_tissue_from_csv`` and ``load_target_statistics_from_coordinates``.

These tests instantiate ``TissueSimulatorMCPServer`` directly (the ``mcp``
package must be installed) and call the async ``_handle_*`` methods.
"""

import asyncio
import csv
import json
import math
import tempfile
import os

import pytest

# Skip the whole module if mcp isn't importable (defensive; it is installed).
mcp = pytest.importorskip("mcp")

from tissue_simulator import TissueSection
from tissue_simulator.mcp.server import TissueSimulatorMCPServer
from tissue_simulator.slicing import TissueSlicer


def _call(coro):
    """Run an async coroutine to completion and return its result."""
    return asyncio.get_event_loop().run_until_complete(coro)


def _result_json(text_contents):
    """Parse the JSON payload from a list of TextContent results."""
    return json.loads(text_contents[0].text)


def _make_tissue_csv(path):
    """Build a small packed tissue, export it to CSV, return cell count."""
    tissue = TissueSection(120, 120, 40, {'a': (5, 8), 'b': (6, 10)})
    tissue.generate_cells(max_attempts=400, seed=7)
    tissue.export_to_csv(path)
    return len(tissue.cells)


def test_tools_registered():
    """Both new tools must be registered and have backing handlers."""
    server = TissueSimulatorMCPServer()

    # Both handler methods must exist.
    assert hasattr(server, "_handle_load_tissue_from_csv")
    assert hasattr(server, "_handle_load_target_statistics_from_coordinates")

    # Try to reach the registered list_tools coroutine. The MCP Server stores
    # request handlers keyed by request type; if we can find the ListTools
    # handler, exercise it and assert both tool names are present.
    list_tools_handler = None
    handlers = getattr(server.server, "request_handlers", None)
    if handlers:
        for req_type, handler in handlers.items():
            name = getattr(req_type, "__name__", str(req_type))
            if "ListTools" in name:
                list_tools_handler = handler
                break

    if list_tools_handler is not None:
        # The handler may expect a request object; call defensively.
        try:
            result = _call(list_tools_handler(object()))
            tools = getattr(result.root, "tools", None) or result.tools
        except TypeError:
            result = _call(list_tools_handler())
            tools = getattr(result.root, "tools", None) or result.tools
        tool_names = {t.name for t in tools}
        assert "load_tissue_from_csv" in tool_names
        assert "load_target_statistics_from_coordinates" in tool_names


def test_load_tissue_from_csv_roundtrip():
    """Export a tissue then load it back via the MCP handler."""
    server = TissueSimulatorMCPServer()
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        n_cells = _make_tissue_csv(path)

        result = _result_json(
            _call(server._handle_load_tissue_from_csv({"filepath": path}))
        )

        assert result.get("status") == "success"
        assert result["num_cells"] == n_cells

        # Side effect: current_tissue is populated with that many cells.
        assert server.current_tissue is not None
        assert len(server.current_tissue.cells) == n_cells

        assert 0 < result["packing_fraction"] < 1
    finally:
        os.remove(path)


def test_load_tissue_missing_filepath():
    """Missing filepath must produce an error JSON."""
    server = TissueSimulatorMCPServer()
    result = _result_json(_call(server._handle_load_tissue_from_csv({})))
    assert "error" in result


def test_load_target_statistics_from_coordinates_full():
    """Coordinate stats must populate proportions, density, and side effects."""
    server = TissueSimulatorMCPServer()
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        _make_tissue_csv(path)

        result = _result_json(
            _call(server._handle_load_target_statistics_from_coordinates(
                {"filepath": path, "network_mode": "radius", "network_radius": 20.0}
            ))
        )

        assert result.get("status") == "success"
        assert result["source"] == "coordinate_csv"

        proportions = result["cell_type_proportions"]
        assert proportions is not None
        assert abs(sum(proportions.values()) - 1.0) < 1e-6

        assert result["target_density"] is not None
        assert 0 < result["target_density"] < 1

        assert result["num_interaction_types"] > 0

        # Side effects on the server.
        assert server.target_stats is not None
        assert server.network_mode == "radius"
        assert server.network_radius == 20.0
    finally:
        os.remove(path)


def test_assign_cell_types_handler_with_seed():
    """
    End-to-end MCP handler test for ``assign_cell_types``:

    1. Build a small packed tissue and set it as ``server.current_tissue``.
    2. Write a graph-coloring-format CSV (the exact shape that
       ``tissue_simulator.graph_coloring.load_target_statistics_from_csv``
       reads — ``nodes_<color>`` + ``edges_<c1>-<c2>`` columns).
    3. Call ``_handle_assign_cell_types`` with an explicit ``seed`` and
       custom annealing/network args.
    4. Assert the result JSON shape (``status``, ``num_nodes_colored``,
       ``color_counts``, ``seed``, ``used_z_position``,
       ``used_network_radius``) matches the handler / docs.
    5. Assert the side-effect: ``server.cell_type_assignment`` is populated
       as a dict of the same length.
    """
    server = TissueSimulatorMCPServer()

    # 1. Build a small packed tissue with three cell types so the slice
    #    yields >5 nodes for the coloring step.
    tissue = TissueSection(
        120, 120, 40,
        {'cancer': (5, 8), 'immune': (6, 9), 'stroma': (5, 7)},
        seed=7,
    )
    tissue.generate_cells(max_attempts=400, seed=7)
    assert len(tissue.cells) > 5
    server.current_tissue = tissue

    # 2. Write a target-statistics CSV in the graph-coloring format.
    colors = ['cancer', 'immune', 'stroma']
    fd, csv_path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        # load_target_statistics_from_csv reads the first row and pulls:
        #   nodes_<color> for each color
        #   edges_<c1>-<c2> for each i<=j pair (in the given color order)
        header_cols = [f"nodes_{c}" for c in colors]
        for i in range(len(colors)):
            for j in range(i, len(colors)):
                header_cols.append(f"edges_{colors[i]}-{colors[j]}")
        row_values = [
            # node_counts
            "5", "4", "3",
            # edge_counts (cancer-cancer, cancer-immune, cancer-stroma,
            #              immune-immune, immune-stroma, stroma-stroma)
            "3", "2", "1", "1", "1", "1",
        ]
        with open(csv_path, "w") as f:
            f.write(",".join(header_cols) + "\n")
            f.write(",".join(row_values) + "\n")

        # 3. Call the handler.
        result = _result_json(
            _call(server._handle_assign_cell_types({
                "target_statistics_csv": csv_path,
                "colors": colors,
                "z_position": 20.0,
                "network_radius": 25.0,
                "seed": 42,
                "initial_temp": 10.0,
                "final_temp": 0.5,
                "cooling_rate": 0.95,
                "max_iterations": 200,
                "verbose": False,
            }))
        )

        # 4. Validate the result shape.
        assert result.get("status") == "success", result
        assert result["num_nodes_colored"] > 0
        color_counts = result["color_counts"]
        assert isinstance(color_counts, dict)
        assert sum(color_counts.values()) == result["num_nodes_colored"]
        # Every color in color_counts must be one of the requested colors.
        assert set(color_counts.keys()).issubset(set(colors))
        assert result["seed"] == 42
        assert result["used_z_position"] == 20.0
        assert result["used_network_radius"] == 25.0

        # 5. Side effects on the server.
        assert server.cell_type_assignment is not None
        assert isinstance(server.cell_type_assignment, dict)
        assert len(server.cell_type_assignment) == result["num_nodes_colored"]
    finally:
        os.remove(csv_path)


def test_coordinates_vs_interaction_table_distinct():
    """
    Coordinate path populates proportions/density; the interaction-table path
    (load_target_statistics with csv_filepath) leaves them None.
    """
    server = TissueSimulatorMCPServer()

    # Tiny interaction-table CSV (precomputed stats, not coordinates).
    fd, table_path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        with open(table_path, "w") as f:
            f.write("type_a,type_b,normalized_interactions\n")
            f.write("a,a,0.1\n")

        result = _result_json(
            _call(server._handle_load_target_statistics({"csv_filepath": table_path}))
        )
        assert result.get("status") == "success"

        # Reaching into target_stats is the cleanest assertion here.
        assert server.target_stats is not None
        assert server.target_stats.cell_type_proportions is None
        assert server.target_stats.target_density is None
    finally:
        os.remove(table_path)


# ---------------------------------------------------------------------------
# PhysiCell export handlers (added in v0.1.13).
# ---------------------------------------------------------------------------


def _build_sliced_server():
    """
    Build an MCP server with a small packed tissue + populated slicer.

    Returns the server, the underlying tissue, and the slicer so callers can
    make per-cell assertions against ``server.current_slicer.slice_cells``.
    """
    server = TissueSimulatorMCPServer()
    tissue = TissueSection(
        120, 120, 40,
        {'cancer': (5, 8), 'immune': (6, 9), 'stroma': (5, 7)},
        seed=7,
    )
    tissue.generate_cells(max_attempts=400, seed=7)
    assert len(tissue.cells) > 5
    server.current_tissue = tissue

    slicer = TissueSlicer(tissue)
    slicer.slice_plane(z_position=20.0)
    assert len(slicer.slice_cells) > 0, "slice must contain at least one cell"
    server.current_slicer = slicer
    return server, tissue, slicer


def test_export_tissue_to_physicell():
    """Tissue PhysiCell export writes a header-prefixed modern CSV."""
    server = TissueSimulatorMCPServer()
    tissue = TissueSection(
        120, 120, 40,
        {'a': (5, 8), 'b': (6, 10)},
        seed=7,
    )
    tissue.generate_cells(max_attempts=400, seed=7)
    assert len(tissue.cells) >= 1
    server.current_tissue = tissue

    result = _result_json(
        _call(server._handle_export_tissue_to_physicell({}))
    )

    assert result.get("status") == "success", result
    assert result["num_cells_exported"] == len(server.current_tissue.cells)
    assert result["format"] == "modern"

    filepath = result["filepath"]
    assert os.path.exists(filepath)

    with open(filepath, newline="") as f:
        header_line = f.readline().rstrip("\r\n")
    assert header_line == "x,y,z,type,volume"


def test_export_slice_to_physicell_2d_default():
    """Default slice export is 2D with z=0 and one row per slice cell."""
    server, _tissue, slicer = _build_sliced_server()

    result = _result_json(
        _call(server._handle_export_slice_to_physicell({}))
    )

    assert result.get("status") == "success", result
    assert result["geometry"] == "2D"
    assert result["used_z"] == 0.0
    assert result["format"] == "modern"
    assert result["num_cells_exported"] == len(slicer.slice_cells)

    filepath = result["filepath"]
    assert os.path.exists(filepath)

    with open(filepath, newline="") as f:
        rows = list(csv.reader(f))

    # First row is header, remaining rows are data; one data row per cell.
    assert len(rows) - 1 == result["num_cells_exported"]
    header, *data_rows = rows
    assert header == ["x", "y", "z", "type", "volume"]
    z_index = header.index("z")
    for row in data_rows:
        z_val = float(row[z_index])
        assert z_val == 0.0


def test_export_slice_to_physicell_3d():
    """3D slice export uses original 3D positions and sphere volumes."""
    server, _tissue, slicer = _build_sliced_server()

    result = _result_json(
        _call(server._handle_export_slice_to_physicell({"geometry": "3D"}))
    )

    assert result.get("status") == "success", result
    assert result["geometry"] == "3D"
    assert result["used_z"] is None
    assert result["format"] == "modern"
    assert result["num_cells_exported"] == len(slicer.slice_cells)

    filepath = result["filepath"]
    assert os.path.exists(filepath)

    with open(filepath, newline="") as f:
        rows = list(csv.reader(f))
    header, *data_rows = rows
    assert header == ["x", "y", "z", "type", "volume"]
    assert len(data_rows) == len(slicer.slice_cells)

    z_index = header.index("z")
    volume_index = header.index("volume")

    # Regression guard: 3D mode preserves the source-cell z positions so at
    # least one row should be non-zero (slice plane is z=20).
    nonzero_z = [float(r[z_index]) for r in data_rows if float(r[z_index]) != 0.0]
    assert nonzero_z, "expected at least one non-zero z in 3D-mode export"

    # Volume column must match the sphere volume of the corresponding
    # slice_cells[i].radius. Row order in the export follows iteration of
    # current_slicer.slice_cells in PhysiCellExporter.export_slice.
    for row, sc in zip(data_rows, slicer.slice_cells):
        expected_volume = (4.0 / 3.0) * math.pi * float(sc.radius) ** 3
        assert math.isclose(float(row[volume_index]), expected_volume, rel_tol=1e-9)


def test_export_slice_to_physicell_no_slicer_errors():
    """Calling the slice export with no slicer must return an error JSON."""
    server = TissueSimulatorMCPServer()
    assert server.current_slicer is None

    result = _result_json(
        _call(server._handle_export_slice_to_physicell({}))
    )

    assert "error" in result, result
    assert "slice" in result["error"].lower()
