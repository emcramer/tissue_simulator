"""
Tests for the MCP server tools, focusing on the two coordinate-CSV tools:
``load_tissue_from_csv`` and ``load_target_statistics_from_coordinates``.

These tests instantiate ``TissueSimulatorMCPServer`` directly (the ``mcp``
package must be installed) and call the async ``_handle_*`` methods.
"""

import asyncio
import json
import tempfile
import os

import pytest

# Skip the whole module if mcp isn't importable (defensive; it is installed).
mcp = pytest.importorskip("mcp")

from tissue_simulator import TissueSection
from tissue_simulator.mcp.server import TissueSimulatorMCPServer


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
