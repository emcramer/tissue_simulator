"""Tests for the size-aware "mechanical" neighbour graph (PhysiCell-style)."""
import asyncio
import json
import warnings

import numpy as np
import pytest

from tissue_simulator import TissueSection
from tissue_simulator.slicing import TissueSlicer
from tissue_simulator.spatial_analysis import (
    SpatialNetworkAnalyzer, MECHANICAL_INTERACTION_FACTOR, CONTACT_TOLERANCE_FACTOR)
from tissue_simulator.replicate_generator import (
    ReplicateGenerator, load_target_statistics_from_tissue)

pytest.importorskip("networkx")

RADII = {"a": (3.5, 5.0), "b": (3.5, 5.0)}


@pytest.fixture(scope="module")
def tissue():
    t = TissueSection(height=80, width=80, thickness=30, cell_radii=RADII, seed=5)
    t.generate_cells(max_attempts=3000, min_spacing=0.0)
    assert len(t.cells) > 60
    return t


def _graph(tissue, **kw):
    an = SpatialNetworkAnalyzer()
    return an.build_network_from_tissue(tissue, **kw)


def _edges(g):
    return {frozenset(e) for e in g.edges()}


def test_default_factor_is_physicell_value():
    assert MECHANICAL_INTERACTION_FACTOR == 1.5


def test_mechanical_at_contact_tolerance_equals_contact(tissue):
    contact = _graph(tissue, mode="contact")
    mech = _graph(tissue, mode="mechanical", interaction_factor=CONTACT_TOLERANCE_FACTOR)
    assert _edges(mech) == _edges(contact)


def test_mechanical_default_is_superset_and_symmetric(tissue):
    contact = _graph(tissue, mode="contact")
    mech = _graph(tissue, mode="mechanical")
    assert _edges(contact) < _edges(mech)
    assert mech.number_of_edges() > contact.number_of_edges()
    assert not any(u == v for u, v in mech.edges())
    for u, v in mech.edges():
        assert mech.has_edge(v, u)
        d = np.linalg.norm(tissue.cells[u].center - tissue.cells[v].center)
        assert d <= 1.5 * (tissue.cells[u].radius + tissue.cells[v].radius) + 1e-9
    assert mech.graph["network_rule"] == {
        "mode": "mechanical", "radius": None, "interaction_factor": 1.5}
    assert contact.graph["network_rule"]["mode"] == "contact"


def test_radius_vs_mechanical_counts(tissue):
    contact = _graph(tissue, mode="contact").number_of_edges()
    mech = _graph(tissue, mode="mechanical").number_of_edges()
    med = float(np.median([c.radius for c in tissue.cells]))
    rad = _graph(tissue, mode="radius", radius=2 * 1.5 * med).number_of_edges()
    print(f"edge counts: contact={contact} mechanical={mech} radius={rad}")
    assert mech > contact and rad > contact
    assert abs(rad - mech) / mech < 0.35  # similar scale ...


def test_unknown_mode_and_bad_factor(tissue):
    with pytest.raises(ValueError):
        _graph(tissue, mode="nope")
    with pytest.raises(ValueError):
        _graph(tissue, mode="mechanical", interaction_factor=0.0)


def test_slice_variant(tissue):
    slicer = TissueSlicer(tissue)
    slicer.slice_plane(z_position=15.0)
    an = SpatialNetworkAnalyzer()
    contact = an.build_network_from_slice(slicer, mode="contact")
    c_edges = _edges(contact)
    mech = an.build_network_from_slice(slicer, mode="mechanical")
    assert c_edges <= _edges(mech)
    assert mech.number_of_edges() > len(c_edges)
    assert mech.graph["network_rule"]["interaction_factor"] == 1.5


def test_from_coordinates_default_rule_and_mismatch_warning(tissue, tmp_path):
    path = str(tmp_path / "t.csv")
    tissue.export_to_csv(path)
    gen = ReplicateGenerator.from_coordinates(path, seed=1)
    expected = {"mode": "mechanical", "radius": None, "interaction_factor": 1.5}
    assert gen.network_mode == "mechanical" and gen.interaction_factor == 1.5
    assert gen.network_rule == expected
    assert gen.target_stats.network_rule == expected

    legacy = ReplicateGenerator.from_coordinates(
        path, network_mode="radius", network_radius=20.0, seed=1)
    assert legacy.network_rule["mode"] == "radius"
    assert legacy.target_stats.network_rule["radius"] == 20.0

    # same rule -> no warning; different rule -> warning
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ReplicateGenerator(gen.target_stats, (80, 80, 30), RADII,
                           network_mode="mechanical", interaction_factor=1.5)
    with pytest.warns(UserWarning, match="network rule"):
        ReplicateGenerator(gen.target_stats, (80, 80, 30), RADII, network_mode="contact")
    with pytest.warns(UserWarning, match="network rule"):
        ReplicateGenerator(gen.target_stats, (80, 80, 30), RADII,
                           network_mode="mechanical", interaction_factor=2.0)


def test_loader_records_rule(tissue):
    t = load_target_statistics_from_tissue(tissue, network_mode="mechanical",
                                           interaction_factor=2.0)
    assert t.network_rule == {"mode": "mechanical", "radius": None,
                              "interaction_factor": 2.0}
    assert load_target_statistics_from_tissue(tissue).network_rule["mode"] == "contact"


def test_mcp_accepts_mechanical(tissue):
    pytest.importorskip("mcp")
    from tissue_simulator.mcp.server import TissueSimulatorMCPServer
    server = TissueSimulatorMCPServer()
    server.current_tissue = tissue
    loop = asyncio.new_event_loop()  # private loop, see test_density_replicates.py
    try:
        data = json.loads(loop.run_until_complete(server._handle_load_target_statistics(
            {"use_current_tissue": True, "network_mode": "mechanical",
             "interaction_factor": 1.25}))[0].text)
        tools = loop.run_until_complete(server.server.request_handlers[
            __import__("mcp.types", fromlist=["ListToolsRequest"]).ListToolsRequest](
            __import__("mcp.types", fromlist=["ListToolsRequest"]).ListToolsRequest(
                method="tools/list")))
    finally:
        loop.close()
    assert data["status"] == "success", data
    assert data["network_rule"]["interaction_factor"] == 1.25
    assert server.interaction_factor == 1.25
    for tool in tools.root.tools:
        props = tool.inputSchema.get("properties", {})
        if "network_mode" in props:
            assert "mechanical" in props["network_mode"]["enum"]
            assert "interaction_factor" in props
