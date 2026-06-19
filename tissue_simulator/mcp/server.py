"""
Model Context Protocol (MCP) server for Tissue Simulator.

This module exposes tissue simulation functionality as MCP tools that can be
called by Large Language Models (LLMs) through the MCP protocol.
"""

import json
import asyncio
from typing import Any, Dict, List, Optional
import tempfile
import os
from pathlib import Path

# MCP server imports
try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import Tool, TextContent, ImageContent, EmbeddedResource
    MCP_AVAILABLE = True
except ImportError:
    MCP_AVAILABLE = False
    print("Warning: MCP library not installed. Install with: pip install mcp")

from ..tissue import TissueSection, load_tissue_from_csv
from ..slicing import TissueSlicer, create_standard_slices
from ..spatial_analysis import SpatialNetworkAnalyzer
from ..graph_coloring import (GraphColorizer,
                               load_target_statistics_from_csv as load_graph_coloring_stats_from_csv)
from ..replicate_generator import (ReplicateGenerator, TargetStatistics,
                                     load_target_statistics_from_csv,
                                     load_target_statistics_from_tissue,
                                     load_target_statistics_from_coordinates)
from ..physicell_export import PhysiCellExporter


class TissueSimulatorMCPServer:
    """
    MCP Server for Tissue Simulator.
    
    Exposes tissue generation, slicing, and analysis capabilities as MCP tools.
    """
    
    def __init__(self):
        """Initialize the MCP server."""
        if not MCP_AVAILABLE:
            raise ImportError("MCP library not installed. Install with: pip install mcp")
        
        self.server = Server("tissue-simulator")
        self.current_tissue: Optional[TissueSection] = None
        self.current_slicer: Optional[TissueSlicer] = None
        self.replicate_generator: Optional[ReplicateGenerator] = None
        self.generated_replicates: List[Tuple[TissueSection, Any]] = []
        self.temp_dir = tempfile.mkdtemp(prefix="tissue_sim_")
        
        # Register tools
        self._register_tools()
        
        # Register tool handlers
        self._register_handlers()
    
    def _register_tools(self):
        """Register available tools with the MCP server."""
        
        @self.server.list_tools()
        async def list_tools() -> List[Tool]:
            """List all available tissue simulator tools."""
            return [
                Tool(
                    name="create_tissue",
                    description=(
                        "Create a 3D tissue section with specified dimensions and cell types. "
                        "This generates the tissue structure but does not populate it with cells yet. "
                        "Use generate_cells after this to populate the tissue."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "height": {
                                "type": "number",
                                "description": "Height of tissue in micrometers (Y-dimension)",
                                "minimum": 50,
                                "maximum": 2000
                            },
                            "width": {
                                "type": "number",
                                "description": "Width of tissue in micrometers (X-dimension)",
                                "minimum": 50,
                                "maximum": 2000
                            },
                            "thickness": {
                                "type": "number",
                                "description": "Thickness of tissue in micrometers (Z-dimension)",
                                "minimum": 20,
                                "maximum": 500
                            },
                            "cell_types": {
                                "type": "object",
                                "description": "Dictionary mapping cell type names to [min_radius, max_radius] in micrometers",
                                "additionalProperties": {
                                    "type": "array",
                                    "items": {"type": "number"},
                                    "minItems": 2,
                                    "maxItems": 2
                                }
                            }
                        },
                        "required": ["height", "width", "thickness", "cell_types"]
                    }
                ),
                
                Tool(
                    name="generate_cells",
                    description=(
                        "Populate the tissue with cells using random sphere packing. "
                        "Must be called after create_tissue. Returns the number of cells generated."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "max_attempts": {
                                "type": "integer",
                                "description": "Maximum failed placement attempts before stopping",
                                "default": 1000,
                                "minimum": 100,
                                "maximum": 10000
                            },
                            "min_spacing": {
                                "type": "number",
                                "description": "Minimum spacing between cell surfaces in micrometers",
                                "default": 0.5,
                                "minimum": 0,
                                "maximum": 10
                            },
                            "allow_boundary_cells": {
                                "type": "boolean",
                                "description": "Allow cells that extend beyond tissue bounds",
                                "default": True
                            }
                        }
                    }
                ),
                
                Tool(
                    name="get_tissue_statistics",
                    description=(
                        "Get comprehensive statistics about the generated tissue including "
                        "cell counts, type distribution, packing efficiency, and average radii."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {}
                    }
                ),
                
                Tool(
                    name="create_slice",
                    description=(
                        "Create a 2D slice through the 3D tissue at any angle. "
                        "Can specify horizontal slice by z-position, angled slice by rotation angles, "
                        "or custom slice by normal vector."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "z_position": {
                                "type": "number",
                                "description": "Z-position for horizontal slice (simple method)"
                            },
                            "angle_x": {
                                "type": "number",
                                "description": "Rotation angle around X-axis in degrees",
                                "default": 0
                            },
                            "angle_y": {
                                "type": "number",
                                "description": "Rotation angle around Y-axis in degrees",
                                "default": 0
                            },
                            "point": {
                                "type": "array",
                                "description": "Point on the slice plane [x, y, z]",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3
                            },
                            "normal": {
                                "type": "array",
                                "description": "Normal vector to slice plane [nx, ny, nz]",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3
                            }
                        }
                    }
                ),
                
                Tool(
                    name="get_slice_statistics",
                    description=(
                        "Get statistics about the current slice including number of cells captured, "
                        "cell type distribution, and distance metrics."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {}
                    }
                ),
                
                Tool(
                    name="create_serial_slices",
                    description=(
                        "Create multiple evenly-spaced parallel slices through the tissue. "
                        "Useful for histology simulation and 3D reconstruction."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "num_slices": {
                                "type": "integer",
                                "description": "Number of parallel slices to create",
                                "default": 5,
                                "minimum": 2,
                                "maximum": 20
                            }
                        },
                        "required": ["num_slices"]
                    }
                ),
                
                Tool(
                    name="export_tissue_csv",
                    description=(
                        "Export the 3D tissue data to a CSV file. "
                        "Returns the file path where data was saved."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filename": {
                                "type": "string",
                                "description": "Name for the CSV file (without path)",
                                "default": "tissue_data.csv"
                            }
                        }
                    }
                ),

                Tool(
                    name="export_tissue_to_physicell",
                    description=(
                        "Export the current 3D tissue as a PhysiCell IC CSV. The current tissue "
                        "must be created first (via `create_tissue` + `generate_cells` or "
                        "`load_tissue_from_csv`). Writes the file to the server's temp directory "
                        "and returns its path. Tissue export is always 3D — for a 2D PhysiCell "
                        "simulation, slice the tissue first and use `export_slice_to_physicell` "
                        "with `geometry=\"2D\"`."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filename": {
                                "type": "string",
                                "description": "Name for the PhysiCell IC CSV file (without path)",
                                "default": "tissue_physicell.csv"
                            },
                            "format": {
                                "type": "string",
                                "description": "PhysiCell CSV schema: 'modern' (header x,y,z,type [,volume]) or 'legacy' (headerless x,y,z,cell_type_id,volume)",
                                "enum": ["modern", "legacy"],
                                "default": "modern"
                            },
                            "include_volume": {
                                "type": "boolean",
                                "description": "Append the volume column. Ignored for 'legacy' format (volume is always required there).",
                                "default": True
                            },
                            "cell_type_mapping": {
                                "type": "object",
                                "description": "Mapping from cell-type name to PhysiCell integer ID. Required when format='legacy'.",
                                "additionalProperties": {"type": "integer"}
                            }
                        }
                    }
                ),

                Tool(
                    name="load_tissue_from_csv",
                    description=(
                        "Load a TissueSection from a coordinate CSV file. "
                        "This is the inverse of export_tissue_csv: the file must have "
                        "columns x,y,z,radius,cell_type,is_boundary. The loaded tissue is "
                        "set as the current tissue so slicing, analysis, and statistics tools "
                        "can run on externally sourced tissue. Note this is DISTINCT from "
                        "load_target_statistics's csv_filepath, which loads a precomputed "
                        "interaction table rather than per-cell coordinates."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filepath": {
                                "type": "string",
                                "description": "Path to the coordinate CSV file (columns x,y,z,radius,cell_type,is_boundary)"
                            },
                            "height": {
                                "type": "number",
                                "description": "Tissue height (Y) in micrometers; if omitted it is inferred from the coordinate bounds"
                            },
                            "width": {
                                "type": "number",
                                "description": "Tissue width (X) in micrometers; if omitted it is inferred from the coordinate bounds"
                            },
                            "thickness": {
                                "type": "number",
                                "description": "Tissue thickness (Z) in micrometers; if omitted it is inferred from the coordinate bounds"
                            },
                            "default_radius": {
                                "type": "number",
                                "description": "Radius to use for cells missing a radius value",
                                "default": 10.0
                            }
                        },
                        "required": ["filepath"]
                    }
                ),

                Tool(
                    name="load_target_statistics_from_coordinates",
                    description=(
                        "Compute FULL target statistics directly from a coordinate CSV "
                        "(columns x,y,z,radius,cell_type,is_boundary). Unlike "
                        "load_target_statistics with csv_filepath -- which reads a precomputed "
                        "interaction table and does NOT populate cell_type_proportions or "
                        "target_density -- this builds the spatial network from the coordinates "
                        "and produces interactions PLUS cell_type_proportions and target_density. "
                        "These statistics will be used to generate replicate tissues."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filepath": {
                                "type": "string",
                                "description": "Path to the coordinate CSV file (columns x,y,z,radius,cell_type,is_boundary)"
                            },
                            "network_mode": {
                                "type": "string",
                                "description": "Network analysis mode: 'contact' or 'radius'",
                                "default": "contact"
                            },
                            "network_radius": {
                                "type": "number",
                                "description": "Distance threshold for 'radius' mode (micrometers)"
                            }
                        },
                        "required": ["filepath"]
                    }
                ),

                Tool(
                    name="assign_cell_types",
                    description=(
                        "Assign cell types to the cells of the current tissue by simulated "
                        "annealing against a target adjacency structure. Slices the current "
                        "tissue at z_position, builds a contact/radius network, loads target "
                        "statistics from a graph-coloring-format CSV, and runs "
                        "GraphColorizer.colorize() with the given seed. Stores the per-node "
                        "coloring on the server (also accessible via subsequent tools)."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "target_statistics_csv": {
                                "type": "string",
                                "description": (
                                    "Path to a graph-coloring-format CSV (with node_counts, "
                                    "edge_counts, neighbor_dist rows — the same shape that "
                                    "tissue_simulator.graph_coloring.load_target_statistics_from_csv "
                                    "reads)."
                                )
                            },
                            "colors": {
                                "type": "array",
                                "description": (
                                    "List of cell-type names (e.g. [\"cancer\", \"immune\", "
                                    "\"stroma\"])."
                                ),
                                "items": {"type": "string"},
                                "minItems": 1
                            },
                            "z_position": {
                                "type": "number",
                                "description": (
                                    "Z-plane for the horizontal slice; defaults to "
                                    "tissue.thickness / 2 if omitted."
                                )
                            },
                            "network_radius": {
                                "type": "number",
                                "description": (
                                    "Distance threshold for the spatial network (micrometers). "
                                    "Defaults to the server's stored network_radius if set, "
                                    "otherwise 50.0."
                                )
                            },
                            "seed": {
                                "type": "integer",
                                "description": (
                                    "RNG seed for the simulated annealing; makes colorize "
                                    "bit-reproducible."
                                )
                            },
                            "initial_temp": {
                                "type": "number",
                                "description": "Starting temperature for simulated annealing",
                                "default": 100.0
                            },
                            "final_temp": {
                                "type": "number",
                                "description": "Stopping temperature for simulated annealing",
                                "default": 0.1
                            },
                            "cooling_rate": {
                                "type": "number",
                                "description": "Temperature decrease rate per step",
                                "default": 0.995
                            },
                            "max_iterations": {
                                "type": "integer",
                                "description": "Maximum number of simulated-annealing iterations",
                                "default": 5000
                            },
                            "verbose": {
                                "type": "boolean",
                                "description": "Print per-iteration progress",
                                "default": False
                            },
                            "warm_start": {
                                "type": "boolean",
                                "description": (
                                    "When true, reuse the server's previously stored "
                                    "cell-type assignment as the initial coloring for this "
                                    "run, IF a prior assignment exists AND its node-id set "
                                    "exactly equals the new graph's node set; otherwise it "
                                    "is ignored."
                                ),
                                "default": False
                            }
                        },
                        "required": ["target_statistics_csv", "colors"]
                    }
                ),

                Tool(
                    name="generate_colored_replicates",
                    description=(
                        "Generate multiple independent cell-type colorings of the current "
                        "tissue slice that all match the same target statistics. Slices the "
                        "tissue at z_position, builds a contact/radius network, loads target "
                        "statistics from a graph-coloring-format CSV, then runs "
                        "GraphColorizer.colorize() num_replicates times. Each replicate "
                        "cold-starts from its own derived seed by default, producing diverse "
                        "but statistically-equivalent labelings. Returns per-replicate color "
                        "counts plus a mean pairwise diversity score; stores the first "
                        "replicate as the active assignment."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "target_statistics_csv": {
                                "type": "string",
                                "description": (
                                    "Path to a graph-coloring-format CSV (same shape that "
                                    "tissue_simulator.graph_coloring.load_target_statistics_from_csv "
                                    "reads)."
                                )
                            },
                            "colors": {
                                "type": "array",
                                "description": "List of cell-type names.",
                                "items": {"type": "string"},
                                "minItems": 1
                            },
                            "num_replicates": {
                                "type": "integer",
                                "description": "Number of colorings to generate (>= 1).",
                                "minimum": 1
                            },
                            "z_position": {
                                "type": "number",
                                "description": (
                                    "Z-plane for the horizontal slice; defaults to "
                                    "tissue.thickness / 2 if omitted."
                                )
                            },
                            "network_radius": {
                                "type": "number",
                                "description": (
                                    "Distance threshold for the spatial network (micrometers). "
                                    "Defaults to the server's stored network_radius if set, "
                                    "otherwise 50.0."
                                )
                            },
                            "seed": {
                                "type": "integer",
                                "description": (
                                    "Base RNG seed; replicate k uses a deterministic seed "
                                    "derived from (seed, k), making the whole set reproducible."
                                )
                            },
                            "warm_start": {
                                "type": "boolean",
                                "description": (
                                    "When true, each replicate after the first warm-starts "
                                    "from the previous replicate's coloring. This speeds "
                                    "convergence but collapses diversity (replicates become "
                                    "near-identical); intended for refinement, not independent "
                                    "replicates."
                                ),
                                "default": False
                            },
                            "initial_temp": {
                                "type": "number",
                                "description": "Starting temperature for simulated annealing",
                                "default": 100.0
                            },
                            "final_temp": {
                                "type": "number",
                                "description": "Stopping temperature for simulated annealing",
                                "default": 0.1
                            },
                            "cooling_rate": {
                                "type": "number",
                                "description": "Temperature decrease rate per step",
                                "default": 0.995
                            },
                            "max_iterations": {
                                "type": "integer",
                                "description": "Maximum simulated-annealing iterations per replicate",
                                "default": 5000
                            },
                            "verbose": {
                                "type": "boolean",
                                "description": "Print per-replicate progress",
                                "default": False
                            }
                        },
                        "required": ["target_statistics_csv", "colors", "num_replicates"]
                    }
                ),

                Tool(
                    name="export_slice_csv",
                    description=(
                        "Export the current 2D slice data to a CSV file. "
                        "Returns the file path where data was saved."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filename": {
                                "type": "string",
                                "description": "Name for the CSV file (without path)",
                                "default": "slice_data.csv"
                            },
                            "include_3d": {
                                "type": "boolean",
                                "description": "Include original 3D coordinates in export",
                                "default": True
                            }
                        }
                    }
                ),

                Tool(
                    name="export_slice_to_physicell",
                    description=(
                        "Export the current slice (created via `create_slice`) as a PhysiCell IC "
                        "CSV. The `geometry` parameter selects between two modes: `\"2D\"` (default) "
                        "writes the slice-plane projection using each cell's 2D coordinates, a "
                        "fixed z (default 0), and disk-area volumes — suited to PhysiCell 2D "
                        "simulations; `\"3D\"` writes each cell's original 3D position from "
                        "`center_3d` and sphere-volume volumes — suited to PhysiCell 3D "
                        "simulations where the slice acted as a spatial filter. The `z` argument "
                        "is ignored when `geometry=\"3D\"`."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filename": {
                                "type": "string",
                                "description": "Name for the PhysiCell IC CSV file (without path)",
                                "default": "slice_physicell.csv"
                            },
                            "geometry": {
                                "type": "string",
                                "description": "'2D' writes slice-plane projection at the supplied z with disk-area volumes; '3D' writes each cell's original 3D position with sphere-volume volumes (z is ignored).",
                                "enum": ["2D", "3D"],
                                "default": "2D"
                            },
                            "z": {
                                "type": "number",
                                "description": "Fixed z written for every row in 2D mode. Ignored in 3D mode.",
                                "default": 0.0
                            },
                            "format": {
                                "type": "string",
                                "description": "PhysiCell CSV schema: 'modern' (header x,y,z,type [,volume]) or 'legacy' (headerless x,y,z,cell_type_id,volume)",
                                "enum": ["modern", "legacy"],
                                "default": "modern"
                            },
                            "include_volume": {
                                "type": "boolean",
                                "description": "Append the volume column. Ignored for 'legacy' format (volume is always required there).",
                                "default": True
                            },
                            "cell_type_mapping": {
                                "type": "object",
                                "description": "Mapping from cell-type name to PhysiCell integer ID. Required when format='legacy'.",
                                "additionalProperties": {"type": "integer"}
                            }
                        }
                    }
                ),

                Tool(
                    name="visualize_tissue",
                    description=(
                        "Create a 3D visualization of the tissue and save as PNG. "
                        "Returns the file path to the saved image."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "elevation": {
                                "type": "number",
                                "description": "Viewing elevation angle in degrees",
                                "default": 20
                            },
                            "azimuth": {
                                "type": "number",
                                "description": "Viewing azimuth angle in degrees",
                                "default": 45
                            },
                            "filename": {
                                "type": "string",
                                "description": "Name for the PNG file (without path)",
                                "default": "tissue_3d.png"
                            }
                        }
                    }
                ),
                
                Tool(
                    name="visualize_slice_2d",
                    description=(
                        "Create a 2D visualization of the current slice and save as PNG. "
                        "Returns the file path to the saved image."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "filename": {
                                "type": "string",
                                "description": "Name for the PNG file (without path)",
                                "default": "slice_2d.png"
                            }
                        }
                    }
                ),
                
                Tool(
                    name="reset_tissue",
                    description=(
                        "Clear the current tissue and start fresh. "
                        "Use this to begin a new simulation."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {}
                    }
                ),
                
                Tool(
                    name="load_target_statistics",
                    description=(
                        "Load target spatial statistics from a CSV file or extract from current tissue. "
                        "These statistics will be used to generate replicate tissues."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "csv_filepath": {
                                "type": "string",
                                "description": "Path to CSV file with target statistics (optional)"
                            },
                            "use_current_tissue": {
                                "type": "boolean",
                                "description": "Extract statistics from current tissue",
                                "default": False
                            },
                            "network_mode": {
                                "type": "string",
                                "description": "Network analysis mode: 'contact' or 'radius'",
                                "default": "contact"
                            },
                            "network_radius": {
                                "type": "number",
                                "description": "Distance threshold for 'radius' mode (micrometers)"
                            }
                        }
                    }
                ),
                
                Tool(
                    name="setup_replicate_generator",
                    description=(
                        "Setup the replicate generator with tissue dimensions and cell parameters. "
                        "Must be called after load_target_statistics."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "height": {
                                "type": "number",
                                "description": "Tissue height in micrometers",
                                "minimum": 50
                            },
                            "width": {
                                "type": "number",
                                "description": "Tissue width in micrometers",
                                "minimum": 50
                            },
                            "thickness": {
                                "type": "number",
                                "description": "Tissue thickness in micrometers",
                                "minimum": 20
                            },
                            "cell_radii": {
                                "type": "object",
                                "description": "Dict mapping cell types to [min_radius, max_radius]",
                                "additionalProperties": {
                                    "type": "array",
                                    "items": {"type": "number"},
                                    "minItems": 2,
                                    "maxItems": 2
                                }
                            },
                            "seed": {
                                "type": "integer",
                                "description": "Random seed for reproducibility"
                            },
                            "method": {
                                "type": "string",
                                "enum": ["radius_tuning", "graph_coloring"],
                                "description": (
                                    "Replicate strategy. 'radius_tuning' (default) repacks and "
                                    "nudges per-type radii to match proportions. 'graph_coloring' "
                                    "packs geometry once and assigns cell types via simulated-"
                                    "annealing graph coloring to match interaction statistics — "
                                    "more consistent and faster-converging for interaction targets."
                                ),
                                "default": "radius_tuning"
                            },
                            "n_restarts": {
                                "type": "integer",
                                "description": (
                                    "For 'graph_coloring': independent SA runs per replicate; "
                                    "the lowest-cost coloring is kept (hardens against bad local "
                                    "minima)."
                                ),
                                "default": 1,
                                "minimum": 1
                            },
                            "radius_optimizer": {
                                "type": "string",
                                "enum": ["heuristic", "differential_evolution"],
                                "description": (
                                    "For 'radius_tuning': proportion tuner. 'heuristic' (default) "
                                    "or gradient-free scipy 'differential_evolution' (slower, more "
                                    "robust)."
                                ),
                                "default": "heuristic"
                            }
                        },
                        "required": ["height", "width", "thickness", "cell_radii"]
                    }
                ),

                Tool(
                    name="generate_replicates",
                    description=(
                        "Generate multiple tissue replicates matching target statistics. "
                        "Requires setup_replicate_generator to be called first."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "num_replicates": {
                                "type": "integer",
                                "description": "Number of replicates to generate",
                                "minimum": 1,
                                "maximum": 100,
                                "default": 1
                            },
                            "max_attempts": {
                                "type": "integer",
                                "description": "Max cell packing attempts per replicate",
                                "default": 1000
                            },
                            "min_spacing": {
                                "type": "number",
                                "description": "Minimum spacing between cells",
                                "default": 0.5
                            },
                            "allow_boundary": {
                                "type": "boolean",
                                "description": "Allow cells extending beyond bounds",
                                "default": True
                            },
                            "max_iterations": {
                                "type": "integer",
                                "description": "Max parameter adjustment iterations",
                                "default": 5,
                                "minimum": 1,
                                "maximum": 20
                            },
                            "tolerance": {
                                "type": "number",
                                "description": "Acceptable divergence threshold",
                                "default": 0.15,
                                "minimum": 0.01,
                                "maximum": 1.0
                            }
                        },
                        "required": ["num_replicates"]
                    }
                ),
                
                Tool(
                    name="export_replicate_statistics",
                    description=(
                        "Export statistics for all generated replicates to CSV files."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "base_filename": {
                                "type": "string",
                                "description": "Base name for output files",
                                "default": "replicates"
                            }
                        }
                    }
                ),
                
                Tool(
                    name="export_replicate_tissues",
                    description=(
                        "Export each replicate tissue to a separate CSV file."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "output_dir": {
                                "type": "string",
                                "description": "Directory for output files",
                                "default": "replicates"
                            }
                        }
                    }
                ),
                
                Tool(
                    name="get_replicate_summary",
                    description=(
                        "Get summary statistics for all generated replicates."
                    ),
                    inputSchema={
                        "type": "object",
                        "properties": {}
                    }
                )
            ]
    
    def _register_handlers(self):
        """Register handlers for each tool."""
        
        @self.server.call_tool()
        async def call_tool(name: str, arguments: Dict[str, Any]) -> List[TextContent]:
            """Handle tool calls from the LLM."""
            
            try:
                if name == "create_tissue":
                    return await self._handle_create_tissue(arguments)
                elif name == "generate_cells":
                    return await self._handle_generate_cells(arguments)
                elif name == "get_tissue_statistics":
                    return await self._handle_get_tissue_statistics(arguments)
                elif name == "create_slice":
                    return await self._handle_create_slice(arguments)
                elif name == "get_slice_statistics":
                    return await self._handle_get_slice_statistics(arguments)
                elif name == "create_serial_slices":
                    return await self._handle_create_serial_slices(arguments)
                elif name == "export_tissue_csv":
                    return await self._handle_export_tissue_csv(arguments)
                elif name == "export_slice_csv":
                    return await self._handle_export_slice_csv(arguments)
                elif name == "export_tissue_to_physicell":
                    return await self._handle_export_tissue_to_physicell(arguments)
                elif name == "export_slice_to_physicell":
                    return await self._handle_export_slice_to_physicell(arguments)
                elif name == "visualize_tissue":
                    return await self._handle_visualize_tissue(arguments)
                elif name == "visualize_slice_2d":
                    return await self._handle_visualize_slice_2d(arguments)
                elif name == "reset_tissue":
                    return await self._handle_reset_tissue(arguments)
                elif name == "load_target_statistics":
                    return await self._handle_load_target_statistics(arguments)
                elif name == "load_tissue_from_csv":
                    return await self._handle_load_tissue_from_csv(arguments)
                elif name == "load_target_statistics_from_coordinates":
                    return await self._handle_load_target_statistics_from_coordinates(arguments)
                elif name == "assign_cell_types":
                    return await self._handle_assign_cell_types(arguments)
                elif name == "generate_colored_replicates":
                    return await self._handle_generate_colored_replicates(arguments)
                elif name == "setup_replicate_generator":
                    return await self._handle_setup_replicate_generator(arguments)
                elif name == "generate_replicates":
                    return await self._handle_generate_replicates(arguments)
                elif name == "export_replicate_statistics":
                    return await self._handle_export_replicate_statistics(arguments)
                elif name == "export_replicate_tissues":
                    return await self._handle_export_replicate_tissues(arguments)
                elif name == "get_replicate_summary":
                    return await self._handle_get_replicate_summary(arguments)
                else:
                    return [TextContent(
                        type="text",
                        text=f"Unknown tool: {name}"
                    )]
            except Exception as e:
                return [TextContent(
                    type="text",
                    text=f"Error executing {name}: {str(e)}"
                )]
    
    async def _handle_create_tissue(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle create_tissue tool call."""
        height = args["height"]
        width = args["width"]
        thickness = args["thickness"]
        cell_types = args["cell_types"]
        
        # Convert cell_types dict to proper format
        cell_radii = {
            name: tuple(radii) 
            for name, radii in cell_types.items()
        }
        
        self.current_tissue = TissueSection(
            height=height,
            width=width,
            thickness=thickness,
            cell_radii=cell_radii
        )
        
        result = {
            "status": "success",
            "message": f"Created tissue: {width}x{height}x{thickness} μm",
            "cell_types": list(cell_radii.keys()),
            "cell_type_radii": cell_types
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_generate_cells(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle generate_cells tool call."""
        if self.current_tissue is None:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue created. Call create_tissue first."})
            )]
        
        max_attempts = args.get("max_attempts", 1000)
        min_spacing = args.get("min_spacing", 0.5)
        allow_boundary = args.get("allow_boundary_cells", True)
        
        num_cells = self.current_tissue.generate_cells(
            max_attempts=max_attempts,
            min_spacing=min_spacing,
            allow_boundary_cells=allow_boundary
        )
        
        stats = self.current_tissue.get_cell_statistics()
        
        result = {
            "status": "success",
            "num_cells_generated": num_cells,
            "interior_cells": stats.get("interior_cells", 0),
            "boundary_cells": stats.get("boundary_cells", 0),
            "packing_fraction": round(stats.get("packing_fraction", 0), 4),
            "cell_type_counts": stats.get("cell_types", {})
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_get_tissue_statistics(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle get_tissue_statistics tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells available."})
            )]
        
        stats = self.current_tissue.get_cell_statistics()
        
        # Make stats JSON serializable
        result = {
            "total_cells": stats["total_cells"],
            "interior_cells": stats["interior_cells"],
            "boundary_cells": stats["boundary_cells"],
            "packing_fraction": round(stats["packing_fraction"], 4),
            "cell_types": stats["cell_types"],
            "average_radii": {k: round(v, 2) for k, v in stats["avg_radii"].items()},
            "tissue_dimensions": {
                "width": self.current_tissue.width,
                "height": self.current_tissue.height,
                "thickness": self.current_tissue.thickness
            }
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_create_slice(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle create_slice tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells available."})
            )]
        
        self.current_slicer = TissueSlicer(self.current_tissue)
        
        # Extract slicing parameters
        kwargs = {}
        if "z_position" in args:
            kwargs["z_position"] = args["z_position"]
        if "angle_x" in args:
            kwargs["angle_x"] = args["angle_x"]
        if "angle_y" in args:
            kwargs["angle_y"] = args["angle_y"]
        if "point" in args:
            kwargs["point"] = tuple(args["point"])
        if "normal" in args:
            kwargs["normal"] = tuple(args["normal"])
        
        slice_cells = self.current_slicer.slice_plane(**kwargs)
        stats = self.current_slicer.get_slice_statistics()
        
        result = {
            "status": "success",
            "num_cells_in_slice": len(slice_cells),
            "plane_point": [round(x, 2) for x in stats["plane_point"]],
            "plane_normal": [round(x, 3) for x in stats["plane_normal"]],
            "cell_type_counts": stats["cell_types"],
            "mean_distance_from_plane": round(stats["mean_distance_from_plane"], 2)
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_get_slice_statistics(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle get_slice_statistics tool call."""
        if self.current_slicer is None:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No slice created. Call create_slice first."})
            )]
        
        stats = self.current_slicer.get_slice_statistics()
        
        result = {
            "num_cells": stats["num_cells"],
            "plane_point": [round(x, 2) for x in stats["plane_point"]],
            "plane_normal": [round(x, 3) for x in stats["plane_normal"]],
            "cell_types": stats["cell_types"],
            "avg_intersection_radii": {k: round(v, 2) for k, v in stats["avg_intersection_radii"].items()},
            "mean_distance_from_plane": round(stats["mean_distance_from_plane"], 2),
            "max_distance_from_plane": round(stats["max_distance_from_plane"], 2)
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_create_serial_slices(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle create_serial_slices tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells available."})
            )]
        
        num_slices = args["num_slices"]
        slicers = create_standard_slices(self.current_tissue, num_slices=num_slices)
        
        # Collect statistics for each slice
        slices_info = []
        for i, slicer in enumerate(slicers, 1):
            stats = slicer.get_slice_statistics()
            slices_info.append({
                "slice_number": i,
                "z_position": round(stats["plane_point"][2], 2),
                "num_cells": stats["num_cells"],
                "cell_types": stats["cell_types"]
            })
        
        # Store the last slicer as current
        self.current_slicer = slicers[-1] if slicers else None
        
        result = {
            "status": "success",
            "num_slices_created": len(slicers),
            "slices": slices_info
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_export_tissue_csv(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle export_tissue_csv tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells to export."})
            )]
        
        filename = args.get("filename", "tissue_data.csv")
        filepath = os.path.join(self.temp_dir, filename)
        
        self.current_tissue.export_to_csv(filepath)
        
        result = {
            "status": "success",
            "filepath": filepath,
            "num_cells_exported": len(self.current_tissue.cells)
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_load_tissue_from_csv(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle load_tissue_from_csv tool call."""
        filepath = args.get("filepath")
        if not filepath:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "filepath is required"})
            )]

        height = args.get("height")
        width = args.get("width")
        thickness = args.get("thickness")
        default_radius = args.get("default_radius", 10.0)

        try:
            self.current_tissue = load_tissue_from_csv(
                filepath,
                height=height,
                width=width,
                thickness=thickness,
                default_radius=default_radius
            )

            stats = self.current_tissue.get_cell_statistics()

            result = {
                "status": "success",
                "filepath": filepath,
                "num_cells": stats['total_cells'],
                "cell_types": stats.get('cell_types', {}),
                "dimensions": {
                    "height": self.current_tissue.height,
                    "width": self.current_tissue.width,
                    "thickness": self.current_tissue.thickness
                },
                "packing_fraction": round(stats['packing_fraction'], 4) if stats.get('packing_fraction') is not None else None
            }

            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]

        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to load tissue: {str(e)}"})
            )]

    async def _handle_load_target_statistics_from_coordinates(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle load_target_statistics_from_coordinates tool call."""
        filepath = args.get("filepath")
        if not filepath:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "filepath is required"})
            )]

        network_mode = args.get("network_mode", "contact")
        network_radius = args.get("network_radius")

        try:
            self.target_stats = load_target_statistics_from_coordinates(
                filepath,
                network_mode=network_mode,
                network_radius=network_radius
            )

            # Store network mode for replicate generator
            self.network_mode = network_mode
            self.network_radius = network_radius

            result = {
                "status": "success",
                "source": "coordinate_csv",
                "filepath": filepath,
                "network_mode": network_mode,
                "num_interaction_types": len(self.target_stats.interaction_stats),
                "cell_type_proportions": self.target_stats.cell_type_proportions,
                "target_cell_count": self.target_stats.target_cell_count,
                "target_density": round(self.target_stats.target_density, 4) if self.target_stats.target_density else None
            }

            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]

        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to load statistics: {str(e)}"})
            )]

    async def _handle_assign_cell_types(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle assign_cell_types tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells available. Call create_tissue/generate_cells or load_tissue_from_csv first."})
            )]

        target_csv = args.get("target_statistics_csv")
        colors = args.get("colors")
        if not target_csv:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "target_statistics_csv is required"})
            )]
        if not colors:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "colors is required (non-empty list of cell-type names)"})
            )]

        z_pos = args.get("z_position")
        if z_pos is None:
            z_pos = self.current_tissue.thickness / 2

        net_radius = args.get("network_radius")
        if net_radius is None:
            net_radius = getattr(self, "network_radius", None) or 50.0

        seed = args.get("seed")
        initial_temp = args.get("initial_temp", 100.0)
        final_temp = args.get("final_temp", 0.1)
        cooling_rate = args.get("cooling_rate", 0.995)
        max_iterations = args.get("max_iterations", 5000)
        verbose = args.get("verbose", False)
        warm_start = args.get("warm_start", False)

        try:
            # 1. Slice the current tissue at z_pos.
            slicer = TissueSlicer(self.current_tissue)
            slicer.slice_plane(z_position=z_pos)

            # 2. Build a radius-mode spatial network from the slice.
            analyzer = SpatialNetworkAnalyzer()
            graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=net_radius)

            # 3. Load target statistics from the graph-coloring-format CSV.
            target_stats = load_graph_coloring_stats_from_csv(target_csv, colors)

            # 4. Run simulated annealing.
            colorizer = GraphColorizer(
                target_graph=graph,
                colors=colors,
                target_statistics=target_stats,
                seed=seed
            )

            # Determine whether a warm-start initial coloring applies: the
            # prior assignment must exist and cover exactly the new graph's nodes.
            initial_coloring = None
            warm_start_applied = False
            if warm_start:
                prior = getattr(self, "cell_type_assignment", None)
                if isinstance(prior, dict) and set(prior.keys()) == set(graph.nodes()):
                    initial_coloring = prior
                    warm_start_applied = True

            assignment = colorizer.colorize(
                initial_temp=initial_temp,
                final_temp=final_temp,
                cooling_rate=cooling_rate,
                max_iterations=max_iterations,
                verbose=verbose,
                initial_coloring=initial_coloring
            )

            # Store on server for downstream tools.
            self.cell_type_assignment = assignment
            self.current_slicer = slicer

            # 5. Build result summary.
            from collections import Counter
            color_counts = dict(Counter(assignment.values()))

            result = {
                "status": "success",
                "num_nodes_colored": len(assignment),
                "color_counts": color_counts,
                "seed": seed,
                "used_z_position": z_pos,
                "used_network_radius": net_radius,
                "warm_start_applied": warm_start_applied
            }

            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]

        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to assign cell types: {str(e)}"})
            )]

    async def _handle_generate_colored_replicates(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle generate_colored_replicates tool call."""
        import numpy as np
        from collections import Counter

        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells available. Call create_tissue/generate_cells or load_tissue_from_csv first."})
            )]

        target_csv = args.get("target_statistics_csv")
        colors = args.get("colors")
        num_replicates = args.get("num_replicates")
        if not target_csv:
            return [TextContent(type="text", text=json.dumps({"error": "target_statistics_csv is required"}))]
        if not colors:
            return [TextContent(type="text", text=json.dumps({"error": "colors is required (non-empty list of cell-type names)"}))]
        if not isinstance(num_replicates, int) or num_replicates < 1:
            return [TextContent(type="text", text=json.dumps({"error": "num_replicates is required and must be an integer >= 1"}))]

        z_pos = args.get("z_position")
        if z_pos is None:
            z_pos = self.current_tissue.thickness / 2

        net_radius = args.get("network_radius")
        if net_radius is None:
            net_radius = getattr(self, "network_radius", None) or 50.0

        seed = args.get("seed")
        warm_start = args.get("warm_start", False)
        initial_temp = args.get("initial_temp", 100.0)
        final_temp = args.get("final_temp", 0.1)
        cooling_rate = args.get("cooling_rate", 0.995)
        max_iterations = args.get("max_iterations", 5000)
        verbose = args.get("verbose", False)

        try:
            # 1. Slice the current tissue and build a radius-mode network.
            slicer = TissueSlicer(self.current_tissue)
            slicer.slice_plane(z_position=z_pos)
            analyzer = SpatialNetworkAnalyzer()
            graph = analyzer.build_network_from_slice(slicer, mode="radius", radius=net_radius)

            # 2. Load target statistics.
            target_stats = load_graph_coloring_stats_from_csv(target_csv, colors)

            # 3. Generate num_replicates colorings.
            colorings = []
            previous = None
            for k in range(num_replicates):
                replicate_seed = (
                    int(np.random.SeedSequence([seed, k]).generate_state(1)[0])
                    if seed is not None else None
                )
                colorizer = GraphColorizer(
                    target_graph=graph,
                    colors=colors,
                    target_statistics=target_stats,
                    seed=replicate_seed,
                )
                initial_coloring = previous if (warm_start and previous is not None) else None
                coloring = colorizer.colorize(
                    initial_temp=initial_temp,
                    final_temp=final_temp,
                    cooling_rate=cooling_rate,
                    max_iterations=max_iterations,
                    verbose=verbose,
                    initial_coloring=initial_coloring,
                )
                colorings.append(coloring)
                previous = coloring

            # Store the first replicate as the active assignment.
            self.cell_type_assignment = colorings[0]
            self.current_slicer = slicer

            # 4. Mean pairwise diversity: fraction of nodes whose color differs,
            #    averaged over all replicate pairs (0.0 when only one replicate).
            nodes = list(graph.nodes())
            pair_diffs = []
            for a in range(len(colorings)):
                for b in range(a + 1, len(colorings)):
                    ca, cb = colorings[a], colorings[b]
                    differing = sum(1 for n in nodes if ca.get(n) != cb.get(n))
                    pair_diffs.append(differing / len(nodes) if nodes else 0.0)
            mean_pairwise_diversity = float(np.mean(pair_diffs)) if pair_diffs else 0.0

            result = {
                "status": "success",
                "num_replicates": num_replicates,
                "num_nodes_colored": len(colorings[0]),
                "replicate_color_counts": [dict(Counter(c.values())) for c in colorings],
                "mean_pairwise_diversity": mean_pairwise_diversity,
                "warm_start": bool(warm_start),
                "seed": seed,
                "used_z_position": z_pos,
                "used_network_radius": net_radius,
            }
            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to generate colored replicates: {str(e)}"})
            )]

    async def _handle_export_slice_csv(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle export_slice_csv tool call."""
        if self.current_slicer is None:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No slice to export. Call create_slice first."})
            )]

        filename = args.get("filename", "slice_data.csv")
        include_3d = args.get("include_3d", True)
        filepath = os.path.join(self.temp_dir, filename)

        self.current_slicer.export_slice_csv(filepath, include_3d=include_3d)

        result = {
            "status": "success",
            "filepath": filepath,
            "num_cells_exported": len(self.current_slicer.slice_cells),
            "include_3d_coordinates": include_3d
        }

        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]

    async def _handle_export_tissue_to_physicell(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle export_tissue_to_physicell tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells to export. Run create_tissue + generate_cells first."})
            )]

        filename = args.get("filename", "tissue_physicell.csv")
        fmt = args.get("format", "modern")
        include_volume = args.get("include_volume", True)
        cell_type_mapping = args.get("cell_type_mapping")

        filepath = os.path.join(self.temp_dir, filename)

        try:
            exporter = PhysiCellExporter()
            exporter.export_tissue(
                self.current_tissue,
                filepath,
                cell_type_mapping=cell_type_mapping,
                include_volume=include_volume,
                format=fmt,
            )
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to export tissue to PhysiCell: {e}"})
            )]

        result = {
            "status": "success",
            "filepath": filepath,
            "num_cells_exported": len(self.current_tissue.cells),
            "format": fmt,
        }
        return [TextContent(type="text", text=json.dumps(result, indent=2))]


    async def _handle_export_slice_to_physicell(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle export_slice_to_physicell tool call."""
        if self.current_slicer is None or not self.current_slicer.slice_cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No slice available to export. Run create_slice first."})
            )]

        filename = args.get("filename", "slice_physicell.csv")
        geometry = args.get("geometry", "2D")
        z = float(args.get("z", 0.0))
        fmt = args.get("format", "modern")
        include_volume = args.get("include_volume", True)
        cell_type_mapping = args.get("cell_type_mapping")

        filepath = os.path.join(self.temp_dir, filename)

        try:
            exporter = PhysiCellExporter()
            exporter.export_slice(
                self.current_slicer.slice_cells,
                filepath,
                geometry=geometry,
                z=z,
                cell_type_mapping=cell_type_mapping,
                include_volume=include_volume,
                format=fmt,
            )
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to export slice to PhysiCell: {e}"})
            )]

        result = {
            "status": "success",
            "filepath": filepath,
            "num_cells_exported": len(self.current_slicer.slice_cells),
            "format": fmt,
            "geometry": geometry,
            "used_z": z if geometry == "2D" else None,
        }
        return [TextContent(type="text", text=json.dumps(result, indent=2))]

    async def _handle_visualize_tissue(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle visualize_tissue tool call."""
        if self.current_tissue is None or not self.current_tissue.cells:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No tissue with cells to visualize."})
            )]
        
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        
        elevation = args.get("elevation", 20)
        azimuth = args.get("azimuth", 45)
        filename = args.get("filename", "tissue_3d.png")
        filepath = os.path.join(self.temp_dir, filename)
        
        # Create visualization
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Simplified visualization for MCP (reduced resolution)
        import numpy as np
        cell_types = list(set(c.cell_type for c in self.current_tissue.cells))
        colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        for cell in self.current_tissue.cells[:100]:  # Limit for performance
            u = np.linspace(0, 2 * np.pi, 10)
            v = np.linspace(0, np.pi, 10)
            x = cell.radius * np.outer(np.cos(u), np.sin(v)) + cell.center[0]
            y = cell.radius * np.outer(np.sin(u), np.sin(v)) + cell.center[1]
            z = cell.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cell.center[2]
            
            color = color_map[cell.cell_type]
            alpha = 0.3 if cell.is_boundary else 0.6
            
            ax.plot_surface(x, y, z, facecolors=np.tile(color, x.shape + (1,)),
                          alpha=alpha, linewidth=0, antialiased=True, shade=False)
        
        ax.set_xlabel('Width (μm)')
        ax.set_ylabel('Height (μm)')
        ax.set_zlabel('Thickness (μm)')
        ax.set_xlim(0, self.current_tissue.width)
        ax.set_ylim(0, self.current_tissue.height)
        ax.set_zlim(0, self.current_tissue.thickness)
        ax.view_init(elev=elevation, azim=azimuth)
        ax.set_title(f'3D Tissue: {len(self.current_tissue.cells)} cells')
        
        plt.tight_layout()
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        result = {
            "status": "success",
            "filepath": filepath,
            "cells_visualized": min(100, len(self.current_tissue.cells))
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_visualize_slice_2d(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle visualize_slice_2d tool call."""
        if self.current_slicer is None:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No slice to visualize. Call create_slice first."})
            )]
        
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        import numpy as np
        
        filename = args.get("filename", "slice_2d.png")
        filepath = os.path.join(self.temp_dir, filename)
        
        # Create 2D visualization
        fig, ax = plt.subplots(figsize=(10, 10))
        
        cell_types = list(set(c.cell_type for c in self.current_slicer.slice_cells))
        colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))
        color_map = dict(zip(cell_types, colors))
        
        stats = self.current_slicer.get_slice_statistics()
        
        for slice_cell in self.current_slicer.slice_cells:
            color = color_map[slice_cell.cell_type]
            max_dist = stats['max_distance_from_plane'] if stats['max_distance_from_plane'] > 0 else 1
            alpha = 1.0 - (slice_cell.distance_from_plane / max_dist) * 0.5
            
            circle = plt.Circle(
                slice_cell.center_2d,
                slice_cell.intersection_radius,
                color=color,
                alpha=alpha,
                edgecolor='black',
                linewidth=0.5
            )
            ax.add_patch(circle)
        
        if self.current_slicer.slice_cells:
            x_coords = [c.center_2d[0] for c in self.current_slicer.slice_cells]
            y_coords = [c.center_2d[1] for c in self.current_slicer.slice_cells]
            margin = 20
            ax.set_xlim(min(x_coords) - margin, max(x_coords) + margin)
            ax.set_ylim(min(y_coords) - margin, max(y_coords) + margin)
        
        ax.set_aspect('equal')
        ax.set_xlabel('U coordinate (μm)')
        ax.set_ylabel('V coordinate (μm)')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'2D Tissue Slice: {len(self.current_slicer.slice_cells)} cells')
        
        plt.tight_layout()
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        result = {
            "status": "success",
            "filepath": filepath,
            "cells_visualized": len(self.current_slicer.slice_cells)
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_reset_tissue(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle reset_tissue tool call."""
        self.current_tissue = None
        self.current_slicer = None
        
        result = {
            "status": "success",
            "message": "Tissue and slice data cleared. Ready for new simulation."
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def _handle_load_target_statistics(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle load_target_statistics tool call."""
        csv_filepath = args.get("csv_filepath")
        use_current = args.get("use_current_tissue", False)
        network_mode = args.get("network_mode", "contact")
        network_radius = args.get("network_radius")
        
        try:
            if csv_filepath:
                # Load from CSV
                self.target_stats = load_target_statistics_from_csv(csv_filepath)
                result = {
                    "status": "success",
                    "source": "csv_file",
                    "filepath": csv_filepath,
                    "num_interaction_types": len(self.target_stats.interaction_stats),
                    "interaction_pairs": [
                        f"{s.type_a}-{s.type_b}" 
                        for s in self.target_stats.interaction_stats
                    ]
                }
            elif use_current:
                # Extract from current tissue
                if self.current_tissue is None or not self.current_tissue.cells:
                    return [TextContent(
                        type="text",
                        text=json.dumps({"error": "No tissue with cells available."})
                    )]
                
                self.target_stats = load_target_statistics_from_tissue(
                    self.current_tissue,
                    network_mode=network_mode,
                    network_radius=network_radius
                )
                
                result = {
                    "status": "success",
                    "source": "current_tissue",
                    "network_mode": network_mode,
                    "num_cells": self.current_tissue.get_cell_statistics()['total_cells'],
                    "num_interaction_types": len(self.target_stats.interaction_stats),
                    "cell_type_proportions": self.target_stats.cell_type_proportions,
                    "target_cell_count": self.target_stats.target_cell_count,
                    "target_density": round(self.target_stats.target_density, 4) if self.target_stats.target_density else None
                }
            else:
                return [TextContent(
                    type="text",
                    text=json.dumps({"error": "Must specify either csv_filepath or use_current_tissue=True"})
                )]
            
            # Store network mode for replicate generator
            self.network_mode = network_mode
            self.network_radius = network_radius
            
            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]
        
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to load statistics: {str(e)}"})
            )]
    
    async def _handle_setup_replicate_generator(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle setup_replicate_generator tool call."""
        if not hasattr(self, 'target_stats') or self.target_stats is None:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "Must load target statistics first using load_target_statistics."})
            )]
        
        height = args["height"]
        width = args["width"]
        thickness = args["thickness"]
        cell_radii_dict = args["cell_radii"]
        seed = args.get("seed")
        method = args.get("method", "radius_tuning")
        n_restarts = args.get("n_restarts", 1)
        radius_optimizer = args.get("radius_optimizer", "heuristic")

        # Convert cell_radii dict to proper format
        cell_radii = {
            name: tuple(radii)
            for name, radii in cell_radii_dict.items()
        }

        try:
            self.replicate_generator = ReplicateGenerator(
                target_stats=self.target_stats,
                tissue_dimensions=(height, width, thickness),
                base_cell_radii=cell_radii,
                network_mode=self.network_mode,
                network_radius=self.network_radius,
                seed=seed,
                method=method,
                n_restarts=n_restarts,
                radius_optimizer=radius_optimizer,
            )

            result = {
                "status": "success",
                "tissue_dimensions": {
                    "width": width,
                    "height": height,
                    "thickness": thickness
                },
                "cell_types": list(cell_radii.keys()),
                "network_mode": self.network_mode,
                "seed": seed,
                "method": method,
                "ready_to_generate": True
            }
            
            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]
        
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to setup generator: {str(e)}"})
            )]
    
    async def _handle_generate_replicates(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle generate_replicates tool call."""
        if self.replicate_generator is None:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "Must setup replicate generator first using setup_replicate_generator."})
            )]
        
        num_replicates = args["num_replicates"]
        max_attempts = args.get("max_attempts", 1000)
        min_spacing = args.get("min_spacing", 0.5)
        allow_boundary = args.get("allow_boundary", True)
        max_iterations = args.get("max_iterations", 5)
        tolerance = args.get("tolerance", 0.15)
        
        try:
            # Generate replicates
            self.generated_replicates = self.replicate_generator.generate_replicates(
                num_replicates=num_replicates,
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary=allow_boundary,
                max_iterations=max_iterations,
                tolerance=tolerance
            )
            
            # Collect summary statistics
            replicate_summaries = []
            for tissue, stats in self.generated_replicates:
                replicate_summaries.append({
                    "replicate_id": stats.replicate_id,
                    "num_cells": stats.num_cells,
                    "cell_types": stats.cell_type_counts,
                    "packing_fraction": round(stats.packing_fraction, 4),
                    "divergence_score": round(stats.divergence_score, 4)
                })
            
            result = {
                "status": "success",
                "num_replicates_generated": len(self.generated_replicates),
                "replicates": replicate_summaries,
                "avg_divergence": round(
                    sum(s.divergence_score for _, s in self.generated_replicates) / len(self.generated_replicates),
                    4
                )
            }
            
            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]
        
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to generate replicates: {str(e)}"})
            )]
    
    async def _handle_export_replicate_statistics(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle export_replicate_statistics tool call."""
        if not self.generated_replicates:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No replicates generated. Call generate_replicates first."})
            )]
        
        base_filename = args.get("base_filename", "replicates")
        filepath = os.path.join(self.temp_dir, base_filename)
        
        try:
            self.replicate_generator.export_replicate_statistics(
                self.generated_replicates,
                filepath
            )
            
            result = {
                "status": "success",
                "summary_file": f"{filepath}_summary.csv",
                "interactions_file": f"{filepath}_interactions.csv",
                "num_replicates": len(self.generated_replicates)
            }
            
            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]
        
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to export statistics: {str(e)}"})
            )]
    
    async def _handle_export_replicate_tissues(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle export_replicate_tissues tool call."""
        if not self.generated_replicates:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No replicates generated. Call generate_replicates first."})
            )]
        
        output_dir = args.get("output_dir", "replicates")
        full_path = os.path.join(self.temp_dir, output_dir)
        
        try:
            self.replicate_generator.export_replicate_tissues(
                self.generated_replicates,
                full_path
            )
            
            result = {
                "status": "success",
                "output_directory": full_path,
                "num_files_created": len(self.generated_replicates)
            }
            
            return [TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]
        
        except Exception as e:
            return [TextContent(
                type="text",
                text=json.dumps({"error": f"Failed to export tissues: {str(e)}"})
            )]
    
    async def _handle_get_replicate_summary(self, args: Dict[str, Any]) -> List[TextContent]:
        """Handle get_replicate_summary tool call."""
        if not self.generated_replicates:
            return [TextContent(
                type="text",
                text=json.dumps({"error": "No replicates generated. Call generate_replicates first."})
            )]
        
        # Compute summary statistics across all replicates
        import numpy as np
        
        num_cells_list = [stats.num_cells for _, stats in self.generated_replicates]
        divergence_list = [stats.divergence_score for _, stats in self.generated_replicates]
        packing_list = [stats.packing_fraction for _, stats in self.generated_replicates]
        
        # Cell type statistics
        all_cell_types = set()
        for _, stats in self.generated_replicates:
            all_cell_types.update(stats.cell_type_counts.keys())
        
        cell_type_stats = {}
        for ct in all_cell_types:
            counts = [stats.cell_type_counts.get(ct, 0) for _, stats in self.generated_replicates]
            total_cells = [stats.num_cells for _, stats in self.generated_replicates]
            proportions = [c/t for c, t in zip(counts, total_cells)]
            
            cell_type_stats[ct] = {
                "mean_count": round(np.mean(counts), 1),
                "std_count": round(np.std(counts), 1),
                "mean_proportion": round(np.mean(proportions), 4),
                "std_proportion": round(np.std(proportions), 4)
            }
        
        result = {
            "num_replicates": len(self.generated_replicates),
            "cell_counts": {
                "mean": round(np.mean(num_cells_list), 1),
                "std": round(np.std(num_cells_list), 1),
                "min": int(np.min(num_cells_list)),
                "max": int(np.max(num_cells_list))
            },
            "divergence_scores": {
                "mean": round(np.mean(divergence_list), 4),
                "std": round(np.std(divergence_list), 4),
                "min": round(np.min(divergence_list), 4),
                "max": round(np.max(divergence_list), 4)
            },
            "packing_fractions": {
                "mean": round(np.mean(packing_list), 4),
                "std": round(np.std(packing_list), 4),
                "min": round(np.min(packing_list), 4),
                "max": round(np.max(packing_list), 4)
            },
            "cell_type_statistics": cell_type_stats
        }
        
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )]
    
    async def run(self):
        """Run the MCP server."""
        async with stdio_server() as (read_stream, write_stream):
            await self.server.run(
                read_stream,
                write_stream,
                self.server.create_initialization_options()
            )


async def main():
    """Main entry point for the MCP server."""
    server = TissueSimulatorMCPServer()
    await server.run()


if __name__ == "__main__":
    asyncio.run(main())
