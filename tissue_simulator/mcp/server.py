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

from ..tissue import TissueSection
from ..slicing import TissueSlicer, create_standard_slices


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
                elif name == "visualize_tissue":
                    return await self._handle_visualize_tissue(arguments)
                elif name == "visualize_slice_2d":
                    return await self._handle_visualize_slice_2d(arguments)
                elif name == "reset_tissue":
                    return await self._handle_reset_tissue(arguments)
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
