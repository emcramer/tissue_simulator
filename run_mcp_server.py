#!/usr/bin/env python3
"""
Tissue Simulator MCP Server Entry Point

Run this script to start the MCP server:
    python run_mcp_server.py

Or make it executable and run directly:
    chmod +x run_mcp_server.py
    ./run_mcp_server.py
"""

import sys
import asyncio
from tissue_simulator.mcp import main

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nMCP server stopped.")
        sys.exit(0)
    except Exception as e:
        print(f"Error running MCP server: {e}", file=sys.stderr)
        sys.exit(1)
