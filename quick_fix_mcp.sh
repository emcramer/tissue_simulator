#!/bin/bash
# Quick fix for MCP configuration

# Find python3
PYTHON=$(which python3)

if [ -z "$PYTHON" ]; then
    echo "Error: python3 not found"
    echo "Please install Python 3 first"
    exit 1
fi

echo "Using Python: $PYTHON"

# Create config
CONFIG_FILE="$HOME/Library/Application Support/Claude/claude_desktop_config.json"
mkdir -p "$HOME/Library/Application Support/Claude"

cat > "$CONFIG_FILE" << EOF
{
  "mcpServers": {
    "tissue-simulator": {
      "command": "$PYTHON",
      "args": [
        "/Users/cramere/tissue_simulator/run_mcp_server.py"
      ]
    }
  }
}
EOF

echo "✓ Configuration created!"
echo "✓ Restart Claude Desktop now"
