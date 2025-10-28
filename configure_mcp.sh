#!/bin/bash
# Automatic MCP configuration script for Tissue Simulator

echo "================================================"
echo "Tissue Simulator MCP Configuration"
echo "================================================"
echo ""

# Find Python
echo "1. Finding Python executable..."
PYTHON_PATH=""

# Try common locations
if command -v python3 &> /dev/null; then
    PYTHON_PATH=$(which python3)
    echo "   Found python3: $PYTHON_PATH"
elif command -v python &> /dev/null; then
    PYTHON_PATH=$(which python)
    echo "   Found python: $PYTHON_PATH"
else
    echo "   ERROR: Python not found!"
    echo "   Please install Python 3.8+ first."
    exit 1
fi

# Verify Python version
PYTHON_VERSION=$($PYTHON_PATH --version 2>&1)
echo "   Version: $PYTHON_VERSION"

# Check if tissue_simulator is installed
echo ""
echo "2. Checking tissue_simulator package..."
if $PYTHON_PATH -c "import tissue_simulator" 2>/dev/null; then
    echo "   ✓ Package installed"
else
    echo "   ✗ Package not installed"
    echo "   Installing now..."
    cd /Users/cramere/tissue_simulator
    $PYTHON_PATH -m pip install -e . --quiet
    if [ $? -eq 0 ]; then
        echo "   ✓ Package installed successfully"
    else
        echo "   ✗ Installation failed"
        exit 1
    fi
fi

# Check if MCP is installed
echo ""
echo "3. Checking MCP library..."
if $PYTHON_PATH -c "import mcp" 2>/dev/null; then
    echo "   ✓ MCP installed"
else
    echo "   ✗ MCP not installed"
    echo "   Installing now..."
    $PYTHON_PATH -m pip install mcp --quiet
    if [ $? -eq 0 ]; then
        echo "   ✓ MCP installed successfully"
    else
        echo "   ✗ Installation failed"
        exit 1
    fi
fi

# Test the server
echo ""
echo "4. Testing MCP server..."
timeout 2 $PYTHON_PATH /Users/cramere/tissue_simulator/run_mcp_server.py 2>&1 | head -1 &
wait $!
if [ $? -eq 124 ]; then
    echo "   ✓ Server starts correctly"
else
    echo "   Testing server startup (should timeout)..."
fi

# Create config file
echo ""
echo "5. Creating Claude Desktop configuration..."
CONFIG_DIR="$HOME/Library/Application Support/Claude"
CONFIG_FILE="$CONFIG_DIR/claude_desktop_config.json"

# Create directory if it doesn't exist
mkdir -p "$CONFIG_DIR"

# Check if config file exists
if [ -f "$CONFIG_FILE" ]; then
    echo "   Config file exists, creating backup..."
    cp "$CONFIG_FILE" "$CONFIG_FILE.backup.$(date +%Y%m%d_%H%M%S)"
    echo "   Backup saved: $CONFIG_FILE.backup.*"
fi

# Write config
cat > "$CONFIG_FILE" << EOF
{
  "mcpServers": {
    "tissue-simulator": {
      "command": "$PYTHON_PATH",
      "args": [
        "/Users/cramere/tissue_simulator/run_mcp_server.py"
      ],
      "env": {}
    }
  }
}
EOF

echo "   ✓ Configuration created"
echo ""
echo "Configuration details:"
echo "   File: $CONFIG_FILE"
echo "   Python: $PYTHON_PATH"
echo "   Script: /Users/cramere/tissue_simulator/run_mcp_server.py"

# Validate JSON
echo ""
echo "6. Validating configuration..."
if $PYTHON_PATH -c "import json; json.load(open('$CONFIG_FILE'))" 2>/dev/null; then
    echo "   ✓ Valid JSON"
else
    echo "   ✗ Invalid JSON"
    exit 1
fi

echo ""
echo "================================================"
echo "✓ Configuration Complete!"
echo "================================================"
echo ""
echo "Next steps:"
echo "1. Restart Claude Desktop completely"
echo "2. The 'tissue-simulator' tools should now be available"
echo ""
echo "To verify:"
echo "   cat '$CONFIG_FILE'"
echo ""
echo "To test manually:"
echo "   $PYTHON_PATH /Users/cramere/tissue_simulator/run_mcp_server.py"
echo "   (Press Ctrl+C to stop)"
echo ""
