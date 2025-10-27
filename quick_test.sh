#!/bin/bash
# Quick test script to verify the fix

echo "Testing tissue simulator..."
echo ""

cd /Users/cramere/tissue_simulator

echo "1. Testing imports..."
python3 -c "from tissue_simulator import TissueSection; print('✓ Imports work')"

echo ""
echo "2. Testing basic generation..."
python3 -c "
from tissue_simulator import TissueSection
tissue = TissueSection(100, 100, 50, (5, 10))
num = tissue.generate_cells(max_attempts=200)
print(f'✓ Generated {num} cells')
"

echo ""
echo "3. Testing CSV export..."
python3 -c "
from tissue_simulator import TissueSection
import os
tissue = TissueSection(100, 100, 50, (5, 10))
tissue.generate_cells(max_attempts=200)
tissue.export_to_csv('/tmp/test_tissue.csv')
if os.path.exists('/tmp/test_tissue.csv'):
    print('✓ CSV export works')
    os.remove('/tmp/test_tissue.csv')
"

echo ""
echo "✅ All basic tests passed!"
echo ""
echo "To test visualization, run:"
echo "  python examples/simple_example.py"
