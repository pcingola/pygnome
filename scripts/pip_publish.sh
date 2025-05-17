#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Clean previous builds
echo "Cleaning previous builds..."
rm -rf build/ dist/ *.egg-info/

# Build the package
echo "Building package..."
python -m build

# Check the package
echo "Checking package with twine..."
python -m twine check dist/*

echo ""
echo "Package built successfully!"
echo ""
echo "To publish to TestPyPI (recommended for testing):"
echo "python -m twine upload --repository testpypi dist/*"
echo ""
echo "To publish to PyPI (production):"
echo "python -m twine upload dist/*"
echo ""
echo "Note: You'll need to have twine installed and be registered on PyPI/TestPyPI"
echo "      pip install twine"
echo "      Register at: https://test.pypi.org/ and https://pypi.org/"