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

# To publish to TestPyPI (recommended for testing):"
# python -m twine upload --repository testpypi dist/*"
# "
# To publish to PyPI (production):"
# python -m twine upload dist/*"
# "
# Note: You'll need to have twine installed and be registered on PyPI/TestPyPI"
#       pip install twine"
#       Register at: https://test.pypi.org/ and https://pypi.org/"

echo "Publishing package to PyPI..."
python -m twine upload dist/*