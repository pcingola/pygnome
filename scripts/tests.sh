#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Run all tests
echo "Running all tests..."
# python -m unittest \
#     discover \
#     -s "${PROJECT_DIR}/pygnome/tests" \
#     -p "test_*.py" \
#     --failfast \
#     -v


# OK
# python -m unittest discover  -s "${PROJECT_DIR}/pygnome/tests/parsers"   -p "test_*.py"  --failfast
# python -m unittest discover  -s "${PROJECT_DIR}/pygnome/tests/sequences" -p "test_*.py"  --failfast

# Testing
# python -m unittest discover  -s "${PROJECT_DIR}/pygnome/tests/genomics"   -p "test_*.py"  -v --failfast
python -m unittest discover  -s "${PROJECT_DIR}/pygnome/tests/search"    -p "test_*.py" -v --failfast