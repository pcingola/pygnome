#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Run all tests
echo "Running all tests..."
python -m unittest discover -s "${PROJECT_DIR}/pygnome/tests" -p "test_*.py"