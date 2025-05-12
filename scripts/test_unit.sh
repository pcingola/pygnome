#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Run unit tests
echo "Running unit tests..."
python -m unittest discover -s "${PROJECT_DIR}/tests" -p "test_*.py"