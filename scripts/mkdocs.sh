#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Run unit tests
echo "Building docs..."
cd "${PROJECT_DIR}"
mkdocs build
