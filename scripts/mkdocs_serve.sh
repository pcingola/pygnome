#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

echo "Running local docs server..."
cd "${PROJECT_DIR}"
pip install -e ".[docs]"
mkdocs serve
