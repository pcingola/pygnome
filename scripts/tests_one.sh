#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# Run one test
cd "${PROJECT_DIR}"
python -m unittest \
    "pygnome.tests.search.test_msi_chromosome_store" \
    --failfast \
    -v
