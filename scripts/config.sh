#!/bin/bash -eu
set -o pipefail

# Define project directory
export PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Add the project root to PYTHONPATH to make the package importable
export PYTHONPATH="${PROJECT_DIR}:${PYTHONPATH:-}"

# Load environment variables from .env if it exists
if [ -f "${PROJECT_DIR}/.env" ]; then
    set -o allexport
    source "${PROJECT_DIR}/.env"
    set +o allexport
fi

# Activate virtual environment if it exists
if [ -d "${PROJECT_DIR}/.venv" ]; then
    source "${PROJECT_DIR}/.venv/bin/activate"
fi

# Check if the package is installed in development mode
if ! python -c "import pygnome" &>/dev/null; then
    echo "Installing pygnome package in development mode..."
    cd "${PROJECT_DIR}" && uv pip install -e .
fi