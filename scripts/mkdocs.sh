#!/bin/bash -eu
set -o pipefail

# Source the config script to set up the environment
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

cd "${PROJECT_DIR}"

# Parse command line arguments
ACTION=${1:-"build"}

# Install documentation dependencies if needed
if ! python -c "import mkdocs" &>/dev/null; then
    echo "Installing documentation dependencies..."
    uv pip install -e ".[docs]"
fi

# Execute the requested action
case "$ACTION" in
    "build")
        echo "Building documentation..."
        mkdocs build
        ;;
    "serve")
        echo "Serving documentation locally..."
        mkdocs serve
        ;;
    "clean")
        echo "Cleaning documentation build..."
        rm -rf docs/
        mkdir -p docs
        touch docs/.nojekyll  # Prevent GitHub Pages from using Jekyll
        ;;
    "publish")
        echo "Pushing documentation..."
        touch docs/.nojekyll
        git add docs/
        git commit -m "Update documentation"
        git push origin main
        ;;
    *)
        echo "Unknown action: $ACTION"
        echo "Usage: $0 [build|serve|clean|publish]"
        exit 1
        ;;
esac
