"""Generate API reference pages."""

import os
from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()

# Define the modules to document
modules = [
    "pygnome.genomics",
    "pygnome.parsers",
    "pygnome.sequences",
    "pygnome.feature_store",
]

# List of manually created files that should not be overwritten
manual_files = [
    Path("api/genomics.md"),
    Path("api/parsers.md"),
    Path("api/sequences.md"),
    Path("api/feature-store.md"),
]

# Generate API reference pages for each module
for module in modules:
    module_path = module.replace(".", "/")
    
    # Handle feature_store module differently
    if module.endswith("feature_store"):
        doc_path = Path("api/feature-store.md")
    else:
        doc_path = Path("api", f"{module.split('.')[-1]}.md")
    
    # Only generate the file if it doesn't exist in manual_files
    if doc_path not in manual_files:
        with mkdocs_gen_files.open(doc_path, "w") as fd:
            fd.write(f"# {module}\n\n")
            fd.write(f"::: {module}\n")
    
    # Always add the file to navigation
    nav[tuple(module.split(".")[1:])] = doc_path.as_posix()

# Generate the navigation
with mkdocs_gen_files.open("api/SUMMARY.md", "w") as nav_file:
    # The issue is that the SUMMARY.md file is in the api/ directory,
    # but the paths in the navigation are also prefixed with api/
    # We'll modify how we add paths to the navigation
    
    # Clear the navigation
    nav = mkdocs_gen_files.Nav()
    
    # Re-add the paths without the api/ prefix for the SUMMARY.md file
    for module in modules:
        if module.endswith("feature_store"):
            doc_path = "feature-store.md"  # Relative to api/ directory
        else:
            doc_path = f"{module.split('.')[-1]}.md"  # Relative to api/ directory
        
        nav[tuple(module.split(".")[1:])] = doc_path
    
    nav_file.writelines(nav.build_literate_nav())
