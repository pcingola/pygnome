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

# Generate API reference pages for each module
for module in modules:
    module_path = module.replace(".", "/")
    doc_path = Path("api", f"{module.split('.')[-1]}.md")
    
    with mkdocs_gen_files.open(doc_path, "w") as fd:
        fd.write(f"# {module}\n\n")
        fd.write(f"::: {module}\n")
    
    nav[module.split(".")[1:]] = doc_path.as_posix()

# Generate the navigation
with mkdocs_gen_files.open("api/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
