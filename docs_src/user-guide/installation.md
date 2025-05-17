# Installation Guide

This guide provides detailed instructions for installing PyGnome and its dependencies.

## Requirements

PyGnome requires:

- Python 3.9 or higher
- NumPy for efficient array operations
- Pathlib for file path handling

## Basic Installation

The simplest way to install PyGnome is using pip:

```bash
pip install pygnome
```

This will install PyGnome and its required dependencies.

## Development Installation

For development or contributing to PyGnome, you can install from the source code:

```bash
# Clone the repository
git clone https://github.com/pcingola/pygnome.git
cd pygnome

# Install in development mode
pip install -e .
```

Installing in development mode (`-e`) allows you to modify the code and have the changes immediately available without reinstalling.

## Installing Optional Dependencies

PyGnome has several optional dependencies for specific functionality:

### For VCF Parsing

For working with VCF files, you may want to install pysam:

```bash
pip install pysam
```

### For Visualization

For visualization features, install matplotlib:

```bash
pip install matplotlib
```

## Verifying Installation

You can verify that PyGnome is installed correctly by importing it in Python:

```python
import pygnome
print(f"PyGnome is installed")
```

## Troubleshooting

### Common Issues

#### Missing NumPy

If you encounter an error about missing NumPy, install it separately:

```bash
pip install numpy
```

#### Compilation Issues

If you encounter compilation issues with C extensions, ensure you have a C compiler installed:

- **Linux**: Install GCC with `sudo apt-get install build-essential` (Ubuntu/Debian) or `sudo yum install gcc` (CentOS/RHEL)
- **macOS**: Install Xcode Command Line Tools with `xcode-select --install`
- **Windows**: Install Visual C++ Build Tools

## Upgrading

To upgrade PyGnome to the latest version:

```bash
pip install --upgrade pygnome
```

## Uninstalling

To uninstall PyGnome:

```bash
pip uninstall pygnome