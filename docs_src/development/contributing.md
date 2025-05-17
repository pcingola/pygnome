# Contributing to PyGnome

Thank you for your interest in contributing to PyGnome! This document provides guidelines and instructions for contributing to the project.

## Getting Started

### Prerequisites

- Python 3.9 or higher
- Git
- A GitHub account

### Setting Up the Development Environment

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/your-username/pygnome.git
   cd pygnome
   ```
3. Set up a virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```
4. Install the package in development mode:
   ```bash
   pip install -e .
   ```
5. Install development dependencies:
   ```bash
   pip install -r requirements-dev.txt
   ```

## Development Workflow

### Branching Strategy

- `main`: The main development branch
- `feature/*`: Feature branches
- `bugfix/*`: Bug fix branches
- `release/*`: Release branches

### Creating a New Feature or Bug Fix

1. Create a new branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b bugfix/issue-number
   ```
2. Make your changes
3. Run tests to ensure your changes don't break existing functionality:
   ```bash
   ./scripts/tests.sh
   ```
4. Commit your changes:
   ```bash
   git add .
   git commit -m "Your descriptive commit message"
   ```
5. Push your branch to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```
6. Create a pull request on GitHub

## Coding Standards

### Python Code Style

- Use type hints in function and method signatures
- Use the new way of defining types (e.g., `dict`, `list`, `| None`, `any`)
- Write short but informative docstrings
- Use `Path` from `pathlib` instead of `str` / `import os` for file operations
- Use Pydantic classes or dataclasses instead of `dict` structures
- Use `Enum` objects instead of hard-coded string values
- Remove unused imports
- Sort functions and methods alphabetically (except underscores)
- Keep files under 500 lines; refactor if they get too long
- Write classes (and enums) in their own files whenever possible

### Testing

- Write unit tests for all new functionality
- Place test files in the `tests` directory
- Use Python's built-in `unittest` framework
- Run tests using the provided script:
  ```bash
  ./scripts/tests.sh
  ```

## Documentation

- Update documentation for any changes to the API
- Write clear and concise documentation
- Include examples where appropriate
- Use markdown for documentation files

## Pull Request Process

1. Ensure your code follows the coding standards
2. Update the documentation as necessary
3. Add or update tests as necessary
4. Ensure all tests pass
5. Submit your pull request
6. Address any feedback from reviewers

## Release Process

1. Update the version number in `pyproject.toml`
2. Update the changelog
3. Create a new release branch: `release/vX.Y.Z`
4. Create a pull request to merge the release branch into `main`
5. Once the pull request is approved and merged, create a new release on GitHub
6. Tag the release with the version number
7. Publish the package to PyPI