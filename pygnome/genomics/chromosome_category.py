"""Chromosome category enumeration."""

from enum import IntEnum


class ChromosomeCategory(IntEnum):
    """Categories for chromosome sorting."""
    NUMERIC = 0      # Standard numbered chromosomes (1-22)
    SEX = 1          # Sex chromosomes (X, Y)
    MITOCHONDRIAL = 2  # Mitochondrial chromosome (M, MT)
    OTHER = 3        # Other chromosomes