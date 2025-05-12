"""Phase enumeration for coding sequences."""

from enum import IntEnum


class Phase(IntEnum):
    """
    Frame phase for coding sequences.
    
    The phase indicates the number of bases that should be skipped to reach
    the first base of the next codon.
    """
    ZERO = 0  # No bases should be skipped
    ONE = 1   # Skip one base to reach the next codon
    TWO = 2   # Skip two bases to reach the next codon