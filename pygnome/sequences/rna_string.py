"""
RnaString class for efficient 2-bit representation of RNA sequences.
"""

from pygnome.sequences.alphabets import RNA_BITS_TO_NT, RNA_NT_TO_BITS
from .base_sequence import BaseSequence


class RnaString(BaseSequence):
    """
    Efficient 2-bit representation of RNA sequences.
    
    This class stores RNA sequences (A, C, G, U) using 2 bits per nucleotide,
    allowing 16 nucleotides to be packed into a single 32-bit integer.
    
    Ambiguous nucleotides (N, R, Y, etc.) are not directly supported in the
    2-bit representation and will be converted to 'A' by default.
    """
    
    @property
    def _NT_TO_BITS(self) -> dict[str, int]:
        """Mapping from nucleotide characters to 2-bit values."""
        return RNA_NT_TO_BITS
    
    @property
    def _BITS_TO_NT(self) -> list[str]:
        return RNA_BITS_TO_NT
    