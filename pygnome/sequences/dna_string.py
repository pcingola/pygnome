"""
DnaString class for efficient 2-bit representation of DNA sequences.
"""

from typing import Dict

from .base_sequence import BaseSequence


class DnaString(BaseSequence):
    """
    Efficient 2-bit representation of DNA sequences.
    
    This class stores DNA sequences (A, C, G, T) using 2 bits per nucleotide,
    allowing 16 nucleotides to be packed into a single 32-bit integer.
    
    Encoding:
    - A (Adenine) = 00 (0)
    - C (Cytosine) = 01 (1)
    - G (Guanine) = 10 (2)
    - T (Thymine) = 11 (3)
    
    Ambiguous nucleotides (N, R, Y, etc.) are not directly supported in the
    2-bit representation and will be converted to 'A' by default.
    """
    
    @property
    def _NT_TO_BITS(self) -> Dict[str, int]:
        """Mapping from nucleotide characters to 2-bit values."""
        return {
            'A': 0,  # 00
            'C': 1,  # 01
            'G': 2,  # 10
            'T': 3,  # 11
        }
    
    @property
    def _BITS_TO_NT(self) -> Dict[int, str]:
        """Mapping from 2-bit values to nucleotide characters."""
        return {
            0: 'A',  # 00
            1: 'C',  # 01
            2: 'G',  # 10
            3: 'T',  # 11
        }
    
    @property
    def _class_name(self) -> str:
        """Name of the class for use in __repr__."""
        return "DnaString"