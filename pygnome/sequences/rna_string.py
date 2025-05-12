"""
RnaString class for efficient 2-bit representation of RNA sequences.
"""

from typing import Dict

from .base_sequence import BaseSequence


class RnaString(BaseSequence):
    """
    Efficient 2-bit representation of RNA sequences.
    
    This class stores RNA sequences (A, C, G, U) using 2 bits per nucleotide,
    allowing 16 nucleotides to be packed into a single 32-bit integer.
    
    Encoding:
    - A (Adenine) = 00 (0)
    - C (Cytosine) = 01 (1)
    - G (Guanine) = 10 (2)
    - U (Uracil) = 11 (3)
    
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
            'U': 3,  # 11
        }
    
    @property
    def _BITS_TO_NT(self) -> Dict[int, str]:
        """Mapping from 2-bit values to nucleotide characters."""
        return {
            0: 'A',  # 00
            1: 'C',  # 01
            2: 'G',  # 10
            3: 'U',  # 11
        }
    
    @property
    def _class_name(self) -> str:
        """Name of the class for use in __repr__."""
        return "RnaString"
    
    def _normalize_nucleotide(self, nt: str) -> str:
        """
        Normalize a nucleotide character, converting 'T' to 'U'.
        
        Args:
            nt: A nucleotide character
            
        Returns:
            Normalized nucleotide character
        """
        nt = super()._normalize_nucleotide(nt)
        # Convert 'T' to 'U' for RNA
        if nt == 'T':
            return 'U'
        return nt