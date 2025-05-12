"""
Strand enumeration for genomic features.
"""

from enum import Enum, auto


class Strand(str, Enum):
    """
    Enumeration of possible strand values in genomic annotation formats.
    
    Values:
        POSITIVE: Forward strand ('+')
        NEGATIVE: Reverse strand ('-')
        UNSTRANDED: Feature is not stranded ('.')
        UNKNOWN: Strandedness is relevant but unknown ('?')
    """
    POSITIVE = "+"
    NEGATIVE = "-"
    UNSTRANDED = "."
    UNKNOWN = "?"
    
    def __str__(self) -> str:
        """Return the string value of the strand."""
        return self.value