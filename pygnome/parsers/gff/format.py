"""
Enumeration of supported genomic annotation file formats.
"""

from enum import Enum, auto


class Format(str, Enum):
    """
    Enumeration of supported genomic annotation file formats.
    
    Values:
        GFF2: Gene Finding Format version 2
        GFF3: Generic Feature Format version 3
        GTF: Gene Transfer Format (GTF2.2)
    """
    GFF2 = "gff2"
    GFF3 = "gff3"
    GTF = "gtf"
    
    def __str__(self) -> str:
        """Return the string value of the format."""
        return self.value