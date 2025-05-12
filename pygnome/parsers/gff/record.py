"""
Record class for representing genomic features from GFF/GTF files.
"""

from dataclasses import dataclass, field
from typing import Any

from ...genomics.strand import Strand


@dataclass
class Record:
    """
    Represents a single feature/annotation line from a GFF/GTF file.
    
    This class provides a unified representation of features across different
    formats (GFF2, GFF3, GTF) with methods to access attributes and convert
    between formats.
    """
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float | None = None
    strand: Strand | None = None
    phase: int | None = None
    attributes: dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate and convert types after initialization."""
        # Convert start and end to integers
        self.start = int(self.start)
        self.end = int(self.end)
        
        # Convert score to float if not None or '.'
        if self.score and self.score != '.':
            self.score = float(self.score)
        else:
            self.score = None
            
        # Convert phase to int if not None or '.'
        if self.phase and self.phase != '.':
            self.phase = int(self.phase)
        else:
            self.phase = None
            
        # Convert strand to Strand enum if not None or '.'
        if isinstance(self.strand, str):
            if self.strand == '.':
                self.strand = None
            else:
                try:
                    self.strand = Strand(self.strand)
                except ValueError:
                    raise ValueError(f"Invalid strand value: {self.strand}")
    
    def get_attribute(self, name: str, default=None) -> Any:
        """Get an attribute value by name."""
        return self.attributes.get(name, default)
    
    def set_attribute(self, name: str, value: Any) -> None:
        """Set an attribute value."""
        self.attributes[name] = value
    
    def __str__(self) -> str:
        """String representation of the record."""
        attrs = "; ".join(f"{k}={v}" for k, v in self.attributes.items())
        return (f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t"
                f"{self.score or '.'}\t{self.strand.value if self.strand else '.'}\t"
                f"{self.phase or '.'}\t{attrs}")