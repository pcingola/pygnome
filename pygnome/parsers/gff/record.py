"""
GffRecord class for representing genomic features from GFF/GTF files.
"""

from dataclasses import dataclass, field
from typing import Any

from ...genomics.genomic_feature import GenomicFeature
from ...genomics.strand import Strand


@dataclass
class GffRecord(GenomicFeature):
    """
    Represents a single feature/annotation line from a GFF/GTF file.
    
    This class provides a unified representation of features across different
    formats (GFF2, GFF3, GTF) with methods to access attributes and convert
    between formats.
    
    Inherits from GenomicFeature to provide standard genomic coordinate functionality.
    """
    source: str
    type: str
    score: float | None = None
    phase: int | None = None
    attributes: dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate and convert types after initialization."""
        # Set a default ID if not provided
        if not self.id:
            self.id = f"{self.type}_{self.chrom}_{self.start}"
        
        # Convert start and end to integers
        self.start = int(self.start)
        self.end = int(self.end)
        
        # Convert score to float if not None or '.'
        if self.score and self.score != '.':
            self.score = float(self.score)
        else:
            self.score = None
            
        # Convert phase to int if not None or '.'
        if self.phase == 0 or (self.phase and self.phase != '.'):
            self.phase = int(self.phase)
        else:
            self.phase = None
            
        # Convert strand to Strand enum if not None or '.'
        if isinstance(self.strand, str):
            if self.strand == '.':
                self.strand = Strand.UNSTRANDED
            else:
                try:
                    self.strand = Strand(self.strand)
                except ValueError:
                    raise ValueError(f"Invalid strand value: {self.strand}")
        elif self.strand is None:
            self.strand = Strand.UNSTRANDED
            
        # Call the parent class's __post_init__ to validate
        super().__post_init__()
    
    def get_attribute(self, name: str, default=None) -> Any:
        """Get an attribute value by name."""
        return self.attributes.get(name, default)
    
    def set_attribute(self, name: str, value: Any) -> None:
        """Set an attribute value."""
        self.attributes[name] = value
    
    def __str__(self) -> str:
        """String representation of the record."""
        attrs = "; ".join(f"{k}={v}" for k, v in self.attributes.items())
        return (f"{self.chrom}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t"
                f"{self.score or '.'}\t{self.strand.value if self.strand != Strand.UNSTRANDED else '.'}\t"
                f"{self.phase or '.'}\t{attrs}")