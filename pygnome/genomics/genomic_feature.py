"""Base class for genomic features."""

from pydantic import BaseModel, Field, validator

from .strand import Strand


class GenomicFeature(BaseModel):
    """
    Base class for all genomic features.
    
    Uses 0-based, half-open interval coordinates [start, end) where start is
    included but end is excluded. Length is calculated as (end - start).
    """
    id: str
    chrom: str
    start: int = Field(..., ge=0, description="Start position (0-based, inclusive)")
    end: int = Field(..., ge=0, description="End position (0-based, exclusive)")
    strand: Strand
    
    @validator('end')
    def end_must_be_after_start(cls, v, values):
        """Validate that end position is not before start position."""
        if 'start' in values and v < values['start']:
            raise ValueError(f"End position ({v}) must be >= start position ({values['start']})")
        return v
    
    def __str__(self) -> str:
        """Return a string representation of the feature."""
        return f"{self.__class__.__name__}({self.id}, {self.chrom}:{self.start}-{self.end}:{self.strand})"
    
    @property
    def length(self) -> int:
        """Return the length of the feature (end - start)."""
        return self.end - self.start