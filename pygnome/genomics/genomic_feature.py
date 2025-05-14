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
    
    def intersects_point(self, position: int) -> bool:
        """Check if the feature intersects with a specific point."""
        # Special case for zero-length features
        if self.start == self.end:
            return position == self.start
        return self.start <= position < self.end

    def intersects_interval(self, start: int, end: int) -> bool:
        """Check if the feature intersects with a specific interval."""
        # Special case for zero-length features
        if self.start == self.end:
            return start <= self.start < end
        return not (self.end <= start or self.start >= end)
    
    def intersects(self, other: 'GenomicFeature') -> bool:
        """Check if the feature intersects with another feature."""
        return self.intersects_interval(other.start, other.end)
    
    def contains(self, other: 'GenomicFeature') -> bool:
        """Check if the feature contains another feature."""
        return self.start <= other.start and other.end < self.end 

    def distance(self, position: int) -> int:
        """Calculate the distance from a point to the feature."""
        if position < self.start:
            return self.start - position
        elif position >= self.end:
            return (position - self.end) + 1
        else:
            return 0