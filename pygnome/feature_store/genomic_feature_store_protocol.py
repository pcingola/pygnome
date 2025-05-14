"""Base interfaces for genomic feature storage and search."""

from typing import Protocol, runtime_checkable

# Maximum distance to search for 'closest' feature
MAX_DISTANCE = 2147483647

from ..genomics.genomic_feature import GenomicFeature


@runtime_checkable
class GenomicFeatureStoreProtocol(Protocol):
    """Protocol for genomic feature stores."""
    
    def add_feature(self, feature: GenomicFeature) -> None:
        """Add a genomic feature to the store."""
        ...
    
    def get_by_position(self, chrom: str, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        ...
    
    def get_by_interval(self, chrom: str, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        ...
    
    def get_nearest(self, chrom: str, position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature | None:
        """Get the nearest feature to the given position."""
        ...
