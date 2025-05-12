"""Base interfaces for genomic feature storage and search."""

from abc import ABC, abstractmethod
from typing import Protocol, runtime_checkable

from ..genomic_feature import GenomicFeature


@runtime_checkable
class FeatureStore(Protocol):
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
    
    def get_nearest(self, chrom: str, position: int, max_distance: int = None) -> GenomicFeature | None:
        """Get the nearest feature to the given position."""
        ...


class ChromosomeFeatureStore(ABC):
    """Base class for chromosome-specific feature storage."""
    
    def __init__(self):
        self.features: list[GenomicFeature] = []
    
    def add_feature(self, feature: GenomicFeature) -> None:
        """Add a feature to this chromosome's store."""
        self.features.append(feature)
    
    @abstractmethod
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        pass
    
    @abstractmethod
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        pass
    
    def get_nearest(self, position: int, max_distance: int = None) -> GenomicFeature | None:
        """Get the nearest feature to the given position."""
        # Start with a small window and expand until we find features
        window = 1
        while max_distance is None or window <= max_distance:
            features = self.get_by_interval(position - window, position + window)
            if features:
                # Find the feature with the smallest distance to the position
                return min(features, key=lambda f: min(abs(f.start - position), abs(f.end - position)))
            window *= 2
        
        return None