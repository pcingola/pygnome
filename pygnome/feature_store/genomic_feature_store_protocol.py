"""Base interfaces for genomic feature storage and search."""

from typing import Protocol, runtime_checkable

# Maximum distance to search for 'closest' feature
MAX_DISTANCE = 2147483647

from ..genomics.genomic_feature import GenomicFeature
from ..genomics.variant import Variant


@runtime_checkable
class GenomicFeatureStoreProtocol(Protocol):
    """
    Protocol for genomic feature stores.

    Feature stores MUST also be context managers (__enter__ and __exit__ methods), which is used during index creation
    """
    
    def add(self, feature: GenomicFeature) -> None:
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
        
    def get_by_variant(self, chrom: str, position: int, ref: str, alt: str) -> list[GenomicFeature]:
        """
        Get all features that match the specific variant (chr, pos, ref, alt).
        
        Args:
            chrom: Chromosome name
            position: Position of the variant (0-based)
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            List of features that match the variant
        """
        ...
        
    def get_variant(self, variant: Variant) -> list[GenomicFeature]:
        """
        Get all features that match the given Variant object.
        
        Args:
            variant: A Variant object with chromosome, position, reference and alternate alleles
            
        Returns:
            List of features that match the variant
        """
        ...

    def __getitem__(self, chrom: str) -> 'ChromosomeFeatureStore': # type: ignore
        """Get a chromosome store by name."""
        ...

    def __iterator__(self):
        """Iterate over all chromosome stores."""
        ...
        
    def trim(self) -> None:
        """
        Trim internal data structures to reduce memory usage.
        
        This method should be called before serialization to reduce
        the size of pickled objects by removing unused allocated memory.
        """
        ...

