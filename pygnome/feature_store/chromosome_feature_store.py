"""Main implementation of the genomic feature store."""

from abc import ABC, abstractmethod

from pygnome.feature_store.genomic_feature_store_protocol import MAX_DISTANCE
from pygnome.genomics.genomic_feature import GenomicFeature


class ChromosomeFeatureStore(ABC):
    """
    Base class for chromosome-specific genomic feature storage.
    It has a list of features and provides methods to add and query them.
    
    Sub-classes typically implement more efficient search mechanisms, by adding an index to the features.

    Usage:
        chrom_store = ChromosomeFeatureStore("chr1")

        # When adding features you can use in a context manager to ensure index build mode is active
        with chrom_store:
            for feature in features:
                chrom_store.add(feature)

        # Alternatively, you can call index_build_start() and index_build_end() manually
        chrom_store.index_build_start()
        for feature in features:
            chrom_store.add(feature)
        chrom_store.index_build_end()
    """
    
    def __init__(self, chromosome: str) -> None:
        self.features: list[GenomicFeature] = []
        self.index_build_mode = False
        self.index_finished = False
        self.chromosome = chromosome
    
    def add(self, feature: GenomicFeature) -> None:
        """Add a feature to this chromosome's store."""
        if not self.index_build_mode:
            raise RuntimeError("Index build mode is not active. Call index_build_start() before adding features.")
        self.features.append(feature)
    
    @abstractmethod
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        pass
    
    @abstractmethod
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        pass
    
    def get_features(self) -> list[GenomicFeature]:
        """Get all features."""
        return self.features

    def __getitem__(self, index: int) -> list[GenomicFeature]:
        """Get all features at a specific index."""
        return self.features[index]


    def get_nearest(self, position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature | None:
        """Get the nearest feature to the given position."""
        # Start with a small window and expand until we find features
        window = 1
        while window <= max_distance:
            features = self.get_by_interval(position - window, position + window)
            if features:
                # Find the feature with the smallest distance to the position
                closest_feature = features[0]
                min_distance = closest_feature.distance(position)
                for feature in features:
                    distance = feature.distance(position)
                    if distance < min_distance:
                        min_distance = distance
                        closest_feature = feature
                return closest_feature
            window *= 2        
        return None
    
    def index_build_start(self) -> None:
        """ This method must be called to build the index before adding features. """
        self.index_build_mode = True
        self.index_finished = False

    def index_build_end(self) -> None:
        """ This method must be called to finish building the index, after adding features, but before quering. """
        self.index_build_mode = False
        self.index_finished = True

    def __enter__(self):
        """Enter the context manager."""
        self.index_build_start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context manager."""
        self.index_build_end()