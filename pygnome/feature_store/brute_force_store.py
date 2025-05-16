from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore, MAX_SAMPLES_TO_SHOW
from pygnome.feature_store.genomic_feature_store_protocol import MAX_DISTANCE
from pygnome.genomics.genomic_feature import GenomicFeature


class BruteForceFeatureStore(ChromosomeFeatureStore):
    """
    A naive brute-force implementation for genomic feature storage.
    This is not memory efficient and is not recommended for large datasets.
    It is primarily for testing purposes.
    """
    
    def __init__(self, chromosome: str):
        super().__init__(chromosome)
    
    def add(self, feature: GenomicFeature) -> None:
        """Add a feature to this chromosome's store."""
        self.features.append(feature)
    
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        return [f for f in self.features if f.intersects_point(position)]
    
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        return [f for f in self.features if f.intersects_interval(start, end)]

    def get_nearest(self, position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature | None:
        """Get the nearest feature to a specific position."""
        nearest_feature = None
        min_distance = max_distance
        
        for feature in self.features:
            distance = feature.distance(position)
            if distance < min_distance:
                min_distance = distance
                nearest_feature = feature
        return nearest_feature
    
    def __str__(self) -> str:
        """Return a string representation of the brute force feature store."""
        status = "building" if self.index_build_mode else "built" if self.index_finished else "uninitialized"
        
        # Sample of features
        sample_size = min(MAX_SAMPLES_TO_SHOW, len(self.features))
        sample = self.features[:sample_size] if self.features else []
        sample_str = ""
        if sample:
            sample_str = f", sample: [{', '.join(str(f.id) for f in sample)}" + (", ...]" if len(self.features) > sample_size else "]")
        
        return (f"BruteForceFeatureStore(chromosome='{self.chromosome}', features={len(self.features)}, "
                f"status={status}{sample_str})")
    
    def __repr__(self) -> str:
        """Return a string representation of the brute force feature store."""
        return self.__str__()
        
        