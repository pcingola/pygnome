"""Main implementation of the genomic feature store."""

import pickle
from enum import Enum
from pathlib import Path

from pygnome.feature_store.binned_store import BinnedGenomicStore
from pygnome.feature_store.brute_force_store import BruteForceFeatureStore
from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore
from pygnome.feature_store.genomic_feature_store_protocol import MAX_DISTANCE, GenomicFeatureStoreProtocol
from pygnome.feature_store.interval_tree_store import IntervalTreeStore
from pygnome.feature_store.msi_chromosome_store import MsiChromosomeStore
from pygnome.genomics.genomic_feature import GenomicFeature



class StoreType(str, Enum):
    """Types of genomic feature stores."""
    INTERVAL_TREE = "interval_tree"
    BINNED = "binned"
    BRUTE_FORCE = "brute_force"
    MSI = "msi"


class GenomicFeatureStore(GenomicFeatureStoreProtocol):
    """Implementation of the genomic feature store."""

    def __init__(self, store_type: StoreType | str = StoreType.INTERVAL_TREE, bin_size: int = 100000):
        """
        Initialize the genomic feature store.
        
        Args:
            store_type: Type of store to use
            bin_size: Size of bins for the binned store
        """
        self.chromosomes: dict[str, ChromosomeFeatureStore] = {}
        
        if isinstance(store_type, str):
            try:
                self.store_type = StoreType(store_type)
            except ValueError:
                raise ValueError(f"Unknown store type: {store_type}")
        else:
            self.store_type = store_type
            
        self.bin_size = bin_size

    def __enter__(self):
        """Enter the context manager when adding features"""
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context manager after adding features, ensures all indeces are built"""
        # No special cleanup needed for this store
        for chrom in self.chromosomes.values():
            chrom.index_build_end()
        return False

    def _get_or_create_chrom_store(self, chrom: str) -> ChromosomeFeatureStore:
        """Get or create a chromosome store."""
        if chrom not in self.chromosomes:
            match self.store_type:
                case StoreType.INTERVAL_TREE:
                    self.chromosomes[chrom] = IntervalTreeStore(chrom)
                case StoreType.BINNED:
                    self.chromosomes[chrom] = BinnedGenomicStore(chrom, self.bin_size)
                case StoreType.BRUTE_FORCE:
                    self.chromosomes[chrom] = BruteForceFeatureStore(chrom)
                case StoreType.MSI:
                    self.chromosomes[chrom] = MsiChromosomeStore(chrom)
                case _:
                    raise ValueError(f"Unknown store type: {self.store_type}")
            # Make sure the chromosome store is in 'index build' mode
            self.chromosomes[chrom].index_build_start()
        return self.chromosomes[chrom]

    def add(self, feature: GenomicFeature) -> None:
        """Add a genomic feature to the store."""
        chrom_store = self._get_or_create_chrom_store(feature.chrom)
        chrom_store.add(feature)

    def add_features(self, features: list[GenomicFeature]) -> None:
        """Add multiple genomic features to the store."""
        for feature in features:
            self.add(feature)

    def get_by_position(self, chrom: str, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        if chrom not in self.chromosomes:
            return []
        return self.chromosomes[chrom].get_by_position(position)

    def get_by_interval(self, chrom: str, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        if chrom not in self.chromosomes:
            return []
        return self.chromosomes[chrom].get_by_interval(start, end)

    def get_nearest(self, chrom: str, position: int, max_distance: int = MAX_DISTANCE) -> GenomicFeature | None:
        """Get the nearest feature to the given position."""
        if chrom not in self.chromosomes:
            return None
        return self.chromosomes[chrom].get_nearest(position, max_distance)
    
    def get_chromosomes(self) -> list[str]:
        """Get all chromosome names in the store."""
        return list(self.chromosomes.keys())
    
    def __getitem__(self, chrom):
        return self.chromosomes.get(chrom)
    
    def get_feature_count(self, chrom: str | None = None) -> int:
        """
        Get the number of features in the store.
        
        Args:
            chrom: If provided, count only features in this chromosome
        """
        if chrom is not None:
            if chrom not in self.chromosomes:
                return 0
            return len(self.chromosomes[chrom].features)
        
        # Count across all chromosomes
        return sum(len(store) for store in self.chromosomes.values())

    def __iterator__(self):
        """Iterate over all features in the store."""
        return self
    
    def __iter__(self):
        """Iterate over all features in the store."""
        return iter(self.chromosomes.values())
    
    def __str__(self) -> str:
        """Return a string representation of the feature store."""
        chrom_counts = [f"{chrom}: {len(store)}" for chrom, store in self.chromosomes.items()]
        keys = sorted(list(self.chromosomes.keys()))
        stores = [f"\t{self.chromosomes[chr]}" for chr in keys]
        return (f"GenomicFeatureStore(type={self.store_type.value}, "
                + f"chromosomes={len(self.chromosomes)}, "
                + f"features={self.get_feature_count()}, "
                + f"{', '.join(chrom_counts)})\n"
                + '\n'.join(stores)
                )

    def __repr__(self) -> str:
        return self.__str__()
    
    def trim(self) -> None:
        """
        Trim internal data structures to reduce memory usage.
        
        This method calls trim() on all chromosome stores
        to reduce memory usage before serialization.
        """
        for chrom_store in self.chromosomes.values():
            chrom_store.trim()
    
    def save(self, filepath: Path) -> None:
        """
        Save the genomic feature store to a file using pickle.
        
        Args:
            filepath: Path to save the store to
        """
        # Ensure the directory exists
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        # Make sure all indices are built before saving
        for chrom in self.chromosomes.values():
            if chrom.index_build_mode:
                chrom.index_build_end()
        
        # Trim all stores to reduce memory usage before serialization
        self.trim()
        
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def load(cls, filepath: Path) -> 'GenomicFeatureStore':
        """
        Load a genomic feature store from a file.
        
        Args:
            filepath: Path to load the store from
            
        Returns:
            The loaded genomic feature store
        """
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        with open(filepath, 'rb') as f:
            store = pickle.load(f)
        
        # Validate the loaded object
        if not isinstance(store, GenomicFeatureStore):
            raise TypeError(f"Loaded object is not a GenomicFeatureStore: {type(store)}")
        
        return store
