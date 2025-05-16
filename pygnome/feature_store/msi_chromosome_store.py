"""
Efficient storage for MSI (Microsatellite Instability) sites.

This module provides classes for efficiently storing and querying millions of MSI sites
per chromosome using NumPy arrays and DnaStringArray for memory efficiency.
"""

import numpy as np

from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore, MAX_SAMPLES_TO_SHOW
from pygnome.genomics.genomic_feature import GenomicFeature
from pygnome.genomics.strand import Strand
from pygnome.sequences.dna_string_array import DnaStringArray
from pygnome.parsers.msi.msi_site_record import MsiSiteRecord


DEFAULT_BIN_SIZE = 100_000  # Default bin size
DEFAULT_FEATURE_COUNT = 1024  # Default number of features to allocate space for
MAX_LOOK_FORWARD = 1000 # Maximum number of features to look forward in the array

class MsiSiteCounter:
    """
    Helper class to count MSI sites and calculate statistics per chromosome.
    
    This class is used in a first pass over the data to count features and
    calculate statistics needed for efficient storage and querying.
    """
    
    def __init__(self, bin_size: int = DEFAULT_BIN_SIZE):
        """
        Initialize an MSI site counter.
        
        Args:
            bin_size: Size of each bin in base pairs
        """
        self.bin_size = bin_size
        self.feature_counts: dict[str, int] = {}  # chromosome -> count
        self.max_lengths_by_bin: dict[str, dict[int, int]] = {}  # chromosome -> {bin_id -> max_length}
        
    def add(self, feature: MsiSiteRecord) -> None:
        """
        Count a feature and track its statistics.
        
        Args:
            feature: The MSI site record to count
        """
        chrom = feature.chrom
        
        # Initialize chromosome data if needed
        if chrom not in self.feature_counts:
            self.feature_counts[chrom] = 0
            self.max_lengths_by_bin[chrom] = {}
            
        # Count the feature
        self.feature_counts[chrom] += 1
        
        # Track max feature length by bin
        bin_id = feature.start // self.bin_size
        current_max = self.max_lengths_by_bin[chrom].get(bin_id, 0)
        self.max_lengths_by_bin[chrom][bin_id] = max(current_max, feature.length)
        
    def get_count(self, chrom: str) -> int:
        """
        Get the feature count for a chromosome.
        
        Args:
            chrom: Chromosome name
            
        Returns:
            Number of features for the chromosome
        """
        return self.feature_counts.get(chrom, 0)
        
    def get_max_lengths(self, chrom: str) -> dict[int, int]:
        """
        Get the max feature lengths by bin for a chromosome.
        
        Args:
            chrom: Chromosome name
            
        Returns:
            dictionary mapping bin IDs to maximum feature length in that bin
        """
        return self.max_lengths_by_bin.get(chrom, {})


class MsiChromosomeStore(ChromosomeFeatureStore):
    """
    Efficient storage for millions of MSI sites in a chromosome.
    
    Uses NumPy arrays and DnaStringArray for memory efficiency.
    Implements binary search for efficient querying.
    """

    def __init__(self, chrom: str, feature_count: int = DEFAULT_FEATURE_COUNT, max_lengths_by_bin: dict[int, int] | None = None, bin_size: int = DEFAULT_BIN_SIZE):
        """
        Initialize an MSI chromosome store with pre-allocated arrays.
        
        Args:
            chrom: Name of the chromosome
            feature_count: Number of features to allocate space for
            max_lengths_by_bin: dictionary mapping bin IDs to maximum feature length in that bin
            bin_size: Size of each bin in base pairs
        """
        super().__init__(chromosome=chrom)  # Parent class uses 'chromosome' parameter
        self.features = None  # type: ignore # We won't use the features field
        self.bin_size = bin_size
        
        # Arrays for efficient storage
        assert feature_count > 0, "Feature count must be positive"
        self._starts = np.zeros(feature_count, dtype=np.uint32)
        self._ends = np.zeros(feature_count, dtype=np.uint32)
        self._repeat_unit_bases = DnaStringArray(initial_strings=feature_count)

        # Binning for efficient lookups
        self._max_feature_length_by_bin = max_lengths_by_bin or {}
        
        # Tracking
        self._feature_count = 0
        self._is_loaded = False
        
    def add(self, feature: MsiSiteRecord) -> None:
        """
        Add an MSI site to the store.
        
        Args:
            feature: The MSI site record to add
        """
        if not self.index_build_mode:
            raise RuntimeError("Index build mode is not active. Call index_build_start() before adding features.")        
        # Store in arrays, do not use the features list
        idx = self._feature_count
        self._starts[idx] = feature.start
        self._ends[idx] = feature.end
        self._repeat_unit_bases.add(feature.repeat_unit_bases)
        self._feature_count += 1
            
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """
        Get all features at a specific position using binary search.
        
        Args:
            position: The position to query
            
        Returns:
            List of features that contain the position
        """
        if not self._is_loaded or self._starts is None or len(self._starts) == 0:
            return []
            
        # Find the first feature that starts at or after the position
        idx = self._binary_search_position(position)
        
        # Look backward for features that contain this position
        bin_id = position // self.bin_size
        max_length = self._max_feature_length_by_bin.get(bin_id, 0)
        
        # Calculate how far back to look based on max feature length
        look_back_idx = idx
        while look_back_idx > 0 and position - self._starts[look_back_idx - 1] <= max_length:
            look_back_idx -= 1
            
        # Collect all features that contain the position
        results = []
        for i in range(look_back_idx, min(len(self._starts), idx + MAX_LOOK_FORWARD)):  # Limit forward search
            if self._starts[i] > position:
                break
            if self._starts[i] <= position < self._ends[i]:
                # Create a lightweight feature object for the result
                results.append(self._create_feature_from_index(i))
                
        return results
        
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """
        Get all features that overlap with the given range using binary search.
        
        Args:
            start: Start position of the interval (inclusive)
            end: End position of the interval (exclusive)
            
        Returns:
            List of features that overlap with the interval
        """
        if not self._is_loaded or self._starts is None or len(self._starts) == 0:
            return []
            
        # Find the first feature that starts at or after the start position
        idx = self._binary_search_position(start)
        
        # Look backward for features that might overlap
        bin_id = start // self.bin_size
        max_length = self._max_feature_length_by_bin.get(bin_id, 0)
        
        # Calculate how far back to look based on max feature length
        look_back_idx = idx
        while look_back_idx > 0 and start - self._starts[look_back_idx - 1] <= max_length:
            look_back_idx -= 1
            
        # Collect all features that overlap with the interval
        results = []
        for i in range(look_back_idx, len(self._starts)):
            if self._starts[i] >= end:
                break
            if self._starts[i] < end and self._ends[i] > start:
                results.append(self._create_feature_from_index(i))
                
        return results
        
    def _binary_search_position(self, position: int) -> int:
        """
        Find the index of the first feature that starts at or after the given position.
        
        Args:
            position: The position to search for
            
        Returns:
            Index of the first feature that starts at or after the position,
            or len(self._starts) if no such feature exists
        """
        left, right = 0, len(self._starts) - 1
        result = len(self._starts)  # Default if no suitable index is found
        
        while left <= right:
            mid = (left + right) // 2
            if self._starts[mid] >= position:
                result = mid
                right = mid - 1
            else:
                left = mid + 1
                
        return result
        
    def _create_feature_from_index(self, idx: int) -> GenomicFeature:
        """
        Create a lightweight feature object from stored data.
        
        Args:
            idx: Index in the arrays
            
        Returns:
            A GenomicFeature object representing the MSI site
        """
        start = int(self._starts[idx])
        end = int(self._ends[idx])
        repeat_unit_bases = self._repeat_unit_bases[idx]
        
        # Create a lightweight feature object (not a full MsiSiteRecord)
        return MsiSiteRecord(
            id=f"MSI_{start}",
            chrom=self.chromosome,
            location=start,
            repeat_unit_bases=repeat_unit_bases,
            repeat_unit_length=len(repeat_unit_bases),
            repeat_times=(end - start) // len(repeat_unit_bases),
        )
        
    def index_build_end(self) -> None:
        """
        Finalize the index build process.
        
        This sorts the arrays by start position for binary search.
        """
        super().index_build_end()
        
        # Sort arrays by start position for binary search
        if self._starts is not None and len(self._starts) > 0 and self._feature_count > 0:
            sort_indices = np.argsort(self._starts[:self._feature_count])
            self._starts = self._starts[sort_indices]
            self._ends = self._ends[sort_indices]
            
            # Reorder the DnaStringArray based on sort_indices
            new_repeat_unit_bases = DnaStringArray(initial_strings=len(self._repeat_unit_bases))
            for new_idx, old_idx in enumerate(sort_indices):
                new_repeat_unit_bases.add(self._repeat_unit_bases[old_idx])
            self._repeat_unit_bases = new_repeat_unit_bases
            
        self._is_loaded = True

    def __getitem__(self, index: int) -> GenomicFeature:
        """Get all features at a specific index."""
        return self._create_feature_from_index(index)
    
    def get_features(self) -> list[GenomicFeature]:
        """Get all features."""
        return [self._create_feature_from_index(i) for i in range(self._feature_count)]
    
    def __len__(self) -> int:
        return self._feature_count
        
    def trim(self) -> None:
        """
        Trim arrays to their actual used size to reduce memory usage.
        
        This is particularly useful before serialization to avoid storing
        large amounts of unused memory in pickle files.
        """
        # Call the parent class implementation first
        super().trim()
        
        # Trim NumPy arrays if they're larger than needed
        if self._starts is not None and len(self._starts) > self._feature_count:
            self._starts = self._starts[:self._feature_count].copy()
            
        if self._ends is not None and len(self._ends) > self._feature_count:
            self._ends = self._ends[:self._feature_count].copy()
            
        # Trim the DnaStringArray if it exists
        if hasattr(self, '_repeat_unit_bases') and self._repeat_unit_bases is not None:
            self._repeat_unit_bases.trim()
            
    def __getstate__(self):
        """
        Prepare the object for pickling.
        
        This method is called by pickle before serialization.
        It trims the arrays to reduce the serialized size.
        
        Returns:
            The object's state dictionary
        """
        # Trim arrays before pickling
        self.trim()
        
        # Return the object's state
        return self.__dict__
        
    def __str__(self) -> str:
        """Return a string representation of the MSI chromosome store."""
        status = "building" if self.index_build_mode else "built" if self.index_finished else "uninitialized"
        loaded_status = "loaded" if self._is_loaded else "not loaded"
        
        # Memory usage statistics
        mem_stats = ""
        if self._starts is not None and self._ends is not None:
            starts_mem = self._starts.nbytes / (1024 * 1024)  # MB
            ends_mem = self._ends.nbytes / (1024 * 1024)  # MB
            mem_stats = f", memory_usage={starts_mem + ends_mem:.2f}MB"
        
        # Bin statistics
        bin_stats = ""
        if self._max_feature_length_by_bin:
            bin_count = len(self._max_feature_length_by_bin)
            avg_max_length = sum(self._max_feature_length_by_bin.values()) / bin_count if bin_count > 0 else 0
            bin_stats = f", bins={bin_count}, avg_max_length={avg_max_length:.1f}"
        
        # Sample of features (create a few sample features)
        sample_str = ""
        if self._feature_count > 0 and self._is_loaded:
            sample_size = min(MAX_SAMPLES_TO_SHOW, self._feature_count)
            sample_indices = range(sample_size)
            sample_features = [self._create_feature_from_index(i) for i in sample_indices]
            sample_str = f", sample: [{', '.join(f.id for f in sample_features)}" + (", ...]" if self._feature_count > sample_size else "]")
        
        return (f"MsiChromosomeStore(chromosome='{self.chromosome}', features={self._feature_count}, "
                f"status={status}, data_status={loaded_status}{mem_stats}{bin_stats}{sample_str})")
    
    def __repr__(self) -> str:
        """Return a string representation of the MSI chromosome store."""
        return self.__str__()