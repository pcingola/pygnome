"""Memory-efficient interval storage for genomic feature queries."""

import numpy as np
from typing import List, Optional

from ..genomic_feature import GenomicFeature
from .base import ChromosomeFeatureStore


class CompactIntervalArray:
    """
    Memory-efficient interval storage using numpy arrays.
    
    This implementation uses numpy arrays for storing intervals and provides
    efficient query operations for overlapping intervals.
    """
    
    def __init__(self, initial_capacity: int = 1000):
        """Initialize with a given capacity."""
        self.capacity = initial_capacity
        self.size = 0
        
        # Use numpy arrays with efficient dtypes
        self.starts = np.zeros(initial_capacity, dtype=np.int32)
        self.ends = np.zeros(initial_capacity, dtype=np.int32)
        self.indices = np.zeros(initial_capacity, dtype=np.int32)
        
        # Track if sorted
        self.is_sorted = True
    
    def _ensure_capacity(self, needed_size: int) -> None:
        """Ensure the arrays have enough capacity."""
        if needed_size > self.capacity:
            new_capacity = max(self.capacity * 2, needed_size)
            
            # Create new arrays with increased capacity
            new_starts = np.zeros(new_capacity, dtype=np.int32)
            new_ends = np.zeros(new_capacity, dtype=np.int32)
            new_indices = np.zeros(new_capacity, dtype=np.int32)
            
            # Copy existing data
            new_starts[:self.size] = self.starts[:self.size]
            new_ends[:self.size] = self.ends[:self.size]
            new_indices[:self.size] = self.indices[:self.size]
            
            # Replace arrays
            self.starts = new_starts
            self.ends = new_ends
            self.indices = new_indices
            self.capacity = new_capacity
    
    def add(self, start: int, end: int, index: int) -> None:
        """Add an interval to the array."""
        self._ensure_capacity(self.size + 1)
        
        self.starts[self.size] = start
        self.ends[self.size] = end
        self.indices[self.size] = index
        self.size += 1
        
        # Mark as unsorted if we add in non-sorted order
        if self.size > 1 and self.is_sorted and start < self.starts[self.size - 2]:
            self.is_sorted = False
    
    def _ensure_sorted(self) -> None:
        """Ensure the intervals are sorted by start position."""
        if not self.is_sorted and self.size > 0:
            # Get the actual data slice
            data_slice = slice(0, self.size)
            
            # Get sort order
            sort_idx = np.argsort(self.starts[data_slice])
            
            # Sort all arrays using this order
            self.starts[data_slice] = self.starts[data_slice][sort_idx]
            self.ends[data_slice] = self.ends[data_slice][sort_idx]
            self.indices[data_slice] = self.indices[data_slice][sort_idx]
            
            self.is_sorted = True
    
    def overlap(self, start: int, end: int) -> np.ndarray:
        """Find all intervals that overlap with the given range."""
        if self.size == 0:
            return np.array([], dtype=np.int32)
        
        self._ensure_sorted()
        
        # Get the actual data slice
        data_slice = slice(0, self.size)
        
        # Find all intervals where end >= start and start <= end
        mask = (self.ends[data_slice] >= start) & (self.starts[data_slice] <= end)
        return self.indices[data_slice][mask]
    
    def at(self, position: int) -> np.ndarray:
        """Find all intervals that contain the given position."""
        return self.overlap(position, position)


class IntervalTreeStore(ChromosomeFeatureStore):
    """Store genomic features using a memory-efficient interval array."""
    
    def __init__(self):
        super().__init__()
        self.interval_array = CompactIntervalArray()
    
    def add_feature(self, feature: GenomicFeature) -> None:
        """Add a feature to the interval array."""
        super().add_feature(feature)
        # Store the index in the features list
        self.interval_array.add(feature.start, feature.end, len(self.features) - 1)
    
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        indices = self.interval_array.at(position)
        return [self.features[idx] for idx in indices]
    
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        indices = self.interval_array.overlap(start, end)
        return [self.features[idx] for idx in indices]