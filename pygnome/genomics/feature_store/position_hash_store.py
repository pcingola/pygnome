"""Position hash implementation for genomic feature storage."""

import numpy as np

from ..genomic_feature import GenomicFeature
from .base import ChromosomeFeatureStore


class PositionHashStore(ChromosomeFeatureStore):
    """
    Store genomic features using a position-based hash table.
    
    This implementation is optimized for exact position queries but can be
    memory-intensive for large features. It's recommended to use this only
    for small features or specific use cases where position queries are dominant.
    """
    
    def __init__(self):
        """Initialize a position hash store."""
        super().__init__()
        # Map from position to array of feature indices
        self.position_map: dict[int, np.ndarray] = {}
        
        # Buffer for collecting indices before creating arrays
        self._position_buffers: dict[int, list[int]] = {}
        
        # Threshold for converting buffer to array
        self._buffer_threshold = 50
    
    def _add_to_position(self, position: int, feature_idx: int) -> None:
        """Add a feature index to a position, using buffer for small counts."""
        # Initialize buffer if needed
        if position not in self._position_buffers:
            self._position_buffers[position] = []
        
        # Add to buffer
        self._position_buffers[position].append(feature_idx)
        
        # Convert buffer to array if it reaches threshold
        if len(self._position_buffers[position]) >= self._buffer_threshold:
            self._convert_buffer_to_array(position)
    
    def _convert_buffer_to_array(self, position: int) -> None:
        """Convert a position buffer to a numpy array for memory efficiency."""
        if position in self._position_buffers and self._position_buffers[position]:
            # Create array from buffer
            if position in self.position_map:
                # Append to existing array
                old_array = self.position_map[position]
                new_array = np.concatenate([
                    old_array,
                    np.array(self._position_buffers[position], dtype=np.int32)
                ])
                self.position_map[position] = new_array
            else:
                # Create new array
                self.position_map[position] = np.array(self._position_buffers[position], dtype=np.int32)
            
            # Clear buffer
            self._position_buffers[position] = []
    
    def _ensure_all_arrays(self) -> None:
        """Convert all buffers to arrays."""
        for position in list(self._position_buffers.keys()):
            if self._position_buffers[position]:
                self._convert_buffer_to_array(position)
    
    def add_feature(self, feature: GenomicFeature) -> None:
        """
        Add a feature to the position hash store.
        
        Note: This can be memory-intensive for large features as it stores
        an entry for every position the feature spans.
        """
        super().add_feature(feature)
        feature_idx = len(self.features) - 1
        
        # For very large features, we might want to skip or use a different approach
        feature_length = feature.end - feature.start + 1
        if feature_length > 1_000_000:  # Arbitrary threshold
            # For extremely large features, just store at start and end positions
            self._add_to_position(feature.start, feature_idx)
            self._add_to_position(feature.end, feature_idx)
            return
        
        # Add to all positions this feature spans
        for pos in range(feature.start, feature.end + 1):
            self._add_to_position(pos, feature_idx)
    
    def _get_position_indices(self, position: int) -> np.ndarray:
        """Get all feature indices for a position."""
        # Convert buffer if needed
        if position in self._position_buffers and self._position_buffers[position]:
            self._convert_buffer_to_array(position)
        
        # Return array or empty array
        if position in self.position_map:
            return self.position_map[position]
        return np.array([], dtype=np.int32)
    
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        indices = self._get_position_indices(position)
        
        # For large features that were only indexed at start/end
        result = [self.features[idx] for idx in indices]
        
        # Check if we need to scan for large features
        for feature in self.features:
            if feature.start <= position <= feature.end and feature.length > 1_000_000:
                if feature not in result:
                    result.append(feature)
        
        return result
    
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """
        Get all features that overlap with the given range.
        
        Note: This is inefficient for large ranges, but included for completeness.
        Consider using a different store type for range queries.
        """
        # This is inefficient for large ranges
        if end - start > 10000:  # Arbitrary threshold
            # For large ranges, scan all features instead
            return [f for f in self.features if f.start <= end and f.end >= start]
        
        # Ensure all buffers are converted to arrays
        self._ensure_all_arrays()
        
        # Get unique feature indices from all positions in the range
        feature_indices = set()
        for pos in range(start, end + 1):
            indices = self._get_position_indices(pos)
            feature_indices.update(indices)
        
        # Return features
        return [self.features[idx] for idx in feature_indices]