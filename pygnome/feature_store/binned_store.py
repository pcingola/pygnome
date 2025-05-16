import numpy as np

from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore, MAX_SAMPLES_TO_SHOW
from pygnome.genomics.genomic_feature import GenomicFeature

DEFAULT_BIN_SIZE = 100_000  # Default bin size in base pairs

class BinnedGenomicStore(ChromosomeFeatureStore):
    """Store genomic features using a memory-efficient binning approach."""

    def __init__(self, chromosome: str, bin_size: int = DEFAULT_BIN_SIZE):
        """
        Initialize a binned genomic store.
        
        Args:
            chromosome: Name of the chromosome
            bin_size: Size of each bin in base pairs
        """
        super().__init__(chromosome=chromosome)
        self.bin_size = bin_size
    
        # Use dictionary of numpy arrays for memory efficiency
        # Each key is a bin_id, and the value is a numpy array of feature indices
        self.bins: dict[int, np.ndarray] = {}
        
        # Track features per bin for efficient array sizing
        self._bin_counts: dict[int, int] = {}
        
        # Buffer for collecting indices before creating arrays
        self._bin_buffers: dict[int, list[int]] = {}
        
        # Threshold for converting buffer to array
        self._buffer_threshold = 100
    
    def _get_bin_ids(self, start: int, end: int) -> set[int]:
        """Get all bin IDs that this range spans."""
        start_bin = start // self.bin_size
        # Special case for zero-length features
        if start == end:
            return {start_bin}
        end_bin = (end - 1) // self.bin_size
        return set(range(start_bin, end_bin + 1))
    
    def _add_to_bin(self, bin_id: int, feature_idx: int) -> None:
        """Add a feature index to a bin, using buffer for small counts."""
        # Initialize buffer if needed
        if bin_id not in self._bin_buffers:
            self._bin_buffers[bin_id] = []
            self._bin_counts[bin_id] = 0
        
        # Add to buffer
        self._bin_buffers[bin_id].append(feature_idx)
        self._bin_counts[bin_id] += 1
        
        # Convert buffer to array if it reaches threshold
        if len(self._bin_buffers[bin_id]) >= self._buffer_threshold:
            self._convert_buffer_to_array(bin_id)
    
    def _convert_buffer_to_array(self, bin_id: int) -> None:
        """Convert a bin buffer to a numpy array for memory efficiency."""
        if bin_id in self._bin_buffers and self._bin_buffers[bin_id]:
            # Create array from buffer
            if bin_id in self.bins:
                # Append to existing array
                old_array = self.bins[bin_id]
                new_array = np.concatenate([
                    old_array,
                    np.array(self._bin_buffers[bin_id], dtype=np.int32)
                ])
                self.bins[bin_id] = new_array
            else:
                # Create new array
                self.bins[bin_id] = np.array(self._bin_buffers[bin_id], dtype=np.int32)
            
            # Clear buffer
            self._bin_buffers[bin_id] = []
    
    def add(self, feature: GenomicFeature) -> None:
        """Add a feature to the binned store."""
        super().add(feature)
        feature_idx = len(self.features) - 1
        
        # Add to all bins this feature spans
        for bin_id in self._get_bin_ids(feature.start, feature.end):
            self._add_to_bin(bin_id, feature_idx)
    
    def _get_bin_indices(self, bin_id: int) -> np.ndarray | None:
        """Get all feature indices for a bin."""
        return self.bins.get(bin_id)
    
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """Get all features at a specific position."""
        bin_id = position // self.bin_size
        
        # Get indices from this bin
        indices = self._get_bin_indices(bin_id)
        if indices is None:
            return []
        
        # Filter features that contain the position
        result = []
        for idx in indices:
            feature: GenomicFeature = self.features[idx]
            if feature.intersects_point(position):
                result.append(feature)
        return result
    
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """Get all features that overlap with the given range."""
        bin_ids = self._get_bin_ids(start, end)
        
        # Get unique feature indices from all relevant bins
        feature_indices = set()
        for bin_id in bin_ids:
            indices = self._get_bin_indices(bin_id)
            if indices is not None:
                feature_indices.update(indices)
        
        # Filter features that actually overlap
        result = []
        for idx in feature_indices:
            feature = self.features[idx]
            if feature.intersects_interval(start, end):
                result.append(feature)
        return result
    
    def index_build_end(self) -> None:
        """Finalize the index build process."""
        super().index_build_end()
        # Convert all remaining buffers to arrays
        for bin_id in list(self._bin_buffers.keys()):
            if self._bin_buffers[bin_id]:
                self._convert_buffer_to_array(bin_id)
                
    def __str__(self) -> str:
        """Return a string representation of the binned genomic store."""
        # Get basic info from parent class
        status = "building" if self.index_build_mode else "built" if self.index_finished else "uninitialized"
        
        # Calculate bin statistics
        bin_count = len(self.bins)
        avg_features_per_bin = len(self.features) / max(1, bin_count) if bin_count > 0 else 0
        
        # Find most populated bins
        bin_sizes = {bin_id: len(indices) for bin_id, indices in self.bins.items()}
        top_bins = sorted(bin_sizes.items(), key=lambda x: x[1], reverse=True)[:3] if bin_sizes else []
        
        # Sample of features
        sample_size = min(MAX_SAMPLES_TO_SHOW, len(self.features))
        sample = self.features[:sample_size] if self.features else []
        sample_str = ""
        if sample:
            sample_str = f", sample: [{', '.join(str(f.id) for f in sample)}" + (", ...]" if len(self.features) > sample_size else "]")
        
        # Bin statistics
        bin_stats = ""
        if top_bins:
            bin_stats = f", top bins: {', '.join(f'bin_{bin_id}({count})' for bin_id, count in top_bins)}"
        
        return (f"BinnedGenomicStore(chromosome='{self.chromosome}', features={len(self.features)}, "
                f"status={status}, bin_size={self.bin_size}, bins={bin_count}, "
                f"avg_features_per_bin={avg_features_per_bin:.1f}{bin_stats}{sample_str})")
    
    def __repr__(self) -> str:
        """Return a string representation of the binned genomic store."""
        return self.__str__()
