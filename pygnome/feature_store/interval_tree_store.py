"""Interval tree implementation for efficient genomic feature queries."""

import numpy as np
from dataclasses import dataclass, field

from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore, MAX_SAMPLES_TO_SHOW
from pygnome.genomics.genomic_feature import GenomicFeature


MAX_INTERVALS_IN_LEAF = 8

@dataclass
class IntervalNode:
    """
    A node in the interval tree.
    
    Each node contains a center point and intervals that overlap with that center point.
    """
    center: int
    # Intervals containing the center point (start, end, index), or 'leaf node
    intervals: list[tuple[int, int, int]] = field(default_factory=list)
    # Child nodes
    left: 'IntervalNode | None' = None
    right: 'IntervalNode | None' = None
    
    def __str__(self) -> str:
        """String representation of the node."""
        return f"IntervalNode(center={self.center}, intervals={len(self.intervals)})"

    def is_center(self, start:int, end: int) -> bool:
        """Check if the interval [start, end) intersects the node's center."""
        return start <= self.center < end

    def is_left(self, start:int, end: int) -> bool:
        """Check if the interval [start, end) is completely to the left."""
        return end <= self.center

    def is_right(self, start:int, end: int) -> bool:
        """Check if the interval [start, end) is completely to the right."""
        return self.center < start


class IntervalTree:
    """
    An interval tree data structure for efficient interval queries.
    
    This implementation uses a centered interval tree approach, which provides
    O(log n + k) query time for finding k intervals that overlap with a given point or range,
    and O(n log n) construction time.
    """
    
    def __init__(self):
        """Initialize an empty interval tree."""
        self.root: IntervalNode | None = None
        self.intervals: list[tuple[int, int, int]] = []  # (start, end, index)
    
    def add(self, start: int, end: int, index: int) -> None:
        """Add an interval to the tree."""
        self.intervals.append((start, end, index))
    
    def build(self) -> None:
        """Build the interval tree from the added intervals."""
        if not self.intervals:
            return
        self.root = self._build_tree(self.intervals)
    
    def _build_tree(self, intervals: list[tuple[int, int, int]]) -> IntervalNode | None:
        """Recursively build the interval tree."""
        if not intervals:
            return None
        
        # Base case: if only one interval, create a leaf node
        if len(intervals) <= MAX_INTERVALS_IN_LEAF:
            start, end, idx = intervals[0]
            center = (start + end) // 2
            node = IntervalNode(center=center)
            node.intervals.extend(intervals)
            return node
        
        # Find the center point (median of all endpoints)
        all_points = []
        for start, end, _ in intervals:
            all_points.append(start)
            all_points.append(end)
        
        all_points.sort()
        if len(all_points) % 2 == 0:
            # If even number of points, take the average median
            center = (all_points[len(all_points) // 2 - 1] + all_points[len(all_points) // 2]) // 2
        else:
            # If odd number of points, take the median
            center = all_points[len(all_points) // 2]
        
        # Create the node for this center
        node = IntervalNode(center=center)
        
        # Divide intervals into left, right, and center
        left_intervals = []
        right_intervals = []
        
        for interval in intervals:
            start, end, idx = interval
            if node.is_center(start, end):
                # If interval contains the center point, add to intervals
                node.intervals.append(interval)
            elif node.is_left(start, end):
                # If interval is completely to the left of center
                left_intervals.append(interval)        
            elif node.is_right(start, end):
                # If interval is completely to the right of center
                right_intervals.append(interval)
            else:
                # This case should not happen due to the way we split intervals
                raise ValueError(f"Interval {interval} does not fit in left/right categories.")
        
        # Recursively build left and right subtrees
        # Only recurse if we're making progress (intervals are being divided)
        if left_intervals:
            node.left = self._build_tree(left_intervals)
        
        if right_intervals:
            node.right = self._build_tree(right_intervals)
        
        return node
    
    def _check_overlapping_intervals(self, node: IntervalNode, start: int, end: int, result: set[int]) -> None:
        """Check if any intervals in the node overlap with the given range."""
        for interval_start, interval_end, idx in node.intervals:
            # Special case for zero-length features
            if interval_start == interval_end and start <= interval_start < end:
                result.add(idx)
            # Normal case
            elif interval_start < end and interval_end > start:
                result.add(idx)
    
    def overlap(self, start: int, end: int) -> set[int]:
        """Find all intervals that overlap with the given range."""
        if not self.root:
            return set()
        result = set()
        self._query_overlap(self.root, start, end, result)
        return result
    
    def _query_overlap(self, node: IntervalNode | None, start: int, end: int, result: set[int]) -> None:
        """Recursively query for overlapping intervals."""
        if not node:
            return
        # Check intervals at this node
        self._check_overlapping_intervals(node, start, end, result)

        if node.is_left(start, end):
            # Query interval is completely to the left of the center
            self._query_overlap(node.left, start, end, result)
        elif node.is_right(start, end):
            # Query interval is completely to the right of the center
            self._query_overlap(node.right, start, end, result)
        else:
            # Query range contains the center point: Check both subtrees
            self._query_overlap(node.left, start, end, result)
            self._query_overlap(node.right, start, end, result)
    
    def at(self, position: int) -> set[int]:
        """Find all intervals that contain the given position."""
        if not self.root:
            return set()
        result = set()
        self._query_at_position(self.root, position, result)
        return result
    
    def _check_node_intervals(self, node: IntervalNode, position: int, result: set[int]) -> None:
        """Check if any intervals in the node contain the position."""
        for interval_start, interval_end, idx in node.intervals:
            # For position queries, we use half-open intervals [start, end)
            # Special case for zero-length features
            if interval_start == interval_end and position == interval_start:
                result.add(idx)
            # Normal case
            elif interval_start <= position < interval_end:
                result.add(idx)
    
    def _query_at_position(self, node: IntervalNode | None, position: int, result: set[int]) -> None:
        """Recursively query for intervals containing a specific position."""
        if not node:
            return
        
        # Check intervals at this node
        self._check_node_intervals(node, position, result)
        
        if position < node.center:
            # If the position is to the left of the center: Only need to check the left subtree
            self._query_at_position(node.left, position, result)    
        elif position > node.center:
            # Check if the position is to the right of the center: Only need to check the right subtree
            self._query_at_position(node.right, position, result)    
        else:
            # Position is exactly at the center: Check both subtrees (position could be at the boundary)
            self._query_at_position(node.left, position, result)
            self._query_at_position(node.right, position, result)


class IntervalTreeStore(ChromosomeFeatureStore):
    """Store genomic features using an efficient interval tree."""
    
    def __init__(self, chromosome: str):
        super().__init__(chromosome=chromosome)
        self.interval_tree = IntervalTree()
        self.tree_built = False
    
    def add(self, feature: GenomicFeature) -> None:
        """Add a feature to the interval tree."""
        super().add(feature)
        # Store the index in the features list
        feature_idx = len(self.features) - 1
        self.interval_tree.add(feature.start, feature.end, feature_idx)
        # Mark tree as needing rebuild
        self.tree_built = False
    
    def index_build_end(self) -> None:
        """Ensure the interval tree is built before querying."""
        super().index_build_end()
        if not self.tree_built:
            self.interval_tree.build()
            self.tree_built = True
    
    def get_by_position(self, position: int) -> list[GenomicFeature]:
        """
        Get all features at a specific position.
        
        Uses half-open intervals [start, end) where start is included but end is excluded.
        """
        indices = self.interval_tree.at(position)
        return [self.features[idx] for idx in indices]
    
    def get_by_interval(self, start: int, end: int) -> list[GenomicFeature]:
        """
        Get all features that overlap with the given range.
        
        Uses half-open intervals [start, end) where start is included but end is excluded.
        """
        # Handle invalid ranges
        if end <= start:
            return []
        # Get indices from interval tree
        indices = self.interval_tree.overlap(start, end)
        return [self.features[idx] for idx in indices]
        
    def __str__(self) -> str:
        """Return a string representation of the interval tree store."""
        status = "building" if self.index_build_mode else "built" if self.index_finished else "uninitialized"
        tree_status = "built" if self.tree_built else "not built"
        
        # Get tree statistics if available
        tree_stats = ""
        if self.tree_built and self.interval_tree.root:
            # Count nodes in the tree (simple BFS)
            node_count = 0
            queue = [self.interval_tree.root]
            while queue:
                node = queue.pop(0)
                node_count += 1
                if node.left:
                    queue.append(node.left)
                if node.right:
                    queue.append(node.right)
            
            # Count intervals at the root level
            root_intervals = len(self.interval_tree.root.intervals) if self.interval_tree.root else 0
            tree_stats = f", tree_nodes={node_count}, root_intervals={root_intervals}"
        
        # Sample of features
        sample_size = min(MAX_SAMPLES_TO_SHOW, len(self.features))
        sample = self.features[:sample_size] if self.features else []
        sample_str = ""
        if sample:
            sample_str = f", sample: [{', '.join(str(f.id) for f in sample)}" + (", ...]" if len(self.features) > sample_size else "]")
        
        return (f"IntervalTreeStore(chromosome='{self.chromosome}', features={len(self.features)}, "
                f"status={status}, tree_status={tree_status}{tree_stats}{sample_str})")
    
    def __repr__(self) -> str:
        """Return a string representation of the interval tree store."""
        return self.__str__()
