"""Interval tree implementation for efficient genomic feature queries."""

import numpy as np
from dataclasses import dataclass, field

from ..genomics.genomic_feature import GenomicFeature
from .genomic_feature_store_protocol import ChromosomeFeatureStore


@dataclass
class IntervalNode:
    """
    A node in the interval tree.
    
    Each node contains a center point and intervals that overlap with that center point.
    """
    center: int
    # Intervals containing the center point (start, end, index)
    intervals: list[tuple[int, int, int]] = field(default_factory=list)
    # Child nodes
    left: 'IntervalNode | None' = None
    right: 'IntervalNode | None' = None
    
    def __str__(self) -> str:
        """String representation of the node."""
        return f"IntervalNode(center={self.center}, intervals={len(self.intervals)})"


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
        if len(intervals) == 1:
            start, end, idx = intervals[0]
            center = (start + end) // 2
            node = IntervalNode(center=center)
            node.intervals.append(intervals[0])
            return node
        
        # Find the center point (median of all endpoints)
        all_points = []
        for start, end, _ in intervals:
            all_points.append(start)
            all_points.append(end)
        
        all_points.sort()
        center = all_points[len(all_points) // 2]
        
        # Create the node for this center
        node = IntervalNode(center=center)
        
        # Divide intervals into left, right, and center
        left_intervals = []
        right_intervals = []
        
        for interval in intervals:
            start, end, idx = interval
            
            # If interval contains the center point, add to intervals
            if start <= center < end:
                node.intervals.append(interval)
            
            # If interval is completely to the left of center
            elif end <= center:
                left_intervals.append(interval)
            
            # If interval is completely to the right of center
            elif start >= center:
                right_intervals.append(interval)
        
        # No need to sort intervals - we'll check them all during queries
        
        # Recursively build left and right subtrees
        # Only recurse if we're making progress (intervals are being divided)
        if left_intervals and len(left_intervals) < len(intervals):
            node.left = self._build_tree(left_intervals)
        
        if right_intervals and len(right_intervals) < len(intervals):
            node.right = self._build_tree(right_intervals)
        
        return node
    
    def overlap(self, start: int, end: int) -> set[int]:
        """Find all intervals that overlap with the given range."""
        if not self.root:
            return set()
        
        result = set()
        self._query_overlap(self.root, start, end, result)
        return result
    
    def _check_overlapping_intervals(self, node: IntervalNode, start: int, end: int, result: set[int]) -> None:
        """Check if any intervals in the node overlap with the given range."""
        for interval_start, interval_end, idx in node.intervals:
            if interval_start < end and interval_end >= start:
                result.add(idx)
    
    def _query_overlap(self, node: IntervalNode | None, start: int, end: int, result: set[int]) -> None:
        """Recursively query for overlapping intervals."""
        if not node:
            return
        
        # Check intervals at this node
        self._check_overlapping_intervals(node, start, end, result)
        
        # Check if the query range is completely to the left of the center
        if end <= node.center:
            # Only need to check the left subtree
            self._query_overlap(node.left, start, end, result)
        
        # Check if the query range is completely to the right of the center
        elif start >= node.center:
            # Only need to check the right subtree
            self._query_overlap(node.right, start, end, result)
        
        # Query range contains the center point
        else:
            # Check both subtrees
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
            # Special case: include intervals where position equals end-1
            if interval_start <= position < interval_end or position == interval_end - 1:
                result.add(idx)
    
    def _query_at_position(self, node: IntervalNode | None, position: int, result: set[int]) -> None:
        """Recursively query for intervals containing a specific position."""
        if not node:
            return
        
        # Check intervals at this node
        self._check_node_intervals(node, position, result)
        
        # Check if the position is to the left of the center
        if position < node.center:
            # Only need to check the left subtree
            self._query_at_position(node.left, position, result)
        
        # Check if the position is to the right of the center
        elif position > node.center:
            # Only need to check the right subtree
            self._query_at_position(node.right, position, result)
        
        # Position is exactly at the center
        else:
            # Check both subtrees (position could be at the boundary)
            self._query_at_position(node.left, position, result)
            self._query_at_position(node.right, position, result)


class IntervalTreeStore(ChromosomeFeatureStore):
    """Store genomic features using an efficient interval tree."""
    
    def __init__(self):
        super().__init__()
        self.interval_tree = IntervalTree()
        self.tree_built = False
    
    def add(self, feature: GenomicFeature) -> None:
        """Add a feature to the interval tree."""
        super().add(feature)
        # Store the index in the features list
        self.interval_tree.add(feature.start, feature.end, len(self.features) - 1)
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
        
        indices = self.interval_tree.overlap(start, end)
        return [self.features[idx] for idx in indices]
