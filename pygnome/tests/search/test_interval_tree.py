"""Benchmark tests for genomic feature interval_tree implementations."""

import unittest
import numpy as np

from pygnome.feature_store.interval_tree_store import IntervalNode, IntervalTree


def get_all_idx(node: IntervalNode) -> set[int]:
    """Get all indices in the interval tree."""
    idxs = set()
    _get_all_idx(node, idxs)
    return idxs


def _get_all_idx(node: IntervalNode, idxs: set[int]) -> None:
    """Recursively collect all indices in the node and its children."""
    if node is None:
        return
    idxs.update([idx for start, end, idx in node.intervals])
    _get_all_idx(node.left, idxs)
    _get_all_idx(node.right, idxs)


class TestIntervalTree(unittest.TestCase):

    def setUp(self):
        """Set up benchmark parameters."""        
        # Parameters
        self.feature_counts = [1, 10, 100, 1000, 10_000, 100_000]
        self.query_counts = 1000
        self.max_position = 1_000_000


    def generate_intervals(self, count):
        """Generate random genomic intervals."""
        intervals = []
        np.random.seed(42)
        
        for i in range(count):
            # Randomly select feature size
            size_class = np.random.choice(["small", "medium", "large"], p=[0.7, 0.2, 0.1])
            
            if size_class == "small":
                # Small intervals (1-100 bp)
                start = np.random.randint(0, self.max_position)
                length = np.random.randint(1, 100)
            elif size_class == "medium":
                # Medium intervals (100-1000 bp)
                start = np.random.randint(0, self.max_position)
                length = np.random.randint(100, 1000)
            else:
                # Large intervals (1000-10000 bp)
                start = np.random.randint(0, self.max_position)
                length = np.random.randint(1000, 10000)
            
            end = start + length
            triplet = (start, end, i)            
            intervals.append(triplet)
        
        return intervals


    def generate_position_queries(self, count):
        """Generate random position queries."""
        np.random.seed(43)  # Different seed from feature generation
        return [np.random.randint(0, self.max_position) for _ in range(count)]


    def generate_interval_queries(self, count):
        """Generate random range queries."""
        np.random.seed(44)  # Different seed
        queries = []
        for _ in range(count):
            start = np.random.randint(0, self.max_position)
            length = np.random.randint(1, 10000)
            end = min(start + length, self.max_position)
            queries.append((start, end))        
        return queries


    def test_insertion(self):
        """Benchmark feature insertion performance."""
        for count in self.feature_counts:
            intervals = self.generate_intervals(count)
            print(f"Testing insertion for {count} intervals...")

            interval_tree = IntervalTree()
            for interval in intervals:
                start, end, idx = interval
                interval_tree.add(start, end, idx)
            interval_tree.build()
        
            # Check that the interval_tree has all the intervals
            expected_idxs = {f[2] for f in intervals}
            idxs = get_all_idx(interval_tree.root)
            self.assertEqual(idxs, expected_idxs, f"IntervalTree stored {len(idxs)} / {len(expected_idxs)} intervals")


    def test_position_query(self):
        """Position query test."""
        for count in self.feature_counts:
            intervals = self.generate_intervals(count)
            queries = self.generate_position_queries(self.query_counts)
            print(f"Testing possition queries for {count} intervals...")

            # Test on different interval_tree types
            interval_tree = IntervalTree()
            for feature in intervals:
                start, end, idx = feature
                interval_tree.add(start, end, idx)
            interval_tree.build()
            
            # Check that the interval_tree can find the intervals
            for pos in queries:
                result_ids = interval_tree.at(pos)
                # Compare with brute-force search
                ref_ids = {idx for start, end, idx in intervals if start <= pos < end}
                self.assertEqual(ref_ids, result_ids, f"Total features found differs for {pos}")


    def test_interval_query(self):
        """Interval query test"""
        for count in self.feature_counts:
            intervals = self.generate_intervals(count)
            queries = self.generate_interval_queries(self.query_counts)
            print(f"Testing interval queries for {count} intervals...")

            # Test on different interval_tree types
            interval_tree = IntervalTree()
            for feature in intervals:
                start, end, idx = feature
                interval_tree.add(start, end, idx)
            interval_tree.build()
            
            # Check that the interval_tree can find the intervals
            for query in queries:
                start, end = query
                result_ids = interval_tree.overlap(start, end)
                # Compare with brute-force search
                ref_ids = {idx for i_start, i_end, idx in intervals if not (i_end <= start or i_start >= end)}
                self.assertEqual(ref_ids, result_ids, f"Total features found differs for {query}")


if __name__ == "__main__":
    unittest.main()