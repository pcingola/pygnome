"""Tests for genomic feature store implementations."""

import unittest
import time
import numpy as np

from pygnome.genomics import GenomicFeature, Strand
from pygnome.genomics.feature_store import (
    BinnedGenomicStore, GenomicFeatureStoreImpl, 
    IntervalTreeStore, PositionHashStore, StoreType
)


class TestFeatureStore(unittest.TestCase):
    """Test cases for genomic feature store implementations."""
    
    def setUp(self):
        """Set up test data."""
        # Create some test features
        self.features = []
        
        # Add some small features
        for i in range(100):
            start = i * 100
            end = start + 50
            self.features.append(
                GenomicFeature(
                    id=f"small_{i}",
                    chrom="chr1",
                    start=start,
                    end=end,
                    strand=Strand.POSITIVE
                )
            )
        
        # Add some medium features
        for i in range(20):
            start = i * 500
            end = start + 300
            self.features.append(
                GenomicFeature(
                    id=f"medium_{i}",
                    chrom="chr1",
                    start=start,
                    end=end,
                    strand=Strand.POSITIVE
                )
            )
        
        # Add some large features
        for i in range(5):
            start = i * 2000
            end = start + 1500
            self.features.append(
                GenomicFeature(
                    id=f"large_{i}",
                    chrom="chr1",
                    start=start,
                    end=end,
                    strand=Strand.POSITIVE
                )
            )
    
    def test_interval_tree_store(self):
        """Test the interval tree store implementation."""
        store = IntervalTreeStore()
        
        # Add features
        for feature in self.features:
            store.add_feature(feature)
        
        # Test position query
        features_at_pos = store.get_by_position(250)
        self.assertEqual(len(features_at_pos), 3)  # small_2, medium_0, large_0
        
        # Test range query
        features_in_range = store.get_by_interval(1000, 1200)
        self.assertEqual(len(features_in_range), 5)  # small_10, small_11, small_12, medium_2, large_0
    
    def test_binned_store(self):
        """Test the binned store implementation."""
        store = BinnedGenomicStore(bin_size=1000)
        
        # Add features
        for feature in self.features:
            store.add_feature(feature)
        
        # Test position query
        features_at_pos = store.get_by_position(250)
        self.assertEqual(len(features_at_pos), 3)  # small_2, medium_0, large_0
        
        # Test range query
        features_in_range = store.get_by_interval(1000, 1200)
        self.assertEqual(len(features_in_range), 5)  # small_10, small_11, small_12, medium_2, large_0
    
    def test_position_hash_store(self):
        """Test the position hash store implementation."""
        store = PositionHashStore()
        
        # Add features
        for feature in self.features:
            store.add_feature(feature)
        
        # Test position query
        features_at_pos = store.get_by_position(250)
        self.assertEqual(len(features_at_pos), 3)  # small_2, medium_0, large_0
        
        # Test range query
        features_in_range = store.get_by_interval(1000, 1200)
        self.assertEqual(len(features_in_range), 5)  # small_10, small_11, small_12, medium_2, large_0
    
    def test_feature_store_impl(self):
        """Test the main feature store implementation."""
        # Test with different store types
        for store_type in [StoreType.INTERVAL_TREE, StoreType.BINNED, StoreType.POSITION_HASH]:
            store = GenomicFeatureStoreImpl(store_type=store_type)
            
            # Add features
            for feature in self.features:
                store.add_feature(feature)
            
            # Test position query
            features_at_pos = store.get_by_position("chr1", 250)
            self.assertEqual(len(features_at_pos), 3)  # small_2, medium_0, large_0
            
            # Test range query
            features_in_range = store.get_by_interval("chr1", 1000, 1200)
            self.assertEqual(len(features_in_range), 5)  # small_10, small_11, small_12, medium_2, large_0
    
    def test_performance_comparison(self):
        """Compare performance of different store implementations."""
        # Skip in normal test runs to avoid slowing down the test suite
        if not hasattr(self, 'run_performance_tests'):
            self.skipTest("Performance tests disabled")
        
        # Create stores
        interval_store = GenomicFeatureStoreImpl(store_type=StoreType.INTERVAL_TREE)
        binned_store = GenomicFeatureStoreImpl(store_type=StoreType.BINNED)
        position_store = GenomicFeatureStoreImpl(store_type=StoreType.POSITION_HASH)
        
        # Create more test data for performance testing
        features = []
        np.random.seed(42)
        for i in range(10000):
            start = np.random.randint(0, 1000000)
            length = np.random.randint(1, 1000)
            end = start + length
            features.append(
                GenomicFeature(
                    id=f"feature_{i}",
                    chrom="chr1",
                    start=start,
                    end=end,
                    strand=Strand.POSITIVE
                )
            )
        
        # Measure insertion time
        stores = {
            "Interval Tree": interval_store,
            "Binned": binned_store,
            "Position Hash": position_store
        }
        
        for name, store in stores.items():
            start_time = time.time()
            for feature in features:
                store.add_feature(feature)
            end_time = time.time()
            print(f"{name} insertion time: {end_time - start_time:.4f} seconds")
        
        # Measure position query time (100 random positions)
        positions = np.random.randint(0, 1000000, size=100)
        
        for name, store in stores.items():
            start_time = time.time()
            for pos in positions:
                store.get_by_position("chr1", pos)
            end_time = time.time()
            print(f"{name} position query time: {end_time - start_time:.4f} seconds")
        
        # Measure range query time (100 random ranges)
        ranges = []
        for _ in range(100):
            start = np.random.randint(0, 1000000)
            end = start + np.random.randint(1, 10000)
            ranges.append((start, end))
        
        for name, store in stores.items():
            start_time = time.time()
            for start, end in ranges:
                store.get_by_interval("chr1", start, end)
            end_time = time.time()
            print(f"{name} range query time: {end_time - start_time:.4f} seconds")


if __name__ == "__main__":
    # To run performance tests:
    # test = TestFeatureStore()
    # test.run_performance_tests = True
    # test.setUp()
    # test.test_performance_comparison()
    unittest.main()