"""Tests for genomic feature store implementations."""

import unittest
import time
import numpy as np

from pygnome.feature_store.brute_force_store import BruteForceFeatureStore
from pygnome.feature_store.chromosome_feature_store import ChromosomeFeatureStore
from pygnome.feature_store.genomic_feature_store_protocol import GenomicFeatureStoreProtocol
from pygnome.genomics import GenomicFeature, Strand
from pygnome.feature_store import (
    BinnedGenomicStore, GenomicFeatureStore, 
    IntervalTreeStore, StoreType
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
            self.features.append(GenomicFeature( id=f"small_{i}", chrom="chr1", start=start, end=end, strand=Strand.POSITIVE))
        
        # Add some medium features
        for i in range(20):
            start = i * 500
            end = start + 300
            self.features.append(GenomicFeature(id=f"medium_{i}", chrom="chr1", start=start, end=end, strand=Strand.POSITIVE))
        
        # Add some large features
        for i in range(5):
            start = i * 2000
            end = start + 1500
            self.features.append(GenomicFeature(id=f"large_{i}", chrom="chr1",start=start,end=end,strand=Strand.POSITIVE))

    def compare_results(self, store_type, expected, got, chr: str, start: int, end: int | None):
        end_str = "" if end is None else f"-{end}"
        # Compare number of features
        exp_str = "\n".join(sorted([f"\t\t{f}" for f in expected]))
        got_str = "\n".join(sorted([f"\t\t{f}" for f in got]))
        if len(expected) != len(got):
            self.assertEqual(len(expected), len(got), f"\nFeature store {store_type} failed to retrieve at {chr}:{start}{end_str}:\n" 
                             + f"\tExpexcted ({len(expected)}):\n{exp_str}\n" 
                             + f"\tGot       ({len(got)}):\n{got_str}\n"
                             )
        # Compare found features
        exp_ids = [f.id for f in expected]
        got_ids = [f.id for f in got]
        if exp_ids != got_ids:
            self.assertEqual(len(expected), len(got), f"\nFeature store {store_type} failed to retrieve at {chr}:{start}{end_str}:\n" 
                             + f"\tExpexcted IDs: : {exp_ids}\n" 
                             + f"\tGot IDs        : {got_ids}\n"
                             )        


    def compare_chr_stores(self, store_reference: ChromosomeFeatureStore, store: ChromosomeFeatureStore):
            """
            Compare two chromosome feature stores.
            One is a reference store (store_reference) and the other is the store being tested (store).
            """
            store_type = store.__class__.__name__
            # Add features
            with store_reference, store:
                for feature in self.features:
                    store.add(feature)
                    store_reference.add(feature)
            
            # Test position query
            chr, pos = 'chr1', 250
            expected = store_reference.get_by_position(pos)
            features_at_pos = store.get_by_position(pos)
            self.compare_results(store_type, expected, features_at_pos, "chr1", pos, None)
            
            # Test range query
            chr, start, end = 'chr1', 1000, 1200
            expected = store_reference.get_by_interval(start, end)
            features_in_range = store.get_by_interval(start, end)
            self.compare_results(store_type, expected, features_in_range, "chr1", start, end)


    def compare_stores(self, store_reference: GenomicFeatureStoreProtocol, store: GenomicFeatureStoreProtocol):
            """
            Compare two feature stores.
            One is a reference store (store_reference) and the other is the store being tested (store).
            """
            store_type = store.store_type
            # Add features
            with store_reference, store:
                for feature in self.features:
                    store.add(feature)
                    store_reference.add(feature)
            
            # Test position query
            chr, pos = 'chr1', 250
            expected = store_reference.get_by_position(chr, pos)
            features_at_pos = store.get_by_position(chr, pos)
            self.compare_results(store_type, expected, features_at_pos, "chr1", pos, None)

            # Test range query
            chr, start, end = 'chr1', 1000, 1200
            expected = store_reference.get_by_interval(chr, start, end)
            features_in_range = store.get_by_interval(chr, start, end)
            self.compare_results(store_type, expected, features_in_range, "chr1", start, end)


    def test_stores(self):
        """Test the interval tree store implementation."""
        for store in [IntervalTreeStore('chr1'), BinnedGenomicStore('chr1', bin_size=1000)]:
            store_naive = BruteForceFeatureStore('chr1')
            self.compare_chr_stores(store_naive, store)
    
    
    def test_feature_store_impl(self):
        """Test the main feature store implementation."""
        # Test with different store types
        for store_type in [StoreType.INTERVAL_TREE, StoreType.BINNED]:
            store = GenomicFeatureStore(store_type=store_type)
            store_naive = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)
            self.compare_stores(store_naive, store)

    
    def test_performance_comparison(self):
        """Compare performance of different store implementations."""        
        # Create stores
        interval_store = GenomicFeatureStore(store_type=StoreType.INTERVAL_TREE)
        binned_store = GenomicFeatureStore(store_type=StoreType.BINNED)
        
        # Create more test data for performance testing
        features = []
        np.random.seed(42)
        for i in range(10000):
            start = np.random.randint(0, 1000000)
            length = np.random.randint(1, 1000)
            end = start + length
            features.append(GenomicFeature( id=f"feature_{i}", chrom="chr1", start=start, end=end, strand=Strand.POSITIVE))
        
        # Measure insertion time
        stores = {
            "Interval Tree": interval_store,
            "Binned": binned_store,
        }
        
        for name, store in stores.items():
            start_time = time.time()
            for feature in features:
                store.add(feature)
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