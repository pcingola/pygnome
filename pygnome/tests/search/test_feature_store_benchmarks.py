"""Benchmark tests for genomic feature store implementations."""

import unittest
import time
import numpy as np

from pygnome.feature_store.binned_store import BinnedGenomicStore
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.feature_store.interval_tree_store import IntervalTreeStore
from pygnome.genomics import GenomicFeature, Strand


class TestFeatureStoreBenchmarks(unittest.TestCase):
    """Benchmark tests for genomic feature store implementations."""
    
    def setUp(self):
        """Set up benchmark parameters."""        
        # Parameters
        self.feature_counts = [1000, 10_000, 100_000]  # Reduced for faster tests
        self.query_counts = 100
        self.chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5"]
        self.max_position = 1_000_000
        
        # Store types
        self.store_types = {
            "Interval Tree": StoreType.INTERVAL_TREE,
            "Binned": StoreType.BINNED,
        }
        
        # Performance thresholds (in seconds)
        # These are generous thresholds that can be adjusted based on hardware
        self.max_insertion_time = {
            1000: 0.01,    # 1 second for 1000 features
            10000: 0.1,  # 10 seconds for 10000 features
            100000: 1.0,  # 10 seconds for 10000 features
        }
        self.max_query_time = {
            1000: 0.01,    # 1 second for querying 1000 features
            10000: 0.1,   # 5 seconds for querying 10000 features
            100000: 1.0,   # 5 seconds for querying 10000 features
        }
    
    def generate_features(self, count):
        """Generate random genomic features."""
        features = []
        np.random.seed(42)
        
        for i in range(count):
            # Randomly select feature size
            size_class = np.random.choice(["small", "medium", "large"], p=[0.7, 0.2, 0.1])
            
            if size_class == "small":
                # Small features (1-100 bp)
                start = np.random.randint(0, self.max_position)
                length = np.random.randint(1, 100)
            elif size_class == "medium":
                # Medium features (100-1000 bp)
                start = np.random.randint(0, self.max_position)
                length = np.random.randint(100, 1000)
            else:
                # Large features (1000-10000 bp)
                start = np.random.randint(0, self.max_position)
                length = np.random.randint(1000, 10000)
            
            end = start + length
            
            # Randomly assign to chromosomes
            chrom = np.random.choice(self.chromosomes)
            
            # Create the feature
            feature = GenomicFeature(id=f"feature_{i}", chrom=chrom, start=start, end=end, strand=Strand.POSITIVE if np.random.random() > 0.5 else Strand.NEGATIVE)            
            features.append(feature)
        
        return features
    
    def generate_position_queries(self, count):
        """Generate random position queries."""
        np.random.seed(43)  # Different seed from feature generation
        queries = []
        
        for _ in range(count):
            chrom = np.random.choice(self.chromosomes)
            pos = np.random.randint(0, self.max_position)
            queries.append((chrom, pos))
        
        return queries
    
    def generate_range_queries(self, count):
        """Generate random range queries."""
        np.random.seed(44)  # Different seed
        queries = []
        
        for _ in range(count):
            chrom = np.random.choice(self.chromosomes)
            start = np.random.randint(0, self.max_position)
            # Range sizes from 100 to 10000
            length = np.random.randint(100, 10000)
            end = min(start + length, self.max_position)
            queries.append((chrom, start, end))
        
        return queries
    
    def test_insertion_benchmark(self):
        """Benchmark feature insertion performance."""
        for count in self.feature_counts:
            features = self.generate_features(count)
            # Test on different store types
            for name, store_type in self.store_types.items():
                print(f"Testing insertion for {count} features, {name}...")
                store = GenomicFeatureStore(store_type=store_type)
                
                # Add fetures, measure insertion time
                start_time = time.time()
                with store:
                    for feature in features:
                        store.add(feature)
                end_time = time.time()
                elapsed = end_time - start_time
                
                # Assert that insertion completes within threshold time (if enabled)
                self.assertLessEqual(elapsed, self.max_insertion_time[count], f"{name} insertion time exceeded threshold for {count} features. Elapsed: {elapsed:.4f} seconds")
            
                # Check that the store has all the features
                for chrom in self.chromosomes:
                    expected_ids = {f.id for f in features if f.chrom == chrom}
                    ids = {f.id for f in store[chrom].get_features()}
                    self.assertEqual(ids, expected_ids, f"{name} did not store all features correctly for {count} features in chromosome '{chrom}'")
                    # Query all features in the store
                    end = max(f.end for f in features if f.chrom == chrom)
                    results = store[chrom].get_by_interval(0, end + 1)
                    result_ids = {f.id for f in results}
                    print(f"LEN: {chrom}, {len(expected_ids)} / {len(result_ids)}")
                    self.assertEqual(result_ids, expected_ids, f"{name} did not return all features correctly for {count} features in chromosome '{chrom}' on query [0, {end + 1}]")


    def test_position_query_benchmark(self):
        """Benchmark position query performance."""
        for count in self.feature_counts:
            features = self.generate_features(count)
            queries = self.generate_position_queries(self.query_counts)
            # Test on different store types
            for name, store_type in self.store_types.items():
                store = GenomicFeatureStore(store_type=store_type)
                store_reference = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)
                
                # Add features
                with store, store_reference:
                    for feature in features:
                        store.add(feature)
                        store_reference.add(feature)
                
                # Measure query time
                start_time = time.time()
                all_results = []
                for chrom, pos in queries:
                    results = store.get_by_position(chrom, pos)
                    result_ids = {f.id for f in results}
                    all_results.append((chrom, pos, result_ids))
                end_time = time.time()
                elapsed = end_time - start_time
                
                # Assert that query completes within threshold time (if enabled)
                self.assertLessEqual(elapsed, self.max_query_time[count], f"{name} position query time exceeded threshold for {count} features")
            
                # Verify store return values
                for chrom, pos, result_ids in all_results:
                    ref_ids = {f.id for f in store_reference.get_by_position(chrom, pos)}
                    self.assertEqual(ref_ids, result_ids, f"Total features found differs for {chrom}:{pos} between {name} and reference")


    def test_range_query_benchmark(self):
        """Benchmark range query performance."""
        for count in self.feature_counts:
            features = self.generate_features(count)
            queries = self.generate_range_queries(self.query_counts)
            
            for name, store_type in self.store_types.items():
                store = GenomicFeatureStore(store_type=store_type)
                store_reference = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)

                # Add features
                with store, store_reference:
                    for feature in features:
                        store.add(feature)
                        store_reference.add(feature)
                
                # Measure query time
                start_time = time.time()
                total_results = 0
                all_results = []
                for chrom, start, end in queries:
                    results = store.get_by_interval(chrom, start, end)
                    if end == 733876:
                        print(f"RESULT: {start}-{end}, {results}")
                    total_results += len(results)
                    all_results.append((chrom, start, end, {f.id for f in results} if results else set()))
                end_time = time.time()
                elapsed = end_time - start_time
                
                # # Assert that query completes within threshold time (if enabled)
                # self.assertLessEqual(elapsed, self.max_query_time[count], f"{name} range query time exceeded threshold for {count} features")
            
                # Verify store return values
                for chrom, start, end, result_ids in all_results:
                    ref_ids = {f.id for f in store_reference.get_by_interval(chrom, start, end)}
                    diff_ids = ref_ids - result_ids if len(ref_ids) > len(result_ids) else result_ids - ref_ids
                    diff = [f for f in features if f.id in diff_ids]
                    self.assertEqual(ref_ids, result_ids, f"\nTotal features found differs for {chrom}:{start}-{end} between {name} and reference\n"
                                     + f"\tExpected ({len(ref_ids)}): {ref_ids}\n"
                                     + f"\tFound    ({len(result_ids)}): {result_ids}\n" 
                                     + f"\tDifferent IDs: {diff_ids}\n\tDiff features: {diff}"
                                     )


    def test_nearest_query_benchmark(self):
        """Benchmark nearest feature query performance."""
        for count in self.feature_counts:
            features = self.generate_features(count)
            queries = self.generate_position_queries(self.query_counts)
            print(f"Testing nearest query for {count} features, {len(queries)} queries...")
            # Test on different store types
            for name, store_type in self.store_types.items():
                print(f"Testing {name} with {count} features...")
                store = GenomicFeatureStore(store_type=store_type)
                store_reference = GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE)
                
                # Add features
                with store, store_reference:
                    for feature in features:
                        store.add(feature)
                        store_reference.add(feature)
                
                # Measure query time
                start_time = time.time()
                all_results = []                
                for chrom, pos in queries:
                    result = store.get_nearest(chrom, pos)
                    all_results.append((chrom, pos, result.id if result else None))
                    print(f"Querying {chrom}:{pos} -> Result: {all_results[-1]}")
                end_time = time.time()
                elapsed = end_time - start_time
                
                # Assert that query completes within threshold time (if enabled)
                self.assertLessEqual(elapsed, self.max_query_time[count], f"{name} nearest query time exceeded threshold for {count} features")

                # Verify store return values
                for chrom, pos, result_id in all_results:
                    ref_result = store_reference.get_nearest(chrom, pos)
                    ref_id = ref_result.id if ref_result else None
                    self.assertEqual(ref_id, result_id, f"Features found differs for {chrom}:{pos} between {name} and reference")


if __name__ == "__main__":
    unittest.main()