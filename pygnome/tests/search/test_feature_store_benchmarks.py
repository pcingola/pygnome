"""Benchmark tests for genomic feature store implementations."""

import unittest
import time
import numpy as np

from pygnome.genomics import GenomicFeature, Strand
from pygnome.genomics.feature_store import (
    GenomicFeatureStoreImpl, StoreType
)


class TestFeatureStoreBenchmarks(unittest.TestCase):
    """Benchmark tests for genomic feature store implementations."""
    
    def setUp(self):
        """Set up benchmark parameters."""
        # Skip benchmarks by default
        if not hasattr(self, 'run_benchmarks'):
            self.skipTest("Benchmarks disabled")
        
        # Parameters
        self.feature_counts = [1000, 10000]  # Reduced for faster tests
        self.query_counts = 100
        self.chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5"]
        self.max_position = 1_000_000
        
        # Store types
        self.store_types = {
            "Interval Tree": StoreType.INTERVAL_TREE,
            "Binned": StoreType.BINNED,
            "Position Hash": StoreType.POSITION_HASH
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
            feature = GenomicFeature(
                id=f"feature_{i}",
                chrom=chrom,
                start=start,
                end=end,
                strand=Strand.POSITIVE if np.random.random() > 0.5 else Strand.NEGATIVE
            )
            
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
        print("\n=== Insertion Performance ===")
        
        for count in self.feature_counts:
            print(f"\nGenerating {count} features...")
            features = self.generate_features(count)
            
            for name, store_type in self.store_types.items():
                print(f"Testing {name} with {count} features...")
                store = GenomicFeatureStoreImpl(store_type=store_type)
                
                # Measure insertion time
                start_time = time.time()
                for feature in features:
                    store.add_feature(feature)
                end_time = time.time()
                
                elapsed = end_time - start_time
                print(f"{name}: {elapsed:.4f} seconds")
    
    def test_position_query_benchmark(self):
        """Benchmark position query performance."""
        print("\n=== Position Query Performance ===")
        
        for count in self.feature_counts:
            print(f"\nGenerating {count} features...")
            features = self.generate_features(count)
            queries = self.generate_position_queries(self.query_counts)
            
            for name, store_type in self.store_types.items():
                print(f"Testing {name} with {count} features...")
                store = GenomicFeatureStoreImpl(store_type=store_type)
                
                # Add features
                for feature in features:
                    store.add_feature(feature)
                
                # Measure query time
                start_time = time.time()
                total_results = 0
                for chrom, pos in queries:
                    results = store.get_by_position(chrom, pos)
                    total_results += len(results)
                end_time = time.time()
                
                elapsed = end_time - start_time
                print(f"{name}: {elapsed:.4f} seconds, found {total_results} features")
    
    def test_range_query_benchmark(self):
        """Benchmark range query performance."""
        print("\n=== Range Query Performance ===")
        
        for count in self.feature_counts:
            print(f"\nGenerating {count} features...")
            features = self.generate_features(count)
            queries = self.generate_range_queries(self.query_counts)
            
            for name, store_type in self.store_types.items():
                print(f"Testing {name} with {count} features...")
                store = GenomicFeatureStoreImpl(store_type=store_type)
                
                # Add features
                for feature in features:
                    store.add_feature(feature)
                
                # Measure query time
                start_time = time.time()
                total_results = 0
                for chrom, start, end in queries:
                    results = store.get_by_interval(chrom, start, end)
                    total_results += len(results)
                end_time = time.time()
                
                elapsed = end_time - start_time
                print(f"{name}: {elapsed:.4f} seconds, found {total_results} features")
    
    def test_nearest_query_benchmark(self):
        """Benchmark nearest feature query performance."""
        print("\n=== Nearest Feature Query Performance ===")
        
        for count in self.feature_counts:
            print(f"\nGenerating {count} features...")
            features = self.generate_features(count)
            queries = self.generate_position_queries(self.query_counts)
            
            for name, store_type in self.store_types.items():
                print(f"Testing {name} with {count} features...")
                store = GenomicFeatureStoreImpl(store_type=store_type)
                
                # Add features
                for feature in features:
                    store.add_feature(feature)
                
                # Measure query time
                start_time = time.time()
                found_count = 0
                for chrom, pos in queries:
                    result = store.get_nearest(chrom, pos)
                    if result is not None:
                        found_count += 1
                end_time = time.time()
                
                elapsed = end_time - start_time
                print(f"{name}: {elapsed:.4f} seconds, found {found_count} features")


if __name__ == "__main__":
    # To run benchmarks:
    # test = TestFeatureStoreBenchmarks()
    # test.run_benchmarks = True
    # test.setUp()
    # test.test_insertion_benchmark()
    # test.test_position_query_benchmark()
    # test.test_range_query_benchmark()
    # test.test_nearest_query_benchmark()
    unittest.main()