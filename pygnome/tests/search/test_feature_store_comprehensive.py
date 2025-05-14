"""Comprehensive tests for genomic feature store implementations."""

import unittest
import numpy as np

from pygnome.feature_store.binned_store import BinnedGenomicStore
from pygnome.feature_store.genomic_feature_store import GenomicFeatureStore, StoreType
from pygnome.feature_store.interval_tree_store import IntervalTreeStore
from pygnome.genomics import GenomicFeature, Strand


class TestFeatureStoreComprehensive(unittest.TestCase):
    """Comprehensive test cases for genomic feature store implementations."""
    
    def setUp(self):
        """Set up test data with various edge cases."""
        # Create test features
        self.features = []
        
        # Features on multiple chromosomes
        self.chromosomes = ["chr1", "chr2", "chr3", "chrX", "chrY"]
        
        # 1. Basic features on different chromosomes
        for chrom in self.chromosomes:
            for i in range(10):
                start = i * 1000
                end = start + 500
                self.features.append(GenomicFeature(id=f"{chrom}_basic_{i}", chrom=chrom, start=start, end=end, strand=Strand.POSITIVE))
        
        # 2. Overlapping features
        for i in range(5):
            start = i * 200
            end = start + 300  # Ensures overlap with next feature
            self.features.append(GenomicFeature(id=f"overlap_{i}", chrom="chr1", start=start, end=end, strand=Strand.POSITIVE))
        
        # 3. Adjacent features (end of one = start of next)
        for i in range(5):
            start = i * 1000
            end = start + 1000  # Changed from 999 to 1000 to include position 999
            self.features.append(GenomicFeature(id=f"adjacent_{i}_a", chrom="chr2", start=start, end=end, strand=Strand.POSITIVE))
            self.features.append(GenomicFeature(id=f"adjacent_{i}_b", chrom="chr2", start=end, end=end + 1000, strand=Strand.POSITIVE))
        
        # 4. Nested features (one completely inside another)
        for i in range(3):
            outer_start = i * 5000
            outer_end = outer_start + 2000
            inner_start = outer_start + 500
            inner_end = outer_end - 500
            
            self.features.append(GenomicFeature(id=f"outer_{i}", chrom="chr3", start=outer_start, end=outer_end, strand=Strand.POSITIVE))
            self.features.append(GenomicFeature(id=f"inner_{i}",chrom="chr3",start=inner_start,end=inner_end,strand=Strand.POSITIVE))
        
        # 5. Zero-length features (start = end)
        for i in range(5):
            pos = i * 1000
            self.features.append(GenomicFeature(id=f"zero_length_{i}", chrom="chrX", start=pos, end=pos, strand=Strand.POSITIVE))
        
        # 6. Very large features
        for i in range(2):
            start = i * 1000000
            end = start + 500000
            self.features.append(GenomicFeature(id=f"large_{i}", chrom="chrY", start=start, end=end, strand=Strand.POSITIVE))
        
        # Create main stores
        self.stores = {
            StoreType.BRUTE_FORCE: GenomicFeatureStore(store_type=StoreType.BRUTE_FORCE),
            StoreType.INTERVAL_TREE: GenomicFeatureStore(store_type=StoreType.INTERVAL_TREE),
            StoreType.BINNED: GenomicFeatureStore(store_type=StoreType.BINNED),
        }
        self.stores_list = list(self.stores.values())

        # Add features to stores        
        for store in self.stores_list:
            with store:
                for feature in self.features:
                        store.add(feature)
    
    def test_empty_queries(self):
        """Test queries that should return empty results."""
        # Test position query on non-existent position
        for store in self.stores_list:
            features = store.get_by_position("chr1", 99999999)
            self.assertEqual(len(features), 0, f"Expected empty result for {store.__class__.__name__}")
        
        # Test range query on non-existent range
        for store in self.stores_list:
            features = store.get_by_interval("chr1", 99999999, 99999999 + 1000)
            self.assertEqual(len(features), 0, f"Expected empty result for {store.__class__.__name__}")
        
        # Test position query on non-existent chromosome
        for store_type, store in self.stores.items():
            features = store.get_by_position("non_existent_chrom", 1000)
            self.assertEqual(len(features), 0, f"Expected empty result for {store_type}")
        
        # Test range query on non-existent chromosome
        for store_type, store in self.stores.items():
            features = store.get_by_interval("non_existent_chrom", 1000, 2000)
            self.assertEqual(len(features), 0, f"Expected empty result for {store_type}")
    
    def test_overlapping_features(self):
        """Test queries with overlapping features."""
        # Position 250 should overlap with overlap_0 and overlap_1
        for store in self.stores_list:
            features = store.get_by_position("chr1", 250)
            ids = {f.id for f in features}
            self.assertIn("overlap_0", ids)
            self.assertIn("overlap_1", ids)
            self.assertIn("chr1_basic_0", ids)
            self.assertEqual(len(features), 3, f"Expected 3 features at position 250 for {store.__class__.__name__}")
        
        # Range 150-350 should overlap with overlap_0, overlap_1, and chr1_basic_0
        for store in self.stores_list:
            features = store.get_by_interval("chr1", 150, 350)
            ids = {f.id for f in features}
            self.assertIn("overlap_0", ids)
            self.assertIn("overlap_1", ids)
            self.assertIn("chr1_basic_0", ids)
            self.assertEqual(len(features), 3, f"Expected 3 features in range 150-350 for {store.__class__.__name__}")
    
    def test_adjacent_features(self):
        """Test queries with adjacent features."""
        # Position at the boundary should only match one feature
        for store_type, store in self.stores.items():
            # Test at the boundary of adjacent_0_a and adjacent_0_b
            features = store.get_by_position("chr2", 999)
            self.assertEqual(len(features), 1, f"Expected 1 feature at boundary for {store_type}")
            self.assertEqual(features[0].id, "adjacent_0_a")
            
            features = store.get_by_position("chr2", 1000)
            # Check that we have the expected features
            feature_ids = [f.id for f in features]
            self.assertIn("adjacent_0_b", feature_ids)
            self.assertIn("adjacent_1_a", feature_ids)
            self.assertIn("chr2_basic_1", feature_ids)
            self.assertEqual(len(features), 3, f"Expected 3 features at boundary for {store_type}")
        
        # Range query spanning the boundary should match both features
        for store_type, store in self.stores.items():
            features = store.get_by_interval("chr2", 999, 1001)
            feature_ids = [f.id for f in features]
            self.assertIn("adjacent_0_a", feature_ids)
            self.assertIn("adjacent_0_b", feature_ids)
            self.assertIn("adjacent_1_a", feature_ids)
            self.assertIn("chr2_basic_1", feature_ids)
            self.assertEqual(len(features), 4, f"Expected 4 features for boundary range for {store_type}")
    
    def test_nested_features(self):
        """Test queries with nested features."""
        # Position in the middle should match both outer and inner
        for store_type, store in self.stores.items():
            features = store.get_by_position("chr3", 1000)
            ids = {f.id for f in features}
            self.assertIn("outer_0", ids)
            self.assertIn("inner_0", ids)
            self.assertIn("chr3_basic_1", ids)
            self.assertEqual(len(features), 3, f"Expected 3 features at position 1000 for {store_type}")
        
        # Range query in the middle should match both outer and inner
        for store_type, store in self.stores.items():
            features = store.get_by_interval("chr3", 1000, 1100)
            ids = {f.id for f in features}
            self.assertIn("outer_0", ids)
            self.assertIn("inner_0", ids)
            self.assertIn("chr3_basic_1", ids)
            self.assertEqual(len(features), 3, f"Expected 3 features in range 1000-1100 for {store_type}")
    
    def test_zero_length_features(self):
        """Test queries with zero-length features."""
        # Position query at exact position should find both features at position 0
        for store_type, store in self.stores.items():
            features = store.get_by_position("chrX", 0)
            feature_ids = [f.id for f in features]
            print(f"Features at position 0 on chrX for {store_type}: {feature_ids}")
            ids = {f.id for f in features}
            self.assertIn("zero_length_0", ids)
            self.assertIn("chrX_basic_0", ids)
            self.assertEqual(len(features), 2, f"Expected 2 features at position 0 for {store_type}")
        
        # Range query spanning the position should find both features
        for store_type, store in self.stores.items():
            features = store.get_by_interval("chrX", 0, 10)
            ids = {f.id for f in features}
            self.assertIn("zero_length_0", ids)
            self.assertIn("chrX_basic_0", ids)
            self.assertEqual(len(features), 2, f"Expected 2 features in range 0-10 for {store_type}")
    
    def test_large_features(self):
        """Test queries with very large features."""
        # Position in the middle of a large feature
        for store_type, store in self.stores.items():
            features = store.get_by_position("chrY", 250000)
            self.assertEqual(len(features), 1, f"Expected 1 large feature for {store_type}")
            self.assertEqual(features[0].id, "large_0")
        
        # Range query within a large feature
        for store_type, store in self.stores.items():
            features = store.get_by_interval("chrY", 250000, 260000)
            self.assertEqual(len(features), 1, f"Expected 1 large feature for {store_type}")
            self.assertEqual(features[0].id, "large_0")
    
    def test_nearest_feature(self):
        """Test nearest feature queries."""
        # Test nearest feature when there's an exact match
        for store_type, store in self.stores.items():
            feature = store.get_nearest("chr1", 500)
            self.assertIsNotNone(feature, f"Expected to find a feature for {store_type}")
            # Should find a feature with start=500 or one that contains 500
            self.assertTrue(
                feature.start == 500 or (feature.start <= 500 and feature.end >= 500),
                f"Expected feature to contain position 500 for {store_type}"
            )
        
        # Test nearest feature when there's no exact match
        for store_type, store in self.stores.items():
            feature = store.get_nearest("chr1", 1550)
            self.assertIsNotNone(feature, f"Expected to find a feature for {store_type}")
            # Should find the closest feature
            self.assertTrue(
                abs(feature.start - 1550) <= 50 or abs(feature.end - 1550) <= 50,
                f"Expected feature to be close to position 1550 for {store_type}"
            )
        
        # Test nearest feature with max_distance
        for store_type, store in self.stores.items():
            # Should find a feature within 100bp
            feature = store.get_nearest("chr1", 1550, max_distance=100)
            self.assertIsNotNone(feature, f"Expected to find a feature within 100bp for {store_type}")
            
            # Should not find a feature within 10bp
            feature = store.get_nearest("chr1", 1550, max_distance=10)
            self.assertIsNone(feature, f"Expected no feature within 10bp for {store_type}")
        
        # Test nearest feature on non-existent chromosome
        for store_type, store in self.stores.items():
            feature = store.get_nearest("non_existent_chrom", 1000)
            self.assertIsNone(feature, f"Expected no feature on non-existent chromosome for {store_type}")
    
    def test_multi_chromosome_queries(self):
        """Test queries across multiple chromosomes."""
        # Count features on each chromosome
        expected_counts = {
            "chr1": 15,  # 10 basic + 5 overlapping
            "chr2": 20,  # 10 basic + 10 adjacent
            "chr3": 16,  # 10 basic + 6 nested
            "chrX": 15,  # 10 basic + 5 zero-length (all zero-length features are now included)
            "chrY": 11,  # 10 basic + 1 large (large_1 is not included because it starts at position 1000000)
        }
        
        # Print out all the zero-length features for debugging
        for i in range(5):
            print(f"zero_length_{i}: start={i * 1000}, end={i * 1000}")

        for store_type, store in self.stores.items():
            for chrom, expected_count in expected_counts.items():
                features = store.get_by_interval(chrom, 0, 1000000)
                if (chrom == "chrX" or chrom == "chrY") and len(features) != expected_count:
                    # Print out the features for debugging
                    feature_ids = [f.id for f in features]
                    print(f"Features on {chrom} for {store_type}: {feature_ids}")
                self.assertEqual(
                    len(features), expected_count,
                    f"Expected {expected_count} features on {chrom} for {store_type}"
                )
    
    def test_add_features_batch(self):
        """Test adding features in batch."""
        # Create new features
        new_features = []
        for i in range(10):
            new_features.append(GenomicFeature(id=f"batch_{i}", chrom="chr4", start=i * 100, end=i * 100 + 50, strand=Strand.POSITIVE))
        
        # Test batch add for GenomicFeatureStore
        for store_type, store in self.stores.items():
            # Add features in batch using context manager
            with store:
                store.add_features(new_features)
            
            # Verify they were added
            for i in range(10):
                features = store.get_by_position("chr4", i * 100)
                self.assertEqual(len(features), 1, f"Expected 1 feature at position {i*100} for {store_type}")
                self.assertEqual(features[0].id, f"batch_{i}")


if __name__ == "__main__":
    unittest.main()