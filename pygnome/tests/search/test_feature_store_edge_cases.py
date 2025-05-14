"""Edge case tests for genomic feature store implementations."""

import unittest
import numpy as np

from pygnome.genomics import GenomicFeature, Strand
from pygnome.genomics.feature_store import (
    BinnedGenomicStore, GenomicFeatureStoreImpl, 
    IntervalTreeStore, PositionHashStore, StoreType
)


class TestFeatureStoreEdgeCases(unittest.TestCase):
    """Edge case tests for genomic feature store implementations."""
    
    def test_empty_store(self):
        """Test operations on empty stores."""
        # Create empty stores
        interval_store = IntervalTreeStore()
        binned_store = BinnedGenomicStore()
        position_store = PositionHashStore()
        
        # Test position query on empty store
        self.assertEqual(len(interval_store.get_by_position(100)), 0)
        self.assertEqual(len(binned_store.get_by_position(100)), 0)
        self.assertEqual(len(position_store.get_by_position(100)), 0)
        
        # Test range query on empty store
        self.assertEqual(len(interval_store.get_by_interval(100, 200)), 0)
        self.assertEqual(len(binned_store.get_by_interval(100, 200)), 0)
        self.assertEqual(len(position_store.get_by_interval(100, 200)), 0)
        
        # Test nearest feature on empty store
        self.assertIsNone(interval_store.get_nearest(100))
        self.assertIsNone(binned_store.get_nearest(100))
        self.assertIsNone(position_store.get_nearest(100))
    
    def test_extreme_coordinates(self):
        """Test features with extreme coordinates."""
        # Create stores
        interval_store = IntervalTreeStore()
        binned_store = BinnedGenomicStore()
        position_store = PositionHashStore()
        
        # Feature with very large coordinates
        large_coord_feature = GenomicFeature(
            id="large_coord",
            chrom="chr1",
            start=1_000_000_000,  # 1 billion
            end=1_000_001_000,
            strand=Strand.POSITIVE
        )
        
        # Feature at position 0
        zero_pos_feature = GenomicFeature(
            id="zero_pos",
            chrom="chr1",
            start=0,
            end=100,
            strand=Strand.POSITIVE
        )
        
        # Add features to stores
        for store in [interval_store, binned_store, position_store]:
            store.add_feature(large_coord_feature)
            store.add_feature(zero_pos_feature)
        
        # Test queries at extreme positions
        for store in [interval_store, binned_store, position_store]:
            # Query at position 0
            features = store.get_by_position(0)
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].id, "zero_pos")
            
            # Query at large position
            features = store.get_by_position(1_000_000_500)
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].id, "large_coord")
            
            # Range query at extreme positions
            features = store.get_by_interval(0, 50)
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].id, "zero_pos")
            
            features = store.get_by_interval(1_000_000_500, 1_000_000_600)
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].id, "large_coord")
    
    def test_invalid_ranges(self):
        """Test queries with invalid ranges."""
        # Create stores
        stores = {
            "interval_tree": GenomicFeatureStoreImpl(store_type=StoreType.INTERVAL_TREE),
            "binned": GenomicFeatureStoreImpl(store_type=StoreType.BINNED),
            "position_hash": GenomicFeatureStoreImpl(store_type=StoreType.POSITION_HASH)
        }
        
        # Add a sample feature
        feature = GenomicFeature(
            id="sample",
            chrom="chr1",
            start=1000,
            end=2000,
            strand=Strand.POSITIVE
        )
        
        for name, store in stores.items():
            store.add_feature(feature)
        
        # Test range where end < start
        for name, store in stores.items():
            features = store.get_by_interval("chr1", 2000, 1000)
            self.assertEqual(len(features), 0, f"Expected empty result for invalid range in {name}")
    
    def test_very_large_features(self):
        """Test handling of extremely large features."""
        # Create a feature spanning 10 million bases
        huge_feature = GenomicFeature(
            id="huge",
            chrom="chr1",
            start=1_000_000,
            end=11_000_000,  # 10 million bp long
            strand=Strand.POSITIVE
        )
        
        # Create stores
        interval_store = IntervalTreeStore()
        binned_store = BinnedGenomicStore()
        position_store = PositionHashStore()
        
        # Add the huge feature
        interval_store.add_feature(huge_feature)
        binned_store.add_feature(huge_feature)
        position_store.add_feature(huge_feature)
        
        # Test position queries at various points
        positions = [1_000_000, 5_000_000, 11_000_000]
        for pos in positions:
            for store in [interval_store, binned_store, position_store]:
                features = store.get_by_position(pos)
                self.assertEqual(len(features), 1,
                                f"Expected 1 feature at position {pos} for {store.__class__.__name__}")
                self.assertEqual(features[0].id, "huge")
        
        # Test position query just outside the range
        for store in [interval_store, binned_store, position_store]:
            features = store.get_by_position(11_000_001)
            self.assertEqual(len(features), 0,
                            f"Expected 0 features at position 11,000,001 for {store.__class__.__name__}")
    
    def test_many_features_same_position(self):
        """Test handling many features at the same position."""
        # Create stores
        interval_store = IntervalTreeStore()
        binned_store = BinnedGenomicStore()
        position_store = PositionHashStore()
        
        # Create 100 features all containing position 1000
        features = []
        for i in range(100):
            features.append(
                GenomicFeature(
                    id=f"overlap_{i}",
                    chrom="chr1",
                    start=1000 - i,  # Each starts at a different position
                    end=1000 + i,    # Each ends at a different position
                    strand=Strand.POSITIVE
                )
            )
        
        # Add features to stores
        for feature in features:
            interval_store.add_feature(feature)
            binned_store.add_feature(feature)
            position_store.add_feature(feature)
        
        # Test position query at the common position
        for store in [interval_store, binned_store, position_store]:
            result = store.get_by_position(1000)
            self.assertEqual(len(result), 100,
                            f"Expected 100 features at position 1000 for {store.__class__.__name__}")
            
            # Verify all features are present
            ids = {f.id for f in result}
            for i in range(100):
                self.assertIn(f"overlap_{i}", ids)
    
    def test_bin_size_variations(self):
        """Test binned store with different bin sizes."""
        # Create features
        features = []
        for i in range(10):
            features.append(
                GenomicFeature(
                    id=f"feature_{i}",
                    chrom="chr1",
                    start=i * 100,
                    end=i * 100 + 50,
                    strand=Strand.POSITIVE
                )
            )
        
        # Test different bin sizes
        bin_sizes = [10, 100, 1000, 10000]
        for bin_size in bin_sizes:
            store = BinnedGenomicStore(bin_size=bin_size)
            
            # Add features
            for feature in features:
                store.add_feature(feature)
            
            # Test position queries
            for i in range(10):
                pos = i * 100 + 25  # Middle of each feature
                result = store.get_by_position(pos)
                self.assertEqual(len(result), 1, 
                                f"Expected 1 feature at position {pos} with bin_size={bin_size}")
                self.assertEqual(result[0].id, f"feature_{i}")
            
            # Test range query spanning multiple features
            result = store.get_by_interval(250, 650)
            self.assertEqual(len(result), 5, 
                            f"Expected 5 features in range 250-650 with bin_size={bin_size}")
    
    def test_feature_removal(self):
        """Test feature removal if implemented."""
        # This is a placeholder for future implementation of feature removal
        # Currently, our stores don't support removal, but this test can be
        # uncommented when that functionality is added
        
        # # Create a store with features
        # store = GenomicFeatureStoreImpl(store_type=StoreType.INTERVAL_TREE)
        # 
        # # Add features
        # for i in range(10):
        #     feature = GenomicFeature(
        #         id=f"feature_{i}",
        #         chrom="chr1",
        #         start=i * 100,
        #         end=i * 100 + 50,
        #         strand=Strand.POSITIVE
        #     )
        #     store.add_feature(feature)
        # 
        # # Verify features exist
        # self.assertEqual(len(store.get_by_position("chr1", 150)), 1)
        #
        # # Remove a feature
        # store.remove_feature("feature_1")
        #
        # # Verify it's gone
        # self.assertEqual(len(store.get_by_position("chr1", 150)), 0)
        pass


if __name__ == "__main__":
    unittest.main()