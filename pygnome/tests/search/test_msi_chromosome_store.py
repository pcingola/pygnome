"""
Unit tests for the MsiChromosomeStore class.
"""

import unittest
import numpy as np
import random
import string
import time

from pygnome.feature_store.msi_chromosome_store import MsiChromosomeStore, MsiSiteCounter
from pygnome.parsers.msi.msi_site_record import MsiSiteRecord
from pygnome.genomics.strand import Strand


class TestMsiChromosomeStore(unittest.TestCase):
    """Test cases for the MsiChromosomeStore class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create some test MSI site records
        self.records = [
            MsiSiteRecord(repeat_unit_length=2, repeat_unit_binary=9, repeat_times=3, left_flank_binary=258, right_flank_binary=409, repeat_unit_bases="GC", left_flank_bases="CAAAG", right_flank_bases="CGCGC", chromosome="chr1", location=10000 ),
            MsiSiteRecord(repeat_unit_length=4, repeat_unit_binary=149, repeat_times=3, left_flank_binary=150, right_flank_binary=685, repeat_unit_bases="GCCC", left_flank_bases="AGCCG", right_flank_bases="GGGTC", chromosome="chr1", location=10100),
            MsiSiteRecord(repeat_unit_length=2, repeat_unit_binary=2, repeat_times=3, left_flank_binary=665, right_flank_binary=614, repeat_unit_bases="AG", left_flank_bases="GGCGC", right_flank_bases="GCGCG", chromosome="chr1", location=10200)
        ]
        
    def test_msi_site_counter(self):
        """Test the MsiSiteCounter class."""
        counter = MsiSiteCounter(bin_size=1000)
        
        # Add records to the counter
        for record in self.records:
            counter.add(record)
            
        # Check counts
        self.assertEqual(counter.get_count("chr1"), 3)
        self.assertEqual(counter.get_count("chr2"), 0)
        
        # Get max lengths from the counter
        max_lengths = counter.get_max_lengths("chr1")
        
        # The length is calculated as end - start
        # For the first record: start=10000, end=10000 + (2*3) = 10006, length = 6
        # For the second record: start=10100, end=10100 + (4*3) = 10112, length = 12
        # For the third record: start=10200, end=10200 + (2*3) = 10206, length = 6
        # All three records are in the same bin (bin 10), so max_lengths[10] = 12
        self.assertEqual(max_lengths[10], 12)  # 10000 // 1000 = 10
        
    def test_msi_chromosome_store(self):
        """Test the MsiChromosomeStore class."""
        # Create a counter and count records
        counter = MsiSiteCounter(bin_size=1000)
        for record in self.records:
            counter.add(record)
            
        # Create a store with pre-allocated arrays
        store = MsiChromosomeStore(
            chromosome="chr1",
            feature_count=counter.get_count("chr1"),
            max_lengths_by_bin=counter.get_max_lengths("chr1"),
            bin_size=1000
        )
        
        # Add records to the store
        store.index_build_start()
        for record in self.records:
            store.add(record)
        store.index_build_end()
        
        # Test get_by_position
        features = store.get_by_position(10000)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0].start, 10000)
        self.assertEqual(features[0].end, 10006)
        
        features = store.get_by_position(10102)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0].start, 10100)
        self.assertEqual(features[0].end, 10112)
        
        features = store.get_by_position(10050)
        self.assertEqual(len(features), 0)
        
        # Test get_by_interval
        features = store.get_by_interval(10000, 10150)
        self.assertEqual(len(features), 2)
        self.assertEqual(features[0].start, 10000)
        self.assertEqual(features[1].start, 10100)
        
        features = store.get_by_interval(10150, 10250)
        self.assertEqual(len(features), 1)
        self.assertEqual(features[0].start, 10200)
        
        features = store.get_by_interval(10050, 10090)
        self.assertEqual(len(features), 0)
        
    def test_empty_store(self):
        """Test an empty MsiChromosomeStore."""
        # Empty store
        store = MsiChromosomeStore(chromosome="chr1")
        store.index_build_start()
        store.index_build_end()
        
        self.assertEqual(len(store.get_by_position(10000)), 0, "Empty store should return no features for position query")
        self.assertEqual(len(store.get_by_interval(10000, 10100)), 0, "Empty store should return no features for interval query")
    
    def _create_test_record(self):
        """Helper method to create a test record."""
        return MsiSiteRecord(
            repeat_unit_length=2,
            repeat_unit_binary=9,
            repeat_times=3,
            left_flank_binary=258,
            right_flank_binary=409,
            repeat_unit_bases="GC",
            left_flank_bases="CAAAG",
            right_flank_bases="CGCGC",
            chromosome="chr1",
            location=10000
        )
    
    def _create_store_with_single_feature(self):
        """Helper method to create a store with a single feature."""
        record = self._create_test_record()
        
        # Create a counter to calculate max lengths
        counter = MsiSiteCounter(bin_size=1000)
        counter.add(record)
        
        # Get max lengths from the counter
        max_lengths = counter.get_max_lengths("chr1")
        
        # Create a store with pre-calculated max lengths
        store = MsiChromosomeStore(
            chromosome="chr1",
            feature_count=1,
            max_lengths_by_bin=max_lengths,
            bin_size=1000
        )
        store.index_build_start()
        store.add(record)
        store.index_build_end()
        
        return store, record
        
    def test_single_feature_exact_match(self):
        """Test finding a feature at its exact start position."""
        store, record = self._create_store_with_single_feature()
        
        features = store.get_by_position(10000)
        self.assertEqual(len(features), 1, "Should find feature at its start position")
        self.assertEqual(features[0].start, 10000, "Feature start position should match")
        self.assertEqual(features[0].end, 10006, "Feature end position should match")
    
    def test_position_within_feature(self):
        """Test finding a feature when querying a position within the feature."""
        store, record = self._create_store_with_single_feature()
        
        features = store.get_by_position(10002)
        self.assertEqual(len(features), 1, "Should find feature when position is within feature")
        self.assertEqual(features[0].start, 10000, "Feature start position should match")
    
    def test_position_outside_feature(self):
        """Test querying a position outside any feature."""
        store, record = self._create_store_with_single_feature()
        
        features = store.get_by_position(10007)
        self.assertEqual(len(features), 0, "Should not find any features when position is outside feature")
    
    def test_interval_overlap(self):
        """Test finding a feature when the query interval overlaps with it."""
        store, record = self._create_store_with_single_feature()
        
        features = store.get_by_interval(9990, 10003)
        self.assertEqual(len(features), 1, "Should find feature when interval overlaps with it")
        self.assertEqual(features[0].start, 10000, "Feature start position should match")
    
    def test_interval_no_overlap(self):
        """Test querying an interval that doesn't overlap with any feature."""
        store, record = self._create_store_with_single_feature()
        
        features = store.get_by_interval(10007, 10100)
        self.assertEqual(len(features), 0, "Should not find any features when interval doesn't overlap")
        
    def test_binary_search(self):
        """Test the binary search implementation."""
        # Create a store with features at specific positions
        store = MsiChromosomeStore(chromosome="chr1", feature_count=5)
        store.index_build_start()
        
        # Add features at positions 10, 20, 30, 40, 50
        for i, pos in enumerate([10, 20, 30, 40, 50]):
            record = MsiSiteRecord(
                repeat_unit_length=2,
                repeat_unit_binary=9,
                repeat_times=3,
                left_flank_binary=258,
                right_flank_binary=409,
                repeat_unit_bases="GC",
                left_flank_bases="CAAAG",
                right_flank_bases="CGCGC",
                chromosome="chr1",
                location=pos
            )
            store.add(record)
            
        store.index_build_end()
        
        # Test binary search at exact positions
        self.assertEqual(store._binary_search_position(10), 0)
        self.assertEqual(store._binary_search_position(20), 1)
        self.assertEqual(store._binary_search_position(30), 2)
        self.assertEqual(store._binary_search_position(40), 3)
        self.assertEqual(store._binary_search_position(50), 4)
        
        # Test binary search between positions
        self.assertEqual(store._binary_search_position(15), 1)
        self.assertEqual(store._binary_search_position(25), 2)
        self.assertEqual(store._binary_search_position(35), 3)
        self.assertEqual(store._binary_search_position(45), 4)
        
        # Test binary search before first position
        self.assertEqual(store._binary_search_position(5), 0)
        
        # Test binary search after last position
        self.assertEqual(store._binary_search_position(55), 5)


    def test_exact_boundary_search(self):
        """
        Test searching for a position where there's a feature with interval [pos-1, pos+1].
        This tests the edge case mentioned in the requirements where we need to find
        features that contain the query position, even if they're centered around it.
        """
        # Create a counter to calculate max lengths
        counter = MsiSiteCounter(bin_size=1000)
        
        # Create three test cases:
        # 1. Feature at [10000-10002] - exactly 1 base before and after position 10001
        # 2. Feature at [20000-20002] - exactly 1 base before and after position 20001
        # 3. Feature at [30000-30002] - exactly 1 base before and after position 30001
        test_records = [
            MsiSiteRecord(
                repeat_unit_length=1,
                repeat_unit_binary=1,
                repeat_times=2,
                left_flank_binary=1,
                right_flank_binary=1,
                repeat_unit_bases="A",
                left_flank_bases="AAAAA",
                right_flank_bases="AAAAA",
                chromosome="chr1",
                location=10000  # Will span 10000-10002
            ),
            MsiSiteRecord(
                repeat_unit_length=1,
                repeat_unit_binary=1,
                repeat_times=2,
                left_flank_binary=1,
                right_flank_binary=1,
                repeat_unit_bases="A",
                left_flank_bases="AAAAA",
                right_flank_bases="AAAAA",
                chromosome="chr1",
                location=20000  # Will span 20000-20002
            ),
            MsiSiteRecord(
                repeat_unit_length=1,
                repeat_unit_binary=1,
                repeat_times=2,
                left_flank_binary=1,
                right_flank_binary=1,
                repeat_unit_bases="A",
                left_flank_bases="AAAAA",
                right_flank_bases="AAAAA",
                chromosome="chr1",
                location=30000  # Will span 30000-30002
            )
        ]
        
        # Add records to the counter to calculate max lengths
        for record in test_records:
            counter.add(record)
        
        # Get max lengths from the counter
        max_lengths = counter.get_max_lengths("chr1")
        
        # Create a store with pre-calculated max lengths
        store = MsiChromosomeStore(chromosome="chr1", feature_count=3, max_lengths_by_bin=max_lengths, bin_size=1000)
        
        # Add the records to the store
        with store:
            for record in test_records:
                store.add(record)
        
        # Test searching for positions exactly in the middle of each feature
        positions_to_test = [10001, 20001, 30001]
        
        for pos in positions_to_test:
            # Get features
            features = store.get_by_position(pos)
            self.assertEqual(len(features), 1, f"Should find exactly one feature at position {pos}")
            self.assertTrue(features[0].start <= pos < features[0].end,
                           f"Feature at position {pos} should contain the position")
            
        # Test searching for positions at the exact boundaries
        boundary_tests = [
            (10000, 1),  # Start boundary - should find the feature
            (10002, 0),  # End boundary - should NOT find the feature (end is exclusive)
            (20000, 1),  # Start boundary - should find the feature
            (20002, 0),  # End boundary - should NOT find the feature (end is exclusive)
        ]
        
        for pos, expected_count in boundary_tests:
            features = store.get_by_position(pos)
            self.assertEqual(len(features), expected_count,
                            f"Should find {expected_count} features at boundary position {pos}")
            
        # Test the specific case mentioned in requirements:
        # When searching for position 'pos' and there is a feature with interval [pos-1, pos+1]
        # Create a feature at [15000-15002] and search for position 15001
        specific_record = MsiSiteRecord(
            repeat_unit_length=1,
            repeat_unit_binary=1,
            repeat_times=2,
            left_flank_binary=1,
            right_flank_binary=1,
            repeat_unit_bases="A",
            left_flank_bases="AAAAA",
            right_flank_bases="AAAAA",
            chromosome="chr1",
            location=15000  # Will span 15000-15002
        )
        
        # Create a counter for this specific test
        specific_counter = MsiSiteCounter(bin_size=1000)
        specific_counter.add(specific_record)
        
        # Get max lengths from the counter
        specific_max_lengths = specific_counter.get_max_lengths("chr1")
        
        # Create a new store for this specific test with pre-calculated max lengths
        specific_store = MsiChromosomeStore(chromosome="chr1", feature_count=1, max_lengths_by_bin=specific_max_lengths, bin_size=1000)
        specific_store.index_build_start()
        specific_store.add(specific_record)
        specific_store.index_build_end()
        
        
        # Search for position 15001 (middle of the feature)
        features = specific_store.get_by_position(15001)
        self.assertEqual(len(features), 1, "Should find the feature that spans [15000-15002] when searching for position 15001")
        self.assertEqual(features[0].start, 15000)
        self.assertEqual(features[0].end, 15002)
        
    def test_large_dataset_search(self):
        """Test search performance with a large dataset of random MSI sites."""
        # Number of random MSI sites to generate
        num_sites = 10000
        
        # Generate random MSI sites
        random.seed(42)  # For reproducibility
        records = []
        
        # Generate sites with positions from 1 to 1,000,000
        positions = sorted([random.randint(1, 1000000) for _ in range(num_sites)])
        
        for pos in positions:
            # Generate random repeat unit (1-6 bases)
            repeat_length = random.randint(1, 6)
            repeat_bases = ''.join(random.choices('ACGT', k=repeat_length))
            
            # Generate random repeat times (2-10)
            repeat_times = random.randint(2, 10)
            
            # Generate random flanking sequences
            left_flank = ''.join(random.choices('ACGT', k=5))
            right_flank = ''.join(random.choices('ACGT', k=5))
            
            # Create the record
            record = MsiSiteRecord(
                repeat_unit_length=repeat_length,
                repeat_unit_binary=random.randint(1, 1000),  # Not important for this test
                repeat_times=repeat_times,
                left_flank_binary=random.randint(1, 1000),  # Not important for this test
                right_flank_binary=random.randint(1, 1000),  # Not important for this test
                repeat_unit_bases=repeat_bases,
                left_flank_bases=left_flank,
                right_flank_bases=right_flank,
                chromosome="chr1",
                location=pos
            )
            records.append(record)
        
        # Create counter and store
        counter = MsiSiteCounter(bin_size=10000)
        for record in records:
            counter.add(record)
            
        store = MsiChromosomeStore(
            chromosome="chr1",
            feature_count=counter.get_count("chr1"),
            max_lengths_by_bin=counter.get_max_lengths("chr1"),
            bin_size=10000
        )
        
        # Add records to the store
        store.index_build_start()
        for record in records:
            store.add(record)
        store.index_build_end()
        
        # Test 1: Search for specific positions
        test_positions = [
            10000,    # Beginning of dataset
            500000,   # Middle of dataset
            990000    # End of dataset
        ]
        
        for pos in test_positions:
            start_time = time.time()
            features = store.get_by_position(pos)
            elapsed = time.time() - start_time
            
            # Verify results with brute force search
            expected = [r for r in records if r.start <= pos < r.end]
            self.assertEqual(len(features), len(expected))
        
        # Test 2: Search for intervals of different sizes
        test_intervals = [
            (10000, 10100),      # Small interval
            (100000, 110000),    # Medium interval
            (400000, 600000)     # Large interval
        ]
        
        for start, end in test_intervals:
            start_time = time.time()
            features = store.get_by_interval(start, end)
            elapsed = time.time() - start_time
            
            # Verify results with brute force search
            expected = [r for r in records if r.start < end and r.end > start]
            self.assertEqual(len(features), len(expected))
            
    def test_worst_case_search(self):
        """Test search performance in worst-case scenarios."""
        # Create a store with features that all have very long spans
        num_sites = 1000
        bin_size = 10000
        
        # Generate sites with positions from 1 to 1,000,000
        # but all with very long repeat units and many repeats
        records = []
        
        for i in range(num_sites):
            pos = i * bin_size  # Spread out evenly
            
            # Long repeat unit (6 bases) with many repeats (50)
            # This creates features that span multiple bins
            repeat_bases = ''.join(random.choices('ACGT', k=6))
            repeat_times = 50
            
            # Create the record
            record = MsiSiteRecord(
                repeat_unit_length=6,
                repeat_unit_binary=random.randint(1, 1000),
                repeat_times=repeat_times,
                left_flank_binary=random.randint(1, 1000),
                right_flank_binary=random.randint(1, 1000),
                repeat_unit_bases=repeat_bases,
                left_flank_bases='AAAAA',
                right_flank_bases='TTTTT',
                chromosome="chr1",
                location=pos
            )
            records.append(record)
        
        # Create counter and store
        counter = MsiSiteCounter(bin_size=bin_size)
        for record in records:
            counter.add(record)
            
        store = MsiChromosomeStore(
            chromosome="chr1",
            feature_count=counter.get_count("chr1"),
            max_lengths_by_bin=counter.get_max_lengths("chr1"),
            bin_size=bin_size
        )
        
        # Add records to the store
        store.index_build_start()
        for record in records:
            store.add(record)
        store.index_build_end()
        
        # Test position search in a bin with many overlapping features
        # This tests the backward lookup mechanism
        test_position = 500000
        
        start_time = time.time()
        features = store.get_by_position(test_position)
        elapsed = time.time() - start_time
        
        # Verify results
        expected = [r for r in records if r.start <= test_position < r.end]
        self.assertEqual(len(features), len(expected))
        
        # Test interval search with a large interval that contains many features
        start_time = time.time()
        features = store.get_by_interval(400000, 600000)
        elapsed = time.time() - start_time
        
        # Verify results
        expected = [r for r in records if r.start < 600000 and r.end > 400000]
        self.assertEqual(len(features), len(expected))
        
    def test_realistic_distribution(self):
        """Test with a more realistic distribution of MSI sites."""
        # Create a store with a realistic distribution of MSI sites
        # MSI sites tend to cluster in certain regions
        num_sites = 10000
        
        # Generate clusters of positions
        clusters = [
            (10000, 20000, 2000),    # (start, end, count) - dense cluster
            (100000, 300000, 3000),  # medium density cluster
            (500000, 900000, 5000)   # sparse cluster
        ]
        
        records = []
        
        for cluster_start, cluster_end, count in clusters:
            # Generate positions within this cluster
            positions = sorted([random.randint(cluster_start, cluster_end) for _ in range(count)])
            
            for pos in positions:
                # Generate random repeat unit (1-4 bases, with 1-2 being more common)
                repeat_length = random.choices([1, 2, 3, 4], weights=[0.4, 0.4, 0.1, 0.1])[0]
                repeat_bases = ''.join(random.choices('ACGT', k=repeat_length))
                
                # Generate random repeat times (3-15, with lower numbers being more common)
                repeat_times = random.choices(range(3, 16), weights=[0.2, 0.2, 0.15, 0.1, 0.1, 0.05, 0.05, 0.05, 0.03, 0.02, 0.02, 0.01, 0.01])[0]
                
                # Create the record
                record = MsiSiteRecord(
                    repeat_unit_length=repeat_length,
                    repeat_unit_binary=random.randint(1, 1000),
                    repeat_times=repeat_times,
                    left_flank_binary=random.randint(1, 1000),
                    right_flank_binary=random.randint(1, 1000),
                    repeat_unit_bases=repeat_bases,
                    left_flank_bases=''.join(random.choices('ACGT', k=5)),
                    right_flank_bases=''.join(random.choices('ACGT', k=5)),
                    chromosome="chr1",
                    location=pos
                )
                records.append(record)
        
        # Create counter and store
        counter = MsiSiteCounter(bin_size=10000)
        for record in records:
            counter.add(record)
            
        store = MsiChromosomeStore(
            chromosome="chr1",
            feature_count=counter.get_count("chr1"),
            max_lengths_by_bin=counter.get_max_lengths("chr1"),
            bin_size=10000
        )
        
        # Add records to the store
        store.index_build_start()
        for record in records:
            store.add(record)
        store.index_build_end()
        
        # Test searches in different density regions
        test_regions = [
            (15000, 15100),    # Dense cluster
            (200000, 200100),  # Medium density
            (700000, 700100)   # Sparse region
        ]
        
        for start, end in test_regions:
            # Position search
            start_time = time.time()
            pos_features = store.get_by_position(start)
            pos_elapsed = time.time() - start_time
            
            # Interval search
            start_time = time.time()
            int_features = store.get_by_interval(start, end)
            int_elapsed = time.time() - start_time
            
            
            # Verify results
            expected_pos = [r for r in records if r.start <= start < r.end]
            expected_int = [r for r in records if r.start < end and r.end > start]
            
            self.assertEqual(len(pos_features), len(expected_pos))
            self.assertEqual(len(int_features), len(expected_int))


if __name__ == "__main__":
    unittest.main()